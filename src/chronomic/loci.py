"""LD-based locus definition and gene mapping (shared by observed and null loci)."""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .common import BIN

R2_LEVELS = (0.5, 0.6, 0.8)


def plink_ld(bfile: Path, snp_ids: list[str], out_prefix: Path, window_kb: int, min_r2: float,
             chrom: str | None = None, threads: int = 2, dprime: bool = False) -> Path:
    """Pairwise LD between each listed SNP and every panel variant within the window."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    lst = out_prefix.with_suffix(".snplist")
    lst.write_text("\n".join(snp_ids) + "\n")
    cmd = [BIN / "plink", "--bfile", bfile, "--r2"] + (["dprime"] if dprime else []) + [
        "--ld-snp-list", lst, "--ld-window-kb", window_kb, "--ld-window", 999999,
        "--ld-window-r2", min_r2, "--threads", threads, "--memory", 4000, "--out", out_prefix]
    if chrom is not None:
        cmd += ["--chr", chrom]
    subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL)
    return out_prefix.with_suffix(".ld")


def read_ld(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", dtype={"CHR_A": str, "CHR_B": str})
    return df


def summarise_ld(ld: pd.DataFrame, index_pos: pd.Series, levels=R2_LEVELS, buddy_r2: float = 0.5) -> pd.DataFrame:
    """Per index SNP: number of LD buddies and interval bounds at each r2 level.

    `index_pos` maps SNP id -> bp (so SNPs without any partner still get a 1-bp interval).
    """
    out = pd.DataFrame(index=index_pos.index)
    out["bp"] = index_pos.values
    ld = ld[ld.SNP_A != ld.SNP_B]
    g = ld[ld.R2 >= buddy_r2].groupby("SNP_A").size()
    out["n_buddies"] = g.reindex(out.index).fillna(0).astype(int)
    for lv in levels:
        sub = ld[ld.R2 >= lv]
        lo = sub.groupby("SNP_A").BP_B.min().reindex(out.index)
        hi = sub.groupby("SNP_A").BP_B.max().reindex(out.index)
        tag = f"{int(round(lv * 100)):02d}"
        out[f"start_r{tag}"] = np.fmin(lo.fillna(out.bp), out.bp).astype(int)
        out[f"end_r{tag}"] = np.fmax(hi.fillna(out.bp), out.bp).astype(int)
    return out


class GeneIndex:
    """Fast overlap / nearest-gene queries against a gene table (one chromosome at a time)."""

    def __init__(self, genes: pd.DataFrame):
        self.genes = genes.reset_index(drop=True)
        self.col = {g: i for i, g in enumerate(self.genes.gene_id)}
        self.by_chrom = {c: sub for c, sub in self.genes.groupby("chrom")}

    def overlap_matrix(self, chrom: pd.Series, start: pd.Series, end: pd.Series, flank: int) -> sparse.csr_matrix:
        """Sparse (n_intervals x n_genes) indicator of gene bodies overlapping [start-flank, end+flank]."""
        n = len(chrom)
        rows, cols = [], []
        chrom = chrom.astype(str).values
        s = start.values.astype(np.int64) - flank
        e = end.values.astype(np.int64) + flank
        idx_all = np.arange(n)
        for c in np.unique(chrom):
            g = self.by_chrom.get(c)
            if g is None:
                continue
            gi = g.index.values
            gs, ge = g.start.values, g.end.values
            sel = idx_all[chrom == c]
            for chunk in np.array_split(sel, max(1, len(sel) // 2000)):
                hit = (gs[None, :] <= e[chunk, None]) & (ge[None, :] >= s[chunk, None])
                r, k = np.nonzero(hit)
                rows.append(chunk[r])
                cols.append(gi[k])
        if rows:
            rows = np.concatenate(rows)
            cols = np.concatenate(cols)
        else:
            rows = cols = np.array([], dtype=int)
        data = np.ones(len(rows), dtype=bool)
        return sparse.csr_matrix((data, (rows, cols)), shape=(n, len(self.genes)), dtype=bool)

    def nearest(self, chrom: pd.Series, pos: pd.Series) -> pd.DataFrame:
        """Nearest gene body (distance 0 if the position lies inside a gene)."""
        out_gene = np.empty(len(chrom), dtype=object)
        out_dist = np.full(len(chrom), np.nan)
        chrom = chrom.astype(str).values
        p = pos.values.astype(np.int64)
        idx_all = np.arange(len(chrom))
        for c in np.unique(chrom):
            g = self.by_chrom.get(c)
            if g is None:
                continue
            gs, ge, gid = g.start.values, g.end.values, g.gene_id.values
            sel = idx_all[chrom == c]
            for chunk in np.array_split(sel, max(1, len(sel) // 2000)):
                d = np.maximum(0, np.maximum(gs[None, :] - p[chunk, None], p[chunk, None] - ge[None, :]))
                j = d.argmin(axis=1)
                out_gene[chunk] = gid[j]
                out_dist[chunk] = d[np.arange(len(chunk)), j]
        return pd.DataFrame({"nearest_gene": out_gene, "nearest_dist": out_dist})


def merge_intervals(df: pd.DataFrame, start_col: str, end_col: str, flank: int) -> pd.Series:
    """Assign a region id to rows whose padded intervals overlap (per chromosome)."""
    region = pd.Series(index=df.index, dtype=object)
    k = 0
    for c, sub in df.sort_values(["chrom", start_col]).groupby("chrom", sort=False):
        cur_end = -1
        for i, r in sub.iterrows():
            s, e = r[start_col] - flank, r[end_col] + flank
            if s > cur_end:
                k += 1
                cur_end = e
            else:
                cur_end = max(cur_end, e)
            region[i] = f"R{k:03d}"
    return region
