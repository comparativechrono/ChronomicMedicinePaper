#!/usr/bin/env python3
"""Score a large random pool of common variants for SNPsnap-style matched-null sampling.

For each pool variant (1000G EUR, MAF >= 1%, > 1 Mb from any chronotype lead) we compute
exactly the same features and gene sets as for the observed leads:
  MAF, number of LD buddies (r2 >= 0.5 within 1 Mb), LD intervals at r2 >= 0.5/0.6/0.8,
  protein-coding genes overlapping each interval (+ flank), nearest gene distance.

Outputs (data/interim/null/):
  pool_features.parquet
  pool_genes_r{50,60,80}.npz      sparse pool x gene indicator matrices (column order = gene_universe.tsv)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, get_logger, load_config, chrom_norm  # noqa: E402
from chronomic import loci as L  # noqa: E402

log = get_logger("03_null_pool")
cfg = load_config()
OUT = INTERIM / "null"
OUT.mkdir(parents=True, exist_ok=True)
BFILE = INTERIM / "ldref" / "eur_maf01"
FLANK = cfg["loci"]["flank_kb"] * 1000
EXCLUDE_BP = 1_000_000
rng = np.random.default_rng(cfg["seed"])

bim = pd.read_csv(BFILE.with_suffix(".bim"), sep="\t", header=None, usecols=[0, 1, 3],
                  names=["chrom", "snp", "bp"], dtype={"chrom": str, "snp": str, "bp": np.int64},
                  engine="pyarrow")
bim["chrom"] = bim.chrom.map(chrom_norm)
frq = pd.read_csv(BFILE.with_suffix(".afreq"), sep="\t", usecols=["ID", "ALT_FREQS"], engine="pyarrow")
bim["maf"] = np.minimum(frq.ALT_FREQS.values, 1 - frq.ALT_FREQS.values)  # same row order as .bim
bim = bim[bim.chrom.isin([str(i) for i in range(1, 23)] + ["X"])]
bim = bim.drop_duplicates("snp", keep=False)

leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t", dtype={"chrom": str})
excl = np.zeros(len(bim), dtype=bool)
for c, sub in leads.groupby("chrom"):
    m = (bim.chrom.values == c)
    pos = bim.bp.values[m]
    lp = np.sort(sub.pos.values)
    j = np.searchsorted(lp, pos)
    d = np.full(pos.shape, np.inf)
    for k in (j - 1, j):
        ok = (k >= 0) & (k < len(lp))
        d[ok] = np.minimum(d[ok], np.abs(pos[ok] - lp[k[ok]]))
    excl[np.where(m)[0][d <= EXCLUDE_BP]] = True
cand = bim[~excl]
n_pool = min(cfg["null_model"]["pool_size"], len(cand))
pool = cand.iloc[np.sort(rng.choice(len(cand), size=n_pool, replace=False))].reset_index(drop=True)
log.info("candidate variants: %d (excluded %d within 1 Mb of leads); pool: %d", len(cand), excl.sum(), len(pool))
del bim, frq, cand

genes = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str})
gi = L.GeneIndex(genes)

feats, mats = [], {lv: [] for lv in L.R2_LEVELS}
for c in [str(i) for i in range(1, 23)] + ["X"]:
    sub = pool[pool.chrom == c].reset_index(drop=True)
    if sub.empty:
        continue
    ldp = L.plink_ld(BFILE, sub.snp.tolist(), OUT / f"tmp_chr{c}", window_kb=cfg["ld"]["window_kb"],
                     min_r2=0.5, chrom=c, dprime=False)
    ld = pd.read_csv(ldp, sep=r"\s+", usecols=["SNP_A", "BP_A", "SNP_B", "BP_B", "R2"])
    s = L.summarise_ld(ld, sub.set_index("snp").bp)
    s = s.reset_index(names="snp")
    s["chrom"] = c
    s["maf"] = sub.maf.values
    near = gi.nearest(s.chrom, s.bp)
    s["nearest_dist"] = near.nearest_dist.values
    for lv in L.R2_LEVELS:
        tag = f"{int(round(lv * 100)):02d}"
        M = gi.overlap_matrix(s.chrom, s[f"start_r{tag}"], s[f"end_r{tag}"], FLANK)
        s[f"n_genes_r{tag}"] = np.asarray(M.sum(axis=1)).ravel()
        mats[lv].append(M)
    feats.append(s)
    for f in OUT.glob(f"tmp_chr{c}.*"):
        f.unlink()
    log.info("chr%s: %d pool variants, %d LD pairs", c, len(sub), len(ld))

F = pd.concat(feats, ignore_index=True)
F.to_parquet(OUT / "pool_features.parquet", index=False)
for lv, ms in mats.items():
    tag = f"{int(round(lv * 100)):02d}"
    sparse.save_npz(OUT / f"pool_genes_r{tag}.npz", sparse.vstack(ms).tocsr())
log.info("pool features written: %d variants; median buddies %.0f; median genes (r2>=0.6) %.1f",
         len(F), F.n_buddies.median(), F.n_genes_r60.median())
