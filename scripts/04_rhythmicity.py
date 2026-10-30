#!/usr/bin/env python3
"""Rhythmic / day-night gene expression annotation.

1. Mouse 12-tissue circadian atlas (Zhang et al. 2014 PNAS; GEO GSE54650; constant darkness,
   CT18-CT64 every 2 h). Per tissue and gene: harmonic regression (24-h cosinor) F-test,
   BH-FDR within tissue; rhythmic = q < 0.05 and peak-to-trough >= 0.2 log2 units.
   Mouse genes mapped to human by one-to-one MGI homology (HGNC IDs).
2. Human GTEx day-night differential expression (Wucher et al. 2023 PLoS Biol, 46 tissues):
   day-night gene = FDR < 0.05 and permutation-robust (proportion of permuted p <= 0.05 with
   |logFC| >= 0.1 at most 5%).

Outputs: data/interim/rhythm/gene_rhythmicity.tsv, results/tables/rhythmicity_gene_level.tsv
         and two extra gene sets appended to data/interim/genesets (family 'rhythmicity').
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic import annot  # noqa: E402
from chronomic.common import INTERIM, RAW, TABLES, bh_fdr, get_logger, save_json  # noqa: E402

log = get_logger("04_rhythmicity")
OUT = INTERIM / "rhythm"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ mouse atlas
with gzip.open(RAW / "rhythm" / "GSE54650_series_matrix.txt.gz", "rt") as fh:
    meta, lines = {}, []
    for line in fh:
        if line.startswith("!Sample_title"):
            meta["title"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!series_matrix_table_begin"):
            break
    expr = pd.read_csv(fh, sep="\t", index_col=0, comment="!")
titles = meta["title"]
tissue = [t.split("_CT")[0] for t in titles]
ct = np.array([float(t.split("_CT")[1]) for t in titles])
expr.columns = titles

ann = pd.read_csv(RAW / "rhythm" / "GPL6246.annot.gz", sep="\t", skiprows=27, dtype=str,
                  usecols=["ID", "Gene symbol", "Gene ID"]).dropna()
ann = ann[~ann["Gene ID"].str.contains("///")]
probe2sym = dict(zip(ann.ID.astype(str), ann["Gene ID"]))  # mouse Entrez Gene ID (robust to symbol renames)
expr.index = expr.index.astype(str)
expr = expr[expr.index.isin(probe2sym)]

w = 2 * np.pi / 24.0
rows = []
for tis in sorted(set(tissue)):
    cols = [i for i, t in enumerate(tissue) if t == tis]
    Y = expr.iloc[:, cols].values  # probes x 24
    t = ct[cols]
    X = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    fit = (X @ beta).T
    rss1 = ((Y - fit) ** 2).sum(axis=1)
    rss0 = ((Y - Y.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    n = len(t)
    F = ((rss0 - rss1) / 2) / (rss1 / (n - 3))
    from scipy.stats import f as fdist
    p = fdist.sf(F, 2, n - 3)
    amp = 2 * np.hypot(beta[1], beta[2])  # peak-to-trough (log2)
    phase = (np.arctan2(beta[2], beta[1]) / w) % 24
    df = pd.DataFrame({"probe": expr.index, "tissue": tis, "p": p, "ptt_log2": amp, "phase_ct": phase,
                       "mean": Y.mean(axis=1)})
    df["q"] = bh_fdr(df.p.values)
    rows.append(df)
M = pd.concat(rows)
M["mouse_symbol"] = M.probe.map(probe2sym)
M = M.sort_values("p").drop_duplicates(["mouse_symbol", "tissue"])  # best probe per gene/tissue
M["rhythmic"] = (M.q < 0.05) & (M.ptt_log2 >= 0.2)
log.info("mouse atlas: %d tissues; rhythmic genes per tissue: %s", M.tissue.nunique(),
         M[M.rhythmic].groupby("tissue").size().to_dict())

hom = pd.read_csv(RAW / "rhythm" / "HOM_MouseHumanSequence.rpt", sep="\t", dtype=str)
mm = hom[hom["NCBI Taxon ID"] == "10090"][["DB Class Key", "EntrezGene ID"]].rename(columns={"EntrezGene ID": "Symbol"})
hs = hom[hom["NCBI Taxon ID"] == "9606"][["DB Class Key", "Symbol", "HGNC ID"]]
one2one = mm.groupby("DB Class Key").filter(lambda d: len(d) == 1).merge(
    hs.groupby("DB Class Key").filter(lambda d: len(d) == 1), on="DB Class Key", suffixes=("_mm", "_hs"))
idm = annot.id_maps()
one2one["gene_id"] = one2one["HGNC ID"].map(idm["hgnc"]).fillna(one2one.Symbol_hs.map(idm["symbol"]))
m2h = dict(zip(one2one.Symbol_mm, one2one.gene_id))
M["gene_id"] = M.mouse_symbol.map(m2h)
# several mouse genes can map to the same human Ensembl gene (e.g. via HGNC symbol fall-back): keep the best
# probe per human gene and tissue so that tissues are counted once
M = M.dropna(subset=["gene_id"]).sort_values("p").drop_duplicates(["gene_id", "tissue"]).reset_index(drop=True)
TIS = {"Adr": "adrenal", "Aor": "aorta", "BFat": "brown fat", "Bstm": "brainstem", "Cer": "cerebellum",
       "Hrt": "heart", "Hyp": "hypothalamus", "Kid": "kidney", "Liv": "liver", "Lun": "lung",
       "Mus": "skeletal muscle", "WFat": "white fat"}
M["tissue"] = M.tissue.map(TIS).fillna(M.tissue)
mouse_gene = M.groupby("gene_id").agg(
    mouse_n_tissues_tested=("tissue", "nunique"),
    mouse_n_rhythmic=("rhythmic", "sum"),
    mouse_rhythmic_tissues=("tissue", lambda s: ";".join(sorted(s[M.loc[s.index, "rhythmic"]]))),
    mouse_min_q=("q", "min"))
M.to_csv(OUT / "mouse_atlas_cosinor.tsv.gz", sep="\t", index=False)

# ------------------------------------------------------------------ human GTEx day-night
d = pd.read_csv(RAW / "rhythm" / "wucher2023_S1_daynight.tab.gz", sep="\t")
d.columns = ["tissue", "gene", "name", "logFC", "AveExpr", "t", "P", "FDR", "B", "perm", "perm_fc"]
d["gene_id"] = d.gene.str.split(".").str[0]
d["daynight"] = (d.FDR < 0.05) & (d.perm_fc <= 0.05)
d = d.sort_values("FDR").drop_duplicates(["gene_id", "tissue"])
human_gene = d.groupby("gene_id").agg(
    gtex_n_tissues_tested=("tissue", "nunique"),
    gtex_n_daynight=("daynight", "sum"),
    gtex_daynight_tissues=("tissue", lambda s: ";".join(sorted(s[d.loc[s.index, "daynight"]].str.replace("_", " ")))),
    gtex_min_fdr=("FDR", "min"))

R = mouse_gene.join(human_gene, how="outer")
R.index.name = "gene_id"
R["symbol"] = R.index.map(annot.symbol_of(R.index))
R.to_csv(OUT / "gene_rhythmicity.tsv", sep="\t")
universe = set(pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t").gene_id)
Ru = R[R.index.isin(universe)]
summ = {
    "mouse_genes_tested": int(Ru.mouse_n_tissues_tested.notna().sum()),
    "mouse_rhythmic_any": int((Ru.mouse_n_rhythmic >= 1).sum()),
    "mouse_rhythmic_ge3": int((Ru.mouse_n_rhythmic >= 3).sum()),
    "gtex_genes_tested": int(Ru.gtex_n_tissues_tested.notna().sum()),
    "gtex_daynight_any": int((Ru.gtex_n_daynight >= 1).sum()),
}
save_json(summ, TABLES / "rhythmicity_summary.json")
log.info("summary %s", summ)

# extra gene sets for the locus-level tests
gs_path = INTERIM / "genesets" / "genesets.tsv"
sm_path = INTERIM / "genesets" / "geneset_summary.tsv"
gs = pd.read_csv(gs_path, sep="\t")
sm = pd.read_csv(sm_path, sep="\t")
gs = gs[gs.family != "rhythmicity"]
sm = sm[sm.family != "rhythmicity"]
new = {
    "RHY_MOUSE_GE3": ("Rhythmic in >= 3 of 12 mouse tissues (Zhang 2014 atlas)",
                      set(Ru.index[Ru.mouse_n_rhythmic >= 3])),
    "RHY_GTEX_DAYNIGHT": ("Day-night differential expression in >= 1 GTEx tissue (Wucher 2023)",
                          set(Ru.index[Ru.gtex_n_daynight >= 1])),
}
add_rows, add_sum = [], []
for sid, (lab, genes) in new.items():
    add_rows += [(sid, "rhythmicity", lab, g) for g in sorted(genes)]
    add_sum.append((sid, "rhythmicity", lab, len(genes), True))
gs = pd.concat([gs, pd.DataFrame(add_rows, columns=gs.columns)])
sm = pd.concat([sm, pd.DataFrame(add_sum, columns=sm.columns)])
gs.to_csv(gs_path, sep="\t", index=False)
sm.to_csv(sm_path, sep="\t", index=False)
with open(INTERIM / "genesets" / "genesets.magma.txt", "w") as fh:
    for sid, sub in gs.groupby("set_id", sort=False):
        fh.write(sid + " " + " ".join(sub.gene_id) + "\n")
log.info("added rhythmicity gene sets: %s", {k: len(v[1]) for k, v in new.items()})
