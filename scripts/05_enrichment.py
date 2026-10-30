#!/usr/bin/env python3
"""Locus-level enrichment of Mendelian / phenotype gene sets at chronotype loci.

Each chronotype locus is compared with matched null loci drawn from a pool of random common
variants (matched on MAF, LD buddies, number of genes in the LD interval and distance to the
nearest gene). Statistic: number of loci containing >= 1 gene of a set.

Analyses: primary (r2 >= 0.6), r2 >= 0.5, r2 >= 0.8, excluding MHC and 17q21.31, nearest gene.
Outputs: results/tables/enrichment_<analysis>.tsv, enrichment_all.tsv, locus_set_hits.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, TABLES, bh_fdr, get_logger, load_config, save_json  # noqa: E402
from chronomic import enrich as E  # noqa: E402
from chronomic import loci as L  # noqa: E402

log = get_logger("05_enrichment")
cfg = load_config()
rng = np.random.default_rng(cfg["seed"])
N_PERM = cfg["null_model"]["n_permutations"]

genes = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str})
gcol = {g: i for i, g in enumerate(genes.gene_id)}
gs = pd.read_csv(INTERIM / "genesets" / "genesets.tsv", sep="\t")
summary = pd.read_csv(INTERIM / "genesets" / "geneset_summary.tsv", sep="\t")
set_ids = summary.loc[summary.tested, "set_id"].tolist()
set_ind = np.zeros((len(genes), len(set_ids)), dtype=bool)
for j, sid in enumerate(set_ids):
    for g in gs.loc[gs.set_id == sid, "gene_id"]:
        set_ind[gcol[g], j] = True
ref_col = set_ids.index("MEND_ANY")

leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t", dtype={"chrom": str})
lg = pd.read_csv(INTERIM / "loci" / "locus_genes.tsv", sep="\t", dtype={"r2_level": str})
pool = pd.read_parquet(INTERIM / "null" / "pool_features.parquet")
pool_mats = {t: sparse.load_npz(INTERIM / "null" / f"pool_genes_r{t}.npz") for t in ("50", "60", "80")}
gi = L.GeneIndex(genes)
pool_near = gi.nearest(pool.chrom, pool.bp)
pool["nearest_gene"] = pool_near.nearest_gene.values


def obs_matrix(tag: str) -> sparse.csr_matrix:
    sub = lg[lg.r2_level == tag]
    row = {ln: i for i, ln in enumerate(leads.locus_no)}
    r = sub.locus_no.map(row).values
    c = sub.gene_id.map(gcol).values
    return sparse.csr_matrix((np.ones(len(r), bool), (r, c)), shape=(len(leads), len(genes)))


def nearest_matrix(gene_ids: pd.Series) -> sparse.csr_matrix:
    c = gene_ids.map(gcol).values
    r = np.arange(len(gene_ids))
    ok = ~pd.isna(c)
    return sparse.csr_matrix((np.ones(ok.sum(), bool), (r[ok], c[ok].astype(int))),
                             shape=(len(gene_ids), len(genes)))


def in_window(chrom, start, end, w):
    return (chrom.astype(str) == w["chrom"]) & (end >= w["start"]) & (start <= w["end"])


excl_windows = [cfg["loci"]["mhc"], cfg["loci"]["inv17q21"]]
ANALYSES = [
    ("primary", "60", "ld", False),
    ("r2_0.5", "50", "ld", False),
    ("r2_0.8", "80", "ld", False),
    ("exclude_MHC_17q21", "60", "ld", True),
    ("nearest_gene", "60", "nearest", False),
]

fam = summary.set_index("set_id").family
lab = summary.set_index("set_id").label
all_res = []
for name, tag, mapping, exclude in ANALYSES:
    obs_df = leads.copy()
    pool_df = pool.copy()
    if mapping == "ld":
        Mo, Mp = obs_matrix(tag), pool_mats[tag]
    else:
        Mo, Mp = nearest_matrix(obs_df.nearest_gene), nearest_matrix(pool_df.nearest_gene)
    keep_o = np.ones(len(obs_df), bool)
    keep_p = np.ones(len(pool_df), bool)
    if exclude:
        for w in excl_windows:
            keep_o &= ~in_window(obs_df.chrom, obs_df[f"start_r{tag}"], obs_df[f"end_r{tag}"], w).values
            keep_p &= ~in_window(pool_df.chrom, pool_df[f"start_r{tag}"], pool_df[f"end_r{tag}"], w).values
    obs_df, Mo = obs_df[keep_o].reset_index(drop=True), Mo[np.where(keep_o)[0]]
    pool_df, Mp = pool_df[keep_p].reset_index(drop=True), Mp[np.where(keep_p)[0]]

    ncol = f"n_genes_r{tag}"
    cands, level = E.match_candidates(obs_df, pool_df, ncol, q=cfg["null_model"]["maf_bins"])
    obs_hits = E.set_hit_matrix(Mo, set_ind).sum(axis=0)
    pool_hits = E.set_hit_matrix(Mp, set_ind)
    null = E.permute(cands, pool_hits, N_PERM, rng)
    res = E.summarise(obs_hits, null, set_ids, ref_col=ref_col)
    probs = E.hit_probabilities(cands, pool_hits)
    res["exp_loci_analytic"] = probs.sum(axis=0)
    res["p_enrich_analytic"] = [E.poisson_binomial_sf(probs[:, j], int(obs_hits[j])) for j in range(len(set_ids))]
    res["fold_analytic"] = np.asarray(obs_hits, float) / res["exp_loci_analytic"].values
    np.save(INTERIM / "null" / f"null_hits_{name}.npy", null)
    np.save(INTERIM / "null" / f"obs_hits_{name}.npy", obs_hits)
    res.insert(0, "analysis", name)
    res.insert(2, "family", res.set_id.map(fam))
    res.insert(3, "label", res.set_id.map(lab))
    res["n_loci"] = len(obs_df)
    res["q_within_family"] = np.nan
    res["q_analytic_within_family"] = np.nan
    res["q_relative_within_family"] = np.nan
    for f, idx in res.groupby("family").groups.items():
        res.loc[idx, "q_within_family"] = bh_fdr(res.loc[idx, "p_enrich"].values)
        res.loc[idx, "q_analytic_within_family"] = bh_fdr(res.loc[idx, "p_enrich_analytic"].values)
        if "p_relative" in res:
            res.loc[idx, "q_relative_within_family"] = bh_fdr(res.loc[idx, "p_relative"].values)
    res.to_csv(TABLES / f"enrichment_{name}.tsv", sep="\t", index=False)
    all_res.append(res)
    lvl = pd.Series(level).value_counts().to_dict()
    log.info("[%s] loci=%d matching=%s MEND_ANY obs=%d exp=%.1f fold=%.2f p=%.4g",
             name, len(obs_df), lvl, obs_hits[ref_col], null[:, ref_col].mean(),
             obs_hits[ref_col] / null[:, ref_col].mean(), res.loc[ref_col, "p_enrich"])
    if name == "primary":
        np.save(INTERIM / "null" / "primary_null_hits.npy", null)
        hits = E.set_hit_matrix(Mo, set_ind)
        rows = []
        Mo_l = Mo.tolil()
        for i in range(len(obs_df)):
            g_idx = Mo_l.rows[i]
            for j, sid in enumerate(set_ids):
                if hits[i, j]:
                    members = [genes.symbol.iloc[k] for k in g_idx if set_ind[k, j]]
                    rows.append((obs_df.locus_no.iloc[i], obs_df.lead_rsid.iloc[i], sid, ";".join(members)))
        pd.DataFrame(rows, columns=["locus_no", "lead_rsid", "set_id", "genes"]).to_csv(
            TABLES / "locus_set_hits.tsv", sep="\t", index=False)
        save_json({"matching_levels": lvl, "n_candidates_median": float(np.median([len(c) for c in cands]))},
                  TABLES / "enrichment_primary_matching.json")

pd.concat(all_res).to_csv(TABLES / "enrichment_all.tsv", sep="\t", index=False)
log.info("done")

# ------------------------------------------------------------------ pre-specified relative contrasts
CONTRASTS = [("HPO_HP0002360", "MEND_ANY"), ("HPO_HP0002360", "NDD"), ("HPO_CIRCADIAN", "MEND_ANY"),
             ("HPO_HP0001657", "MEND_ANY"), ("HPO_HP0001639", "MEND_ANY"), ("EXP_INTD_MONOAMINE", "MEND_ANY"),
             ("INH_AD", "MEND_ANY"), ("HPO_HP0001263", "MEND_ANY"), ("HPO_HP0002360", "INH_AD")]
rel = []
for name, *_ in ANALYSES:
    null = np.load(INTERIM / "null" / f"null_hits_{name}.npy")
    obs = np.load(INTERIM / "null" / f"obs_hits_{name}.npy")
    for num, den in CONTRASTS:
        i, j = set_ids.index(num), set_ids.index(den)
        fo = obs[i] / obs[j] if obs[j] else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            fn = null[:, i] / null[:, j]
        fn = fn[np.isfinite(fn)]
        rel.append({"analysis": name, "numerator": num, "denominator": den, "obs_numerator": int(obs[i]),
                    "obs_denominator": int(obs[j]), "obs_fraction": fo, "null_mean_fraction": fn.mean(),
                    "p_relative": (1 + (fn >= fo).sum()) / (1 + fn.size)})
pd.DataFrame(rel).to_csv(TABLES / "enrichment_relative.tsv", sep="\t", index=False)
log.info("relative contrasts written")
