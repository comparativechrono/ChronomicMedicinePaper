#!/usr/bin/env python3
"""Disease-group screen: enrichment of PanelApp disease gene panels at chronotype loci.

Complements the HPO phenotype sets of 05_enrichment.py with clinically defined disease groups
(Genomics England PanelApp green genes, retrieved by 01b_fetch_panelapp.py). Each panel is tested with
the same matched-null framework as 05_enrichment.py (primary null plus the interval and MHC/17q21.31
sensitivity analyses; exact Poisson-binomial p-values) and Benjamini-Hochberg FDR across panels.

Panels: public panels excluding super-panels and "Additional findings" lists; panels whose green-gene set
is identical to another panel's are tested once; 10-3,000 genes in the gene universe.
Outputs: results/tables/enrichment_panelapp.tsv, results/tables/panelapp_panels_tested.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic import annot  # noqa: E402
from chronomic import enrich as E  # noqa: E402
from chronomic.common import INTERIM, RAW, TABLES, bh_fdr, get_logger, load_config  # noqa: E402

log = get_logger("05b_disease_panels")
cfg = load_config()

genes = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str})
gcol = {g: i for i, g in enumerate(genes.gene_id)}
universe = set(genes.gene_id)

# ------------------------------------------------------------------ panels -> gene sets
P = pd.read_csv(RAW / "panelapp" / "panels.tsv", sep="\t")
G = pd.read_csv(RAW / "panelapp" / "green_genes.tsv", sep="\t", dtype=str)
idm = annot.id_maps()
G["gene_id"] = G.ensembl_id_grch37.where(G.ensembl_id_grch37.isin(universe))
G["gene_id"] = G.gene_id.fillna(G.hgnc_id.map(idm["hgnc"])).fillna(G.gene_symbol.map(idm["symbol"]))
G = G[G.gene_id.isin(universe)]
_types = P.types.fillna("").str.split(";")
excl = _types.map(lambda ts: "Super Panel" in ts or "Additional Findings" in ts)   # components are kept
P = P[~excl]
sets = {int(pid): frozenset(sub.gene_id) for pid, sub in G.groupby("panel_id")}
P = P[P.panel_id.isin(sets)].copy()
P["n_genes_universe"] = P.panel_id.map(lambda i: len(sets[i]))
P = P[P.n_genes_universe.between(10, 3000)].sort_values("panel_id")
P = P[~P.panel_id.map(sets).duplicated()]          # identical gene sets tested once
P.to_csv(TABLES / "panelapp_panels_tested.tsv", sep="\t", index=False)
pids = P.panel_id.tolist()
# members of the tested panels, and a MAGMA set file for the gene-level test in 07_magma_genesets.py
# (short set names: MAGMA truncates long ones)
pd.DataFrame([(pid, g) for pid in pids for g in sorted(sets[pid])], columns=["panel_id", "gene_id"]).to_csv(
    INTERIM / "genesets" / "panelapp_members.tsv", sep="\t", index=False)
with open(INTERIM / "genesets" / "panelapp.magma.txt", "w") as fh:
    for pid in pids:
        fh.write(f"PANEL_{pid} " + " ".join(sorted(sets[pid])) + "\n")
set_ind = np.zeros((len(genes), len(pids)), dtype=bool)
for j, pid in enumerate(pids):
    set_ind[[gcol[g] for g in sets[pid]], j] = True
log.info("panels tested: %d (of %d public panels)", len(pids), len(excl))

# ------------------------------------------------------------------ loci and null pool (as in 05_enrichment)
leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t", dtype={"chrom": str})
lg = pd.read_csv(INTERIM / "loci" / "locus_genes.tsv", sep="\t", dtype={"r2_level": str})
pool = pd.read_parquet(INTERIM / "null" / "pool_features.parquet")
row = {ln: i for i, ln in enumerate(leads.locus_no)}


def obs_matrix(tag: str) -> sparse.csr_matrix:
    sub = lg[lg.r2_level == tag]
    return sparse.csr_matrix((np.ones(len(sub), bool), (sub.locus_no.map(row).values, sub.gene_id.map(gcol).values)),
                             shape=(len(leads), len(genes)))


def in_window(chrom, start, end, w):
    return (chrom.astype(str) == w["chrom"]) & (end >= w["start"]) & (start <= w["end"])


Mo = obs_matrix("60")
obs_hit_mat = E.set_hit_matrix(Mo, set_ind)
obs_hits = obs_hit_mat.sum(axis=0)
res = P[["panel_id", "panel_name", "disease_group", "disease_sub_group", "types", "n_genes_universe"]].reset_index(drop=True)
res["obs_loci"] = np.asarray(obs_hits, int)
excl_windows = [cfg["loci"]["mhc"], cfg["loci"]["inv17q21"]]
for name, tag, exclude in (("primary", "60", False), ("r2_0.5", "50", False), ("r2_0.8", "80", False),
                           ("exclude_MHC_17q21", "60", True)):
    mo, mp = obs_matrix(tag), sparse.load_npz(INTERIM / "null" / f"pool_genes_r{tag}.npz")
    od, pdf = leads, pool
    if exclude:
        ko = np.ones(len(leads), bool)
        kp = np.ones(len(pool), bool)
        for w in excl_windows:
            ko &= ~in_window(leads.chrom, leads[f"start_r{tag}"], leads[f"end_r{tag}"], w).values
            kp &= ~in_window(pool.chrom, pool[f"start_r{tag}"], pool[f"end_r{tag}"], w).values
        od, mo = leads[ko].reset_index(drop=True), mo[np.where(ko)[0]]
        pdf, mp = pool[kp].reset_index(drop=True), mp[np.where(kp)[0]]
    cands, _ = E.match_candidates(od, pdf, f"n_genes_r{tag}", q=cfg["null_model"]["maf_bins"])
    ph = E.set_hit_matrix(mp, set_ind)
    oh = E.set_hit_matrix(mo, set_ind).sum(axis=0)
    probs = E.hit_probabilities(cands, ph)
    exp = probs.sum(axis=0)
    p = np.array([E.poisson_binomial_sf(probs[:, j], int(oh[j])) for j in range(len(pids))])
    if name != "primary":
        res[f"obs_{name}"] = np.asarray(oh, int)
    res[f"exp_{name}"] = exp
    res[f"fold_{name}"] = np.asarray(oh, float) / exp
    res[f"p_{name}"] = p
    res[f"q_{name}"] = bh_fdr(p)
sym = genes.set_index("gene_id").symbol
res["genes_in_loci"] = [";".join(sorted({sym[genes.gene_id.iloc[k]] for i in np.where(obs_hit_mat[:, j])[0]
                                         for k in Mo[i].indices if set_ind[k, j]})) for j in range(len(pids))]
res = res.sort_values("p_primary")
res.to_csv(TABLES / "enrichment_panelapp.tsv", sep="\t", index=False)
log.info("\n%s", res.head(25)[["panel_name", "n_genes_universe", "obs_loci", "exp_primary", "p_primary", "q_primary",
                               "p_r2_0.8", "p_exclude_MHC_17q21"]].to_string())
