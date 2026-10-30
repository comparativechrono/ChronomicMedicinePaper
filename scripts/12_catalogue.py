#!/usr/bin/env python3
"""Catalogue of Mendelian disease genes at chronotype loci with all evidence layers.

One row per (locus, Mendelian gene) for genes in the r2 >= 0.6 interval (+10 kb):
diseases, inheritance, ClinGen class, phenotype flags, nearest/L2G status, MAGMA p, exome burden
p and direction, QTL colocalisation, LoF constraint, mouse/human rhythmicity.
An evidence count summarises five lines of gene-level support, which are not independent
(L2G >= 0.5 for the gene in a credible set at the locus, nearest gene to the lead, MAGMA Bonferroni,
exome burden p < 0.01, QTL colocalisation with eCAVIAR CLPP >= 0.1 for that gene).
Outputs: results/tables/chronomic_gene_catalogue.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic import annot  # noqa: E402
from chronomic.common import INTERIM, TABLES, get_logger  # noqa: E402

log = get_logger("12_catalogue")
leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t")
lg = pd.read_csv(INTERIM / "loci" / "locus_genes.tsv", sep="\t", dtype={"r2_level": str})
lg = lg[lg.r2_level == "60"]
gs = pd.read_csv(INTERIM / "genesets" / "genesets.tsv", sep="\t")
S = {s: set(sub.gene_id) for s, sub in gs.groupby("set_id")}
links = pd.read_csv(INTERIM / "genesets" / "gene_disease_links.tsv", sep="\t")
systems = pd.read_csv(INTERIM / "genesets" / "geneset_summary.tsv", sep="\t")
sys_ids = systems.loc[systems.family == "hpo_system", ["set_id", "label"]].values

cat = lg[lg.gene_id.isin(S["MEND_ANY"])].merge(
    leads[["locus_no", "lead_rsid", "chrom", "pos", "meta_p_from_z", "gene_context", "nearest_gene",
           "nearest_gene_symbol"]], on="locus_no")


def diseases(g):
    d = links[(links.gene_id == g) & links.source.isin(["OMIM", "Orphanet"])].drop_duplicates("disease_name")
    return "; ".join(str(x) for x in d.disease_name.dropna().head(6))


def disease_ids(g):
    d = links[(links.gene_id == g) & links.source.isin(["OMIM", "Orphanet"])]
    return ";".join(sorted(str(x) for x in set(d.disease_id.dropna()))[:10])


clingen = annot.clingen_validity()
cg = clingen[clingen.classification.isin(["Definitive", "Strong", "Moderate"])].groupby("gene_id").classification.apply(
    lambda s: ";".join(sorted(set(s))))
cat["diseases"] = cat.gene_id.map(diseases)
cat["disease_ids"] = cat.gene_id.map(disease_ids)
cat["inheritance"] = cat.gene_id.map(lambda g: ";".join(x for x, s in [("AD", "INH_AD"), ("AR", "INH_AR"), ("XL", "INH_XL")]
                                                    if g in S[s]))
cat["clingen"] = cat.gene_id.map(cg)
cat["organ_systems"] = cat.gene_id.map(lambda g: "; ".join(l for s, l in sys_ids if g in S[s]))
for flag, sid in [("sleep_disturbance", "HPO_HP0002360"), ("circadian_hpo", "HPO_CIRCADIAN"),
                  ("monoamine_intd", "EXP_INTD_MONOAMINE"), ("neurodevelopmental", "NDD"),
                  ("prolonged_QT", "HPO_HP0001657"), ("hypertrophic_cm", "HPO_HP0001639"),
                  ("go_circadian", "GO_CIRCADIAN")]:
    cat[flag] = cat.gene_id.isin(S[sid])
cat["is_nearest_gene"] = cat.gene_id == cat.nearest_gene

ot = pd.read_csv(TABLES / "locus_opentargets.tsv", sep="\t")
cat = cat.merge(ot[["locus_no", "l2g_gene", "l2g_score", "l2g_gene_id", "l2g_genes_ge05"]], on="locus_no", how="left")
cat["is_l2g_top"] = cat.gene_id == cat.l2g_gene_id
l2g_gene = pd.read_csv(TABLES / "locus_l2g_genes.tsv", sep="\t").set_index(["locus_no", "gene_id"]).l2g_score
cat["l2g_gene_score"] = [l2g_gene.get((L, g), np.nan) for L, g in zip(cat.locus_no, cat.gene_id)]
coloc = pd.read_csv(TABLES / "locus_qtl_coloc.tsv", sep="\t")
qtl = coloc[coloc.right_type.astype(str).str.contains("qtl", case=False)]
qmax = qtl.groupby(["locus_no", "qtl_gene_id"]).clpp.max()
qtl = qtl[qtl.clpp.fillna(0) >= 0.1]
qg = qtl.groupby(["locus_no", "qtl_gene_id"]).biosample.apply(lambda s: ";".join(sorted(set(s.dropna())))[:300])
cat["qtl_max_clpp"] = [qmax.get((L, g)) for L, g in zip(cat.locus_no, cat.gene_id)]
cat["qtl_coloc_tissues"] = [qg.get((L, g)) for L, g in zip(cat.locus_no, cat.gene_id)]

mg = pd.read_csv(INTERIM / "magma" / "chronotype.genes.out", sep=r"\s+").set_index("GENE")
cat["magma_p"] = cat.gene_id.map(mg.P)
bonf = 0.05 / len(mg)
bur = pd.read_csv(TABLES / "burden_gene_level.tsv", sep="\t", index_col=0)
cat["burden_p"] = cat.gene_id.map(bur.p_burden)
cat["burden_p_ultrarare"] = cat.gene_id.map(bur["p_burden_maf_le_0.1pct"])
cat["burden_direction"] = cat.gene_id.map(bur.direction)
cat["burden_best_mask"] = cat.gene_id.map(bur.best_mask)
con = annot.gnomad_constraint().set_index("gene_id")   # descriptive annotation only
cat["LOEUF_v2"] = cat.gene_id.map(con.LOEUF)
rh = pd.read_csv(INTERIM / "rhythm" / "gene_rhythmicity.tsv", sep="\t", index_col=0)
cat["mouse_rhythmic_n"] = cat.gene_id.map(rh.mouse_n_rhythmic)
cat["mouse_rhythmic_tissues"] = cat.gene_id.map(rh.mouse_rhythmic_tissues)
cat["gtex_daynight_tissues"] = cat.gene_id.map(rh.gtex_daynight_tissues)

cat["ev_l2g"] = cat.l2g_gene_score >= 0.5
cat["ev_nearest"] = cat.is_nearest_gene
cat["ev_magma"] = cat.magma_p < bonf
cat["ev_burden"] = cat.burden_p < 0.01
cat["ev_qtl"] = cat.qtl_coloc_tissues.notna()
cat["evidence_count"] = cat[["ev_l2g", "ev_nearest", "ev_magma", "ev_burden", "ev_qtl"]].sum(axis=1)
cat = cat.sort_values(["evidence_count", "meta_p_from_z"], ascending=[False, True])
cols = ["locus_no", "lead_rsid", "chrom", "pos", "meta_p_from_z", "gene_context", "symbol", "gene_id", "diseases",
        "disease_ids", "inheritance", "clingen", "sleep_disturbance", "circadian_hpo", "monoamine_intd",
        "neurodevelopmental", "prolonged_QT", "hypertrophic_cm", "go_circadian", "organ_systems",
        "nearest_gene_symbol", "is_nearest_gene", "l2g_gene", "l2g_score", "is_l2g_top", "l2g_gene_score", "magma_p",
        "burden_p",
        "burden_p_ultrarare", "burden_direction", "burden_best_mask", "qtl_max_clpp", "qtl_coloc_tissues", "LOEUF_v2",
        "mouse_rhythmic_n", "mouse_rhythmic_tissues", "gtex_daynight_tissues",
        "ev_l2g", "ev_nearest", "ev_magma", "ev_burden", "ev_qtl", "evidence_count"]
cat[cols].to_csv(TABLES / "chronomic_gene_catalogue.tsv", sep="\t", index=False)
log.info("catalogue: %d locus-gene rows, %d genes, %d loci; evidence>=3: %s", len(cat), cat.gene_id.nunique(),
         cat.locus_no.nunique(), ", ".join(cat[cat.evidence_count >= 3].symbol))
