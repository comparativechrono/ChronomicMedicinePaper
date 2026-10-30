#!/usr/bin/env python3
"""Collect every number quoted in the manuscript into results/key_numbers.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, TABLES, ROOT, env_snapshot, save_json  # noqa: E402

K: dict = {}
K["loci"] = json.load(open(INTERIM / "loci" / "loci_summary.json"))
K["genesets"] = json.load(open(INTERIM / "genesets" / "provenance.json"))
K["genesets"]["hpo_release"] = K["genesets"]["hpo_version"].split("/")[-1]
K["genesets"]["sizes"] = pd.read_csv(INTERIM / "genesets" / "geneset_summary.tsv", sep="\t").set_index("set_id").n_genes.to_dict()
K["matching"] = json.load(open(TABLES / "enrichment_primary_matching.json"))
E = pd.read_csv(TABLES / "enrichment_all.tsv", sep="\t")
cols = ["obs_loci", "exp_loci_analytic", "fold", "fold_analytic", "p_enrich", "p_enrich_analytic",
        "q_analytic_within_family", "null_q025", "null_q975", "n_loci"]
K["enrichment"] = {a: {s: {c: float(r[c]) for c in cols} for s, r in sub.set_index("set_id").iterrows()}
                   for a, sub in E.groupby("analysis")}
R = pd.read_csv(TABLES / "enrichment_relative.tsv", sep="\t")
K["relative"] = {f"{r.analysis}|{r.numerator}/{r.denominator}": {"obs": r.obs_fraction, "null": r.null_mean_fraction,
                                                                 "p": r.p_relative, "num": int(r.obs_numerator),
                                                                 "den": int(r.obs_denominator)} for r in R.itertuples()}
M = pd.read_csv(TABLES / "magma_genesets.tsv", sep="\t")
K["magma_sets"] = {f"{r.model}|{r.set_id}": {"p": r.P, "beta": r.BETA, "ngenes": int(r.NGENES)} for r in M.itertuples()}
mg = pd.read_csv(INTERIM / "magma" / "chronotype.genes.out", sep=r"\s+")
K["magma_genes"] = {"n_genes": int(len(mg)), "bonferroni": 0.05 / len(mg), "n_sig": int((mg.P < 0.05 / len(mg)).sum())}
uni = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t").set_index("gene_id").symbol
mg["symbol"] = mg.GENE.map(uni)
K["magma_genes"]["selected"] = {s: float(mg.loc[mg.symbol == s, "P"].iloc[0]) for s in
                                ["PER2", "PER3", "FBXL3", "DDC", "GCH1", "SLC18A2", "DNAJC12", "KCNH2", "MTNR1B",
                                 "CAV3", "TCF4"] if (mg.symbol == s).any()}
B = pd.read_csv(TABLES / "burden_convergence.tsv", sep="\t")
K["burden_convergence"] = {f"{r.masks}|{r.gene_group}": {"n": int(r.n_genes), "mean": r.mean_neglog10p,
                                                         "null_mean": r.null_mean_neglog10p, "p_mean": r.p_mean,
                                                         "n01": int(r["n_p_lt_0.01"]), "exp01": r["exp_p_lt_0.01"],
                                                         "p_count": r.p_count, "genes": r["genes_p_lt_0.01"],
                                                         "trimmed_genes": r.trimmed_genes,
                                                         "mean_trim": r.mean_neglog10p_trimmed,
                                                         "null_mean_trim": r.null_mean_neglog10p_trimmed,
                                                         "p_mean_trim": r.p_mean_trimmed,
                                                         "l3o_null_mean": r.l3o_null_mean_neglog10p,
                                                         "l3o_p_mean": r.l3o_p_mean,
                                                         "l3o_n01": r["l3o_n_p_lt_0.01"],
                                                         "l3o_exp01": r["l3o_exp_p_lt_0.01"],
                                                         "l3o_p_count": r.l3o_p_count}
                           for _, r in B.iterrows()}
K["burden_summary"] = json.load(open(TABLES / "burden_summary.json"))
G = pd.read_csv(TABLES / "burden_gene_level.tsv", sep="\t", index_col=0)
K["burden_genes"] = {r.best_symbol: {"p": r.p_burden, "p_ultrarare": r["p_burden_maf_le_0.1pct"], "mask": r.best_mask,
                                     "trait": r.best_trait, "OR": r.best_odds_ratio, "direction": r.direction,
                                     "best_p": r.best_p_value, "pLOF_best_p": r.pLOF_best_p_value,
                                     "pLOF_best_mask": r.pLOF_best_mask}
                     for _, r in G[G.best_symbol.isin(["PER2", "PER3", "FBXL3", "MTNR1B", "GCH1", "DDC", "KCNH2",
                                                       "DNAJC12", "SLC18A2"])].iterrows()}
mono = G[G.best_symbol.isin(["GCH1", "DDC", "DNAJC12", "SLC18A2"])].p_burden
K["burden_monoamine_range"] = {"min": float(mono.min()), "max": float(mono.max())}
import yaml as _yaml  # noqa: E402
K["burden_sample_size"] = sorted({s["sample_size"] for acc in ("GCST90081559", "GCST90081560", "GCST90081561",
                                                               "GCST90081562")
                                  for s in _yaml.safe_load(open(ROOT / "data" / "raw" / "burden" /
                                                                f"{acc}_buildGRCh38.tsv.gz-meta.yaml"))["samples"]})[0]
K["opentargets"] = json.load(open(TABLES / "opentargets_summary.json"))
C = pd.read_csv(TABLES / "chronomic_gene_catalogue.tsv", sep="\t")
K["catalogue"] = {"rows": int(len(C)), "genes": int(C.gene_id.nunique()), "loci": int(C.locus_no.nunique()),
                  "evidence_ge3": int((C.evidence_count >= 3).sum()), "evidence_ge4": int((C.evidence_count >= 4).sum()),
                  "sleep_genes": sorted(C[C.sleep_disturbance].symbol.unique().tolist()),
                  "monoamine_genes": sorted(C[C.monoamine_intd].symbol.unique().tolist()),
                  "qt_genes": sorted(C[C.prolonged_QT].symbol.unique().tolist()),
                  "hcm_genes": sorted(C[C.hypertrophic_cm].symbol.unique().tolist()),
                  "per_gene": {r.symbol: {"evidence": int(r.evidence_count), "l2g": r.l2g_gene, "l2g_score": r.l2g_score,
                                          "nearest": bool(r.is_nearest_gene), "qtl": r.qtl_coloc_tissues,
                                          "mouse_rhythmic_n": r.mouse_rhythmic_n, "lead": r.lead_rsid,
                                          "LOEUF": r.LOEUF_v2, "qtl_max_clpp": r.qtl_max_clpp,
                                          "ev_l2g": bool(r.ev_l2g), "magma_p": r.magma_p, "burden_p": r.burden_p}
                               for r in C.itertuples()}}
K["rhythm"] = json.load(open(TABLES / "rhythmicity_summary.json"))
K["coloc_kcnh2"] = pd.read_csv(TABLES / "coloc_kcnh2.tsv", sep="\t").to_dict("records")
CC = pd.read_csv(TABLES / "coloc_kcnh2_conditional.tsv", sep="\t")
CS = pd.read_csv(TABLES / "coloc_kcnh2_conditional_signals.tsv", sep="\t")


LEAD_RS = "rs2072413"


def _pick(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).values if not isinstance(v, float) else np.isclose(df[k].values, v)
    return df[m].iloc[0]


joint = _pick(CC, window_kb=250, ridge=1e-3, conditioning="all_other_signals", p12=5e-6)
top1 = _pick(CC, window_kb=250, ridge=1e-3, conditioning="top_qt_signal_only", p12=5e-6)
T1 = CC[CC.conditioning == "top_qt_signal_only"]
sig = CS[(CS.window_kb == 250) & np.isclose(CS.ridge, 1e-3)]
dom = joint.qt_signals.split(";")[0]
domz = CS[(CS.window_kb == 250) & (CS.qt_signal == dom)].set_index("ridge").qt_z_conditional
tagrow = sig[sig.qt_signal == joint.tagging_signal].iloc[0]
weak = sig.loc[sig.qt_z_marginal.abs().idxmin()]
from scipy.stats import norm as _norm  # noqa: E402
DIR = pd.read_csv(TABLES / "coloc_kcnh2_direction.tsv", sep="\t").set_index("rsid")
K["coloc_kcnh2_conditional"] = {
    # primary: QT statistics conditioned on the dominant QT signal only (fewest LD-dependent assumptions)
    "PP_H4": float(top1["PP.H4"]), "PP_H3": float(top1["PP.H3"]), "dominant_signal": dom,
    "PP_H4_top_min": float(T1["PP.H4"].min()), "PP_H4_top_max": float(T1["PP.H4"].max()), "n_top": int(len(T1)),
    "qt_z_lead": float(top1.qt_z_lead_marginal), "qt_z_lead_cond_top": float(top1.qt_z_lead_conditional),
    # joint stepwise model (sensitivity)
    "n_signals": int(joint.n_qt_signals), "signals": joint.qt_signals,
    "n_signals_range": [int(CC.n_qt_signals.min()), int(CC.n_qt_signals.max())],
    "tagging_signal": joint.tagging_signal, "r2_tag_lead": float(tagrow.r2_with_chronotype_lead),
    "PP_H4_joint": float(joint["PP.H4"]), "PP_H3_joint": float(joint["PP.H3"]),
    "PP_H4_min": float(CC["PP.H4"].min()), "PP_H4_max": float(CC["PP.H4"].max()), "n_analyses": int(len(CC)),
    "all_top_shared_snps": ";".join(sorted(CC.top_shared_snp.unique())),
    "max_PP_H4_other_signals_all_settings": float(CS[CS.qt_signal != joint.tagging_signal]["PP.H4"].max()),
    "dominant_z_marginal": float(sig[sig.qt_signal == dom].qt_z_marginal.iloc[0]),
    "dominant_z_conditional_by_ridge": {str(k): float(v) for k, v in domz.items()},
    "weakest_signal": weak.qt_signal, "weakest_signal_marginal_p": float(2 * _norm.sf(abs(weak.qt_z_marginal))),
    # direction of effect for the chronotype lead allele
    "lead_allele": DIR.loc[LEAD_RS, "allele"], "lead_meta_lnOR_morning": float(DIR.loc[LEAD_RS, "meta_lnOR_morning_person_lead"]),
    "lead_z_chronotype": float(DIR.loc[LEAD_RS, "z_chronotype_ukb"]), "lead_beta_qt_ms": float(DIR.loc[LEAD_RS, "beta_qt_ms_marginal"]),
    "lead_z_qt_cond_top": float(DIR.loc[LEAD_RS, "z_qt_cond_top_signal"])}
qtw = pd.read_csv(INTERIM / "coloc_qt_chr7.tsv", sep="\t")
K["qt_gwas_n_max"] = int(qtw.n.max())
K["kcnh2_ld"] = pd.read_csv(TABLES / "coloc_kcnh2_ld.tsv", sep="\t").iloc[0].to_dict()
genes_tab = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str}).set_index("symbol")
leads_tab = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t", dtype={"chrom": str})
g = genes_tab.loc["MTNR1B"]
ll = leads_tab[leads_tab.chrom == g.chrom]
dist = np.minimum((ll.pos - g.end).abs(), (ll.pos - g.start).abs())
K["mtnr1b_nearest_lead"] = {"lead": ll.loc[dist.idxmin(), "lead_rsid"], "distance_kb": float(dist.min() / 1000)}
cfgp = __import__("yaml").safe_load(open(ROOT / "config" / "config.yaml"))
K["loci_pool"] = f"{cfgp['null_model']['pool_size']:,}"
ot = json.load(open(ROOT / "data" / "raw" / "opentargets" / "meta.json"))["meta"]["apiVersion"]
K["otapi"] = f"{ot['x']}.{ot['y']}.{ot['z']}"
# usage of HPO temporal ('chronomic') modifiers and timing-related phenotype terms in disease annotations
hp = pd.read_csv(ROOT / "data" / "raw" / "hpo" / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
mods = {"HP:0025302": "Diurnal", "HP:0410401": "Worse in evening", "HP:0031914": "Fluctuating",
        "HP:0025301": "Nocturnal"}
K["hpo_temporal_modifiers"] = {
    lab: {"annotations_as_modifier": int(hp.modifier.fillna("").str.contains(code).sum()),
          "diseases_as_modifier": int(hp.loc[hp.modifier.fillna("").str.contains(code), "database_id"].nunique()),
          "annotations_as_phenotype": int((hp.hpo_id == code).sum())}
    for code, lab in mods.items()}
K["hpo_annotation_rows"] = int(len(hp))
sys.path.insert(0, str(ROOT / "src"))
from chronomic import annot as _annot  # noqa: E402
import networkx as _nx  # noqa: E402
_G = _annot.hpo_graph()
_tp = _nx.ancestors(_G, "HP:0011008") | {"HP:0011008"}
_mod = hp.modifier.fillna("")
K["hpo_temporal_pattern_modifier_rows"] = int(_mod.apply(lambda s: any(x in _tp for x in s.split(";") if x)).sum())
K["hpo_temporal_pattern_modifier_counts"] = {
    _G.nodes[c]["name"]: int(n) for c, n in _mod.str.split(";").explode().loc[lambda s: s.isin(_tp)].value_counts().items()}
_fb = hp[(hp.database_id == "OMIM:606220") & (hp.hpo_id == "HP:0006979")]
K["fbxl3_sleep_annotation"] = {"frequency": _fb.frequency.iloc[0], "reference": _fb.reference.iloc[0],
                               "disease": _fb.disease_name.iloc[0]}
K["hpo_diseases"] = int(hp.database_id.nunique())
K["constraint_cuts"] = json.load(open(TABLES / "constraint_class_cuts.json"))
hits = pd.read_csv(TABLES / "locus_set_hits.tsv", sep="\t")
K["set_hit_genes"] = {s: ";".join(sorted({g for x in hits[hits.set_id == s].genes for g in x.split(";")}))
                      for s in ["EXP_INTD_MONOAMINE", "EXP_CSF_HVA", "CUR_COMPLEMENT_IUIS", "HPO_CIRCADIAN",
                                "HPO_HP0001657", "HPO_SLEEP_noNDD", "GO_CIRCADIAN"]}
K["magma_floor"] = 1e-15
_sl = C[C.sleep_disturbance].drop_duplicates("gene_id")
K["sleep_genes_in_loci"] = {"n": int(len(_sl)), "n_AD": int(_sl.inheritance.fillna("").str.contains("AD").sum()),
                            "n_loci": int(_sl.locus_no.nunique())}
_px = leads_tab[leads_tab.ld_index_type.astype(str).str.startswith("best_ukb")].copy()
_px["p"] = _px.ld_index_type.str.extract(r"P=([0-9.eE+-]+)").astype(float)
_px["dist_kb"] = (_px.index_bp - _px.pos).abs() / 1000
K["proxy_leads"] = {"n": int(len(_px)), "max_p": float(_px.p.max()), "max_dist_kb": float(_px.dist_kb.max())}
K["hpo_temporal_unused_total"] = int(sum(K["hpo_temporal_modifiers"][k]["annotations_as_modifier"]
                                         for k in ("Diurnal", "Fluctuating", "Worse in evening")))
_loo = pd.read_csv(TABLES / "magma_monoamine_leave_one_out.tsv", sep="\t")
K["magma_monoamine_loo"] = {r.set_id.replace("MONO_minus_", "minus_"): float(r.P) for r in _loo.itertuples()}
K["environment"] = env_snapshot()
K["n_manifest_files"] = int(len(pd.read_csv(ROOT / "data" / "manifest.tsv", sep="\t")))
save_json(K, ROOT / "results" / "key_numbers.json")
print("key numbers written")
