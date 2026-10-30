#!/usr/bin/env python3
"""Generate the manuscript's main tables (pandoc Markdown) from the results files."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, ROOT, TABLES  # noqa: E402

OUT = ROOT / "manuscript" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def fp(p: float, magma: bool = False) -> str:
    """p-value as text with pandoc superscript exponent (MAGMA values below its numerical floor as <1e-15)."""
    if p != p:
        return "–"
    if magma and p < 1e-15:
        return "<1×10^−15^"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10^−{int(e[1:])}^" if e.startswith("-") else f"{m}×10^{int(e)}^"
    if p >= 0.995:
        return "1.00"
    if p >= 0.1:
        return f"{p:.2f}"
    return f"{p:.{1 - int(np.floor(np.log10(p)))}f}"   # two significant digits, trailing zeros kept


def fexp(x: float) -> str:
    return f"{x:.2f}" if x < 1 else f"{x:.1f}"


def ffold(x: float) -> str:
    return f"{x:.0f}" if x >= 10 else f"{x:.2f}"


E = pd.read_csv(TABLES / "enrichment_all.tsv", sep="\t")
P = E[E.analysis == "primary"].set_index("set_id")
M = pd.read_csv(TABLES / "magma_genesets.tsv", sep="\t")
Mm = M[M.model == "marginal"].set_index("set_id").P
sizes = pd.read_csv(INTERIM / "genesets" / "geneset_summary.tsv", sep="\t").set_index("set_id").n_genes

ROWS = [
    ("Mendelian disease genes", [("MEND_ANY", "Any Mendelian (OMIM or Orphanet)"), ("MEND_CLINGEN", "ClinGen moderate–definitive"),
                                 ("INH_AD", "Autosomal dominant disorders"), ("INH_AR", "Autosomal recessive disorders")]),
    ("Positive control", [("GO_CIRCADIAN", "GO: circadian rhythm")]),
    ("Pre-specified phenotypes", [("HPO_HP0002360", "Sleep disturbance"),
                                 ("HPO_SLEEP_noNDD", "Sleep disturbance, excluding neurodevelopmental genes"),
                                 ("HPO_CIRCADIAN", "Sleep–wake cycle disturbance or abnormal melatonin"),
                                 ("HPO_HP0001657", "Prolonged QT interval"), ("HPO_HP0001639", "Hypertrophic cardiomyopathy"),
                                 ("HPO_HP0004431", "Reduced circulating complement (HPO)"),
                                 ("HPO_HP0001263", "Global developmental delay"),
                                 ("HPO_HP0001249", "Intellectual disability"), ("HPO_HP0001250", "Seizure")]),
    ("Post hoc", [("EXP_INTD_MONOAMINE", "Monoamine neurotransmitter disorders (12 genes)"),
                  ("EXP_CSF_HVA", "Decreased CSF homovanillic acid"),
                  ("CUR_COMPLEMENT_IUIS", "Complement deficiency genes (IUIS)")]),
    ("Rhythmic expression", [("RHY_MOUSE_GE3", "Rhythmic in ≥3 mouse tissues"), ("RHY_GTEX_DAYNIGHT", "Human day–night (GTEx)")]),
]
lines = ["| Gene set | Genes | Loci obs/exp (fold) | p | q | MAGMA p |",
         "|:" + "-" * 30 + "|" + "-" * 7 + ":|:" + "-" * 16 + "|:" + "-" * 10 + "|:" + "-" * 10 + "|:" + "-" * 9 + "|"]
for grp, sets in ROWS:
    lines.append(f"| ***{grp}*** | | | | | |")
    for sid, lab in sets:
        p = P.loc[sid]
        q = "–" if grp == "Post hoc" else fp(p.q_analytic_within_family)   # q is not meaningful after selection
        lines.append(
            f"| {lab} | {int(sizes[sid]):,} | {int(p.obs_loci)}/{fexp(p.exp_loci_analytic)} ({ffold(p.fold_analytic)}) "
            f"| {fp(p.p_enrich_analytic)} | {q} | {fp(Mm.get(sid, np.nan), True)} |")
(OUT / "table1.md").write_text("\n".join(lines) + "\n")

# ---------------------------------------------------------------- Table 3 (candidate genes)
cat = pd.read_csv(TABLES / "chronomic_gene_catalogue.tsv", sep="\t").set_index("symbol")
DIS = {"PER3": "Familial advanced sleep phase syndrome 3", "PER2": "Familial advanced sleep phase syndrome 1",
       "FBXL3": "Intellectual developmental disorder with short stature, facial anomalies and speech defects",
       "DDC": "AADC deficiency", "GCH1": "Dopa-responsive dystonia; GTPCH deficiency",
       "DNAJC12": "DNAJC12-deficient hyperphenylalaninaemia", "SLC18A2": "Brain dopamine–serotonin vesicular transport disease",
       "TCF4": "Pitt–Hopkins syndrome", "YWHAZ": "Popov–Chang syndrome", "GNAO1": "GNAO1 encephalopathy (DEE17/NEDIM)",
       "KCNH2": "Long QT syndrome 2; short QT syndrome 1", "CAV3": "Long QT syndrome 9; HCM; rippling muscle disease",
       "MYBPC3": "Hypertrophic cardiomyopathy 4"}
lines = ["| Gene | Lead variant | Mendelian disorder (inheritance) | L2G score (top gene if different) | MAGMA p | Exome burden p (direction) | QTL CLPP | Rhythmic (mouse) | Evid. |",
         "|:" + "-" * 11 + "|:" + "-" * 13 + "|:" + "-" * 19 + "|:" + "-" * 12 + "|:" + "-" * 10 + "|:" + "-" * 13
         + "|:" + "-" * 7 + "|:" + "-" * 11 + "|:" + "-" * 5 + ":|"]
for g in DIS:
    r = cat.loc[g]
    own = f"{r.l2g_gene_score:.2f}" if r.l2g_gene_score == r.l2g_gene_score else "–"
    l2g = own if (not isinstance(r.l2g_gene, str) or r.l2g_gene == g) else f"{own} ({r.l2g_gene} {r.l2g_score:.2f})"
    short = {"morningness": "morning", "eveningness": "evening"}
    bur = (f"{fp(r.burden_p)} ({short.get(r.burden_direction, r.burden_direction)})" if r.burden_p < 0.05
           else fp(r.burden_p)) \
        if r.burden_p == r.burden_p else "–"
    qtl = (f"{r.qtl_max_clpp:.2f}" if r.qtl_max_clpp >= 0.1 else f"{r.qtl_max_clpp:.3f}") \
        if r.qtl_max_clpp == r.qtl_max_clpp else "–"
    rh = f"{int(r.mouse_rhythmic_n)}/12" if r.mouse_rhythmic_n == r.mouse_rhythmic_n else "–"
    inh = str(r.inheritance).replace(";", "/") if isinstance(r.inheritance, str) else "–"
    lines.append(f"| *{g}* | {r.lead_rsid} | {DIS[g]} ({inh}) | {l2g} | {fp(r.magma_p, True)} | {bur} | {qtl} | {rh} | {int(r.evidence_count)} |")
(OUT / "table3.md").write_text("\n".join(lines) + "\n")
# ---------------------------------------------------------------- Table 2 (PanelApp screen, FDR < 0.05)
pr = pd.read_csv(TABLES / "enrichment_panelapp.tsv", sep="\t")
pmg = pd.read_csv(TABLES / "magma_panelapp.tsv", sep="\t").set_index("panel_id").P
NAMES = {"DDG2P": "DDG2P (developmental disorders)"}
lines = ["| PanelApp panel | Genes | Loci obs/exp (fold) | p | q | MAGMA p | Panel genes in loci |",
         "|:" + "-" * 26 + "|" + "-" * 8 + ":|:" + "-" * 14 + "|:" + "-" * 10 + "|:" + "-" * 8 + "|:" + "-" * 10 + "|:"
         + "-" * 27 + "|"]
for r in pr[pr.q_primary < 0.05].sort_values("p_primary").itertuples():
    g = r.genes_in_loci.split(";")
    gtxt = ", ".join(f"*{x}*" for x in g) if len(g) <= 10 else f"{len(g)} genes (Supplementary Table S3)"
    lines.append(f"| {NAMES.get(r.panel_name, r.panel_name)} | {int(r.n_genes_universe):,} | {int(r.obs_loci)}/"
                 f"{fexp(r.exp_primary)} ({ffold(r.fold_primary)}) | {fp(r.p_primary)} | {fp(r.q_primary)} | "
                 f"{fp(pmg.get(r.panel_id, np.nan), True)} | {gtxt} |")
(OUT / "table2.md").write_text("\n".join(lines) + "\n")
print((OUT / "table1.md").read_text())
print((OUT / "table2.md").read_text())
print((OUT / "table3.md").read_text())
