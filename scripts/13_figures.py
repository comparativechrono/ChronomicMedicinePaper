#!/usr/bin/env python3
"""Publication figures (PDF + 300-dpi PNG) in results/figures/.

Fig 1  Study design and the chronomic medicine framework
Fig 2  Locus-level enrichment against matched null loci, with gene-level (MAGMA) p-values: (A) Mendelian,
       pre-specified, post hoc and control gene sets; (B) HPO organ systems
Fig 3  Chronomic candidate genes: evidence grid
Fig 4  Monoamine synthesis and storage genes at chronotype loci (post hoc)
Fig 5  Rare/common/Mendelian convergence: (A) exome burden signal by gene group;
       (B) common-variant (MAGMA) vs rare-variant (burden) gene evidence
Supp S1  KCNH2 chronotype-QT colocalisation, unconditional and conditional
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import FIGURES, INTERIM, TABLES  # noqa: E402

# ---- design tokens (reference palette, light mode; print figures) ----
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID, SURF = "#e6e5e1", "#ffffff"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3 (validated all-pairs)
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none",
    "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "figure.facecolor": SURF, "axes.facecolor": SURF})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor=SURF)
    plt.close(fig)


def fmt_p(p, perm=False):
    if perm and p <= 1 / 10001 + 1e-12:
        return "≤1.0×10$^{-4}$"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10$^{{{int(e)}}}$"
    if p >= 0.1:
        return f"{p:.2f}"
    return f"{p:.{1 - int(np.floor(np.log10(p)))}f}"   # two significant digits, trailing zeros kept


E = pd.read_csv(TABLES / "enrichment_all.tsv", sep="\t")
K_per_gene = json.load(open(TABLES.parent / "key_numbers.json"))["catalogue"]["per_gene"]


# ======================================================================== Fig 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 63); ax.axis("off")

    def box(x, y, w, h, title, lines, fc="#f4f8fd", ec=S1):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.0", fc=fc, ec=ec, lw=0.8))
        ax.text(x + 1.3, y + h - 1.5, title, ha="left", va="top", fontsize=7.4, weight="bold", color=INK)
        ax.text(x + 1.3, y + h - 5.3, "\n".join(lines), ha="left", va="top", fontsize=6.1, color=INK2, linespacing=1.4)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=0.8, color=MUTED))

    ax.text(0.5, 61.5, "Common variation", fontsize=7.6, weight="bold", color=S1)
    box(0.5, 41, 31, 18.5, "351 chronotype loci",
        ["Jones et al. 2019 (n = 697,828)", "LD intervals: 1000G EUR, r² ≥ 0.6", "UKB summary statistics → MAGMA",
         "Open Targets fine-mapping, L2G"])
    ax.text(0.5, 38.3, "Mendelian disease knowledge", fontsize=7.6, weight="bold", color=S2)
    box(0.5, 20, 31, 16.5, "Disease genes & phenotypes",
        ["OMIM (via HPO), Orphanet, ClinGen", "HPO annotations (2026-09 release)",
         "PanelApp clinical gene panels (screen)", "pre-specified: sleep/circadian, cardiac,",
         "immune, neurodevelopmental, growth"], fc="#fdf3ee", ec=S2)
    ax.text(0.5, 17.3, "Rare variation & chronobiology", fontsize=7.6, weight="bold", color=S3)
    box(0.5, 0.5, 31, 15, "Exomes and rhythms",
        ["UKB exome burden (n = 388,330)", "mouse 12-tissue circadian atlas", "GTEx day–night expression"],
        fc="#eefaf5", ec=S3)

    box(35.8, 30.5, 30, 29, "Locus-level enrichment",
        ["null loci matched on MAF, LD,", "gene count and distance", "10,000 permutations; exact",
         "Poisson-binomial p-values", "relative tests vs Mendelian, NDD", "and dominant disease genes",
         "sensitivity: r² thresholds,", "MHC/17q21 removed, nearest gene"],
        fc="#ffffff", ec=MUTED)
    box(35.8, 14.5, 30, 13, "Gene-level tests",
        ["MAGMA competitive gene sets,", "marginal and conditioned on", "NDD or Mendelian genes"], fc="#ffffff", ec=MUTED)
    box(35.8, 0.5, 30, 11, "Convergence",
        ["exome burden vs matched genes", "conditional coloc at KCNH2 (QT)"], fc="#ffffff", ec=MUTED)

    box(70.3, 26, 29.4, 33.5, "Chronomic medicine",
        ["integrating chronobiology with", "genomic medicine: genes where", "common variation is associated",
         "with circadian timing and rare", "variation causes disease with sleep,", "circadian or time-of-day features", "",
         "Question: do chronotype loci", "converge on Mendelian disease", "genes, and for which phenotypes?"],
        fc="#f6f6f4", ec=INK2)
    box(70.3, 0.5, 29.4, 22, "Outputs",
        ["locus- and gene-level enrichment", "candidate-gene catalogue", "(Mendelian genes at loci,", "with gene-level evidence)",
         "open, version-pinned pipeline"], fc="#f6f6f4", ec=INK2)
    arrow(31.8, 50, 35.7, 50)
    arrow(31.8, 27, 35.7, 40)
    arrow(31.8, 27, 35.7, 21)
    arrow(31.8, 8, 35.7, 6)
    arrow(66.2, 45, 70.0, 42)
    arrow(66.2, 21, 70.0, 16)
    arrow(66.2, 6, 70.0, 10)
    save(fig, "Fig1_study_design")


# ======================================================================== Fig 2
ROWS = [
    ("Mendelian disease genes", [("MEND_ANY", "Any Mendelian (OMIM or Orphanet)"), ("MEND_OMIM", "OMIM"),
                                 ("MEND_ORPHA", "Orphanet (assessed)"), ("MEND_CLINGEN", "ClinGen moderate–definitive"),
                                 ("INH_AD", "Autosomal dominant"), ("INH_AR", "Autosomal recessive")]),
    ("Positive control", [("GO_CIRCADIAN", "GO: circadian rhythm")]),
    ("Pre-specified phenotypes", [("HPO_HP0002360", "Sleep disturbance"),
                                 ("HPO_SLEEP_noNDD", "  without neurodevelopmental features"),
                                 ("HPO_CIRCADIAN", "Sleep–wake cycle / melatonin"),
                                 ("HPO_HP0001657", "Prolonged QT interval"), ("HPO_HP0011675", "Arrhythmia"),
                                 ("HPO_HP0001639", "Hypertrophic cardiomyopathy"),
                                 ("HPO_HP0001638", "Cardiomyopathy"),
                                 ("HPO_HP0004431", "Reduced circulating complement"),
                                 ("HPO_HP0002721", "Immunodeficiency"),
                                       ("HPO_HP0001263", "Global developmental delay"),
                                       ("HPO_HP0001249", "Intellectual disability"), ("HPO_HP0001250", "Seizure"),
                                       ("HPO_HP0004322", "Short stature")]),
    ("Post hoc", [("EXP_INTD_MONOAMINE", "Monoamine neurotransmitter disorders"),
                  ("EXP_CSF_HVA", "Decreased CSF homovanillic acid"),
                  ("CUR_COMPLEMENT_IUIS", "Complement deficiency genes (IUIS)")]),
    ("Rhythmic expression", [("RHY_MOUSE_GE3", "Rhythmic in ≥3 mouse tissues"),
                             ("RHY_GTEX_DAYNIGHT", "Human day–night (GTEx)")]),
]


def fig2():
    P = E[E.analysis == "primary"].set_index("set_id")
    MG = pd.read_csv(TABLES / "magma_genesets.tsv", sep="\t")
    Mm = MG[MG.model == "marginal"].set_index("set_id").P

    def fmt_mp(sid):
        p = Mm.get(sid, np.nan)
        return "–" if p != p else fmt_p(p)

    COLS = ((1.04, "obs / exp\nloci"), (1.32, "p, matched\nnull loci"), (1.62, "MAGMA p\n(gene level)"))
    POSTHOC = {sid for grp, sets in ROWS if grp == "Post hoc" for sid, _ in sets}
    items = []
    for grp, sets in ROWS:
        items.append(("hdr", grp, None))
        items += [("row", lab, sid) for sid, lab in sets]
    n = len(items)
    fig = plt.figure(figsize=(7.2, 9.8))
    ax = fig.add_axes([0.33, 0.405, 0.30, 0.57])
    lo, hi = 0.25, 250
    ax.set_xscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(-0.7, n - 0.3)
    tr = ax.get_yaxis_transform()  # x in axes fraction, y in data
    y = np.arange(n)[::-1]
    for yi, (kind, lab, sid) in zip(y, items):
        if kind == "hdr":
            ax.text(-0.98, yi, lab, transform=tr, ha="left", va="center", fontsize=7.2, weight="bold", color=INK)
            continue
        ax.text(-0.94, yi, lab, transform=tr, ha="left", va="center", fontsize=6.8, color=INK)
        r = P.loc[sid]
        f = r.fold_analytic if r.obs_loci > 0 else lo * 1.08
        if sid in POSTHOC:   # q is not meaningful after selection: no significance fill
            ax.plot([f], [yi], marker="D", ms=3.8, mfc=INK2, mec=INK2, ls="none", zorder=3, clip_on=False)
        else:
            sig = r.q_analytic_within_family < 0.05
            ax.plot([f], [yi], marker="o", ms=4.4, mfc=S1 if sig else SURF, mec=S1, mew=1.1, ls="none", zorder=3,
                    clip_on=False)
        if r.obs_loci == 0:
            ax.text(lo * 1.35, yi, "0 observed", fontsize=6, color=INK2, va="center")
        ev = r.exp_loci_analytic
        ax.text(1.04, yi, f"{int(r.obs_loci)} / {ev:.2f}" if ev < 1 else f"{int(r.obs_loci)} / {ev:.1f}", transform=tr,
                fontsize=6.5, color=INK2, va="center")
        ax.text(1.32, yi, fmt_p(r.p_enrich_analytic), transform=tr, fontsize=6.5, color=INK2, va="center")
        ax.text(1.62, yi, fmt_mp(sid), transform=tr, fontsize=6.5, color=INK2, va="center")
    for x, s in COLS:
        ax.text(x, n - 0.2, s, transform=tr, fontsize=6.5, color=INK, weight="bold", va="bottom")
    ax.text(-0.98, n + 1.6, "A", transform=tr, fontsize=10, weight="bold", va="bottom")
    ax.axvline(1, color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks([0.5, 1, 2, 5, 10, 50, 200]); ax.set_xticklabels(["0.5", "1", "2", "5", "10", "50", "200"])
    ax.set_xlabel("Fold enrichment (observed / expected loci)")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    h1 = ax.plot([], [], "o", mfc=S1, mec=S1, ms=4.4, label="FDR q < 0.05 within set family")[0]
    h3 = ax.plot([], [], "o", mfc=SURF, mec=S1, ms=4.4, label="FDR q ≥ 0.05")[0]
    h4 = ax.plot([], [], "D", mfc=INK2, mec=INK2, ms=3.8, label="post hoc set (q not shown)")[0]
    ax.legend(handles=[h1, h3, h4], loc="upper left", bbox_to_anchor=(-1.0, -0.055), fontsize=6.5, ncol=3,
              handletextpad=0.3, title="Null loci matched on MAF, LD, gene count and distance",
              title_fontsize=6.5, alignment="left")

    axb = fig.add_axes([0.33, 0.04, 0.30, 0.25])
    Ps = E[(E.analysis == "primary") & (E.family == "hpo_system")].set_index("set_id")
    Ps = Ps[Ps.exp_loci_analytic >= 5].sort_values("fold_analytic")
    yy = np.arange(len(Ps))
    axb.set_xlim(0.75, 1.9); axb.set_ylim(-0.7, len(Ps) - 0.3)
    trb = axb.get_yaxis_transform()
    for yi, sid in zip(yy, Ps.index):
        r = Ps.loc[sid]
        axb.plot([r.fold_analytic], [yi], "o", ms=3.8, mfc=S1 if r.q_analytic_within_family < 0.05 else SURF,
                 mec=S1, mew=1, zorder=3)
        lab = r.label.replace("Abnormality of the ", "").replace("Abnormality of ", "")
        axb.text(-0.94, yi, lab[0].upper() + lab[1:], transform=trb, fontsize=6.3, color=INK, va="center")
        axb.text(1.04, yi, f"{int(r.obs_loci)} / {r.exp_loci_analytic:.1f}", transform=trb, fontsize=6.3, color=INK2,
                 va="center")
        axb.text(1.32, yi, fmt_p(r.p_enrich_analytic), transform=trb, fontsize=6.3, color=INK2, va="center")
        axb.text(1.62, yi, fmt_mp(sid), transform=trb, fontsize=6.3, color=INK2, va="center")
    for x, s in COLS:
        axb.text(x, len(Ps) - 0.2, s, transform=trb, fontsize=6.5, color=INK, weight="bold", va="bottom")
    axb.axvline(1, color=MUTED, lw=0.7, ls=(0, (3, 2)))
    axb.set_yticks([]); axb.spines["left"].set_visible(False)
    axb.set_xlabel("Fold enrichment — HPO organ systems (expected ≥ 5 loci)")
    axb.grid(axis="x", color=GRID, lw=0.5); axb.set_axisbelow(True)
    axb.text(-0.98, len(Ps) + 0.3, "B", transform=trb, fontsize=10, weight="bold", va="bottom")
    save(fig, "Fig2_locus_enrichment")


# ======================================================================== Fig 5
def fig5():
    B = pd.read_csv(TABLES / "burden_convergence.tsv", sep="\t")
    G = pd.read_csv(TABLES / "burden_gene_level.tsv", sep="\t", index_col=0)
    mg = pd.read_csv(INTERIM / "magma" / "chronotype.genes.out", sep=r"\s+").set_index("GENE")
    fig = plt.figure(figsize=(7.2, 3.3))
    ax = fig.add_axes([0.25, 0.16, 0.25, 0.74])
    sub = B[B.masks == "all masks"].iloc[::-1].reset_index(drop=True)
    for i, r in sub.iterrows():
        ax.plot([r.null_q025_neglog10p, r.null_q975_neglog10p], [i, i], color=GRID, lw=3.2, solid_capstyle="round")
        ax.plot([r.null_mean_neglog10p], [i], "|", color=MUTED, ms=7, mew=1.2)
        ax.plot([r.mean_neglog10p], [i], "o", color=INK2 if "outside" in r.gene_group else S1, ms=5)
        ptxt = fmt_p(r.p_mean, perm=True)
        ptxt = ptxt if ptxt.startswith("≤") else f"= {ptxt}"
        trim = f" ({fmt_p(r.p_mean_trimmed, perm=True)})" if r.p_mean_trimmed == r.p_mean_trimmed else ""
        ax.text(1.06, i, f"p {ptxt}{trim}", va="center", fontsize=6.4, color=INK2, transform=ax.get_yaxis_transform())
        lab = r.gene_group.replace(" in chronotype loci", " in loci").replace("Sleep-disturbance", "Sleep-disturbance") + f" (n={r.n_genes})"
        ax.text(-0.04, i, lab, va="center", ha="right", fontsize=6.6, color=INK, transform=ax.get_yaxis_transform())
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.6, len(sub) - 0.4)
    ax.set_xlabel("Mean −log$_{10}$ p, exome burden (chronotype)")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.text(-1.05, len(sub) - 0.2, "A", fontsize=10, weight="bold", transform=ax.get_yaxis_transform())
    ax.text(-0.85, -0.22, "Grey band: 95% range of CAF-matched random gene sets; tick: null mean. In brackets: p with the\n"
            "three most significant genes removed from the group and from every matched set (PER2, PER3, FBXL3 for the\n"
            "Mendelian and sleep groups).", fontsize=5.8, color=INK2, transform=ax.transAxes, va="top")

    axb = fig.add_axes([0.72, 0.16, 0.26, 0.74])
    d = G.join(mg.P.rename("magma_p"), how="inner")
    d["x"] = -np.log10(d.magma_p).clip(upper=16)
    d["y"] = -np.log10(d.p_burden).clip(upper=16)
    rest = d[~d.in_locus]
    axb.scatter(rest.x, rest.y, s=3, c="#cfcec9", lw=0, rasterized=True, label="Other genes")
    loc = d[d.in_locus & ~d.sleep_mendelian]
    axb.scatter(loc.x, loc.y, s=6, c=S1, lw=0, alpha=0.85, label="Genes in chronotype loci")
    sl = d[d.in_locus & d.sleep_mendelian]
    axb.scatter(sl.x, sl.y, s=16, c=S2, edgecolors=SURF, lw=0.6, label="Mendelian sleep-disturbance genes in loci", zorder=3)
    for sym, dx, dy in (("PER3", -0.4, 0.4), ("PER2", 0.3, 0.3), ("FBXL3", 0.3, 0.2), ("MTNR1B", 0.3, 0.1)):
        r = d[d.best_symbol == sym]
        if len(r):
            axb.annotate(sym, (r.x.iloc[0], r.y.iloc[0]), xytext=(r.x.iloc[0] + dx, r.y.iloc[0] + dy), fontsize=6.6,
                         color=INK, style="italic")
    axb.axhline(-np.log10(0.05 / len(G)), color=MUTED, lw=0.6, ls=(0, (3, 2)))
    axb.axvline(-np.log10(0.05 / len(mg)), color=MUTED, lw=0.6, ls=(0, (3, 2)))
    axb.set_xlabel("Common variants: MAGMA −log$_{10}$ p")
    axb.set_ylabel("Rare variants: burden −log$_{10}$ p")
    axb.set_xlim(0, 16.8); axb.set_ylim(0, 15)
    axb.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), fontsize=6, handletextpad=0.2, borderaxespad=0.0,
               markerscale=1.2)
    axb.text(-0.25, 1.2, "B", fontsize=10, weight="bold", transform=axb.transAxes)
    save(fig, "Fig5_rare_common_convergence")


# ======================================================================== Fig 3
def fig3():
    C = pd.read_csv(TABLES / "chronomic_gene_catalogue.tsv", sep="\t")
    keep = C[C.sleep_disturbance | C.monoamine_intd | C.prolonged_QT | C.circadian_hpo | C.hypertrophic_cm].copy()
    keep = keep.drop_duplicates("symbol")
    keep["grp"] = np.select([keep.monoamine_intd, keep.circadian_hpo, keep.sleep_disturbance,
                             keep.prolonged_QT | keep.hypertrophic_cm], [0, 1, 2, 3], 4)
    keep = keep.sort_values(["grp", "evidence_count", "meta_p_from_z"], ascending=[True, False, True])
    cols = [("ev_l2g", "L2G ≥ 0.5"), ("ev_nearest", "Nearest gene"), ("ev_magma", "MAGMA (Bonferroni)"),
            ("ev_burden", "Exome burden p<0.01"), ("ev_qtl", "QTL coloc. (CLPP ≥ 0.1)"),
            ("sleep_disturbance", "Sleep disturbance"), ("circadian_hpo", "Sleep–wake / melatonin"),
            ("monoamine_intd", "Monoamine disorder"), ("neurodevelopmental", "Neurodevelopmental"),
            ("prolonged_QT", "Prolonged QT"), ("hypertrophic_cm", "Hypertrophic CM"),
            ("go_circadian", "GO circadian"), ("rhy_mouse", "Rhythmic ≥3 mouse tissues"),
            ("rhy_gtex", "Human day–night"), ("constrained", "LOEUF < 0.35")]
    keep = keep.assign(rhy_mouse=keep.mouse_rhythmic_n >= 3, rhy_gtex=keep.gtex_daynight_tissues.notna(),
                       constrained=keep.LOEUF_v2 < 0.35)
    M = keep[[c for c, _ in cols]].fillna(False).astype(bool).values
    n, m = M.shape
    fig = plt.figure(figsize=(7.2, 0.155 * n + 1.5))
    ax = fig.add_axes([0.24, 0.02, 0.62, 0.84 * (0.155 * n) / (0.155 * n + 1.5) + 0.06])
    groups = [(0, 5, S1), (5, 11, S2), (11, 15, S3)]
    for gi, (a, b, col) in enumerate(groups):
        for j in range(a, b):
            for i in range(n):
                if M[i, j]:
                    ax.add_patch(plt.Rectangle((j + 0.08, n - 1 - i + 0.08), 0.84, 0.84, fc=col, ec="none"))
                else:
                    ax.add_patch(plt.Rectangle((j + 0.08, n - 1 - i + 0.08), 0.84, 0.84, fc="#f2f1ee", ec="none"))
    ax.set_xlim(0, m); ax.set_ylim(0, n); ax.axis("off")
    for j, (_, lab) in enumerate(cols):
        ax.text(j + 0.5, n + 0.25, lab, rotation=55, ha="left", va="bottom", fontsize=6.4, color=INK)
    for i, r in enumerate(keep.itertuples()):
        ax.text(-0.15, n - 1 - i + 0.5, f"{r.symbol}", ha="right", va="center", fontsize=6.6, color=INK, style="italic")
        ax.text(-2.4, n - 1 - i + 0.5, f"{r.lead_rsid}", ha="right", va="center", fontsize=6, color=INK2)
        ax.text(m + 0.2, n - 1 - i + 0.5, f"{int(r.evidence_count)}", ha="left", va="center", fontsize=6.4, color=INK2)
    ax.text(m + 0.2, n + 0.25, "Gene\nevidence", ha="left", va="bottom", fontsize=6.2, color=INK, weight="bold")
    for a, b, col, lab in ((0, 5, S1, "Gene-level evidence"), (5, 11, S2, "Mendelian phenotype"),
                           (11, 15, S3, "Chronobiology / constraint")):
        ax.plot([a + 0.1, b - 0.1], [-0.35, -0.35], color=col, lw=2.2, solid_capstyle="round", clip_on=False)
        ax.text((a + b) / 2, -0.75, lab, ha="center", va="top", fontsize=6.4, color=INK)
    save(fig, "Fig3_candidate_genes")


# ======================================================================== Fig 4
def fig4():
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.set_xlim(0, 104); ax.set_ylim(-11.5, 40); ax.axis("off")
    HIT = {"GCH1", "DDC", "DNAJC12", "SLC18A2"}

    def node(x, y, txt, w=11, h=5.2, kind="metab"):
        fc = {"metab": "#ffffff", "context": "#f2f1ee"}.get(kind, "#fdf3ee" if txt in HIT else "#f2f1ee")
        ec = MUTED if kind in ("metab", "context") else (S2 if txt in HIT else MUTED)
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.25,rounding_size=1.0",
                                    fc=fc, ec=ec, lw=1.2 if txt in HIT else 0.7))
        ax.text(x, y, txt, ha="center", va="center", fontsize=6.8,
                weight="bold" if txt in HIT else "normal", style="italic" if kind == "gene" else "normal",
                color=INK2 if kind == "context" else INK)

    def arr(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=INK2))

    # serotonin -> melatonin row
    y1, y2 = 30, 12
    for x, t in ((8, "Tryptophan"), (32, "5-HTP"), (56, "Serotonin"), (80, "N-acetyl-\nserotonin")):
        node(x, y1, t, w=12)
    node(96, y1 - 10, "Melatonin", w=12, kind="context")
    for x, t in ((8, "Tyrosine"), (32, "L-DOPA"), (56, "Dopamine"), (80, "Vesicular\nmonoamines")):
        node(x, y2, t, w=12)
    for xa, xb, y in ((14, 26, y1), (38, 50, y1), (62, 74, y1), (14, 26, y2), (38, 50, y2), (62, 74, y2)):
        arr(xa, y, xb, y)
    arr(86, y1 - 1, 93, y1 - 7.2)
    node(20, y1 + 6.2, "TPH1/2", w=9, h=4.2, kind="gene")
    node(20, y2 - 6.2, "TH", w=9, h=4.2, kind="gene")
    node(44, y1 + 6.2, "DDC", w=9, h=4.2, kind="gene")
    node(44, y2 - 6.2, "DDC", w=9, h=4.2, kind="gene")
    node(68, y1 + 6.2, "AANAT", w=9, h=4.2, kind="gene")
    node(68, y2 - 6.2, "SLC18A2", w=10, h=4.2, kind="gene")
    node(96, y1 + 1, "ASMT", w=8, h=4.2, kind="gene")
    node(20, 21, "GCH1", w=9, h=4.2, kind="gene")
    node(31, 21, "DNAJC12", w=10.5, h=4.2, kind="gene")
    ax.text(42, 21, "BH4 cofactor & chaperone\nfor TPH and TH", fontsize=6.2, color=INK2, va="center")
    ax.add_patch(FancyArrowPatch((20, 23.2), (20, y1 - 3), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=S2))
    ax.add_patch(FancyArrowPatch((20, 18.8), (20, y2 + 3), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=S2))
    ax.add_patch(plt.Rectangle((2, -3.9), 2.4, 2.4, fc="#fdf3ee", ec=S2, lw=1.2))
    ax.text(5.5, -2.7, "gene inside a chronotype GWAS locus (r² ≥ 0.6 interval); all four are inherited monoamine "
            "neurotransmitter disorder genes\n(GTPCH deficiency / dopa-responsive dystonia, AADC deficiency, "
            "DNAJC12 deficiency, VMAT2 deficiency). Grey: context only; no melatonin-pathway gene lies in a locus.",
            fontsize=6.0, color=INK2, va="center")
    ev = {g: K_per_gene.get(g, {}) for g in ("DDC", "GCH1", "DNAJC12", "SLC18A2")}

    def l2g(g):
        s = ev[g].get("l2g_score")
        return f"top L2G gene ({s:.2f})" if ev[g].get("l2g") == g and s == s else "not the L2G gene"
    lines = [f"DDC: {l2g('DDC')}; MAGMA p < 10$^{{-15}}$; plasma pQTL (CLPP {ev['DDC'].get('qtl_max_clpp', float('nan')):.2f})",
             f"GCH1: nearest gene; {l2g('GCH1')}",
             f"DNAJC12: {l2g('DNAJC12')}; nearest gene SIRT1 (itself a circadian regulator)",
             "SLC18A2: no gene-level support beyond locus overlap",
             "Exome burden: none of the four (p ≥ 0.12)"]
    ax.text(1.0, -9.0, "Gene-level support:  " + "   ·   ".join(lines[:2]) + "\n" + "   ·   ".join(lines[2:]),
            fontsize=6.0, color=INK2, va="center")
    save(fig, "Fig4_monoamine_pathway")


# ======================================================================== Supplementary figure
def figS1():
    p = TABLES / "coloc_kcnh2.tsv"
    if not p.exists():
        return
    from scipy.stats import norm
    R = pd.read_csv(p, sep="\t")
    Z = pd.read_csv(TABLES / "coloc_kcnh2_conditional_snps.tsv", sep="\t")
    KC = json.load(open(TABLES.parent / "key_numbers.json"))["coloc_kcnh2_conditional"]
    lead = 150647969
    fig, axs = plt.subplots(3, 1, figsize=(7.2, 5.0), sharex=True, gridspec_kw={"hspace": 0.18})
    panels = ((axs[0], "z_chronotype", S1, "Chronotype (UK Biobank, n = 449,734)"),
              (axs[1], "z_qt", S2, "QT interval (GERA, UK Biobank, PAGE; n ≤ 117,532)"),
              (axs[2], "z_qt_cond_top_signal", S2, f"QT interval conditioned on the dominant QT signal "
                                                 f"{KC['dominant_signal']} (approximate; 1000G EUR LD)"))
    sig = Z[Z.is_qt_signal]
    for ax, col, colr, lab in panels:
        lp = -norm.logsf(Z[col].abs()) / np.log(10) - np.log10(2)
        ax.scatter(Z.pos / 1e6, lp, s=3.5, c=colr, lw=0, rasterized=True)
        ax.axvline(lead / 1e6, color=MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.set_ylabel("−log$_{10}$ p")
        ax.text(0.01, 0.95, lab, transform=ax.transAxes, fontsize=6.8, color=INK, va="top")
        ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
        if col == "z_qt":
            s_lp = -norm.logsf(sig[col].abs()) / np.log(10) - np.log10(2)
            ax.scatter(sig.pos / 1e6, s_lp, s=16, facecolors="none", edgecolors=INK, lw=0.7, zorder=3)
    tag = Z[Z.rsid == KC["tagging_signal"]].iloc[0]
    ytag = -norm.logsf(abs(tag.z_qt_cond_top_signal)) / np.log(10) - np.log10(2)
    axs[2].set_ylim(-4, ytag * 1.45)
    axs[2].annotate(f"{KC['tagging_signal']} / rs2072413 (adjacent; r² = {KC['r2_tag_lead']:.2f})", (tag.pos / 1e6, ytag),
                    xytext=(40, -6), textcoords="offset points", fontsize=6.2, color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    axs[2].axvspan(150.642044, 150.675403, color="#f2f1ee", zorder=0)
    axs[2].text(150.6587, ytag * 0.45, "KCNH2", ha="center", fontsize=6.6, style="italic")
    axs[2].set_xlabel("Chromosome 7 position (Mb, GRCh37)")
    r = R[R.window_kb == 250].iloc[0]
    axs[0].set_title(f"Unconditional coloc (±250 kb): PP.H4 = {r['PP.H4']:.2f}, PP.H3 = {r['PP.H3']:.2f}.   "
                     f"Conditional on {KC['dominant_signal']}: PP.H4 = {KC['PP_H4']:.2f} "
                     f"({KC['PP_H4_min']:.2f}–{KC['PP_H4_max']:.2f} over {KC['n_analyses']} sensitivity analyses).\n"
                     "Open circles (middle): QT signals selected by stepwise conditional analysis. Dashed line: "
                     "chronotype lead rs2072413.", loc="left", fontsize=6.2, color=INK2)
    save(fig, "FigS1_kcnh2_coloc")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, figS1):
        f()
        print("done", f.__name__)
