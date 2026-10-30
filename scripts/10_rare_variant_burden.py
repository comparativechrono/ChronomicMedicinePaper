#!/usr/bin/env python3
"""Rare-variant (exome) burden evidence for chronotype at genes in chronotype GWAS loci.

Data: Backman et al. 2021 (Nature 599:628) UK Biobank exome gene-burden results for the four
categories of the chronotype question (UKB field 1180), each analysed as a binary trait
(GWAS Catalog GCST90081559-62). Masks: M1 = pLOF; M3 = pLOF + deleterious missense; five
allele-frequency ceilings each.

Gene-level evidence: Cauchy (ACAT) combination of the 10 masks within each trait, then of the
four traits (valid under arbitrary dependence). Direction: sign of log(OR) for 'definitely a
morning person' minus that for 'definitely an evening person' in the most significant mask.

Tests: genes in chronotype loci (r2 >= 0.6 intervals) and subsets vs. random gene sets matched
on cumulative allele frequency (CAF) of qualifying variants (10,000 draws). Statistics: mean -log10 p,
the number of genes with p < 0.01. Sensitivity: (i) trimmed mean, removing the three most significant
genes from the observed set and from every draw; (ii) leave-three-out, comparing the remaining genes with
untrimmed matched draws of the same size (conservative). Permutation p-values have a floor of 1/(1 + 10,000).
Outputs: results/tables/burden_gene_level.tsv, burden_convergence.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import cauchy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, RAW, TABLES, get_logger, load_config, save_json  # noqa: E402

log = get_logger("10_rare_variant_burden")
cfg = load_config()
rng = np.random.default_rng(cfg["seed"])
rng_l3o = np.random.default_rng(cfg["seed"] + 1)
N_DRAW = 10000
TRIM = 3   # leave-top-3-out sensitivity analysis

TRAITS = {"GCST90081559": "definitely_morning", "GCST90081560": "more_morning",
          "GCST90081561": "more_evening", "GCST90081562": "definitely_evening"}


def acat(p: np.ndarray) -> float:
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-16)
    p = p[~np.isnan(p)]
    if p.size == 0:
        return np.nan
    t = np.mean(np.tan((0.5 - p) * np.pi))
    return float(cauchy.sf(t))


frames = []
for acc, trait in TRAITS.items():
    d = pd.read_csv(RAW / "burden" / f"{acc}_buildGRCh38.tsv.gz", sep="\t",
                    usecols=["Name", "odds_ratio", "p_value", "effect_allele_frequency"])
    d["gene_id"] = d.Name.str.extract(r"\((ENSG\d+)\)", expand=False)
    d["symbol"] = d.Name.str.extract(r"^([^(]+)\(", expand=False)
    d["mask"] = d.Name.str.split(".GENE.").str[1]
    d["trait"] = trait
    frames.append(d)
B = pd.concat(frames, ignore_index=True)
B["logOR"] = np.log(B.odds_ratio)

per_trait = B.groupby(["gene_id", "trait"]).p_value.apply(lambda s: acat(s.values)).unstack()
omni = per_trait.apply(lambda r: acat(r.values), axis=1).rename("p_burden")
# sensitivity: only masks with MAF <= 0.1% (drop the '.1' = MAF <= 1% masks) to avoid tagging common variants
Br = B[~B["mask"].str.endswith(".1")]
omni_rare = (Br.groupby(["gene_id", "trait"]).p_value.apply(lambda s: acat(s.values)).unstack()
             .apply(lambda r: acat(r.values), axis=1).rename("p_burden_maf_le_0.1pct"))
best = B.loc[B.groupby("gene_id").p_value.idxmin(), ["gene_id", "symbol", "mask", "trait", "odds_ratio", "p_value"]]
best = best.set_index("gene_id").add_prefix("best_")
plof = B[B["mask"].str.startswith("M1.")]
best_plof = (plof.loc[plof.groupby("gene_id").p_value.idxmin(), ["gene_id", "mask", "trait", "p_value"]]
             .set_index("gene_id").add_prefix("pLOF_best_"))
# direction: logOR(definitely morning) - logOR(definitely evening) in the best mask
piv = B.pivot_table(index=["gene_id", "mask"], columns="trait", values="logOR")
dirn = []
for g, m in zip(best.index, best.best_mask):
    try:
        r = piv.loc[(g, m)]
        dirn.append(np.sign(r["definitely_morning"] - r["definitely_evening"]))
    except KeyError:
        dirn.append(np.nan)
caf = B[B["mask"] == "M3.1"].groupby("gene_id").effect_allele_frequency.max().rename("caf_M3_1")
G = pd.concat([per_trait.add_prefix("p_"), omni, omni_rare, best, best_plof, caf], axis=1)
G["direction"] = dirn
G["direction"] = G.direction.map({1.0: "morningness", -1.0: "eveningness"})
n_genes = len(G)
G["exome_wide_sig"] = G.p_burden < 0.05 / n_genes
log.info("genes with burden results: %d; exome-wide significant (p<%.2e): %s", n_genes, 0.05 / n_genes,
         ", ".join(G[G.exome_wide_sig].best_symbol))

# ------------------------------------------------------------------ convergence with GWAS loci
lg = pd.read_csv(INTERIM / "loci" / "locus_genes.tsv", sep="\t", dtype={"r2_level": str})
locus_genes = set(lg[lg.r2_level == "60"].gene_id)
leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t")
nearest = set(leads.nearest_gene.dropna())
gs = pd.read_csv(INTERIM / "genesets" / "genesets.tsv", sep="\t")
sets = {s: set(sub.gene_id) for s, sub in gs.groupby("set_id")}
G["in_locus"] = G.index.isin(locus_genes)
G["nearest_to_lead"] = G.index.isin(nearest)
G["mendelian"] = G.index.isin(sets["MEND_ANY"])
G["sleep_mendelian"] = G.index.isin(sets["HPO_HP0002360"])
G["intd_monoamine"] = G.index.isin(sets["EXP_INTD_MONOAMINE"])
G = G.dropna(subset=["p_burden", "p_burden_maf_le_0.1pct", "caf_M3_1"])
G["caf_bin"] = pd.qcut(np.log10(G.caf_M3_1.clip(lower=1e-7)), 10, labels=False, duplicates="drop")
G["score"] = -np.log10(G.p_burden)

tests = {
    "All genes in chronotype loci": G.in_locus,
    "Nearest gene to lead variant": G.nearest_to_lead,
    "Mendelian genes in chronotype loci": G.in_locus & G.mendelian,
    "Non-Mendelian genes in chronotype loci": G.in_locus & ~G.mendelian,
    "Sleep-disturbance Mendelian genes in loci": G.in_locus & G.sleep_mendelian,
    "Mendelian genes outside loci (reference)": ~G.in_locus & G.mendelian,
}
bins = G.groupby("caf_bin").indices
rows = []
for pcol, label in [("p_burden", "all masks"), ("p_burden_maf_le_0.1pct", "masks with MAF <= 0.1%")]:
    G["score"] = -np.log10(G[pcol])
    for name, mask in tests.items():
        sel = G[mask.values]
        obs_mean = sel.score.mean()
        obs_n01 = int((sel[pcol] < 0.01).sum())
        counts = sel.caf_bin.value_counts()
        draws_mean = np.empty(N_DRAW)
        draws_n01 = np.empty(N_DRAW)
        draws_trim = np.empty(N_DRAW)
        score = G.score.values
        p01 = (G[pcol].values < 0.01)
        idx_lists = [(bins[b], n) for b, n in counts.items()]
        # leave-top-3-out: the same trimming is applied to the observed set and to every matched draw
        K = TRIM if len(sel) > TRIM + 1 else 0
        for k in range(N_DRAW):
            pick = np.concatenate([rng.choice(ix, size=n, replace=False) for ix, n in idx_lists])
            sc = score[pick]
            draws_mean[k] = sc.mean()
            draws_n01[k] = p01[pick].sum()
            if K:
                draws_trim[k] = np.sort(sc)[:-K].mean()
        top = sel.nsmallest(K, pcol) if K else sel.iloc[:0]
        rest = sel.drop(top.index)
        obs_trim = rest.score.mean()
        # leave-three-out (conservative): the remaining genes vs untrimmed CAF-matched draws of the same size,
        # using a separate random stream so that the main results are unchanged
        l3o_mean, l3o_n01 = np.empty(N_DRAW), np.empty(N_DRAW)
        if K:
            idx_rest = [(bins[b], n) for b, n in rest.caf_bin.value_counts().items()]
            for k in range(N_DRAW):
                pick = np.concatenate([rng_l3o.choice(ix, size=n, replace=False) for ix, n in idx_rest])
                l3o_mean[k] = score[pick].mean()
                l3o_n01[k] = p01[pick].sum()
        obs_rest_n01 = int((rest[pcol] < 0.01).sum())
        rows.append({"masks": label, "gene_group": name, "n_genes": len(sel), "mean_neglog10p": obs_mean,
                     "null_mean_neglog10p": draws_mean.mean(),
                     "null_q025_neglog10p": np.quantile(draws_mean, 0.025),
                     "null_q975_neglog10p": np.quantile(draws_mean, 0.975),
                     "p_mean": (1 + (draws_mean >= obs_mean).sum()) / (1 + N_DRAW),
                     "n_p_lt_0.01": obs_n01, "exp_p_lt_0.01": draws_n01.mean(),
                     "p_count": (1 + (draws_n01 >= obs_n01).sum()) / (1 + N_DRAW),
                     "genes_p_lt_0.01": ";".join(sel[sel[pcol] < 0.01].best_symbol.astype(str)),
                     "trimmed_genes": ";".join(top.best_symbol.astype(str)),
                     "mean_neglog10p_trimmed": obs_trim if K else np.nan,
                     "null_mean_neglog10p_trimmed": draws_trim.mean() if K else np.nan,
                     "p_mean_trimmed": (1 + (draws_trim >= obs_trim).sum()) / (1 + N_DRAW) if K else np.nan,
                     "l3o_null_mean_neglog10p": l3o_mean.mean() if K else np.nan,
                     "l3o_p_mean": (1 + (l3o_mean >= obs_trim).sum()) / (1 + N_DRAW) if K else np.nan,
                     "l3o_n_p_lt_0.01": obs_rest_n01 if K else np.nan,
                     "l3o_exp_p_lt_0.01": l3o_n01.mean() if K else np.nan,
                     "l3o_p_count": (1 + (l3o_n01 >= obs_rest_n01).sum()) / (1 + N_DRAW) if K else np.nan,
                     "n_draws": N_DRAW})
G["score"] = -np.log10(G.p_burden)
conv = pd.DataFrame(rows)
conv.to_csv(TABLES / "burden_convergence.tsv", sep="\t", index=False)
G.sort_values("p_burden").to_csv(TABLES / "burden_gene_level.tsv", sep="\t")
log.info("\n%s", conv.drop(columns="genes_p_lt_0.01").to_string())
lam = np.median(-2 * np.log(G.p_burden)) / 1.386  # chi2(2) median for Fisher-like scale
save_json({"n_genes": n_genes, "bonferroni": 0.05 / n_genes,
           "exome_wide_significant": G[G.exome_wide_sig].best_symbol.tolist(),
           "median_p_burden": float(G.p_burden.median())}, TABLES / "burden_summary.json")
