#!/usr/bin/env python3
"""Colocalisation of the chronotype and QT-interval signals at KCNH2 (long QT syndrome 2 gene).

Chronotype: UK Biobank continuous chronotype (Jones et al. 2019; GCST007576), GRCh37.
QT interval: GERA + UK Biobank + PAGE meta-analysis (Hoffmann et al. 2022; GCST90165291), GRCh37.
Method: coloc ABF (single causal variant), priors p1 = p2 = 1e-4, p12 = 5e-6, windows of
+/-100, 250 and 500 kb around the chronotype lead rs2072413. NB the two GWAS share UK Biobank
participants, which coloc tolerates (it compares association patterns, not independent effects).
Conditional analysis: independent QT signals by stepwise selection on approximate conditional z-scores (1000G EUR
signed LD), then coloc of chronotype with each QT signal conditioned on the others, over a grid of windows, LD
regularisation, conditioning sets and priors.
Outputs: results/tables/coloc_kcnh2.tsv, coloc_kcnh2_snps.tsv, coloc_kcnh2_ld.tsv, coloc_kcnh2_conditional.tsv,
         coloc_kcnh2_conditional_signals.tsv, coloc_kcnh2_conditional_snps.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.coloc import coloc_abf, labf_quant, sdy_est, snp_pp_h4  # noqa: E402
from chronomic.common import RAW, ROOT, TABLES, get_logger, load_config  # noqa: E402

log = get_logger("11_coloc_kcnh2")
cfg = load_config()
CHR, LEAD_BP, LEAD = "7", 150647969, "rs2072413"
WMAX = 500_000


def region(path, chrom_col, pos_col, sep="\t", **kw):
    keep = []
    for ch in pd.read_csv(path, sep=sep, chunksize=2_000_000, dtype={chrom_col: str}, **kw):
        m = (ch[chrom_col].astype(str) == CHR) & ch[pos_col].between(LEAD_BP - WMAX, LEAD_BP + WMAX)
        keep.append(ch[m])
    return pd.concat(keep)


cache1 = ROOT / "data" / "interim" / "coloc_chronotype_chr7.tsv"
if not cache1.exists():
    region(ROOT / cfg["gwas"]["sumstats_primary"], "CHR", "BP").to_csv(cache1, sep="\t", index=False)
chrono = pd.read_csv(cache1, sep="\t")
cache2 = ROOT / "data" / "interim" / "coloc_qt_chr7.tsv"
if not cache2.exists():
    region(RAW / "qt" / "GCST90165291_QT_buildGRCh37.tsv.gz", "chromosome", "base_pair_location").to_csv(
        cache2, sep="\t", index=False)
qt = pd.read_csv(cache2, sep="\t")

chrono = chrono.rename(columns={"SNP": "rsid", "BP": "pos", "ALLELE1": "ea", "ALLELE0": "oa", "A1FREQ": "eaf",
                                "BETA": "beta", "SE": "se", "P_BOLT_LMM": "p"})
# In GCST90165291_buildGRCh37 the beta and effect_allele_frequency columns refer to the allele listed in the
# other_allele column. This was checked at variants with well-established QT-prolonging alleles (NOS1AP rs12143842 T
# and rs10494366 G; KCNH2 rs1805123 G, which shortens QT) and by allele-frequency concordance with UK Biobank
# (see the log). The columns are therefore swapped here so that "ea" is the allele to which beta refers.
qt = qt.rename(columns={"variant_id_when_present": "rsid", "base_pair_location": "pos", "other_allele": "ea",
                        "effect_allele": "oa", "effect_allele_frequency": "eaf", "standard_error": "se",
                        "p_value": "p"})
m = chrono.merge(qt, on="pos", suffixes=("_c", "_q"))
same = (m.ea_c == m.ea_q) & (m.oa_c == m.oa_q)
flip = (m.ea_c == m.oa_q) & (m.oa_c == m.ea_q)
m = m[same | flip].copy()
m.loc[flip[same | flip].values, "beta_q"] *= -1
m = m[(m.eaf_c.between(0.01, 0.99)) & m.se_q.notna() & (m.se_q > 0)]
m["maf"] = np.minimum(m.eaf_c, 1 - m.eaf_c)
_fq = np.where(m.ea_c == m.ea_q, m.eaf_q, 1 - m.eaf_q)
log.info("allele-frequency concordance with UK Biobank (QT n >= 100,000): mean |diff| = %.3f",
         float(np.abs(m.eaf_c - _fq)[m.n >= 100_000].mean()))

n_c = cfg["gwas"]["sumstats_primary_n"]
rows, snp_rows = [], None
for w in (100_000, 250_000, 500_000):
    s = m[(m.pos - LEAD_BP).abs() <= w]
    sd1 = sdy_est(s.se_c.values ** 2, s.maf.values, n_c)
    sd2 = sdy_est(s.se_q.values ** 2, s.maf.values, s.n.values)
    l1 = labf_quant(s.beta_c.values, s.se_c.values, sd1)
    l2 = labf_quant(s.beta_q.values, s.se_q.values, sd2)
    res = coloc_abf(l1, l2)
    pp = snp_pp_h4(l1, l2)
    top = s.iloc[int(np.argmax(pp))]
    res |= {"window_kb": w // 1000, "sdY_chronotype": sd1, "sdY_qt": sd2,
            "min_p_chronotype": float(s.p_c.min()), "min_p_qt": float(s.p_q.min()),
            "top_shared_snp": top.rsid_c, "top_shared_snp_pp": float(pp.max()),
            "lead_in_data": bool((s.rsid_c == LEAD).any())}
    rows.append(res)
    if w == 250_000:
        snp_rows = s.assign(snp_pp_h4=pp)[["rsid_c", "pos", "ea_c", "beta_c", "se_c", "p_c", "beta_q", "se_q", "p_q",
                                          "snp_pp_h4"]].sort_values("snp_pp_h4", ascending=False)
R = pd.DataFrame(rows)
R.to_csv(TABLES / "coloc_kcnh2.tsv", sep="\t", index=False)
snp_rows.head(50).to_csv(TABLES / "coloc_kcnh2_snps.tsv", sep="\t", index=False)
lead = m[m.rsid_c == LEAD]
# LD (1000G EUR) between the chronotype lead, the top chronotype SNP and the top QT SNP (+/-250 kb)
import subprocess  # noqa: E402
w250 = m[(m.pos - LEAD_BP).abs() <= 250_000]
top_c = w250.sort_values("p_c").iloc[0].rsid_c
top_q = w250.sort_values("p_q").iloc[0].rsid_c
snps = [LEAD, top_c, top_q]
lst = ROOT / "data" / "interim" / "coloc_kcnh2_snps.txt"
lst.write_text("\n".join(dict.fromkeys(snps)) + "\n")
subprocess.run([str(ROOT / "bin" / "plink"), "--bfile", str(ROOT / "data" / "interim" / "ldref" / "eur_maf01"),
                "--extract", str(lst), "--r2", "inter-chr", "--ld-window-r2", "0", "--out",
                str(ROOT / "data" / "interim" / "coloc_kcnh2_ld")], check=True, stdout=subprocess.DEVNULL)
ldt = pd.read_csv(ROOT / "data" / "interim" / "coloc_kcnh2_ld.ld", sep=r"\s+")
r2 = {frozenset((a, b)): r for a, b, r in zip(ldt.SNP_A, ldt.SNP_B, ldt.R2)}
pd.DataFrame([{"lead": LEAD, "top_chronotype_snp": top_c, "top_qt_snp": top_q,
               "lead_p_chronotype": float(lead.p_c.iloc[0]), "lead_p_qt": float(lead.p_q.iloc[0]),
               "lead_beta_qt": float(lead.beta_q.iloc[0]),
               "top_chronotype_p": float(w250.p_c.min()), "top_qt_p": float(w250.p_q.min()),
               "r2_lead_topchrono": r2.get(frozenset((LEAD, top_c)), 1.0),
               "r2_lead_topqt": r2.get(frozenset((LEAD, top_q)), 1.0),
               "r2_topchrono_topqt": r2.get(frozenset((top_c, top_q)), 1.0)}]).to_csv(
    TABLES / "coloc_kcnh2_ld.tsv", sep="\t", index=False)
log.info("lead %s: chronotype beta=%.4f p=%.2g; QT beta=%.3f p=%.2g", LEAD,
         lead.beta_c.iloc[0], lead.p_c.iloc[0], lead.beta_q.iloc[0], lead.p_q.iloc[0])
log.info("\n%s", R.to_string())

# ------------------------------------------------------------------ conditional colocalisation
# The QT region contains several independent signals (Arking et al. 2014; Hoffmann et al. 2022), which violates
# coloc's single-causal-variant assumption. We identify independent QT signals by forward stepwise selection on
# approximate conditional z-scores using signed LD from the 1000G EUR panel
# (z_j|S = z_j - R_jS R_SS^-1 z_S, scaled by sqrt(1 - R_jS R_SS^-1 R_Sj)); see Yang et al. 2012 (GCTA-COJO).
# Each QT signal, conditioned on all others, is then tested for colocalisation with chronotype.
# Caveats: the QT meta-analysis is multi-ancestry (GERA, UK Biobank, PAGE) and the reference panel has 503
# Europeans, so the conditional statistics are approximate; we restrict to variants genotyped in all cohorts
# (n >= 95% of the maximum) and report a grid of sensitivity analyses.
import itertools  # noqa: E402


def load_window(W):
    s = m[(m.pos - LEAD_BP).abs() <= W].copy()
    s = s[s.n >= 0.95 * m.n.max()]
    win = ROOT / "data" / "interim" / f"coloc_kcnh2_window{W // 1000}"
    win.with_suffix(".snps").write_text("\n".join(s.rsid_c) + "\n")
    subprocess.run([str(ROOT / "bin" / "plink"), "--bfile", str(ROOT / "data" / "interim" / "ldref" / "eur_maf01"),
                    "--extract", str(win.with_suffix(".snps")), "--r", "square", "--make-just-bim", "--out", str(win)],
                   check=True, stdout=subprocess.DEVNULL)
    bim = pd.read_csv(win.with_suffix(".bim"), sep="\t", header=None, names=["chr", "snp", "cm", "bp", "a1", "a2"])
    LDm = np.loadtxt(win.with_suffix(".ld"))
    pos_in = {r: i for i, r in enumerate(bim.snp)}
    s = s[s.rsid_c.isin(pos_in)].copy()
    b = bim.set_index("snp").loc[s.rsid_c]
    s["sign"] = np.where(s.ea_c.values == b.a1.values, 1.0, np.where(s.ea_c.values == b.a2.values, -1.0, np.nan))
    s = s[~np.isnan(s.sign)].reset_index(drop=True)
    ii = np.array([pos_in[r] for r in s.rsid_c])
    return s, LDm[np.ix_(ii, ii)]


def conditional(z, R, S):
    if not S:
        return np.array(z, copy=True), np.ones_like(z)
    RS = R[:, S]
    inv = np.linalg.inv(R[np.ix_(S, S)])
    zc_ = z - RS @ inv @ z[S]
    var = np.clip(1 - np.einsum("ij,jk,ik->i", RS, inv, RS), 1e-6, None)
    return zc_ / np.sqrt(var), var


def stepwise(z, R, thr=5.45, max_k=10):
    S = []
    for _ in range(max_k):
        zz, var = conditional(z, R, S)
        zz[S] = 0
        zz[var < 0.05] = 0
        j = int(np.argmax(np.abs(zz)))
        if abs(zz[j]) < thr:
            break
        S.append(j)
    return S


def coloc_conditioned(s, R, zc, zq, others, p12=5e-6):
    zq_c, var = conditional(zq, R, others)
    keep = var > 0.05
    se_q_c = s.se_q.values / np.sqrt(var)
    sd1 = sdy_est(s.se_c.values[keep] ** 2, s.maf.values[keep], n_c)
    sd2 = sdy_est(se_q_c[keep] ** 2, s.maf.values[keep], s.n.values[keep])
    l1 = labf_quant(zc[keep] * s.se_c.values[keep], s.se_c.values[keep], sd1)
    l2 = labf_quant(zq_c[keep] * se_q_c[keep], se_q_c[keep], sd2)
    res = coloc_abf(l1, l2, p12=p12)
    pp = snp_pp_h4(l1, l2)
    return res, s.rsid_c[keep].iloc[int(np.argmax(pp))], float(pp.max()), zq_c


rows, signal_rows = [], []
for W, ridge, thr in itertools.product((100_000, 250_000), (1e-3, 1e-2, 5e-2), (5.45,)):
    s, R0 = load_window(W)
    R = (R0 + np.eye(len(R0)) * ridge) / (1 + ridge)
    zc = (s.beta_c / s.se_c * s.sign).to_numpy(dtype=float, copy=True)
    zq = (s.beta_q / s.se_q * s.sign).to_numpy(dtype=float, copy=True)
    S = stepwise(zq, R, thr)
    names = s.rsid_c.tolist()
    li = names.index(LEAD) if LEAD in names else None
    # the QT signal that tags the chronotype lead (highest LD with rs2072413)
    tag = max(S, key=lambda j: R0[j, li] ** 2) if li is not None and S else None
    for j in S:
        others = [x for x in S if x != j]
        res, top, top_pp, zq_c = coloc_conditioned(s, R, zc, zq, others)
        signal_rows.append({"window_kb": W // 1000, "ridge": ridge, "qt_signal": names[j],
                            "qt_z_marginal": float(zq[j]), "qt_z_conditional": float(zq_c[j]),
                            "r2_with_chronotype_lead": float(R0[j, li] ** 2) if li is not None else np.nan,
                            **res, "top_shared_snp": top, "top_shared_pp": top_pp})
    variants = {"all_other_signals": [x for x in S if x != tag],
                "top_qt_signal_only": [S[0]] if S and S[0] != tag else [],
                "all_other_signals_drop_weakest": [x for x in S if x != tag][:-1]}
    if W == 250_000 and ridge == 1e-3:
        _, _, _, zq_c_all = coloc_conditioned(s, R, zc, zq, variants["all_other_signals"])
        _, _, _, zq_c_top = coloc_conditioned(s, R, zc, zq, variants["top_qt_signal_only"])
        # z-scores are for the chronotype GWAS effect allele (ea_c)
        snp_tab = pd.DataFrame({"rsid": s.rsid_c, "pos": s.pos, "allele": s.ea_c, "z_chronotype": zc * s.sign,
                                "z_qt": zq * s.sign, "z_qt_cond_other_signals": zq_c_all * s.sign,
                                "z_qt_cond_top_signal": zq_c_top * s.sign,
                                "is_qt_signal": s.rsid_c.isin([names[x] for x in S])})
        snp_tab.to_csv(TABLES / "coloc_kcnh2_conditional_snps.tsv", sep="\t", index=False)
    for label, others in variants.items():
        for p12 in (5e-6, 1e-5):
            res, top, top_pp, zq_c = coloc_conditioned(s, R, zc, zq, others, p12=p12)
            rows.append({"window_kb": W // 1000, "ridge": ridge, "conditioning": label, "p12": p12,
                         "n_qt_signals": len(S), "qt_signals": ";".join(names[x] for x in S),
                         "tagging_signal": names[tag] if tag is not None else None,
                         "conditioned_on": ";".join(names[x] for x in others),
                         "qt_z_lead_marginal": float(zq[li]), "qt_z_lead_conditional": float(zq_c[li]),
                         **res, "top_shared_snp": top, "top_shared_pp": top_pp, "nsnps_window": len(s)})
C = pd.DataFrame(rows)
C.to_csv(TABLES / "coloc_kcnh2_conditional.tsv", sep="\t", index=False)
pd.DataFrame(signal_rows).to_csv(TABLES / "coloc_kcnh2_conditional_signals.tsv", sep="\t", index=False)
prim = C[(C.window_kb == 250) & (C.ridge == 1e-3) & (C.conditioning == "all_other_signals") & (C.p12 == 5e-6)].iloc[0]
log.info("primary conditional coloc: tagging signal %s, PP.H4=%.3f (range over %d sensitivity analyses %.2f-%.2f)",
         prim.tagging_signal, prim["PP.H4"], len(C), C["PP.H4"].min(), C["PP.H4"].max())
log.info("\n%s", C[["window_kb", "ridge", "conditioning", "p12", "n_qt_signals", "tagging_signal", "qt_z_lead_marginal",
                   "qt_z_lead_conditional", "PP.H3", "PP.H4", "top_shared_snp"]].to_string())

# ------------------------------------------------------------------ direction of effect at the shared signal
# Chronotype direction from the Jones et al. meta-analysis (ln OR for being a morning person, allele1), aligned to
# the UK Biobank effect allele; QT direction from the marginal and conditional QT z-scores for the same allele.
leads_tab = pd.read_csv(ROOT / "data" / "interim" / "loci" / "leads.tsv", sep="\t").set_index("lead_rsid")
lr = leads_tab.loc[LEAD]
st = pd.read_csv(TABLES / "coloc_kcnh2_conditional_snps.tsv", sep="\t").set_index("rsid")
dr = []
for rs in [LEAD, prim.tagging_signal]:
    r = st.loc[rs]
    lnor = float(lr.meta_lnOR) if r.allele == lr.allele1 else (-float(lr.meta_lnOR) if r.allele == lr.allele2 else np.nan)
    dr.append({"rsid": rs, "allele": r.allele, "z_chronotype_ukb": r.z_chronotype,
               "meta_lnOR_morning_person_lead": lnor if rs == LEAD else np.nan,
               "beta_qt_ms_marginal": float(m.set_index("rsid_c").loc[rs, "beta_q"]),
               "z_qt_marginal": r.z_qt, "z_qt_cond_top_signal": r.z_qt_cond_top_signal,
               "z_qt_cond_other_signals": r.z_qt_cond_other_signals})
D = pd.DataFrame(dr)
D.to_csv(TABLES / "coloc_kcnh2_direction.tsv", sep="\t", index=False)
log.info("direction at the shared signal:\n%s", D.to_string())

