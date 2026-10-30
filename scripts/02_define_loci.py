#!/usr/bin/env python3
"""Define LD-based loci around the 351 chronotype lead variants and map genes.

Lead variants: Jones et al. 2019 Supplementary Data 1 (UK Biobank + 23andMe meta-analysis).
LD reference: 1000 Genomes phase 3 European panel (n = 503), variants with MAF >= 1%.

Leads absent from the panel are represented by (i) a panel SNP at the same position, else
(ii) the most strongly associated panel variant within 100 kb in the UK Biobank-only
summary statistics.

Outputs (data/interim/loci/):
  leads.tsv           one row per lead with LD-index variant, features, intervals, genes
  regions.tsv         lead loci merged where intervals overlap (primary r2 >= 0.6)
  locus_genes.tsv     long table of region x gene (protein-coding) at each r2 level
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, ROOT, get_logger, load_config, save_json, chrom_norm  # noqa: E402
from chronomic import loci as L  # noqa: E402

log = get_logger("02_define_loci")
cfg = load_config()
OUT = INTERIM / "loci"
OUT.mkdir(parents=True, exist_ok=True)
BFILE = INTERIM / "ldref" / "eur_maf01"
FLANK = cfg["loci"]["flank_kb"] * 1000

# ------------------------------------------------------------------ leads
sd1 = pd.read_excel(ROOT / cfg["gwas"]["lead_snps_source"], header=2)
sd1 = sd1[pd.to_numeric(sd1["Locus No."], errors="coerce").notna()].copy()
leads = pd.DataFrame({
    "locus_no": sd1["Locus No."].astype(int),
    "lead_rsid": sd1["Lead Variant"].astype(str),
    "chrom": sd1["CHR"].astype(str).map(chrom_norm),
    "pos": sd1["POS"].astype(int),
    "allele1": sd1["Allele1"], "allele2": sd1["Allele2"], "freq1": sd1["Freq1"],
    "meta_lnOR": sd1["Meta ln OR"], "meta_z": sd1["Zscore (Corrected)"],
    "meta_p": pd.to_numeric(sd1["Meta P (Corrected)"], errors="coerce"),
    "gene_context": sd1["Gene Context"],
    "jones_locus_start": sd1["Locus Start Pos"].astype(int),
    "jones_locus_end": sd1["Locus End Pos"].astype(int),
    "novel": sd1["Novel Locus?"],
})
# P values below double precision were written as 0 in the spreadsheet; recompute from Z
from scipy.stats import norm  # noqa: E402
leads["meta_p_from_z"] = 2 * norm.sf(leads.meta_z.abs())
log.info("lead variants: %d", len(leads))

bim = pd.read_csv(BFILE.with_suffix(".bim"), sep="\t", header=None, usecols=[0, 1, 3],
                  names=["chrom", "snp", "bp"], dtype={"chrom": str, "snp": str, "bp": np.int64},
                  engine="pyarrow")
bim["chrom"] = bim.chrom.map(chrom_norm).astype("category")

leads["in_panel"] = leads.lead_rsid.isin(bim.snp)
need = leads[~leads.in_panel]
idx, how = [], []
ss = pd.DataFrame(columns=["SNP", "CHR", "BP", "P_BOLT_LMM"])
same_pos = pd.DataFrame(columns=["chrom", "bp", "snp"])
if len(need):
    log.info("%d leads absent from panel: %s", len(need), ", ".join(need.lead_rsid))
    need_keys = set(zip(need.chrom.astype(str), need.pos.astype(int)))
    cand = bim[bim.bp.isin(need.pos.astype(int))]
    cand = cand[[(str(c), int(b)) in need_keys for c, b in zip(cand.chrom, cand.bp)]]
    same_pos = pd.DataFrame({"chrom": cand.chrom.astype(str).values, "bp": cand.bp.values,
                             "snp": cand.snp.values})
    # UK Biobank-only statistics near the missing leads (chunked to keep memory low)
    keep = []
    for ch in pd.read_csv(ROOT / cfg["gwas"]["sumstats_primary"], sep="\t",
                          usecols=["SNP", "CHR", "BP", "P_BOLT_LMM"], dtype={"CHR": str},
                          chunksize=2_000_000):
        ch["CHR"] = ch.CHR.map(chrom_norm)
        m = np.zeros(len(ch), dtype=bool)
        for _, r in need.iterrows():
            m |= (ch.CHR.values == r.chrom) & (np.abs(ch.BP.values - r.pos) <= 100_000)
        keep.append(ch[m])
    ss = pd.concat(keep)
    ss = ss[ss.SNP.isin(bim.snp)]
for _, r in leads.iterrows():
    if r.in_panel:
        idx.append(r.lead_rsid); how.append("lead")
        continue
    same = same_pos[(same_pos.chrom.astype(str) == r.chrom) & (same_pos.bp == r.pos)]
    if len(same):
        idx.append(same.snp.iloc[0]); how.append("same_position")
        continue
    win = ss[(ss.CHR == r.chrom) & (ss.BP.between(r.pos - 100_000, r.pos + 100_000))]
    if len(win):
        best = win.sort_values("P_BOLT_LMM").iloc[0]
        idx.append(best.SNP); how.append(f"best_ukb_within_100kb(P={best.P_BOLT_LMM:.1e})")
    else:
        idx.append(None); how.append("none")
leads["ld_index"] = idx
leads["ld_index_type"] = how
log.info("LD index assignment: %s", leads.ld_index_type.str.split("(").str[0].value_counts().to_dict())
leads = leads[leads.ld_index.notna()].copy()
sub = bim[bim.snp.isin(leads.ld_index)].drop_duplicates("snp").set_index("snp")
leads["index_bp"] = leads.ld_index.map(sub.bp)
frq = pd.read_csv(BFILE.with_suffix(".afreq"), sep="\t", usecols=["ID", "ALT_FREQS"], engine="pyarrow")
frq = frq[frq.ID.isin(leads.ld_index)].set_index("ID").ALT_FREQS
leads["maf"] = leads.ld_index.map(np.minimum(frq, 1 - frq))
del bim

# ------------------------------------------------------------------ LD (r2 >= 0.5 within 1 Mb)
ld_path = L.plink_ld(BFILE, leads.ld_index.unique().tolist(), OUT / "leads_ld",
                     window_kb=cfg["ld"]["window_kb"], min_r2=0.5, dprime=True)
ld = L.read_ld(ld_path)
log.info("LD pairs (r2>=0.5): %d", len(ld))
summ = L.summarise_ld(ld, leads.drop_duplicates("ld_index").set_index("ld_index").index_bp)
leads = leads.merge(summ.drop(columns="bp"), left_on="ld_index", right_index=True, how="left")

# ------------------------------------------------------------------ genes
genes = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str})
gi = L.GeneIndex(genes)
near = gi.nearest(leads.chrom, leads.pos)
leads["nearest_gene"] = near.nearest_gene.values
leads["nearest_gene_symbol"] = leads.nearest_gene.map(genes.set_index("gene_id").symbol)
leads["nearest_dist"] = near.nearest_dist.values
long_rows = []
for lv in L.R2_LEVELS:
    tag = f"{int(round(lv * 100)):02d}"
    M = gi.overlap_matrix(leads.chrom, leads[f"start_r{tag}"], leads[f"end_r{tag}"], FLANK)
    leads[f"n_genes_r{tag}"] = np.asarray(M.sum(axis=1)).ravel()
    for i, row in enumerate(M.tolil().rows):
        for j in row:
            long_rows.append((leads.locus_no.iloc[i], tag, genes.gene_id.iloc[j], genes.symbol.iloc[j]))
lg = pd.DataFrame(long_rows, columns=["locus_no", "r2_level", "gene_id", "symbol"])

# ------------------------------------------------------------------ merge overlapping loci into regions
leads["region"] = L.merge_intervals(leads, "start_r60", "end_r60", FLANK)
reg = (leads.sort_values("meta_p_from_z")
       .groupby("region")
       .agg(chrom=("chrom", "first"), start=("start_r60", "min"), end=("end_r60", "max"),
            n_leads=("locus_no", "size"), index_locus=("locus_no", "first"),
            index_snp=("ld_index", "first"), maf=("maf", "first"), n_buddies=("n_buddies", "first"),
            nearest_dist=("nearest_dist", "first"), leads=("lead_rsid", lambda s: ";".join(s)))
       .reset_index())
Mreg = gi.overlap_matrix(reg.chrom, reg.start, reg.end, FLANK)
reg["n_genes"] = np.asarray(Mreg.sum(axis=1)).ravel()
reg["genes"] = [";".join(genes.symbol.iloc[r]) for r in Mreg.tolil().rows]
lg = lg.merge(leads[["locus_no", "region"]], on="locus_no")

leads.to_csv(OUT / "leads.tsv", sep="\t", index=False)
reg.to_csv(OUT / "regions.tsv", sep="\t", index=False)
lg.to_csv(OUT / "locus_genes.tsv", sep="\t", index=False)
stats = {
    "n_leads_input": int(len(sd1)), "n_leads_mapped": int(len(leads)),
    "index_types": leads.ld_index_type.str.split("(").str[0].value_counts().to_dict(),
    "n_regions": int(len(reg)), "n_multi_lead_regions": int((reg.n_leads > 1).sum()),
    "median_region_kb": float(((reg.end - reg.start) / 1000).median()),
    "regions_with_no_gene": int((reg.n_genes == 0).sum()),
    "median_genes_per_region": float(reg.n_genes.median()),
    "n_unique_genes_r60": int(lg[lg.r2_level == "60"].gene_id.nunique()),
}
save_json(stats, OUT / "loci_summary.json")
log.info("summary: %s", stats)
