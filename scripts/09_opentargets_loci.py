#!/usr/bin/env python3
"""Link the 351 chronotype loci to Open Targets fine-mapped credible sets, L2G genes and
colocalisations.

A credible set is linked to a locus when any of its 95% credible-set variants is the locus lead
(or its LD index) or is in r2 >= 0.6 (1000G EUR) with it.

Outputs: results/tables/locus_opentargets.tsv (one row per locus),
         results/tables/locus_qtl_coloc.tsv (colocalisations with eCAVIAR CLPP >= 0.01, or coloc H4 >= 0.8 for
         standard coloc; gene-level QTL evidence requires CLPP >= 0.1),
         results/tables/opentargets_summary.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import INTERIM, RAW, TABLES, get_logger, save_json  # noqa: E402

log = get_logger("09_opentargets_loci")
D = json.load(open(RAW / "opentargets" / "credible_sets.json"))
meta = D["meta"]
leads = pd.read_csv(INTERIM / "loci" / "leads.tsv", sep="\t")
ld = pd.read_csv(INTERIM / "loci" / "leads_ld.ld", sep=r"\s+")
ld = ld[ld.R2 >= 0.6]
tag2locus = defaultdict(set)
idx2locus = dict(zip(leads.ld_index, leads.locus_no))
for a, b in zip(ld.SNP_A, ld.SNP_B):
    if a in idx2locus:
        tag2locus[b].add(idx2locus[a])
for r in leads.itertuples():
    tag2locus[r.lead_rsid].add(r.locus_no)
    tag2locus[r.ld_index].add(r.locus_no)

gs = pd.read_csv(INTERIM / "genesets" / "genesets.tsv", sep="\t")
mend = set(gs.loc[gs.set_id == "MEND_ANY", "gene_id"])

links, coloc_rows, l2g_rows = [], [], []
for cs in D["credible_sets"]:
    rs = set(cs["variant"]["rsIds"] or [])
    for row in cs["locus"]["rows"]:
        if row.get("is95CredibleSet") in (True, None):
            rs |= set(row["variant"]["rsIds"] or [])
    loci = set().union(*[tag2locus[r] for r in rs if r in tag2locus]) if rs else set()
    l2g = sorted(cs["l2GPredictions"]["rows"], key=lambda x: -x["score"])
    top = l2g[0] if l2g else None
    for L in loci:
        l2g_rows += [{"locus_no": L, "gene_id": x["target"]["id"], "symbol": x["target"]["approvedSymbol"],
                      "l2g_score": x["score"]} for x in l2g]
        links.append({"locus_no": L, "studyLocusId": cs["studyLocusId"], "studyId": cs["studyId"],
                      "finemapping": cs["finemappingMethod"], "cs_lead": ";".join(cs["variant"]["rsIds"] or []),
                      "l2g_top_gene_id": top["target"]["id"] if top else None,
                      "l2g_top_symbol": top["target"]["approvedSymbol"] if top else None,
                      "l2g_top_score": top["score"] if top else np.nan,
                      "l2g_genes_ge0.5": ";".join(x["target"]["approvedSymbol"] for x in l2g if x["score"] >= 0.5)})
        for c in cs["colocalisation"]["rows"]:
            o = c["otherStudyLocus"]
            st = o["study"] or {}
            strong = (c["h4"] or 0) >= 0.8 if c["colocalisationMethod"].upper().startswith("COLOC") and \
                "ECAVIAR" not in c["colocalisationMethod"].upper() else (c["clpp"] or 0) >= 0.01
            if not strong:
                continue
            coloc_rows.append({"locus_no": L, "chronotype_study": cs["studyId"], "method": c["colocalisationMethod"],
                               "h4": c["h4"], "clpp": c["clpp"], "right_type": c["rightStudyType"],
                               "right_study": o["studyId"], "trait_or_project": st.get("traitFromSource") or st.get("projectId"),
                               "qtl_gene": (st.get("target") or {}).get("approvedSymbol"),
                               "qtl_gene_id": (st.get("target") or {}).get("id"),
                               "biosample": (st.get("biosample") or {}).get("biosampleName")})
links = pd.DataFrame(links)
l2g_genes = (pd.DataFrame(l2g_rows, columns=["locus_no", "gene_id", "symbol", "l2g_score"])
             .groupby(["locus_no", "gene_id", "symbol"], as_index=False).l2g_score.max())
l2g_genes.to_csv(TABLES / "locus_l2g_genes.tsv", sep="\t", index=False)
coloc = pd.DataFrame(coloc_rows).drop_duplicates()
coloc.to_csv(TABLES / "locus_qtl_coloc.tsv", sep="\t", index=False)

best = (links.sort_values("l2g_top_score", ascending=False).groupby("locus_no")
        .agg(n_credible_sets=("studyLocusId", "nunique"), studies=("studyId", lambda s: ";".join(sorted(set(s)))),
             l2g_gene_id=("l2g_top_gene_id", "first"), l2g_gene=("l2g_top_symbol", "first"),
             l2g_score=("l2g_top_score", "first"),
             l2g_genes_ge05=("l2g_genes_ge0.5", lambda s: ";".join(sorted({g for x in s.dropna() for g in x.split(";") if g})))))
CLPP_STRONG = 0.1
# Open Targets 26.09 reports COLOC_PIP_ECAVIAR for these credible sets: a PIP-based coloc H4 normalised so that
# H3 + H4 = 1 (less stringent than standard coloc PP.H4) and the eCAVIAR CLPP. Gene-level QTL support therefore
# uses CLPP alone.
coloc["strong"] = coloc.clpp.fillna(0) >= CLPP_STRONG
qtl_all = coloc[coloc.right_type.str.contains("qtl", case=False, na=False)]
qtl = qtl_all[qtl_all.strong]
qtl_genes = qtl.groupby("locus_no").qtl_gene.apply(lambda s: ";".join(sorted(set(s.dropna()))))
qtl_genes_weak = qtl_all.groupby("locus_no").qtl_gene.apply(lambda s: ";".join(sorted(set(s.dropna()))))
gw = coloc[coloc.right_type == "gwas"].groupby("locus_no").trait_or_project.apply(
    lambda s: ";".join(sorted(set(s.dropna()))[:40]))
out = leads[["locus_no", "lead_rsid", "chrom", "pos", "nearest_gene_symbol"]].merge(best, on="locus_no", how="left")
out["qtl_coloc_genes"] = out.locus_no.map(qtl_genes)
out["qtl_coloc_genes_clpp001"] = out.locus_no.map(qtl_genes_weak)
out["gwas_coloc_traits"] = out.locus_no.map(gw)
out["l2g_gene_is_mendelian"] = out.l2g_gene_id.isin(mend)
out.to_csv(TABLES / "locus_opentargets.tsv", sep="\t", index=False)

linked = out.n_credible_sets.notna()
summ = {"open_targets_api": meta["meta"]["apiVersion"], "open_targets_data": meta["meta"]["dataVersion"],
        "retrieved_utc": meta["retrieved_utc"], "n_credible_sets": len(D["credible_sets"]),
        "loci_linked": int(linked.sum()), "loci_with_l2g_gene": int(out.l2g_gene.notna().sum()),
        "loci_l2g_ge_0.5": int((out.l2g_score >= 0.5).sum()),
        "l2g_gene_equals_nearest": int((out.l2g_gene == out.nearest_gene_symbol).sum()),
        "l2g_genes_mendelian": int(out.l2g_gene_is_mendelian.sum()),
        "loci_with_qtl_coloc": int(out.qtl_coloc_genes.notna().sum()),
        "loci_with_qtl_coloc_clpp001": int(out.qtl_coloc_genes_clpp001.notna().sum()),
        "qtl_genes_clpp01": int(qtl.qtl_gene.nunique()), "qtl_genes_clpp001": int(qtl_all.qtl_gene.nunique()),
        "coloc_methods": coloc.method.value_counts().to_dict(),
        "coloc_rows": int(len(coloc)), "coloc_rows_qtl": int(len(qtl_all)),
        "coloc_min_h4_pip": float(coloc.h4.min()), "coloc_min_clpp": float(coloc.clpp.min()),
        "credible_sets_by_study_method": {f"{s}|{m}": int(n) for (s, m), n in
                                          pd.DataFrame([(c["studyId"], c["finemappingMethod"]) for c in D["credible_sets"]],
                                                       columns=["s", "m"]).value_counts().items()}}
save_json(summ, TABLES / "opentargets_summary.json")
log.info("%s", summ)
