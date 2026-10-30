#!/usr/bin/env python3
"""Build the gene universe and all gene sets used for enrichment testing.

Outputs (data/interim/genesets/):
  gene_universe.tsv          protein-coding genes (GENCODE v48, GRCh37 lift) with HGNC symbols
  gene_disease_links.tsv     every gene-disease pair used (OMIM via HPO, Orphanet, ClinGen)
  genesets.tsv               long format: set_id, family, label, gene_id
  geneset_summary.tsv        set sizes and definitions
  genesets.magma.txt         MAGMA set-annotation file (one set per line)
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx
import obonet
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic import annot  # noqa: E402
from chronomic.common import INTERIM, RAW, get_logger, load_config, save_json  # noqa: E402

log = get_logger("01_build_gene_sets")
cfg = load_config()
OUT = INTERIM / "genesets"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ universe
pc = annot.protein_coding_genes(include_x=True)
sym = annot.symbol_of(pc.gene_id)
pc["symbol"] = pc.gene_id.map(sym)
universe = set(pc.gene_id)
pc.to_csv(OUT / "gene_universe.tsv", sep="\t", index=False)
log.info("protein-coding universe: %d genes (autosomes + X)", len(universe))

# ------------------------------------------------------------------ gene-disease links
omim = annot.omim_gene_disease()
omim = omim[omim.assoc == "MENDELIAN"].copy()
omim["source"] = "OMIM"
hpoa = annot.hpo_annotations()
names = hpoa.drop_duplicates("database_id").set_index("database_id").disease_name.to_dict()
omim["disease_name"] = omim.disease_id.map(names)

orpha = annot.orphanet_gene_disease()
orpha_causal = orpha[orpha.assoc_type.isin(annot.ORPHA_CAUSAL) & (orpha.status == "Assessed")].copy()
orpha_causal["source"] = "Orphanet"

clingen = annot.clingen_validity()
clingen_ok = clingen[clingen.classification.isin(cfg["mendelian"]["clingen_levels"])].copy()
clingen_ok["source"] = "ClinGen"
clingen_ok = clingen_ok.rename(columns={"mondo_id": "disease_id"})

links = pd.concat([
    omim[["gene_id", "disease_id", "disease_name", "source"]],
    orpha_causal[["gene_id", "disease_id", "disease_name", "source"]],
    clingen_ok[["gene_id", "disease_id", "disease_name", "source"]].assign(
        disease_name=lambda d: d.disease_name + " [" + clingen_ok.classification + "]"),
], ignore_index=True).drop_duplicates()
links = links[links.gene_id.isin(universe)]
links["symbol"] = links.gene_id.map(sym)
links.to_csv(OUT / "gene_disease_links.tsv", sep="\t", index=False)
log.info("gene-disease links: %d pairs, %d genes", len(links), links.gene_id.nunique())

# ------------------------------------------------------------------ disease -> HPO terms (with ancestors)
G = annot.hpo_graph()
term_anc_cache: dict[str, set[str]] = {}


def with_ancestors(t: str) -> set[str]:
    if t not in term_anc_cache:
        term_anc_cache[t] = (nx.descendants(G, t) | {t}) if t in G else {t}  # child->parent edges
    return term_anc_cache[t]


dis_terms: dict[str, set[str]] = defaultdict(set)
for d, t in zip(hpoa.database_id, hpoa.hpo_id):
    dis_terms[d] |= with_ancestors(t)

# Mendelian diseases with phenotype annotations, linked to genes (OMIM MENDELIAN + Orphanet causal)
gd = pd.concat([omim[["gene_id", "disease_id"]], orpha_causal[["gene_id", "disease_id"]]]).drop_duplicates()
gd = gd[gd.gene_id.isin(universe)]
gene_terms: dict[str, set[str]] = defaultdict(set)
for g, d in zip(gd.gene_id, gd.disease_id):
    gene_terms[g] |= dis_terms.get(d, set())

inh = annot.hpo_inheritance()
INH = {"HP:0000006": "AD", "HP:0000007": "AR", "HP:0001417": "XL"}
gene_inh: dict[str, set[str]] = defaultdict(set)
dis_inh = defaultdict(set)
for d, t in zip(inh.database_id, inh.hpo_id):
    for code, lab in INH.items():
        if code in with_ancestors(t):
            dis_inh[d].add(lab)
for g, d in zip(gd.gene_id, gd.disease_id):
    gene_inh[g] |= dis_inh.get(d, set())

# ------------------------------------------------------------------ assemble sets
sets: list[tuple[str, str, str, set[str]]] = []  # id, family, label, genes


def add(set_id, family, label, genes):
    genes = set(genes) & universe
    sets.append((set_id, family, label, genes))


add("MEND_ANY", "mendelian_source", "Any Mendelian disease gene (OMIM or Orphanet)",
    set(omim.gene_id) | set(orpha_causal.gene_id))
add("MEND_OMIM", "mendelian_source", "OMIM Mendelian disease gene", set(omim.gene_id))
add("MEND_ORPHA", "mendelian_source", "Orphanet disease-causing gene (assessed)", set(orpha_causal.gene_id))
add("MEND_CLINGEN", "mendelian_source", "ClinGen gene-disease validity Moderate/Strong/Definitive",
    set(clingen_ok.gene_id))
add("MEND_HICONF", "mendelian_source", "High-confidence Mendelian (ClinGen M/S/D or in both OMIM and Orphanet)",
    set(clingen_ok.gene_id) | (set(omim.gene_id) & set(orpha_causal.gene_id)))

for code, lab in [("AD", "autosomal dominant"), ("AR", "autosomal recessive"), ("XL", "X-linked")]:
    add(f"INH_{code}", "inheritance", f"Mendelian gene with {lab} disease",
        {g for g, s in gene_inh.items() if code in s})

# HPO organ systems = children of 'Phenotypic abnormality'
systems = sorted(G.predecessors("HP:0000118"), key=lambda t: G.nodes[t]["name"])
for t in systems:
    genes = {g for g, ts in gene_terms.items() if t in ts}
    add(f"SYS_{t.replace(':', '')}", "hpo_system", G.nodes[t]["name"], genes)

# A priori phenotype sets, specified before the enrichment analyses (not pre-registered): phenotypes with
# established links to circadian physiology (sleep; cardiac electrophysiology; complement and immunity;
# neurodevelopment; growth)
APRIORI = [
    ("HP:0002360", "Sleep disturbance"),
    ("HP:0011675", "Arrhythmia"),
    ("HP:0001657", "Prolonged QT interval"),
    ("HP:0001638", "Cardiomyopathy"),
    ("HP:0001639", "Hypertrophic cardiomyopathy"),
    ("HP:0004431", "Reduced circulating complement concentration"),
    ("HP:0002721", "Immunodeficiency"),
    ("HP:0001263", "Global developmental delay"),
    ("HP:0001249", "Intellectual disability"),
    ("HP:0001250", "Seizure"),
    ("HP:0004322", "Short stature"),
]
for t, lab in APRIORI:
    genes = {g for g, ts in gene_terms.items() if t in ts}
    add(f"HPO_{t.replace(':', '')}", "hpo_apriori", f"{G.nodes[t]['name']}", genes)

# Circadian Mendelian set: sleep-wake cycle disturbance OR abnormal pineal melatonin secretion
circ_terms = {"HP:0006979", "HP:0012689"}
add("HPO_CIRCADIAN", "hpo_apriori", "Circadian phenotype (sleep-wake cycle disturbance or abnormal melatonin)",
    {g for g, ts in gene_terms.items() if ts & circ_terms})

# Curated complement deficiency genes (IUIS 2022 classification of inborn errors of immunity, complement
# deficiencies). HPO annotations of several complement deficiencies (e.g. C2 deficiency, OMIM 217000) carry no
# complement term, so the HPO-based set above misses them.
IUIS_COMPLEMENT = ["C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B",
                   "C8G", "C9", "MASP2", "FCN3", "SERPING1", "CFB", "CFD", "CFP", "CFI", "CFH", "CFHR1", "CFHR2",
                   "CFHR3", "CFHR4", "CFHR5", "THBD", "CD46", "CD55", "CD59"]
add("CUR_COMPLEMENT_IUIS", "curated_apriori", "Complement deficiency genes (IUIS 2022 classification)",
    {annot.id_maps()["symbol"][s] for s in IUIS_COMPLEMENT if s in annot.id_maps()["symbol"]})

# Auxiliary sets for confounding checks
gdd = {g for g, ts in gene_terms.items() if "HP:0001263" in ts}
idis = {g for g, ts in gene_terms.items() if "HP:0001249" in ts}
seiz = {g for g, ts in gene_terms.items() if "HP:0001250" in ts}
sleep = {g for g, ts in gene_terms.items() if "HP:0002360" in ts}
add("NDD", "hpo_auxiliary", "Neurodevelopmental (global developmental delay or intellectual disability)", gdd | idis)
add("HPO_SLEEP_noNDD", "hpo_auxiliary", "Sleep disturbance, excluding genes with developmental delay, ID or seizures",
    sleep - (gdd | idis | seiz))

# Exploratory (defined after inspecting the genes driving the sleep-disturbance signal)
INTD_MONOAMINE = ["TH", "DDC", "GCH1", "PTS", "SPR", "QDPR", "PCBD1", "DNAJC12", "SLC6A3", "SLC18A2", "DBH", "MAOA"]
add("EXP_INTD_MONOAMINE", "exploratory", "Inherited monoamine neurotransmitter disorder genes (iNTD; curated)",
    {annot.id_maps()["symbol"][s] for s in INTD_MONOAMINE})
add("EXP_CSF_HVA", "exploratory", "Decreased CSF homovanillic acid (HP:0003785)",
    {g for g, ts in gene_terms.items() if "HP:0003785" in ts})

# Gene Ontology 'circadian rhythm' (GO:0007623) incl. descendants -- non-Mendelian positive control
go = obonet.read_obo(RAW / "go" / "go-basic.obo")
circ_go = nx.ancestors(go, "GO:0007623") | {"GO:0007623"}
m = annot.id_maps()
go_genes = set()
with gzip.open(RAW / "go" / "goa_human.gaf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        f = line.split("\t")
        if f[4] in circ_go and "NOT" not in f[3]:
            gid = m["symbol"].get(f[2])
            if gid:
                go_genes.add(gid)
add("GO_CIRCADIAN", "positive_control", "GO: circadian rhythm (GO:0007623 and descendants)", go_genes)

# constraint (general GWAS-locus property; used as a covariate/contrast)
thr = {}
for ver, sid in [("v2.1.1", "LOEUF_D1"), ("v4.1", "LOEUF4_D1")]:
    con = annot.gnomad_constraint(ver)
    con = con[con.gene_id.isin(universe)].dropna(subset=["LOEUF"])
    thr[ver] = float(con.LOEUF.quantile(0.1))
    add(sid, "constraint", f"Most LoF-constrained decile (gnomAD {ver} LOEUF < {thr[ver]:.3f})",
        set(con.loc[con.LOEUF < thr[ver], "gene_id"]))

# ------------------------------------------------------------------ filter & write
lo, hi = cfg["mendelian"]["min_genes_per_set"], cfg["mendelian"]["max_genes_per_set"]
rows, summ = [], []
for sid, fam, lab, genes in sets:
    keep = lo <= len(genes) <= max(hi, len(universe))
    if fam in ("mendelian_source",) or keep:
        rows += [(sid, fam, lab, g) for g in sorted(genes)]
    summ.append((sid, fam, lab, len(genes), keep or fam == "mendelian_source"))
gs = pd.DataFrame(rows, columns=["set_id", "family", "label", "gene_id"])
gs.to_csv(OUT / "genesets.tsv", sep="\t", index=False)
summary = pd.DataFrame(summ, columns=["set_id", "family", "label", "n_genes", "tested"])
summary.to_csv(OUT / "geneset_summary.tsv", sep="\t", index=False)
with open(OUT / "genesets.magma.txt", "w") as fh:
    for sid, sub in gs.groupby("set_id", sort=False):
        fh.write(sid + " " + " ".join(sub.gene_id) + "\n")

save_json({"hpo_version": annot.hpo_version(), "n_universe": len(universe),
           "n_links": int(len(links)), "loeuf_decile_thresholds": thr}, OUT / "provenance.json")
log.info("\n%s", summary.to_string())
