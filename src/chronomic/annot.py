"""Readers for gene models and Mendelian disease knowledge bases.

All gene-level outputs are keyed on Ensembl gene IDs (unversioned) resolved
through the HGNC complete set, so OMIM/HPO (NCBI Gene IDs), Orphanet
(HGNC/Ensembl cross-references) and ClinGen (HGNC IDs) can be combined.
"""
from __future__ import annotations

import gzip
import re
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import networkx as nx
import numpy as np
import obonet
import pandas as pd

from .common import RAW, chrom_norm

GTF = RAW / "genes" / "gencode.v48lift37.basic.annotation.gtf.gz"
HGNC = RAW / "genes" / "hgnc_complete_set.txt"
HPO_OBO = RAW / "hpo" / "hp.obo"
HPOA = RAW / "hpo" / "phenotype.hpoa"
G2D = RAW / "hpo" / "genes_to_disease.txt"
ORPHA6 = RAW / "orphanet" / "en_product6.xml"
CLINGEN = RAW / "clingen" / "clingen_gene_validity.csv"
GNOMAD2 = RAW / "gnomad" / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

ORPHA_CAUSAL = {
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in",
}

# ---------------------------------------------------------------- genes


@lru_cache(maxsize=1)
def gencode_genes() -> pd.DataFrame:
    """Gene records (GRCh37 coordinates) from GENCODE v48 lift37 basic."""
    rows = []
    attr_re = re.compile(r'(\S+) "([^"]*)"')
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            a = dict(attr_re.findall(f[8]))
            rows.append((a["gene_id"].split(".")[0], a.get("gene_name"), a.get("gene_type"),
                         chrom_norm(f[0]), int(f[3]), int(f[4]), f[6]))
    df = pd.DataFrame(rows, columns=["gene_id", "symbol", "biotype", "chrom", "start", "end", "strand"])
    # GENCODE lift37 can contain PAR_Y duplicates / multiple mappings: keep the longest record per ID
    df["len"] = df.end - df.start
    df = (df.sort_values("len", ascending=False).drop_duplicates("gene_id")
            .drop(columns="len").sort_values(["chrom", "start"]).reset_index(drop=True))
    return df


def protein_coding_genes(include_x: bool = True) -> pd.DataFrame:
    g = gencode_genes()
    keep = [str(i) for i in range(1, 23)] + (["X"] if include_x else [])
    return g[(g.biotype == "protein_coding") & g.chrom.isin(keep)].reset_index(drop=True)


@lru_cache(maxsize=1)
def hgnc() -> pd.DataFrame:
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    h = h[h.status == "Approved"].copy()
    h["entrez_id"] = h["entrez_id"].fillna("")
    return h[["hgnc_id", "symbol", "alias_symbol", "prev_symbol", "entrez_id",
              "ensembl_gene_id", "locus_group", "locus_type"]]


@lru_cache(maxsize=1)
def id_maps() -> dict[str, dict[str, str]]:
    """Maps from HGNC id / Entrez id / symbol (incl. previous & alias) to Ensembl gene id."""
    h = hgnc().dropna(subset=["ensembl_gene_id"])
    by_hgnc = dict(zip(h.hgnc_id, h.ensembl_gene_id))
    by_entrez = {e: g for e, g in zip(h.entrez_id, h.ensembl_gene_id) if e}
    by_symbol = dict(zip(h.symbol, h.ensembl_gene_id))
    # previous symbols then aliases, never overriding an approved symbol
    for col in ("prev_symbol", "alias_symbol"):
        for syms, g in zip(h[col].fillna(""), h.ensembl_gene_id):
            for s in syms.split("|"):
                if s and s not in by_symbol:
                    by_symbol[s] = g
    gc = gencode_genes()
    for gid, sym in zip(gc.gene_id, gc.symbol):  # fall back on GENCODE names
        by_symbol.setdefault(sym, gid)
    return {"hgnc": by_hgnc, "entrez": by_entrez, "symbol": by_symbol}


def symbol_of(gene_ids) -> dict[str, str]:
    gc = gencode_genes().set_index("gene_id").symbol.to_dict()
    h = dict(zip(hgnc().ensembl_gene_id, hgnc().symbol))
    return {g: h.get(g, gc.get(g, g)) for g in gene_ids}


# ---------------------------------------------------------------- HPO


@lru_cache(maxsize=1)
def hpo_graph() -> nx.MultiDiGraph:
    """HPO as a graph with edges child -> parent (obonet convention)."""
    return obonet.read_obo(HPO_OBO)


def hpo_descendants(term: str) -> set[str]:
    g = hpo_graph()
    return nx.ancestors(g, term) | {term}  # edges point child->parent


def hpo_name(term: str) -> str:
    return hpo_graph().nodes[term].get("name", term)


def hpo_version() -> str:
    with open(HPO_OBO) as fh:
        for line in fh:
            if line.startswith("data-version"):
                return line.split(":", 1)[1].strip()
    return "unknown"


@lru_cache(maxsize=1)
def hpo_annotations() -> pd.DataFrame:
    """Positive phenotype annotations (aspect P) for OMIM/ORPHA diseases."""
    a = pd.read_csv(HPOA, sep="\t", comment="#", dtype=str)
    a = a[(a.aspect == "P") & (a.qualifier.fillna("") != "NOT")]
    return a[["database_id", "disease_name", "hpo_id", "frequency"]].drop_duplicates()


@lru_cache(maxsize=1)
def hpo_inheritance() -> pd.DataFrame:
    a = pd.read_csv(HPOA, sep="\t", comment="#", dtype=str)
    a = a[(a.aspect == "I") & (a.qualifier.fillna("") != "NOT")]
    return a[["database_id", "hpo_id"]].drop_duplicates()


# ---------------------------------------------------------------- gene-disease links


@lru_cache(maxsize=1)
def omim_gene_disease() -> pd.DataFrame:
    """OMIM gene-disease pairs (HPO genes_to_disease, MENDELIAN, from mim2gene_medgen)."""
    g = pd.read_csv(G2D, sep="\t", dtype=str)
    g = g[g.disease_id.str.startswith("OMIM:")].copy()
    m = id_maps()
    g["entrez"] = g.ncbi_gene_id.str.replace("NCBIGene:", "", regex=False)
    g["gene_id"] = g.entrez.map(m["entrez"]).fillna(g.gene_symbol.map(m["symbol"]))
    g = g.rename(columns={"association_type": "assoc"})
    return g[["gene_id", "gene_symbol", "disease_id", "assoc"]].dropna(subset=["gene_id"]).drop_duplicates()


@lru_cache(maxsize=1)
def orphanet_gene_disease() -> pd.DataFrame:
    """Orphanet disorder-gene associations with type and validation status."""
    m = id_maps()
    rows = []
    for _, el in ET.iterparse(ORPHA6, events=("end",)):
        if el.tag != "Disorder":
            continue
        code = el.findtext("OrphaCode")
        name = el.findtext("Name")
        dtype = el.findtext("DisorderType/Name")
        for assoc in el.iter("DisorderGeneAssociation"):
            gene = assoc.find("Gene")
            sym = gene.findtext("Symbol")
            refs = {r.findtext("Source"): r.findtext("Reference")
                    for r in gene.iter("ExternalReference")}
            gid = None
            if refs.get("HGNC"):
                gid = m["hgnc"].get("HGNC:" + refs["HGNC"])
            if gid is None and refs.get("Ensembl"):
                gid = refs["Ensembl"]
            if gid is None:
                gid = m["symbol"].get(sym)
            rows.append((gid, sym, f"ORPHA:{code}", name, dtype,
                         assoc.findtext("DisorderGeneAssociationType/Name"),
                         assoc.findtext("DisorderGeneAssociationStatus/Name")))
        el.clear()
    df = pd.DataFrame(rows, columns=["gene_id", "gene_symbol", "disease_id", "disease_name",
                                     "disorder_type", "assoc_type", "status"])
    return df.dropna(subset=["gene_id"]).drop_duplicates()


@lru_cache(maxsize=1)
def clingen_validity() -> pd.DataFrame:
    raw = pd.read_csv(CLINGEN, skiprows=[0, 1, 2, 3, 5], dtype=str)
    raw.columns = [c.strip() for c in raw.columns]
    m = id_maps()
    raw["gene_id"] = raw["GENE ID (HGNC)"].map(m["hgnc"]).fillna(raw["GENE SYMBOL"].map(m["symbol"]))
    return raw.rename(columns={"GENE SYMBOL": "gene_symbol", "DISEASE LABEL": "disease_name",
                               "DISEASE ID (MONDO)": "mondo_id", "MOI": "moi",
                               "CLASSIFICATION": "classification"})[
        ["gene_id", "gene_symbol", "disease_name", "mondo_id", "moi", "classification"]].dropna(subset=["gene_id"])


@lru_cache(maxsize=1)
def gnomad_constraint() -> pd.DataFrame:
    """Per-gene LoF constraint from gnomAD v2.1.1 (Karczewski et al. 2020). Used only as a descriptive
    annotation of candidate genes (catalogue, Fig. 3); no test is adjusted for constraint."""
    c = pd.read_csv(GNOMAD2, sep="\t", compression="gzip",
                    usecols=["gene", "gene_id", "oe_lof_upper", "pLI", "mis_z", "obs_lof", "exp_lof"])
    c = c.dropna(subset=["oe_lof_upper"]).sort_values("oe_lof_upper").drop_duplicates("gene_id")
    return c.rename(columns={"oe_lof_upper": "LOEUF", "mis_z": "missense_z",
                             "obs_lof": "lof.obs", "exp_lof": "lof.exp"})
