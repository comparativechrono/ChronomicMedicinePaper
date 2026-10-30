#!/usr/bin/env python3
"""MAGMA v1.10 gene analysis of UK Biobank chronotype summary statistics.

Gene windows: 35 kb upstream / 10 kb downstream of protein-coding genes (GENCODE v48, GRCh37).
LD reference: 1000 Genomes phase 3 European panel (MAGMA build, dbSNP 151 synonyms).
Runs the primary (continuous chronotype, GCST007576) and secondary (morning person,
GCST007565) GWAS in two parallel batches each.

Outputs (data/interim/magma/): <trait>.genes.out / .genes.raw
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import BIN, INTERIM, RAW, ROOT, get_logger, load_config, chrom_norm  # noqa: E402

log = get_logger("06_magma_genes")
cfg = load_config()
OUT = INTERIM / "magma"
OUT.mkdir(parents=True, exist_ok=True)
REF = RAW / "ldref" / "g1000_eur"
MAGMA = BIN / "magma"

genes = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t", dtype={"chrom": str})
loc = OUT / "genes.loc"
genes[["gene_id", "chrom", "start", "end", "strand", "symbol"]].to_csv(loc, sep="\t", header=False, index=False)

up, down = cfg["magma"]["window_up_kb"], cfg["magma"]["window_down_kb"]
annot = OUT / "genes"
if not annot.with_suffix(".genes.annot").exists():
    subprocess.run([MAGMA, "--annotate", f"window={up},{down}", "--snp-loc", f"{REF}.bim",
                    "--gene-loc", loc, "--out", annot], check=True)

traits = {
    "chronotype": (ROOT / cfg["gwas"]["sumstats_primary"], cfg["gwas"]["sumstats_primary_n"]),
    "morning_person": (ROOT / cfg["gwas"]["sumstats_secondary"], cfg["gwas"]["sumstats_secondary_n"]),
}
only = sys.argv[1:] or list(traits)
for trait in only:
    path, n = traits[trait]
    pv = OUT / f"{trait}.pval"
    if not pv.exists():
        with open(pv, "w") as fh:
            fh.write("SNP\tP\tN\n")
            for ch in pd.read_csv(path, sep="\t", usecols=["SNP", "P_BOLT_LMM"], chunksize=2_000_000):
                ch["N"] = n
                ch.to_csv(fh, sep="\t", header=False, index=False)
    prefix = OUT / trait
    if not prefix.with_suffix(".genes.out").exists():
        procs = []
        for b in (1, 2):
            cmd = [MAGMA, "--bfile", REF, "--pval", pv, "ncol=N", "--gene-annot",
                   annot.with_suffix(".genes.annot"), "--gene-model", "snp-wise=mean",
                   "--batch", str(b), "2", "--out", prefix]
            procs.append(subprocess.Popen([str(c) for c in cmd], stdout=subprocess.DEVNULL))
        for p in procs:
            if p.wait() != 0:
                raise SystemExit(f"MAGMA batch failed for {trait}")
        subprocess.run([MAGMA, "--merge", prefix, "--out", prefix], check=True, stdout=subprocess.DEVNULL)
    log.info("%s gene analysis complete", trait)
