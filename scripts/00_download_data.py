#!/usr/bin/env python3
"""Download every public input used in the analysis and record provenance.

Each file is written to data/raw/<subdir>/<name>. A manifest (data/manifest.tsv)
records the source URL, the download date (UTC), the file size and its MD5
checksum so that the exact inputs can be checked or re-fetched later.

Usage:
    python scripts/00_download_data.py            # fetch anything missing
    python scripts/00_download_data.py --only hpo # fetch one group
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import shutil
import sys
import zipfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "data" / "manifest.tsv"

GWAS_FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000"
JONES_ESM = ("https://media.springernature.com/original/springer-static/esm/"
             "art%3A10.1038%2Fs41467-018-08259-7/MediaObjects/41467_2018_8259_MOESM{n}_ESM.{ext}")
HPO_REL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
SURF = "https://vu.data.surfsara.nl/index.php/s/{token}/download"

# group, sub-directory, local file name, URL
FILES = [
    # --- software (not redistributed; MAGMA licence forbids redistribution) ---
    ("tools", "tools", "plink_linux_x86_64_20231211.zip",
     "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"),
    ("tools", "tools", "plink2_linux_x86_64_20250129.zip",
     "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip"),
    ("tools", "tools", "magma_v1.10_static.zip", SURF.format(token="lxDgt2dNdNr6DYt")),
    # --- LD reference: 1000 Genomes phase 3 European panel (GRCh37), MAGMA build ---
    ("ldref", "ldref", "g1000_eur.zip", SURF.format(token="VZNByNwpD8qqINe")),
    # --- chronotype GWAS (Jones et al. 2019, Nat Commun 10:343) ---
    ("gwas", "gwas", "GCST007576_chronotype_UKB.txt.gz",
     f"{GWAS_FTP}/GCST007576/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz"),
    ("gwas", "gwas", "GCST007576_README.txt", f"{GWAS_FTP}/GCST007576/chronotype_README.txt"),
    ("gwas", "gwas", "GCST007565_morning_person_UKB.txt.gz",
     f"{GWAS_FTP}/GCST007565/morning_person_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3_logORs.txt.gz"),
    ("gwas", "gwas", "GCST007565_README.txt", f"{GWAS_FTP}/GCST007565/morning_person_README.txt"),
] + [
    ("jones", "gwas/jones2019_esm", f"jones2019_SupplementaryData{n - 3}.{ext}",
     JONES_ESM.format(n=n, ext=ext))
    for n, ext in [(4, "xlsx"), (5, "xlsx"), (6, "xlsx"), (7, "xlsx"), (8, "xlsx"), (9, "xlsx"),
                   (10, "xlsx"), (11, "xlsx"), (12, "xlsx"), (13, "xlsx"), (14, "xlsx"),
                   (15, "xlsx"), (16, "xlsx"), (17, "txt"), (18, "xlsx")]
] + [
    ("jones", "gwas/jones2019_esm", "jones2019_SupplementaryInformation.pdf",
     JONES_ESM.format(n=1, ext="pdf")),
    ("jones", "gwas/jones2019_esm", "jones2019_DescriptionOfAdditionalSupplementaryFiles.pdf",
     JONES_ESM.format(n=3, ext="pdf")),
    # --- gene models: GENCODE v48 (current gene set) lifted to GRCh37 ---
    ("genes", "genes", "gencode.v48lift37.basic.annotation.gtf.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh37_mapping/"
     "gencode.v48lift37.basic.annotation.gtf.gz"),
    ("genes", "genes", "hgnc_complete_set.txt",
     "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"),
    # --- Mendelian disease knowledge ---
    ("hpo", "hpo", "hp.obo", f"{HPO_REL}/hp.obo"),
    ("hpo", "hpo", "phenotype.hpoa", f"{HPO_REL}/phenotype.hpoa"),
    ("hpo", "hpo", "genes_to_disease.txt", f"{HPO_REL}/genes_to_disease.txt"),
    ("hpo", "hpo", "genes_to_phenotype.txt", f"{HPO_REL}/genes_to_phenotype.txt"),
    ("orphanet", "orphanet", "en_product6.xml", "https://www.orphadata.com/data/xml/en_product6.xml"),
    ("clingen", "clingen", "clingen_gene_validity.csv",
     "https://search.clinicalgenome.org/kb/gene-validity/download"),
    # --- Gene Ontology (positive-control 'circadian rhythm' gene set) ---
    ("go", "go", "go-basic.obo", "https://current.geneontology.org/ontology/go-basic.obo"),
    ("go", "go", "goa_human.gaf.gz", "https://current.geneontology.org/annotations/goa_human.gaf.gz"),
    # --- rare-variant gene burden results for chronotype (Backman et al. 2021; European-ancestry subset of the
    #     UKB 454,787 exomes; sample size read from the GWAS Catalog metadata files) ---
    ("burden", "burden", "GCST90081559_buildGRCh38.tsv.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081559/GCST90081559_buildGRCh38.tsv.gz"),
    ("burden", "burden", "GCST90081560_buildGRCh38.tsv.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081560/GCST90081560_buildGRCh38.tsv.gz"),
    ("burden", "burden", "GCST90081561_buildGRCh38.tsv.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081561/GCST90081561_buildGRCh38.tsv.gz"),
    ("burden", "burden", "GCST90081562_buildGRCh38.tsv.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081562/GCST90081562_buildGRCh38.tsv.gz"),
    ("burden", "burden", "GCST90081559_buildGRCh38.tsv.gz-meta.yaml",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081559/GCST90081559_buildGRCh38.tsv.gz-meta.yaml"),
    ("burden", "burden", "GCST90081560_buildGRCh38.tsv.gz-meta.yaml",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081560/GCST90081560_buildGRCh38.tsv.gz-meta.yaml"),
    ("burden", "burden", "GCST90081561_buildGRCh38.tsv.gz-meta.yaml",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081561/GCST90081561_buildGRCh38.tsv.gz-meta.yaml"),
    ("burden", "burden", "GCST90081562_buildGRCh38.tsv.gz-meta.yaml",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90081001-GCST90082000/GCST90081562/GCST90081562_buildGRCh38.tsv.gz-meta.yaml"),
    # --- day-night gene expression across GTEx tissues (Wucher et al. 2023 PLoS Biol; CC BY 4.0) ---
    ("rhythm", "rhythm", "wucher2023_S1_daynight.tab.gz", "https://ndownloader.figshare.com/files/38858706"),
    # --- mouse 12-tissue circadian atlas (Zhang et al. 2014 PNAS; GEO GSE54650) + annotation + orthologs ---
    ("rhythm", "rhythm", "GSE54650_series_matrix.txt.gz",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54650/matrix/GSE54650_series_matrix.txt.gz"),
    ("rhythm", "rhythm", "GPL6246.annot.gz", "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6246/annot/GPL6246.annot.gz"),
    ("rhythm", "rhythm", "HOM_MouseHumanSequence.rpt",
     "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"),
    # --- LoF constraint (gnomAD v2.1.1, Karczewski et al. 2020): descriptive annotation of candidate genes ---
    ("gnomad", "gnomad", "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
     "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"),
    # --- QT interval GWAS for colocalisation at KCNH2 (Hoffmann et al. 2022 Genetics; GERA+UKB+PAGE) ---
    ("qt", "qt", "GCST90165291_QT_buildGRCh37.tsv.gz",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165291/"
     "GCST90165291_buildGRCh37.tsv.gz"),
]


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    headers = {"User-Agent": "Mozilla/5.0 (chronomic-medicine reproducibility pipeline)"}
    with requests.get(url, stream=True, timeout=120, headers=headers) as r:
        r.raise_for_status()
        r.raw.decode_content = True  # undo HTTP-level gzip transfer encoding
        with open(tmp, "wb") as fh:
            shutil.copyfileobj(r.raw, fh, length=1 << 20)
    tmp.rename(dest)


def read_manifest() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    if MANIFEST.exists():
        for line in MANIFEST.read_text().splitlines()[1:]:
            parts = line.split("\t")
            rows[parts[0]] = parts
    return rows


def write_manifest(rows: dict[str, list[str]]) -> None:
    header = ["path", "url", "downloaded_utc", "bytes", "md5"]
    with open(MANIFEST, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for key in sorted(rows):
            fh.write("\t".join(rows[key]) + "\n")


def unpack_tools() -> None:
    tools = RAW / "tools"
    bindir = ROOT / "bin"
    bindir.mkdir(exist_ok=True)
    for zname, members in [("plink_linux_x86_64_20231211.zip", ["plink"]),
                           ("plink2_linux_x86_64_20250129.zip", ["plink2"]),
                           ("magma_v1.10_static.zip", ["magma"])]:
        zpath = tools / zname
        if not zpath.exists():
            continue
        with zipfile.ZipFile(zpath) as zf:
            for info in zf.infolist():
                base = os.path.basename(info.filename)
                if base in members:
                    target = bindir / base
                    if target.exists():
                        continue
                    with zf.open(info) as src, open(target, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                    target.chmod(0o755)


def unpack_ldref() -> None:
    zpath = RAW / "ldref" / "g1000_eur.zip"
    out = RAW / "ldref"
    if zpath.exists() and not (out / "g1000_eur.bed").exists():
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="restrict to these groups")
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    args = ap.parse_args()

    rows = read_manifest()
    failures = []
    for group, sub, name, url in FILES:
        if args.only and group not in args.only:
            continue
        dest = RAW / sub / name
        rel = str(dest.relative_to(ROOT))
        if dest.exists() and not args.force:
            if rel not in rows or rows[rel][2] == "unknown":
                stamp = dt.datetime.fromtimestamp(dest.stat().st_mtime, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                rows[rel] = [rel, url, stamp, str(dest.stat().st_size), md5sum(dest)]
            continue
        print(f"[download] {name}", flush=True)
        try:
            fetch(url, dest)
        except Exception as exc:  # keep going; report at the end
            failures.append((name, url, repr(exc)))
            print(f"  FAILED: {exc}", file=sys.stderr)
            continue
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        rows = read_manifest() | {k: v for k, v in rows.items() if k not in read_manifest()}
        rows[rel] = [rel, url, stamp, str(dest.stat().st_size), md5sum(dest)]
        write_manifest(rows)
    disk = read_manifest()
    disk.update({k: v for k, v in rows.items() if k not in disk or disk[k][2] == "unknown"})
    write_manifest(disk)
    unpack_tools()
    unpack_ldref()
    if failures:
        print("\nFailed downloads:", file=sys.stderr)
        for f in failures:
            print("  ", *f, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
