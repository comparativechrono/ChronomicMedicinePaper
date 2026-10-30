# Towards chronomic medicine

**Mendelian disease genes and rare-variant convergence at chronotype loci.**
Malik BT & Hearn TJ. Code, results and manuscript.

*Chronomic medicine* is the integration of chronobiology with genomic medicine. It uses genetic variation in
timing pathways, together with the time-of-day dependence of gene expression, phenotypes and treatments,
to understand, diagnose and manage disease.

## Main findings

* **No general enrichment.** Mendelian disease genes are not enriched at the 351 chronotype loci of Jones et al.
  (2019) once null loci are matched for allele frequency, LD, gene count and distance. The loci are, however,
  enriched for loss-of-function-constrained genes, and constraint accounts for most apparent organ-system and
  neurodevelopmental enrichments.
* **Phenotype-specific convergence is modest.**
  * Genes for Mendelian disorders with sleep disturbance are over-represented, but the signal is attenuated by
    constraint matching, is carried by genes that also cause neurodevelopmental disorders, and is not separable
    from an excess of autosomal dominant disease genes.
  * In a post hoc analysis, 4 of 12 monoamine neurotransmitter disorder genes (*GCH1*, *DDC*, *DNAJC12*,
    *SLC18A2*) lie in chronotype loci; gene-level support is strongest for *DDC*.
* **Rare variants converge on clock genes.** In UK Biobank exomes, rare variants in *PER2* and *PER3* are
  associated with morningness (exome-wide significant), and deleterious *FBXL3* variants nominally with
  eveningness. These three genes account for the excess burden signal among Mendelian genes at chronotype loci.
* **A candidate shared signal at *KCNH2* (post hoc).** After conditioning on the dominant QT signal, a
  sub-genome-wide chronotype signal colocalises with a secondary QT-interval signal (PP.H4 0.87; 0.87–0.97 across
  36 sensitivity analyses). The analysis is approximate, using an external LD reference, and needs replication.
* **Catalogue.** `results/supplementary/Supplementary_Table_S5_chronomic_gene_catalogue.xlsx` catalogues 236
  Mendelian disease genes at 123 chronotype loci, with fine-mapping, colocalisation, burden, constraint and
  circadian expression evidence.

Every number quoted in the manuscript is rendered from `results/key_numbers.json` at build time
(`manuscript/build.py`), so the text cannot drift from the analysis.

## Repository layout

```
config/config.yaml          all analysis parameters (thresholds, windows, seeds)
scripts/00-16_*.py          numbered pipeline stages (see Makefile for order)
src/chronomic/              shared code: annotation readers, LD/locus mapping, matched-null enrichment, coloc
tests/                      unit tests for the statistical helpers
data/manifest.tsv           every input: URL or API endpoint, retrieval time (UTC), size, MD5
results/tables/             analysis outputs (TSV/JSON)
results/figures/            figures (PDF + 300 dpi PNG)
results/supplementary/      Supplementary Tables S1-S6 (xlsx)
results/key_numbers.json    every statistic quoted in the manuscript
manuscript/                 manuscript template, generated tables, bibliography, Word output
```

`data/raw/`, `data/interim/` and `bin/` are created by the pipeline and are not tracked.

## Reproducing the analysis

Requirements: Linux, Python 3.13 (`pip install -r requirements.txt`), pandoc ≥ 3.1, ~40 GB of free disk
space and internet access. PLINK 1.9/2.0 and MAGMA 1.10 are downloaded automatically. MAGMA's licence does
not allow redistribution, so it is fetched from its authors' site.

```bash
make all        # download → gene sets → loci → null pool → rhythmicity → enrichment → MAGMA → Open Targets
                # → exome burden → coloc → catalogue → figures → tables → key numbers → manuscript → tests
```

Individual stages can be re-run with `make <stage>` (see `Makefile`). The full run takes about 1.5 hours on two
cores. The longest steps are the downloads, MAGMA, and the Open Targets API queries (Platform release 26.09).

## Data sources

| Resource | Version / accession | Use |
|---|---|---|
| Chronotype GWAS (Jones et al. 2019, *Nat Commun* 10:343) | Supplementary Data 1; GWAS Catalog GCST007576, GCST007565 | lead variants; UKB summary statistics |
| 1000 Genomes phase 3 EUR (MAGMA build) | 503 individuals, GRCh37 | LD, MAGMA reference |
| GENCODE / HGNC | v48 (GRCh37 lift) / Oct 2026 | genes and identifiers |
| HPO | 2026-09-01 | OMIM gene–disease links, phenotype annotations |
| Orphanet | product 6, June 2026 | disease-causing gene associations |
| ClinGen gene–disease validity | downloaded 2026-10-30 | validated Mendelian genes |
| IUIS classification of inborn errors of immunity (Tangye et al. 2022) | complement deficiencies | curated complement gene set |
| gnomAD constraint | v2.1.1 (primary), v4.1 (sensitivity) | LoF constraint (LOEUF) |
| Gene Ontology | downloaded 2026-10-30 | positive-control circadian gene set |
| Open Targets Platform | 26.09 (GraphQL API) | credible sets, L2G, colocalisation |
| UK Biobank exome burden (Backman et al. 2021) | GCST90081559–62 (European ancestry, n = 388,330) | rare-variant convergence |
| QT interval GWAS (Hoffmann et al. 2022) | GCST90165291 | colocalisation at *KCNH2* |
| Mouse circadian atlas (Zhang et al. 2014) | GEO GSE54650 | rhythmic expression |
| GTEx day–night expression (Wucher et al. 2023) | figshare 10.6084/m9.figshare.21906252 (CC BY 4.0) | day–night expression |

Each resource remains subject to its own terms of use.

## History of this repository

The repository previously held an unrelated earlier draft (Saoulidis, Seed, Malik & Hearn). That content
is preserved in the history of `main` (commit `7e8411f`).

## Licence and citation

Code: MIT (see `LICENSE`). Manuscript text and figures: © the authors pending publication.
Please cite the article when published (see `CITATION.cff`); until then cite this repository.

The analysis code and manuscript were developed with the assistance of an AI system (Claude, Anthropic).
The authors are responsible for reviewing and verifying all analyses and text.
