# Towards chronomic medicine

**Mendelian disease genes and rare-variant convergence at chronotype loci.**
Malik BT & Hearn TJ. Code, results and manuscript.

*Chronomic medicine* is the integration of chronobiology with genomic medicine. It uses genetic variation in
timing pathways, together with the time-of-day dependence of gene expression, phenotypes and treatments,
to understand, diagnose and manage disease.

## Main findings

* **No general enrichment.** Mendelian disease genes are not enriched at the 351 chronotype loci of Jones et al.
  (2019) once null loci are matched for allele frequency, LD, gene count and distance. Genes for dominant and
  neurodevelopmental disorders, and for five HPO organ systems, are modestly over-represented within loci, with
  little gene-level (MAGMA) support.
* **Clinical disease panels.** Of 248 Genomics England PanelApp panels, six are enriched at chronotype loci
  (FDR < 0.05): the neurotransmitter disorders panel, driven by four monoamine genes (*GCH1*, *DDC*, *DNAJC12*,
  *SLC18A2*), and five broad developmental panels with little gene-level support. No cardiac, immune or
  sleep-related panel is enriched.
* **Phenotype-specific convergence is modest.**
  * Genes for Mendelian disorders with sleep disturbance are over-represented, but specificity for sleep is
    uncertain. Most of these genes also cause neurodevelopmental disorders, genes for other syndromic features
    (such as craniosynostosis and clefting) show similar excesses, and gene-level support is nominal.
  * In a post hoc analysis, 4 of 12 monoamine neurotransmitter disorder genes (*GCH1*, *DDC*, *DNAJC12*,
    *SLC18A2*) lie in chronotype loci; gene-level support is strongest for *DDC*. The PanelApp result rests on
    the same four genes, so it supports rather than replicates this finding.
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

## Manuscript and display items

| Item | File |
|---|---|
| Manuscript (Word) | `manuscript/Malik_Hearn_Towards_chronomic_medicine_2026.docx`, built from `manuscript/manuscript.template.md` |
| Figure 1. Study design and the chronomic medicine framework | `results/figures/Fig1_study_design.pdf` / `.png` |
| Figure 2. Enrichment of gene sets at chronotype loci | `results/figures/Fig2_locus_enrichment.pdf` / `.png` |
| Figure 3. Evidence for chronomic candidate genes | `results/figures/Fig3_candidate_genes.pdf` / `.png` |
| Figure 4. Monoamine synthesis and storage genes (post hoc) | `results/figures/Fig4_monoamine_pathway.pdf` / `.png` |
| Figure 5. Convergence of rare-variant burden with chronotype loci | `results/figures/Fig5_rare_common_convergence.pdf` / `.png` |
| Supplementary Figure S1. Chronotype and QT-interval association at *KCNH2* | `results/figures/FigS1_kcnh2_coloc.pdf` / `.png` |
| Table 1. Locus- and gene-level enrichment of selected gene sets | in the manuscript; source `manuscript/tables/table1.md` |
| Table 2. PanelApp disease panels enriched at chronotype loci | in the manuscript; source `manuscript/tables/table2.md` |
| Table 3. Selected chronomic candidate genes | in the manuscript; source `manuscript/tables/table3.md` |
| Supplementary Table S1. Chronotype loci, LD intervals, genes and Open Targets annotations | `results/supplementary/Supplementary_Table_S1_chronotype_loci.xlsx` |
| Supplementary Table S2. Gene universe and gene sets, including PanelApp panels | `results/supplementary/Supplementary_Table_S2_gene_sets.xlsx` |
| Supplementary Table S3. Locus- and gene-level enrichment, relative tests and PanelApp screen | `results/supplementary/Supplementary_Table_S3_enrichment.xlsx` |
| Supplementary Table S4. Exome burden analyses | `results/supplementary/Supplementary_Table_S4_rare_variant_burden.xlsx` |
| Supplementary Table S5. Catalogue of Mendelian disease genes at chronotype loci | `results/supplementary/Supplementary_Table_S5_chronomic_gene_catalogue.xlsx` |
| Supplementary Table S6. *KCNH2* colocalisation and data provenance | `results/supplementary/Supplementary_Table_S6_colocalisation_and_provenance.xlsx` |

Each workbook opens with a "Key" sheet that describes its sheets and columns. The Word build uses copies of the
figures in `manuscript/figures/`.

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
| gnomAD constraint | v2.1.1 | per-gene LoF constraint (LOEUF), descriptive annotation of candidate genes only |
| Gene Ontology | downloaded 2026-10-30 | positive-control circadian gene set |
| Open Targets Platform | 26.09 (GraphQL API) | credible sets, L2G, colocalisation |
| UK Biobank exome burden (Backman et al. 2021) | GCST90081559–62 (European ancestry, n = 388,330) | rare-variant convergence |
| QT interval GWAS (Hoffmann et al. 2022) | GCST90165291 | colocalisation at *KCNH2* |
| Genomics England PanelApp | REST API v1, green genes (retrieved 2026-10-30) | clinical disease-panel screen |
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
