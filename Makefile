# Reproduce every result, figure, table and the manuscript from public data.
# Usage: make all        (full pipeline; ~1.5 h on 2 cores, ~5 GB of downloads)
#        make <target>   (single stage; stages assume earlier ones have run)
PY ?= python3

.PHONY: all download genesets loci null rhythm enrichment magma opentargets burden coloc \
        catalogue figures supptables numbers manuscript test clean

all: download genesets loci null rhythm enrichment magma opentargets burden coloc \
     catalogue figures supptables numbers manuscript test

download:    ; $(PY) scripts/00_download_data.py
genesets:    ; $(PY) scripts/01_build_gene_sets.py
loci:        ; $(PY) scripts/02_define_loci.py
null:        ; $(PY) scripts/03_null_pool.py
rhythm:      ; $(PY) scripts/04_rhythmicity.py          # appends rhythmicity gene sets; run after genesets
enrichment:  ; $(PY) scripts/05_enrichment.py
magma:
	$(PY) scripts/06_magma_genes.py chronotype
	$(PY) scripts/07_magma_genesets.py
opentargets:
	$(PY) scripts/08_fetch_opentargets.py
	$(PY) scripts/08b_fetch_opentargets_qtl.py
	$(PY) scripts/09_opentargets_loci.py
burden:      ; $(PY) scripts/10_rare_variant_burden.py
coloc:       ; $(PY) scripts/11_coloc_kcnh2.py
catalogue:   ; $(PY) scripts/12_catalogue.py
figures:     ; $(PY) scripts/13_figures.py && mkdir -p manuscript/figures && cp results/figures/* manuscript/figures/
supptables:  ; $(PY) scripts/14_supplementary_tables.py
numbers:     ; $(PY) scripts/15_key_numbers.py
manuscript:
	$(PY) scripts/16_manuscript_tables.py
	cd manuscript && $(PY) make_reference_docx.py && $(PY) build.py
test:        ; $(PY) -m pytest -q tests/

clean:
	rm -rf data/interim results/tables/* results/figures/* results/supplementary/* results/key_numbers.json
