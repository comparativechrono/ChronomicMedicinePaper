# Submission notes (for the authors)

## To do before submission

- [ ] Update affiliations, corresponding e-mail and ORCID iDs (title page, marked in brackets).
- [ ] Confirm the CRediT statement, funding and competing interests (marked in brackets).
- [ ] Read the "Use of generative AI" statement and adapt it to the target journal's policy.
- [ ] Merge the GitHub pull request, create a release and archive it on Zenodo (GitHub–Zenodo integration), then
      add the DOI to "Code availability" and to `CITATION.cff`.
- [ ] Check the abstract and reference style against the chosen journal. The manuscript uses AMA numbered style;
      to change it, replace `manuscript/ama.csl` and run `make manuscript`.
- [ ] Optional replication:
  - Test the post hoc monoamine finding in an independent chronotype GWAS when one becomes available, or with
    objective sleep-timing GWAS.
  - Confirm the *KCNH2* conditional colocalisation with in-sample LD (for example, COJO or SuSiE-based
    colocalisation in UK Biobank).

## Journal shortlist (check current limits and fees)

1. **npj Biological Timing and Sleep** (Nature Portfolio, OA). Best conceptual fit. Abstract ≤ 150 words, no
   main-text limit; the abstract would need shortening.
2. **HGG Advances** (ASHG/Cell Press, OA). Best fit for the genetics and methods. Check limits.
3. **Journal of Biological Rhythms** (SAGE/SRBR). No publication fees. Structured abstract ≤ 300 words; no stated
   length limits.
4. **European Journal of Human Genetics**: rare plus common genetics.
5. **SLEEP Advances / SLEEP**: sleep-genetics audience.
6. **Orphanet Journal of Rare Diseases**: if the rare-disease angle is emphasised.
7. **Communications Biology**: ~5,000 words, ≤ 10 display items; the main text would need shortening.
8. **Genome Medicine** (stretch).

The structured abstract is about 300 words. The main text (Introduction to Discussion) is about 7,000 words, with
5 figures, 3 tables, 1 supplementary figure and 6 supplementary tables. Journals with a ~5,000-word limit would need
the Methods moved partly to the Supplement.

## Rebuilding the Word file

```bash
make numbers manuscript      # re-renders every statistic from results/key_numbers.json
```

Edit `manuscript/manuscript.template.md`, not `manuscript.md`. Statistics are written as `{{path@format}}`
placeholders, and the build stops if any placeholder cannot be resolved.
