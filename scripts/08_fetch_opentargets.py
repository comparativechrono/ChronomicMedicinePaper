#!/usr/bin/env python3
"""Retrieve fine-mapped chronotype credible sets, locus-to-gene (L2G) scores and
colocalisations from the Open Targets Platform GraphQL API.

Studies: GCST007576 (Jones 2019, UKB chronotype, PICS), GCST007565 (Jones 2019, UKB morning
person, PICS), GCST90691736 (Karczewski 2025 Pan-UKB EUR chronotype, SuSiE-inf).

Outputs (data/raw/opentargets/): meta.json, credible_sets.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import RAW, get_logger, save_json  # noqa: E402

log = get_logger("08_fetch_opentargets")
URL = "https://api.platform.opentargets.org/api/v4/graphql"
OUT = RAW / "opentargets"
OUT.mkdir(parents=True, exist_ok=True)
STUDIES = ["GCST007576", "GCST007565", "GCST90691736"]


def q(query: str, variables: dict, tries: int = 5) -> dict:
    for k in range(tries):
        try:
            r = requests.post(URL, json={"query": query, "variables": variables}, timeout=300)
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                raise RuntimeError(d["errors"])
            return d["data"]
        except Exception as exc:  # transient API errors: back off and retry
            log.warning("retry %d: %s", k + 1, exc)
            time.sleep(5 * (k + 1))
    raise RuntimeError("Open Targets query failed")


meta = q("{ meta { apiVersion { x y z } dataVersion { year month iteration } } }", {})
meta["retrieved_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
save_json(meta, OUT / "meta.json")

STUDY_Q = """query($id:String!){ study(studyId:$id){ id traitFromSource nSamples publicationFirstAuthor
  publicationDate pubmedId credibleSets(page:{index:0,size:1}){ count } } }"""
CS_Q = """query($id:String!,$i:Int!,$n:Int!){ study(studyId:$id){ credibleSets(page:{index:$i,size:$n}){ rows {
  studyLocusId studyId chromosome position region finemappingMethod confidence credibleSetIndex
  pValueMantissa pValueExponent beta
  variant { id rsIds mostSevereConsequence { label } }
  l2GPredictions(page:{index:0,size:20}){ rows { score target { id approvedSymbol } } }
} } } }"""
LOCUS_Q = """query($sl:String!){ credibleSet(studyLocusId:$sl){
  locus(page:{index:0,size:300}){ count rows { posteriorProbability is95CredibleSet r2Overall variant { id rsIds } } } } }"""
COLOC_Q = """query($sl:String!,$i:Int!){ credibleSet(studyLocusId:$sl){
  colocalisation(page:{index:$i,size:200}){ count rows { h4 clpp colocalisationMethod rightStudyType
       betaRatioSignAverage numberColocalisingVariants
       otherStudyLocus { studyLocusId studyId qtlGeneId
          study { traitFromSource studyType projectId biosample { biosampleName } target { id approvedSymbol } } } } } } }"""

records, studies = [], {}
for sid in STUDIES:
    s = q(STUDY_Q, {"id": sid})["study"]
    studies[sid] = s
    n = s["credibleSets"]["count"]
    log.info("%s %s n=%s: %d credible sets", sid, s["traitFromSource"], s["nSamples"], n)
    size = 10
    for i in range(0, (n + size - 1) // size):
        rows = q(CS_Q, {"id": sid, "i": i, "n": size})["study"]["credibleSets"]["rows"]
        records.extend(rows)
        time.sleep(0.1)
for k, rec in enumerate(records):
    sl = rec["studyLocusId"]
    rec["locus"] = q(LOCUS_Q, {"sl": sl})["credibleSet"]["locus"]
    co = q(COLOC_Q, {"sl": sl, "i": 0})["credibleSet"]["colocalisation"]
    rows = co["rows"]
    for i in range(1, (co["count"] + 199) // 200):
        rows += q(COLOC_Q, {"sl": sl, "i": i})["credibleSet"]["colocalisation"]["rows"]
    rec["colocalisation"] = {"count": co["count"], "rows": rows}
    if k % 50 == 0:
        log.info("details fetched for %d/%d credible sets", k, len(records))
save_json({"meta": meta, "studies": studies, "credible_sets": records}, OUT / "credible_sets.json")
from chronomic.common import record_in_manifest  # noqa: E402
for f in ("meta.json", "credible_sets.json"):
    record_in_manifest(OUT / f, "https://api.platform.opentargets.org/api/v4/graphql (Open Targets Platform 26.09; queries in scripts/08_fetch_opentargets.py and 08b_fetch_opentargets_qtl.py)")
log.info("saved %d credible sets", len(records))
