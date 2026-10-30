#!/usr/bin/env python3
"""Add molecular-QTL colocalisations (eQTL, pQTL, sQTL, tuQTL and single-cell QTLs) to the
Open Targets credible sets fetched by 08_fetch_opentargets.py (the API returns GWAS-GWAS
colocalisations by default)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import RAW, get_logger, save_json  # noqa: E402

log = get_logger("08b_fetch_opentargets_qtl")
URL = "https://api.platform.opentargets.org/api/v4/graphql"
QTL_TYPES = ["eqtl", "pqtl", "sqtl", "tuqtl", "sceqtl", "scpqtl", "scsqtl", "sctuqtl"]
Q = """query($sl:String!,$i:Int!,$t:[StudyTypeEnum!]){ credibleSet(studyLocusId:$sl){
  colocalisation(studyTypes:$t, page:{index:$i,size:200}){ count rows { h4 clpp colocalisationMethod rightStudyType
       betaRatioSignAverage numberColocalisingVariants
       otherStudyLocus { studyLocusId studyId qtlGeneId
          study { traitFromSource studyType projectId biosample { biosampleName } target { id approvedSymbol } } } } } } }"""


def q(variables, tries=5):
    for k in range(tries):
        try:
            r = requests.post(URL, json={"query": Q, "variables": variables}, timeout=300)
            r.raise_for_status()
            d = r.json()
            if "errors" in d:
                raise RuntimeError(d["errors"])
            return d["data"]["credibleSet"]["colocalisation"]
        except Exception as exc:
            log.warning("retry %d: %s", k + 1, exc)
            time.sleep(5 * (k + 1))
    raise RuntimeError("query failed")


path = RAW / "opentargets" / "credible_sets.json"
D = json.load(open(path))
for k, cs in enumerate(D["credible_sets"]):
    co = q({"sl": cs["studyLocusId"], "i": 0, "t": QTL_TYPES})
    rows = co["rows"]
    for i in range(1, (co["count"] + 199) // 200):
        rows += q({"sl": cs["studyLocusId"], "i": i, "t": QTL_TYPES})["rows"]
    cs["colocalisation"]["rows"] = [r for r in cs["colocalisation"]["rows"] if r["rightStudyType"] == "gwas"] + rows
    cs["colocalisation"]["count_qtl"] = co["count"]
    if k % 100 == 0:
        log.info("%d/%d", k, len(D["credible_sets"]))
save_json(D, path)
from chronomic.common import record_in_manifest  # noqa: E402
record_in_manifest(path, "https://api.platform.opentargets.org/api/v4/graphql (Open Targets Platform 26.09; queries in scripts/08_fetch_opentargets.py and 08b_fetch_opentargets_qtl.py)")
log.info("done")
