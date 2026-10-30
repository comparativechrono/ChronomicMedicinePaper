#!/usr/bin/env python3
"""Retrieve Genomics England PanelApp disease gene panels (green, diagnostic-grade genes).

PanelApp panels are expert-curated gene lists for clinically defined disease groups (for example
"Long QT syndrome" or "Genetic epilepsy syndromes"). We retrieve the list of public panels and every
gene rated green (confidence level 3) on any panel, through the PanelApp REST API (v1).

Outputs (data/raw/panelapp/):
  panels.tsv        one row per public panel: id, name, version, disease group, panel types
  green_genes.tsv   one row per (panel, green gene): symbol, HGNC ID, GRCh37 Ensembl ID, inheritance
  green_genes_raw.json.gz  raw API records (cache)
  meta.json         retrieval time (UTC) and counts
Both files are recorded in data/manifest.tsv. PanelApp changes continuously, so an existing retrieval is
reused; pass --refresh to fetch the current panels (results may then differ from the manuscript).
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import RAW, get_logger, record_in_manifest  # noqa: E402

log = get_logger("01b_fetch_panelapp")
API = "https://panelapp.genomicsengland.co.uk/api/v1"
OUT = RAW / "panelapp"
OUT.mkdir(parents=True, exist_ok=True)


def get(url: str) -> dict:
    for k in range(5):
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # transient API errors: back off and retry
            log.warning("retry %d for %s: %s", k + 1, url, exc)
            time.sleep(5 * (k + 1))
    raise RuntimeError(url)


def paged(url: str) -> list[dict]:
    rows = []
    while url:
        d = get(url)
        rows += d["results"]
        url = d["next"]
        time.sleep(0.1)
    return rows


ap = argparse.ArgumentParser()
ap.add_argument("--refresh", action="store_true", help="fetch the current PanelApp instead of reusing the recorded retrieval")
args = ap.parse_args()
if all((OUT / f).exists() for f in ("panels.tsv", "green_genes.tsv", "meta.json")) and not args.refresh:
    # PanelApp changes continuously: reuse the recorded retrieval so that re-runs reproduce the published results
    log.info("reusing the PanelApp retrieval of %s (pass --refresh to fetch the current panels)",
             json.load(open(OUT / "meta.json"))["retrieved_utc"])
    sys.exit(0)
if args.refresh:
    (OUT / "green_genes_raw.json.gz").unlink(missing_ok=True)

retrieved = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
panels = paged(f"{API}/panels/")
P = pd.DataFrame([{"panel_id": p["id"], "panel_name": p["name"], "version": p["version"],
                   "disease_group": p.get("disease_group"), "disease_sub_group": p.get("disease_sub_group"),
                   "n_genes": p["stats"]["number_of_genes"],
                   "types": ";".join(t["name"] for t in p.get("types", []))} for p in panels])
P.to_csv(OUT / "panels.tsv", sep="\t", index=False)
log.info("panels: %d", len(P))

cache = OUT / "green_genes_raw.json.gz"
if cache.exists():
    genes = json.load(gzip.open(cache, "rt"))
else:
    genes = paged(f"{API}/genes/?confidence_level=3")
    json.dump(genes, gzip.open(cache, "wt"))


def grch37_ensembl(gd: dict) -> str | None:
    eg = gd.get("ensembl_genes") or {}
    if isinstance(eg, str):
        try:
            eg = json.loads(eg)
        except ValueError:
            return None
    b37 = eg.get("GRch37") or {} if isinstance(eg, dict) else {}
    for v in (b37.values() if isinstance(b37, dict) else []):
        if isinstance(v, dict) and v.get("ensembl_id"):
            return v["ensembl_id"]
    return None


rows = []
for g in genes:
    gd = g.get("gene_data") or {}
    ens_id = grch37_ensembl(gd)
    rows.append({"panel_id": g["panel"]["id"], "panel_name": g["panel"]["name"],
                 "gene_symbol": gd.get("hgnc_symbol") or gd.get("gene_symbol") or g.get("entity_name"),
                 "hgnc_id": gd.get("hgnc_id"), "ensembl_id_grch37": ens_id,
                 "mode_of_inheritance": g.get("mode_of_inheritance"),
                 "confidence_level": g.get("confidence_level")})
G = pd.DataFrame(rows)
G.to_csv(OUT / "green_genes.tsv", sep="\t", index=False)
json.dump({"api": API, "retrieved_utc": retrieved, "n_panels": len(P), "n_green_gene_entries": len(G)},
          open(OUT / "meta.json", "w"), indent=1)
for f in ("panels.tsv", "green_genes.tsv"):
    record_in_manifest(OUT / f, f"{API} (Genomics England PanelApp REST API; scripts/01b_fetch_panelapp.py)")
log.info("green gene entries: %d across %d panels", len(G), G.panel_id.nunique())
