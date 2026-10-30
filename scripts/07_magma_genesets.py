#!/usr/bin/env python3
"""MAGMA competitive gene-set analyses (complements the locus-level permutation tests).

Models (all with MAGMA's default gene-level covariates: gene size, gene density, inverse MAC):
  marginal              each set vs all other genes
  cond_NDD              + neurodevelopmental set (GDD or ID)
  cond_MEND             + any Mendelian disease gene set
Outputs: results/tables/magma_genesets.tsv
"""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import BIN, INTERIM, TABLES, bh_fdr, get_logger  # noqa: E402

log = get_logger("07_magma_genesets")
MAG = INTERIM / "magma"
MAGMA = BIN / "magma"
sets_file = INTERIM / "genesets" / "genesets.magma.txt"

runs = {
    "marginal": [],
    "cond_NDD": ["--model", "condition-hide=NDD"],
    "cond_MEND": ["--model", "condition-hide=MEND_ANY"],
}
frames = []
for trait in ["chronotype"]:
    raw = MAG / f"{trait}.genes.raw"
    for name, extra in runs.items():
        out = MAG / f"{trait}.gs_{name}"
        cmd = [MAGMA, "--gene-results", raw, "--set-annot", sets_file] + extra + ["--out", out]
        subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL)
        res = pd.read_csv(Path(str(out) + ".gsa.out"), sep=r"\s+", comment="#")
        res = res[res.TYPE == "SET"] if "TYPE" in res.columns else res
        res["model"] = name
        res["trait"] = trait
        frames.append(res)
# leave-one-gene-out for the post hoc monoamine set (marginal model; separate file so that FDR is unaffected)
uni_sym = pd.read_csv(INTERIM / "genesets" / "gene_universe.tsv", sep="\t").set_index("symbol").gene_id
mono = next(line.split() for line in open(sets_file) if line.startswith("EXP_INTD_MONOAMINE "))
loo_file = MAG / "monoamine_leave_one_out.txt"
with open(loo_file, "w") as fh:
    fh.write(" ".join(mono) + "\n")
    for s in ["DDC", "GCH1", "DNAJC12", "SLC18A2"]:
        fh.write(" ".join([f"MONO_minus_{s}"] + [g for g in mono[1:] if g != uni_sym[s]]) + "\n")  # short: MAGMA truncates
out = MAG / "chronotype.gs_monoamine_loo"
subprocess.run([str(c) for c in [MAGMA, "--gene-results", MAG / "chronotype.genes.raw", "--set-annot", loo_file,
                                 "--out", out]], check=True, stdout=subprocess.DEVNULL)
loo = pd.read_csv(Path(str(out) + ".gsa.out"), sep=r"\s+", comment="#")
loo = loo.rename(columns={"VARIABLE": "set_id"}) if "VARIABLE" in loo.columns else loo.rename(columns={"FULL_NAME": "set_id"})
loo.to_csv(TABLES / "magma_monoamine_leave_one_out.tsv", sep="\t", index=False)
log.info("monoamine leave-one-out:\n%s", loo[["set_id", "NGENES", "BETA", "P"]].to_string())

# PanelApp clinical disease panels (written by 05b_disease_panels.py): marginal model, FDR across panels
pfile = INTERIM / "genesets" / "panelapp.magma.txt"
if pfile.exists():
    out = MAG / "chronotype.gs_panelapp"
    subprocess.run([str(c) for c in [MAGMA, "--gene-results", MAG / "chronotype.genes.raw", "--set-annot", pfile,
                                     "--out", out]], check=True, stdout=subprocess.DEVNULL)
    pm = pd.read_csv(Path(str(out) + ".gsa.out"), sep=r"\s+", comment="#")
    pm = pm.rename(columns={"VARIABLE": "set_id"}) if "VARIABLE" in pm.columns else pm.rename(columns={"FULL_NAME": "set_id"})
    pm.insert(0, "panel_id", pm.set_id.str.replace("PANEL_", "", regex=False).astype(int))
    pm["q_across_panels"] = bh_fdr(pm.P.values)
    pm.to_csv(TABLES / "magma_panelapp.tsv", sep="\t", index=False)
    log.info("PanelApp panels in MAGMA: %d; p < 0.05: %d", len(pm), int((pm.P < 0.05).sum()))
else:
    log.warning("no PanelApp gene sets (run scripts/05b_disease_panels.py first); panel MAGMA skipped")

R = pd.concat(frames, ignore_index=True)
R = R.rename(columns={"VARIABLE": "set_id"}) if "VARIABLE" in R.columns else R.rename(columns={"FULL_NAME": "set_id"})
fam = pd.read_csv(INTERIM / "genesets" / "geneset_summary.tsv", sep="\t").set_index("set_id")
R["family"] = R.set_id.map(fam.family)
R["label"] = R.set_id.map(fam.label)
R["q_within_family"] = np.nan
for (m, f), idx in R.groupby(["model", "family"]).groups.items():
    R.loc[idx, "q_within_family"] = bh_fdr(R.loc[idx, "P"].values)
R.to_csv(TABLES / "magma_genesets.tsv", sep="\t", index=False)
key = ["MEND_ANY", "INH_AD", "GO_CIRCADIAN", "HPO_HP0002360", "HPO_CIRCADIAN", "NDD",
       "HPO_HP0001657", "HPO_HP0001639", "HPO_HP0004431", "EXP_INTD_MONOAMINE", "EXP_CSF_HVA", "RHY_MOUSE_GE3",
       "RHY_GTEX_DAYNIGHT"]
piv = R[R.set_id.isin(key)].pivot_table(index="set_id", columns="model", values="P").reindex(key)
log.info("\n%s", piv.to_string())
