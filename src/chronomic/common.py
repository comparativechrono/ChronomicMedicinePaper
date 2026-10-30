"""Paths, configuration and small helpers used by every script."""
from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
LOGS = ROOT / "results" / "logs"
BIN = ROOT / "bin"

for _d in (INTERIM, TABLES, FIGURES, LOGS):
    _d.mkdir(parents=True, exist_ok=True)


def load_config() -> dict:
    with open(ROOT / "config" / "config.yaml") as fh:
        return yaml.safe_load(fh)


def get_logger(name: str) -> logging.Logger:
    log = logging.getLogger(name)
    if not log.handlers:
        log.setLevel(logging.INFO)
        fmt = logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        fh = logging.FileHandler(LOGS / f"{name}.log", mode="w")
        fh.setFormatter(fmt)
        log.addHandler(sh)
        log.addHandler(fh)
    return log


def run(cmd: list[str], log: logging.Logger | None = None) -> None:
    """Run an external command, raising on failure."""
    if log:
        log.info("$ %s", " ".join(map(str, cmd)))
    subprocess.run([str(c) for c in cmd], check=True)


def bh_fdr(p: np.ndarray | pd.Series) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs preserved)."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def save_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, (set,)):
        return sorted(o)
    return str(o)


def chrom_norm(c) -> str:
    c = str(c).replace("chr", "")
    return {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(c, c)


def env_snapshot() -> dict:
    import platform
    import importlib.metadata as md
    pkgs = {}
    for p in ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "obonet", "networkx",
              "pyyaml", "openpyxl", "requests", "pyarrow"]:
        try:
            pkgs[p] = md.version(p)
        except md.PackageNotFoundError:
            pass
    return {"python": platform.python_version(), "packages": pkgs,
            "plink": _tool_version([BIN / "plink", "--version"]),
            "plink2": _tool_version([BIN / "plink2", "--version"]),
            "magma": _tool_version([BIN / "magma", "--version"])}


def _tool_version(cmd) -> str | None:
    try:
        out = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=30)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception:
        return None


def record_in_manifest(path: Path, url: str) -> None:
    """Add or update a manifest row (path, url, time, bytes, MD5) for data retrieved by API rather than download."""
    import datetime as _dt
    import hashlib as _hl
    man = ROOT / "data" / "manifest.tsv"
    rows = {}
    lines = man.read_text().splitlines() if man.exists() else ["path\turl\tdownloaded_utc\tbytes\tmd5"]
    for line in lines[1:]:
        rows[line.split("\t")[0]] = line
    h = _hl.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    when = _dt.datetime.fromtimestamp(path.stat().st_mtime, _dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    rel = str(path.resolve().relative_to(ROOT))
    rows[rel] = "\t".join([rel, url, when, str(path.stat().st_size), h.hexdigest()])
    man.write_text("\n".join([lines[0]] + [rows[k] for k in sorted(rows)]) + "\n")
