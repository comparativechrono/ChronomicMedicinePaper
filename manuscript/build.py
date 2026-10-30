#!/usr/bin/env python3
"""Render the manuscript: fill every {{path@format}} placeholder from results/key_numbers.json,
insert generated tables, then build the Word document with pandoc + citeproc.

Placeholders: {{a::b::c@fmt}} walks key_numbers.json (list indices allowed);
formats: p (p-value, superscript exponent), pm (MAGMA p with relation sign and floor),
pr (permutation p with relation sign and floor), genes (semicolon list as italic gene names),
i (integer), f0/f1/f2/f3 (fixed decimals),
pct (fraction as percentage, 1 dp), raw (as is). {{TABLE:name}} inserts tables/name.md.
Unresolved placeholders abort the build, so no number in the text can silently go stale.

Usage: python manuscript/build.py
"""
from __future__ import annotations

import json
import math
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
K = json.load(open(ROOT / "results" / "key_numbers.json"))
OUTNAME = "Malik_Hearn_Towards_chronomic_medicine_2026.docx"
PAGEBREAK = '\n```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n'


def lookup(path: str):
    node = K
    for part in path.split("::"):
        if isinstance(node, list):
            node = node[int(part)]
        else:
            node = node[part]
    return node


def fmt(v, f: str) -> str:
    if f in ("", "raw"):
        return str(v)
    if f == "genes":  # "A;B;C" -> "*A*, *B* and *C*"
        g = [f"*{x}*" for x in str(v).split(";") if x]
        return g[0] if len(g) == 1 else ", ".join(g[:-1]) + " and " + g[-1]
    v = float(v)
    if f == "pm":   # MAGMA p-value with relation sign; MAGMA's numerical floor is ~1e-15
        return "< 1×10^−15^" if v < 1e-15 else "= " + fmt(v, "p")
    if f == "pr":   # permutation p-value with relation sign; floor 1/(1 + 10,000)
        return "≤ 1.0×10^−4^" if v <= 1 / 10001 + 1e-12 else "= " + fmt(v, "p")
    if f == "p":
        if v < 1e-3:
            m, e = f"{v:.1e}".split("e")
            return f"{m}×10^−{int(e[1:])}^"
        if v >= 0.995:
            return "1.00"
        if v >= 0.1:
            return f"{v:.2f}"
        return f"{v:.{1 - int(math.floor(math.log10(v)))}f}"   # two significant digits, trailing zeros kept
    if f == "i":
        return f"{int(round(v)):,}"
    if f == "pct":
        return f"{100 * v:.1f}%"
    if f.startswith("f"):
        return f"{v:.{int(f[1:])}f}"
    raise ValueError(f)


def render(text: str) -> str:
    text = re.sub(r"\{\{TABLE:([a-z0-9_]+)\}\}", lambda m: (HERE / "tables" / f"{m.group(1)}.md").read_text(), text)

    def sub(m):
        path, _, f = m.group(1).partition("@")
        try:
            return fmt(lookup(path.strip()), f.strip())
        except Exception as exc:  # fail loudly
            raise SystemExit(f"Unresolved placeholder {{{{{m.group(1)}}}}}: {exc!r}")

    text = re.sub(r"\{\{([^{}]+)\}\}", sub, text)
    return text.replace("\\newpage", PAGEBREAK)


def main() -> int:
    src = (HERE / "manuscript.template.md").read_text()
    md = render(src)
    (HERE / "manuscript.md").write_text(md)
    cmd = ["pandoc", "manuscript.md", "--citeproc", "--bibliography", "references.bib", "--csl", "ama.csl",
           "--reference-doc", "reference.docx", "--resource-path", str(HERE), "-o", OUTNAME]
    subprocess.run(cmd, check=True, cwd=HERE)
    print("built", HERE / OUTNAME)
    return 0


if __name__ == "__main__":
    sys.exit(main())
