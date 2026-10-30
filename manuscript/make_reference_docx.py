#!/usr/bin/env python3
"""Create manuscript/reference.docx (pandoc style template): Arial 11 pt, 1.5 line spacing,
black headings, A4 with 2.54 cm margins, continuous line numbers and a centred page-number footer."""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
tmp = Path(tempfile.mkdtemp())
default = tmp / "default.docx"
with open(default, "wb") as fh:
    fh.write(subprocess.run(["pandoc", "--print-default-data-file", "reference.docx"], check=True,
                            capture_output=True).stdout)
src = tmp / "src"
with zipfile.ZipFile(default) as z:
    z.extractall(src)

# ---------------- styles
st = (src / "word/styles.xml").read_text()
st = re.sub(r'<w:rFonts [^>]*/>', '<w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Arial" w:cs="Arial" />', st)
st = st.replace('<w:sz w:val="24" />\n        <w:szCs w:val="24" />', '<w:sz w:val="22" />\n        <w:szCs w:val="22" />', 1)
st = st.replace('<w:spacing w:after="200" />', '<w:spacing w:after="120" w:line="360" w:lineRule="auto" />', 1)
st = st.replace('<w:spacing w:before="180" w:after="180" />', '<w:spacing w:before="0" w:after="160" />')
st = re.sub(r'<w:color w:val="[0-9A-Fa-f]{6}" ?/>', '<w:color w:val="000000" />', st)
st = re.sub(r'<w:color w:val="[0-9A-Fa-f]{6}" w:themeColor="[^"]+"[^>]*/>', '<w:color w:val="000000" />', st)


def set_size(style_id: str, half_points: int, bold: bool = True) -> None:
    global st
    m = re.search(rf'(<w:style [^>]*w:styleId="{style_id}">)(.*?)(</w:style>)', st, re.S)
    if not m:
        return
    body = m.group(2)
    body = re.sub(r'<w:sz w:val="\d+" ?/>', '', body)
    body = re.sub(r'<w:szCs w:val="\d+" ?/>', '', body)
    if "<w:rPr>" in body:
        body = body.replace("<w:rPr>", f'<w:rPr><w:sz w:val="{half_points}" /><w:szCs w:val="{half_points}" />'
                            + ("<w:b />" if bold and "<w:b />" not in body else ""), 1)
    else:
        body += f'<w:rPr><w:sz w:val="{half_points}" /><w:szCs w:val="{half_points}" />{"<w:b />" if bold else ""}</w:rPr>'
    st = st[:m.start()] + m.group(1) + body + m.group(3) + st[m.end():]


for sid, size in (("Title", 32), ("Heading1", 28), ("Heading2", 24), ("Heading3", 22)):
    set_size(sid, size)
set_size("Bibliography", 20, bold=False)
# tables: thin top/bottom rule, header rule, compact 9 pt text
st = st.replace('<w:style w:type="table" w:default="1" w:styleId="Table">',
                '<w:style w:type="table" w:default="1" w:styleId="Table">')
m = re.search(r'(<w:style w:type="table"[^>]*w:styleId="Table">)(.*?)(</w:style>)', st, re.S)
if m:
    tbl = ('<w:name w:val="Table" /><w:basedOn w:val="TableNormal" /><w:uiPriority w:val="99" /><w:qFormat />'
           '<w:pPr><w:spacing w:before="20" w:after="20" w:line="240" w:lineRule="auto" /></w:pPr>'
           '<w:rPr><w:sz w:val="17" /><w:szCs w:val="17" /></w:rPr>'
           '<w:tblPr><w:tblBorders><w:top w:val="single" w:sz="8" w:space="0" w:color="000000" />'
           '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000" />'
           '<w:insideH w:val="single" w:sz="2" w:space="0" w:color="BFBFBF" /></w:tblBorders>'
           '<w:tblCellMar><w:left w:w="60" w:type="dxa" /><w:right w:w="60" w:type="dxa" /></w:tblCellMar></w:tblPr>'
           '<w:tblStylePr w:type="firstRow"><w:rPr><w:b /></w:rPr><w:tcPr><w:tcBorders>'
           '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000" /></w:tcBorders></w:tcPr></w:tblStylePr>')
    st = st[:m.start()] + m.group(1) + tbl + m.group(3) + st[m.end():]
(src / "word/styles.xml").write_text(st)

# ---------------- footer with page number
(src / "word/footer1.xml").write_text(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<w:p><w:pPr><w:suppressLineNumbers /><w:jc w:val="center" /></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin" /></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate" /></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end" /></w:r>'
    '</w:p></w:ftr>')
rels = (src / "word/_rels/document.xml.rels").read_text()
if "footer1.xml" not in rels:
    rels = rels.replace("</Relationships>",
                        '<Relationship Id="rIdFooter1" '
                        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
                        'Target="footer1.xml"/></Relationships>')
    (src / "word/_rels/document.xml.rels").write_text(rels)
ct = (src / "[Content_Types].xml").read_text()
if "footer1.xml" not in ct:
    ct = ct.replace("</Types>", '<Override PartName="/word/footer1.xml" '
                    'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>')
    (src / "[Content_Types].xml").write_text(ct)

# ---------------- section properties: A4, margins, line numbers, footer
doc = (src / "word/document.xml").read_text()
sect = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter1" />'
        '<w:pgSz w:w="11906" w:h="16838" />'
        '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0" />'
        '<w:lnNumType w:countBy="1" w:distance="283" w:restart="continuous" />'
        '<w:cols w:space="708" /></w:sectPr>')
if "<w:sectPr" in doc:
    doc = re.sub(r"<w:sectPr\s*/>|<w:sectPr\b.*?</w:sectPr>", lambda _m: sect, doc, count=1, flags=re.S)
else:
    doc = doc.replace("</w:body>", sect + "</w:body>")
(src / "word/document.xml").write_text(doc)

out = HERE / "reference.docx"
if out.exists():
    out.unlink()
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(src.rglob("*")):
        if f.is_file():
            z.write(f, f.relative_to(src).as_posix())
shutil.rmtree(tmp)
print("wrote", out)
