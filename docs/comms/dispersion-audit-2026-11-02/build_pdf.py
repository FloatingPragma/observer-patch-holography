"""Render the four-page IPI Markdown note with ReportLab and mathtext.

Dependencies: reportlab, matplotlib, pillow, pypdf. The Markdown is the source;
page breaks are explicit so the short technical note stays under five pages.
"""
from __future__ import annotations

import hashlib
import html
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, Image, KeepTogether
from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = HERE.parent/"OPH_FROZEN_DISPERSION_AUDIT_2026-11-02.md"
OUTPUT = SOURCE.with_suffix(".pdf")
SCRATCH = ROOT/"tmp/pdfs/oph-dispersion-math"
SCRATCH.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

FONTS = Path("/System/Library/Fonts/Supplemental")
for name, file in (("NoteSerif", "Times New Roman.ttf"),
                   ("NoteSerif-Bold", "Times New Roman Bold.ttf"),
                   ("NoteSerif-Italic", "Times New Roman Italic.ttf"),
                   ("NoteSerif-BoldItalic", "Times New Roman Bold Italic.ttf")):
    pdfmetrics.registerFont(TTFont(name, str(FONTS/file)))
pdfmetrics.registerFontFamily("NoteSerif", normal="NoteSerif", bold="NoteSerif-Bold",
                              italic="NoteSerif-Italic", boldItalic="NoteSerif-BoldItalic")
matplotlib.rcParams.update({"mathtext.fontset":"stix", "font.family":"STIXGeneral",
                           "savefig.transparent": True})
INK = colors.HexColor("#172532")
BLUE = colors.HexColor("#235676")
WIDTH = A4[0]-88
styles = {
    "body": ParagraphStyle("body", fontName="NoteSerif", fontSize=10.5, leading=13.4,
                           spaceAfter=6, textColor=INK, splitLongWords=False),
    "h1": ParagraphStyle("h1", fontName="NoteSerif-Bold", fontSize=20, leading=22.5,
                         spaceAfter=9, textColor=INK),
    "h2": ParagraphStyle("h2", fontName="NoteSerif-Bold", fontSize=13, leading=15.5,
                         spaceBefore=5, spaceAfter=7, textColor=BLUE, keepWithNext=True),
    "h3": ParagraphStyle("h3", fontName="NoteSerif-Bold", fontSize=10.5, leading=12,
                         spaceBefore=3, spaceAfter=4, textColor=BLUE, keepWithNext=True),
    "subtitle": ParagraphStyle("subtitle", fontName="NoteSerif", fontSize=9.3, leading=12,
                               spaceAfter=10, textColor=BLUE),
    "refs": ParagraphStyle("refs", fontName="NoteSerif", fontSize=8.5, leading=10.4,
                           spaceAfter=3, textColor=INK, splitLongWords=False),
    "table": ParagraphStyle("table", fontName="NoteSerif", fontSize=9.3, leading=11.7,
                            textColor=INK),
    "table_center": ParagraphStyle("table_center", fontName="NoteSerif", fontSize=9.3,
                                   leading=11.7, textColor=INK, alignment=TA_CENTER),
}


def math_image(expr, size):
    expr = re.sub(r"\\(mathbf|mathcal|mathrm|mathbb)\s+([A-Za-z])", r"\\\1{\2}", expr)
    expr = re.sub(r"\\frac([0-9])([0-9])", r"\\frac{\1}{\2}", expr)
    expr = re.sub(r"\\le(?![A-Za-z])", r"\\leq", expr)
    expr = re.sub(r"\\ge(?![A-Za-z])", r"\\geq", expr)
    expr = re.sub(r"\\sqrt([0-9])", r"\\sqrt{\1}", expr)
    expr = re.sub(r"\\sqrt\s+([A-Za-z])", r"\\sqrt{\1}", expr)
    key = hashlib.sha256(f"transparent-v2:{size}:{expr}".encode()).hexdigest()[:18]
    path = SCRATCH/f"{key}.png"
    dpi = 300
    if not path.exists():
        mathtext.math_to_image("$"+expr+"$", str(path), dpi=dpi,
                              prop=FontProperties(size=size), format="png", color="#172532")
    with PILImage.open(path) as im:
        width, height = im.size
    return path, width*72/dpi, height*72/dpi


def inline(text, size=10.5):
    # Isolate formula and link placeholders before escaping source text.
    placeholders = []
    def slot(value):
        placeholders.append(value)
        return f"@@PLACEHOLDER{len(placeholders)-1}@@"
    def equation(match):
        path, w, h = math_image(match.group(1), size)
        return slot(f'<img src="{path}" width="{w:.3f}" height="{h:.3f}" valign="middle"/>')
    text = re.sub(r"\$([^$]+)\$", equation, text)
    def link(match):
        target = match.group(2)
        if not target.startswith("http"):
            target = "https://github.com/FloatingPragma/observer-patch-holography/blob/main/docs/comms/" + target
        return slot(f'<a href="{html.escape(target, quote=True)}" color="#235676">{html.escape(match.group(1))}</a>')
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
    text = html.escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`([^`]+)`", r'<font size="8.2">\1</font>', text)
    for index, value in enumerate(placeholders):
        text = text.replace(f"@@PLACEHOLDER{index}@@", value)
    return text


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#becbd3"))
    canvas.setLineWidth(.45)
    canvas.line(44, 35, A4[0]-44, 35)
    canvas.setFont("NoteSerif", 8)
    canvas.setFillColor(BLUE)
    canvas.drawString(44, 23, "OPH | Dispersion, direct tests and field reconstruction | 2 November 2026")
    canvas.drawRightString(A4[0]-44, 23, str(doc.page)+" / 4")
    canvas.restoreState()


def main():
    text = SOURCE.read_text()
    text = text.replace("\u2013", "-").replace("\u2014", "-").replace("\u2011", "-")
    # Keep the portable Markdown source's prose punctuation consistent too.
    SOURCE.write_text(text)
    story = []
    refs = False
    for page_index, page in enumerate(text.split("<!-- PAGEBREAK -->")):
        if page_index:
            story.append(PageBreak())
        for block in re.split(r"\n\s*\n", page.strip()):
            if not block.strip():
                continue
            if block.startswith("$$"):
                expr = block.strip()[2:-2].strip()
                path, w, h = math_image(expr, 11.4)
                scale = min(1., WIDTH/w)
                img = Image(str(path), width=w*scale, height=h*scale)
                img.hAlign = "CENTER"
                story.extend((Spacer(1, 3), img, Spacer(1, 7)))
            elif block.startswith("#"):
                prefix, heading = block.split(" ", 1)
                key = {"#":"h1", "##":"h2", "###":"h3"}[prefix]
                if heading == "Sources and reproducibility":
                    refs = True
                story.append(Paragraph(inline(heading, styles[key].fontSize), styles[key]))
            elif block.startswith("|"):
                rows = [line.strip().strip("|").split("|") for line in block.splitlines()]
                rows = [row for row in rows if not all(re.match(r"^[\s:\-]+$", cell) for cell in row)]
                parsed = []
                for row_index, row in enumerate(rows):
                    parsed.append([Paragraph(inline(("**"+cell.strip()+"**") if row_index==0 else cell.strip(),9.3),
                                              styles["table" if i==0 else "table_center"])
                                   for i,cell in enumerate(row)])
                table = Table(parsed, colWidths=[173,70,72,94,WIDTH-409])
                table.setStyle(TableStyle([
                    ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf0f4")),
                    ("LINEBELOW",(0,0),(-1,0),.6,BLUE),
                    ("LINEBELOW",(0,-1),(-1,-1),.45,colors.HexColor("#a8bac5")),
                    ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
                    ("TOPPADDING",(0,0),(-1,-1),5),
                    ("BOTTOMPADDING",(0,0),(-1,-1),5),
                    ("LEFTPADDING",(0,0),(-1,-1),5),
                    ("RIGHTPADDING",(0,0),(-1,-1),5),
                ]))
                story.extend((table,Spacer(1,8)))
            else:
                key = "refs" if refs else ("subtitle" if page_index==0 and block.startswith("**Technical reply") else "body")
                story.append(Paragraph(inline(" ".join(block.splitlines()),styles[key].fontSize),styles[key]))
    doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, rightMargin=44, leftMargin=44,
                            topMargin=38, bottomMargin=46,
                            title="OPH: how records become space, clocks and testable physics",
                            author="Bernhard Mueller - technical audit prepared with Codex")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    reader = PdfReader(str(OUTPUT))
    counts = [len((p.extract_text() or "").split()) for p in reader.pages]
    print({"output":str(OUTPUT),"pages":len(reader.pages),"words_per_page":counts})
    if len(reader.pages)!=4:
        raise SystemExit("Page budget failed: revise layout or prose; expected four pages.")


if __name__ == "__main__":
    main()
