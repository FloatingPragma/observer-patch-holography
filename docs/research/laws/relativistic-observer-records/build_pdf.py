"""Typeset report.md and its vector figures with Pandoc 3 and Typst 0.14."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
STYLE = r'''
#let ink = rgb("#172830")
#let teal = rgb("#086b70")
#let cyan = rgb("#4af1f2")
#let gray = rgb("#566870")
#let rule = rgb("#c6d7d8")
#let sans = ("Helvetica Neue", "Arial", "Libertinus Serif")
#let mono = ("Menlo", "DejaVu Sans Mono")
#let caret(col, unit) = box(width: 15*unit, height: 14*unit, {
  for (x,y,w) in ((6,0,3),(5,2,5),(4,4,3),(8,4,3),(3,6,3),(9,6,3),
    (2,8,3),(10,8,3),(1,10,3),(11,10,3),(0,12,3),(12,12,3)) {
    place(top+left, dx: x*unit, dy: y*unit, rect(width:w*unit, height:2*unit, fill:col, stroke:none))
  }
})
#set document(title: "Clocks and spacetime from observer records", author: "Bernhard Mueller, Pragma Research",
  description: "Observer Patch Holography (OPH) experiments on traveling clocks, field measurements and candidate repair laws.")
#set text(font: ("Libertinus Serif", "Palatino"), size: 11.7pt, fill: ink, lang:"en")
#show math.equation: set text(font: ("STIX Two Math", "New Computer Modern Math"))
#set par(leading: .61em, spacing: .92em, justify: false)
#show link: it => text(fill: teal)[#it]
#set image(width: 100%)
#show figure.where(kind:image): it => align(center, it.body)
#show figure.where(kind:table): set text(font: sans, size: 9.75pt)
#show table.cell.where(y:0): set text(weight:"bold", fill:teal)
#set table(inset:(x:6pt,y:5.4pt), fill:(x,y)=>if y==0 {rgb("#e3f1ef")},
  stroke:(x,y)=>(bottom: .45pt+rule))
#show raw.where(block:true): it => block(width:100%,inset:9pt,fill:rgb("#edf4f3"),text(size:8.2pt,it))
#show heading.where(level:2): it => block(width: 100%, above:0pt, below:13pt)[
  #set text(font:sans,size:27pt,weight:"bold",tracking:-.55pt)
  #set par(leading:.22em)
  #it.body
  #v(6pt)
  #line(length:100%,stroke:1pt+teal)
]
#set page(paper:"a4",margin:(x:21mm,y:23mm),fill:rgb("#0b181e"))
'''

COVER_START = r'''
#set text(fill:rgb("#f2f6f4"))
#grid(columns:(1fr, auto), text(font:mono,fill:cyan,size:10pt,tracking:1pt)[PRAGMA RESEARCH],
  text(font:mono,size:9pt,fill:rgb("#adc9cf"))[OPH / RESEARCH REPORT])
#v(4mm)
#line(length:100%,stroke:.55pt+teal)
#v(11mm)
#caret(cyan,1.55mm)
#v(9mm)
#set par(leading:.36em)
'''

BODY_START = r'''
#set page(fill:white, margin:(x:21mm,top:24mm,bottom:21mm),
  header: grid(columns:(auto,1fr,auto),column-gutter:6pt,caret(teal,.22mm),
    text(font:mono,size:8.3pt,fill:teal)[PRAGMA RESEARCH],
    text(font:sans,size:8.5pt,fill:gray)[Observer Patch Holography / Observer records]),
  footer: context [
    #line(length:100%,stroke:.45pt+rule)
    #v(1.5mm)
    #grid(columns:(1fr,auto),text(font:sans,fill:gray,size:8.2pt)[Clocks and spacetime from observer records],
      text(font:mono,fill:teal,size:8.2pt)[#counter(page).display("1")])
  ])
#set text(fill:ink,size:11.7pt)
#set par(leading:.61em,spacing:.92em)
'''


def typeset():
    sections = subprocess.run(
        ["pandoc", "--from=markdown-implicit_figures", "--to=typst", "--wrap=none", "report.md"],
        cwd=HERE, check=True, capture_output=True, text=True
    ).stdout.split("\n== ")
    opening, cover_intro, cover_scope = sections[0].strip().split("\n\n")
    title, _, cover_tag = opening.splitlines()
    title = title.removeprefix("= ")
    cover = STYLE + COVER_START
    cover += f'#text(font:sans,size:40pt,weight:"bold",tracking:-1pt)[{title}]\n'
    cover += '#v(6mm)\n#text(font:sans,size:11pt,fill:cyan)[' + cover_tag + ']\n'
    cover += '#v(9mm)\n#set text(size:15pt)\n#set par(leading:.53em)\n' + cover_intro + '\n\n'
    cover += '#v(3mm)\n#text(size:12.2pt,fill:rgb("#c9dadb"))[' + cover_scope + ']\n'
    cover += '#v(1fr)\n#line(length:100%,stroke:.55pt+teal)\n'
    cover += r'''
#let chapter(mark, name) = link(mark)[#text(font:sans,size:10.5pt,fill:rgb("#d3e4e4"))[#name]
  #h(1fr)#text(font:mono,size:9pt,fill:cyan)[#context counter(page).at(mark).first()]]
#grid(columns:(1fr,1fr),column-gutter:14mm,row-gutter:4mm,
 chapter(<what-special-relativity-asks-of-a-clock>,"Proper time and clocks"),
 chapter(<what-the-continuum-proof-establishes>,"The field proof"),
 chapter(<which-repair-laws-do-these-experiments-select>,"Which laws survive"),
 chapter(<how-to-enumerate-repair-laws-systematically>,"The search program"))
'''
    parts = [cover, BODY_START]
    for section in sections[1:]:
        # The default equal-width text columns waste space on short labels.
        section = section.replace("columns: (50%, 50%),", "columns: (30%, 70%),")
        section = section.replace("columns: (20%, 26.67%, 26.67%, 26.67%),", "columns: (35%, 22%, 22%, 21%),")
        section = re.sub(r"align: \(([^)]+)\),", lambda match: "align: (" + match[1].replace("auto", "left") + "),", section)
        if section.startswith("The rest of"):
            section = re.sub(r"columns: \([^)]+\),", "columns: (9%, 45%, 46%),", section, count=1)
        parts.append("== " + section.strip() + "\n")
        parts.append("#pagebreak(weak:true)\n")
    with tempfile.TemporaryDirectory() as temp:
        scratch = Path(temp)
        shutil.copytree(HERE / "figures", scratch / "figures")
        (scratch / "report.typ").write_text("\n".join(parts[:-1]))
        subprocess.run(["typst", "compile", "report.typ", "report.pdf"], cwd=scratch, check=True)
        (HERE / "report.pdf").write_bytes((scratch / "report.pdf").read_bytes())


if __name__ == "__main__":
    typeset()
