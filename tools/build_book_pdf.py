#!/usr/bin/env python3
"""Build a print-ready PDF for the OPH book Markdown sources."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPO_ROOT.parent
BOOK_DIR = REPO_ROOT / "book"
HEADER_FILE = REPO_ROOT / "tools" / "book_pdf_header.tex"
DEFAULT_OUTPUT = REPO_ROOT / "book" / "reverse-engineering-reality-book.pdf"
BUILD_DIR = WORKSPACE_ROOT / "temp" / "book_pdf_build"
BOOK_COVER_ASSET = REPO_ROOT / "assets" / "book-cover.svg"

TITLE = "Reverse Engineering Reality"
SUBTITLE = "Observer Patch Holography as a Theory of Everything"
AUTHOR = "Bernhard Mueller"

SUPERSCRIPT_TRANSLATION = str.maketrans(
    {
        "⁰": "0",
        "¹": "1",
        "²": "2",
        "³": "3",
        "⁴": "4",
        "⁵": "5",
        "⁶": "6",
        "⁷": "7",
        "⁸": "8",
        "⁹": "9",
        "⁺": "+",
        "⁻": "-",
    }
)

SUBSCRIPT_TRANSLATION = str.maketrans(
    {
        "₀": "0",
        "₁": "1",
        "₂": "2",
        "₃": "3",
        "₄": "4",
        "₅": "5",
        "₆": "6",
        "₇": "7",
        "₈": "8",
        "₉": "9",
    }
)

TEX_SYMBOL_REPLACEMENTS = [
    ("ℏ", r"\ensuremath{\hbar}"),
    ("ℓ", r"\ensuremath{\ell}"),
    ("ℂ", r"\ensuremath{\mathbb{C}}"),
    ("ℤ", r"\ensuremath{\mathbb{Z}}"),
    ("±", r"\ensuremath{\pm}"),
    ("×", r"\ensuremath{\times}"),
    ("Δ", r"\ensuremath{\Delta}"),
    ("Λ", r"\ensuremath{\Lambda}"),
    ("Ψ", r"\ensuremath{\Psi}"),
    ("Ω", r"\ensuremath{\Omega}"),
    ("α", r"\ensuremath{\alpha}"),
    ("β", r"\ensuremath{\beta}"),
    ("γ", "gamma"),
    ("δ", r"\ensuremath{\delta}"),
    ("ε", r"\ensuremath{\epsilon}"),
    ("ζ", r"\ensuremath{\zeta}"),
    ("θ", r"\ensuremath{\theta}"),
    ("π", r"\ensuremath{\pi}"),
    ("ρ", r"\ensuremath{\rho}"),
    ("σ", r"\ensuremath{\sigma}"),
    ("ψ", r"\ensuremath{\psi}"),
    ("ω", r"\ensuremath{\omega}"),
    ("φ", r"\ensuremath{\phi}"),
    ("ϕ", r"\ensuremath{\varphi}"),
    ("𝜙", r"\ensuremath{\phi}"),
    ("†", r"\ensuremath{\dagger}"),
    ("•", "*"),
    ("→", r"\ensuremath{\rightarrow}"),
    ("↔", r"\ensuremath{\leftrightarrow}"),
    ("∂", r"\ensuremath{\partial}"),
    ("√", r"\ensuremath{\sqrt{}}"),
    ("∝", r"\ensuremath{\propto}"),
    ("∩", r"\ensuremath{\cap}"),
    ("≅", r"\ensuremath{\cong}"),
    ("≈", r"\ensuremath{\approx}"),
    ("≤", r"\ensuremath{\leq}"),
    ("≥", r"\ensuremath{\geq}"),
    ("⊆", r"\ensuremath{\subseteq}"),
    ("⊗", r"\ensuremath{\otimes}"),
    ("◇", r"\ensuremath{\Diamond}"),
    ("⟨", r"\ensuremath{\langle}"),
    ("⟩", r"\ensuremath{\rangle}"),
    ("‑", "-"),
]


def chapter_key(path: Path) -> int:
    match = re.match(r"chapter-(\d+)-", path.name)
    if not match:
        raise ValueError(f"unexpected chapter filename: {path.name}")
    return int(match.group(1))


def source_files() -> list[Path]:
    chapters = sorted(BOOK_DIR.glob("chapter-*.md"), key=chapter_key)
    appendices = sorted(BOOK_DIR.glob("appendix-*.md"))
    return [BOOK_DIR / "prologue.md", *chapters, *appendices, BOOK_DIR / "epilogue.md"]


def ensure_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"required tool not found in PATH: {name}")


def run(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def convert_svg_assets(out_dir: Path) -> dict[str, Path]:
    assets = [
        BOOK_COVER_ASSET,
        REPO_ROOT / "assets" / "pixel-constant.svg",
        REPO_ROOT / "assets" / "OPH_Unification_Diagram.svg",
        *sorted((REPO_ROOT / "assets" / "book_diagrams").glob("*.svg")),
    ]
    out_dir.mkdir(parents=True, exist_ok=True)
    converted: dict[str, Path] = {}
    for source_path in assets:
        relative_path = source_path.relative_to(REPO_ROOT)
        original_ref = f"../{relative_path.as_posix()}"
        output_path = out_dir / relative_path.with_suffix(".pdf").relative_to("assets")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        raw_output_path = output_path.with_name(f"{output_path.stem}.raw.pdf")
        run(
            [
                "rsvg-convert",
                "-f",
                "pdf",
                "-o",
                str(raw_output_path),
                str(source_path),
            ]
        )
        run(
            [
                "gs",
                "-q",
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=pdfwrite",
                "-dCompatibilityLevel=1.5",
                f"-sOutputFile={output_path}",
                str(raw_output_path),
            ]
        )
        raw_output_path.unlink()
        converted[original_ref] = output_path.resolve()
    return converted


def rewrite_heading(line: str, unnumbered: bool) -> str:
    match = re.match(r"^(#{1,6})\s+(.*)$", line)
    if not match:
        return line

    marks, raw_title = match.groups()
    title = raw_title.strip()
    title = re.sub(r"^Chapter\s+\d+:\s*", "", title)
    title = re.sub(r"^\d+(?:\.\d+)*\s+", "", title)

    if unnumbered:
        return f"{marks} {title} {{.unnumbered}}\n"
    return f"{marks} {title}\n"


def rewrite_assets(text: str, asset_map: dict[str, Path]) -> str:
    for source_ref, built_path in asset_map.items():
        text = text.replace(f"({source_ref})", f"({built_path.as_posix()})")
    return text


def normalize_tex_symbols(text: str) -> str:
    text = re.sub(
        r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻]+",
        lambda match: rf"\textsuperscript{{{match.group(0).translate(SUPERSCRIPT_TRANSLATION)}}}",
        text,
    )
    text = re.sub(
        r"[₀₁₂₃₄₅₆₇₈₉]+",
        lambda match: rf"\textsubscript{{{match.group(0).translate(SUBSCRIPT_TRANSLATION)}}}",
        text,
    )
    for source, replacement in TEX_SYMBOL_REPLACEMENTS:
        text = text.replace(source, replacement)

    text = re.sub(
        r"\\\[\s*\n\s*\(m_u,m_d,m_s,m_c,m_b,m_t\)=\s*\n\s*\(0\.00216,\\ 0\.00470,\\ 0\.0935,\\ 1\.273,\\ 4\.183,\\ 172\.3523553288311\)\\,\\mathrm\{GeV\}\.\s*\n\s*\\\]",
        lambda _match: (
            "\\[\n\\begin{aligned}\n"
            "(m_u,m_d,m_s,m_c,m_b,m_t)={}&(0.00216,\\ 0.00470,\\ 0.0935,\\\\\n"
            "&1.273,\\ 4.183,\\ 172.3523553288311)\\,\\mathrm{GeV}.\n"
            "\\end{aligned}\n\\]"
        ),
        text,
    )
    text = text.replace(
        "\\[\n  (m_u,m_d,m_s,m_c,m_b,m_t)=\n  (0.00216,\\ 0.00470,\\ 0.0935,\\ 1.273,\\ 4.183,\\ 172.3523553288311)\\,\\mathrm{GeV}.\n  \\]",
        "\\[\n\\begin{aligned}\n(m_u,m_d,m_s,m_c,m_b,m_t)={}&(0.00216,\\ 0.00470,\\ 0.0935,\\\\\n&1.273,\\ 4.183,\\ 172.3523553288311)\\,\\mathrm{GeV}.\n\\end{aligned}\n\\]",
    )
    text = text.replace(
        "\\[\n\\alpha^{-1}(0)=137.035999177(21), \\qquad\nP=1.630968209403959324879279847782648941\\ldots.\n\\]",
        "\\[\n\\begin{aligned}\n\\alpha^{-1}(0)&=137.035999177(21),\\\\\nP&=1.630968209403959324879279847782648941\\ldots.\n\\end{aligned}\n\\]",
    )
    return text


def polish_tex_layout(text: str) -> str:
    def mark_unnumbered_chapter(match: re.Match[str]) -> str:
        chapter_command = match.group(1)
        raw_title = " ".join(match.group(2).split())
        running_head = raw_title
        if raw_title.startswith("Prologue:"):
            running_head = "Prologue"
        elif raw_title.startswith("Epilogue:"):
            running_head = "Epilogue"
        return f"{chapter_command}\n\\chaptermark{{{running_head}}}"

    text = re.sub(
        r"(\\chapter\*\{(.*?)\})",
        mark_unnumbered_chapter,
        text,
        flags=re.DOTALL,
    )

    return text.replace(
        "{\\def\\LTcaptype{none} % do not increment counter\n\\begin{longtable}",
        "{\\small\\setlength{\\tabcolsep}{4.5pt}\\def\\LTcaptype{none} % do not increment counter\n\\begin{longtable}",
    )


def render_manuscript(asset_map: dict[str, Path]) -> str:
    metadata = f"""---
title: "{TITLE}"
subtitle: "{SUBTITLE}"
author: "{AUTHOR}"
date: ""
documentclass: book
classoption:
  - twoside
  - openany
fontsize: 10pt
lang: en-US
colorlinks: true
linkcolor: Accent
urlcolor: Accent
toc: true
toc-depth: 0
numbersections: true
secnumdepth: 1
geometry:
  - paperwidth=6in
  - paperheight=9in
  - top=0.72in
  - bottom=0.9in
  - inner=0.82in
  - outer=0.7in
---

"""
    sections = [metadata]

    for path in source_files():
        unnumbered = path.name in {"prologue.md", "epilogue.md"}
        rewritten_lines = [rewrite_heading(line, unnumbered) for line in path.read_text().splitlines(keepends=True)]
        body = "".join(rewritten_lines).strip()
        body = rewrite_assets(body, asset_map)
        sections.append(body)
        sections.append("")

    return "\n\n".join(sections).rstrip() + "\n"


def build(output: Path) -> None:
    ensure_tool("pandoc")
    ensure_tool("tectonic")
    ensure_tool("rsvg-convert")
    ensure_tool("gs")

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    asset_dir = BUILD_DIR / "assets"
    asset_map = convert_svg_assets(asset_dir)

    manuscript_path = BUILD_DIR / "book_manuscript.md"
    tex_path = BUILD_DIR / "book_manuscript.tex"
    manuscript_path.write_text(render_manuscript(asset_map))

    run(
        [
            "pandoc",
            str(manuscript_path),
            "--standalone",
            "--from=markdown+raw_tex+tex_math_dollars+tex_math_single_backslash",
            "--to=latex",
            "--top-level-division=chapter",
            "--number-sections",
            f"--include-in-header={HEADER_FILE}",
            f"--output={tex_path}",
        ],
        cwd=REPO_ROOT,
    )

    tex = tex_path.read_text()
    tex = normalize_tex_symbols(tex)
    tex = polish_tex_layout(tex)
    tex_path.write_text(tex)

    run(
        [
            "tectonic",
            "-X",
            "compile",
            tex_path.name,
            "--keep-logs",
            f"--outdir={BUILD_DIR}",
        ],
        cwd=BUILD_DIR,
    )

    built_pdf = BUILD_DIR / f"{tex_path.stem}.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built_pdf, output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the OPH book PDF.")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output PDF path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build(args.output.resolve())
    print(f"Built book PDF at {args.output.resolve()}")


if __name__ == "__main__":
    main()
