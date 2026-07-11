#!/usr/bin/env python3
"""Export the canonical OPH paper set to local Markdown copies."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPO_ROOT.parent
PAPER_DIR = REPO_ROOT / "paper"
EXTRA_DIR = REPO_ROOT / "extra"
COSMOLOGY_DIR = REPO_ROOT / "cosmology"
DEFAULT_OUT = WORKSPACE_ROOT / "markdown"
NON_PAPER_TEX = {
    "appendix_B_bft_qecc_extensions.tex",
    "release_info.tex",
}
DEFAULT_CORE_PAPERS = sorted(
    path for path in PAPER_DIR.glob("*.tex") if path.name not in NON_PAPER_TEX
)
DEFAULT_SUPPLEMENTAL_PAPERS = [
]
DEFAULT_EXTRA_PAPERS = sorted(EXTRA_DIR.glob("*.tex"))
DEFAULT_RELEASED_COSMOLOGY_PAPERS = [
    COSMOLOGY_DIR / "oph_dark_matter_paper.tex",
]
DEFAULT_SOURCES = [
    *DEFAULT_CORE_PAPERS,
    *DEFAULT_SUPPLEMENTAL_PAPERS,
    *DEFAULT_EXTRA_PAPERS,
    *DEFAULT_RELEASED_COSMOLOGY_PAPERS,
]
BUILD_INFO_NAME = "_build_info.json"
COMPACT_PAPER_SOURCE = (
    PAPER_DIR / "recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex"
)
COMPACT_MARKDOWN_SOURCE = PAPER_DIR / "tex_fragments" / "PAPER.tex"
MARKDOWN_SOURCE_OVERRIDES = {
    COMPACT_PAPER_SOURCE: COMPACT_MARKDOWN_SOURCE,
}


def postprocess_markdown(text: str) -> str:
    def protect_math_table_pipes(match: re.Match[str]) -> str:
        body = match.group(1).replace("|", r"\vert{}")
        return f"$`{body}`$"

    def protect_table_row_pipes(line: str) -> str:
        if not line.lstrip().startswith("|"):
            return line
        return re.sub(r"\$`([^`]*\|[^`]*)`\$", protect_math_table_pipes, line)

    def table_matches(header: str, required: tuple[str, ...]) -> bool:
        normalized = header.lower()
        normalized = re.sub(r"<sub>\s*([^<]+?)\s*</sub>", r"_\1", normalized)
        normalized = re.sub(r"\$`_\{\\mathrm\{([^}]+)\}\}`\$", r"_\1", normalized)
        normalized = re.sub(r"\s+", "", normalized)
        return all(term.lower().replace(" ", "") in normalized for term in required)

    table_repairs: list[tuple[tuple[str, ...], list[str]]] = [
        (
            ("Delta f_k", "Relative weight"),
            [
                "| 2 | 41.6 | 0.500 |",
                "| 3 | 65.9 | 0.667 |",
                "| 4 | 83.2 | 0.750 |",
                "| 5 | 96.5 | 0.800 |",
                "| 6 | 107.5 | 0.833 |",
            ],
        ),
        (
            ("p_0", "p_1", "p_2", "m_plaq"),
            [
                "| 0.2 | 0.4395 | 0.2803 | 0.2803 | 0.1500 | 0.1500 | 0.154 | 2.22 |",
                "| 0.5 | 0.7509 | 0.1245 | 0.1245 | 0.5989 | 0.5989 | 0.309 | 1.75 |",
                "| 1.0 | 0.9606 | 0.0197 | 0.0197 | 1.2956 | 1.2956 | 0.454 | 4.07 |",
                "| 1.5 | 0.9851 | 0.0074 | 0.0074 | 1.6288 | 1.6288 | 0.509 | 7.06 |",
                "| 2.0 | 0.9921 | 0.0039 | 0.0039 | 1.8440 | 1.8440 | 0.542 | 10.10 |",
            ],
        ),
        (
            ("p_0", "p_1", "g_ent"),
            [
                "| 0.5 | 0.8266 | 0.1734 | 0.391 | 0.249 |",
                "| 1.0 | 0.9612 | 0.0388 | 0.803 | 0.357 |",
                "| 2.0 | 0.9917 | 0.0083 | 1.194 | 0.436 |",
            ],
        ),
        (
            ("Measured ratio", "Deviation from"),
            [
                "| 0.5 | 2.25 | 14% |",
                "| 1.0 | 2.51 | 4% |",
                "| 2.0 | 2.619 | &lt; 0.1% |",
            ],
        ),
        (
            ("p_triv", "p_sign", "p_std", "log-ratio"),
            [
                "| 0.5 | 0.909 | 0.0013 | 0.089 | 1.09 | 1.01 | 8.4% | 2.17 |",
                "| 1.0 | 0.980 | 7.5e-5 | 0.020 | 1.58 | 1.54 | 2.8% | 2.06 |",
                "| 2.0 | 0.996 | 4.3e-6 | 0.004 | 2.06 | 2.04 | 1.0% | 2.02 |",
                "| 5.0 | 0.9993 | 1.0e-7 | 0.00066 | 2.68 | 2.67 | 0.3% | 2.006 |",
                "| 12 | 0.9999 | 3.0e-9 | 0.00011 | 3.27 | 3.27 | 0.1% | 2.002 |",
                "| 100 | 1.0000 | 6.1e-13 | 2.0e-6 | 4.69 | 4.69 | 0.009% | 2.0002 |",
            ],
        ),
        (
            ("extracted t", "mean", "g_ent", "gap"),
            [
                "| 0.3 | 0.314 &plusmn; 0.0005 | 0.224 | 1.92 |",
                "| 0.5 | 0.539 &plusmn; 0.0025 | 0.293 | 1.83 |",
                "| 0.8 | 0.896 &plusmn; 0.012 | 0.378 | 1.72 |",
                "| 1.0 | 1.144 &plusmn; 0.025 | 0.427 | 1.64 |",
            ],
        ),
        (
            ("bare g", "extracted t", "A_eff"),
            [
                "| 0.3 | 0.314 | 2.093 |",
                "| 0.5 | 0.539 | 2.156 |",
                "| 0.8 | 0.896 | 2.240 |",
                "| 1.0 | 1.144 | 2.288 |",
            ],
        ),
    ]

    def repair_known_tables(markdown: str) -> str:
        def is_damaged_table_body(rows: list[str]) -> bool:
            if not rows:
                return False
            damaged_rows = 0
            for row in rows:
                cells = [cell.strip() for cell in row.strip().strip("|").split("|")]
                if not cells:
                    continue
                nonempty = sum(bool(cell) for cell in cells)
                if nonempty == 0 or (len(cells) >= 3 and nonempty <= max(1, len(cells) // 4)):
                    damaged_rows += 1
            return damaged_rows > 0

        lines = markdown.splitlines()
        repaired: list[str] = []
        index = 0
        while index < len(lines):
            line = lines[index]
            replacement_rows = None
            if line.lstrip().startswith("|") and index + 1 < len(lines):
                for required, rows in table_repairs:
                    if table_matches(line, required):
                        replacement_rows = rows
                        break
            if replacement_rows is None:
                repaired.append(line)
                index += 1
                continue

            repaired.append(line)
            repaired.append(lines[index + 1])
            index += 2
            body_rows: list[str] = []
            while index < len(lines) and lines[index].lstrip().startswith("|"):
                body_rows.append(lines[index])
                index += 1
            repaired.extend(replacement_rows if is_damaged_table_body(body_rows) else body_rows)
        return "\n".join(repaired)

    text = re.sub(
        r"(\*\*Paper release:\*\* `[^`]+`)(?=\*\*Released:\*\*)",
        r"\1  \n",
        text,
        count=1,
    )
    text = "\n".join(protect_table_row_pipes(line) for line in text.splitlines())
    text = repair_known_tables(text)
    return text


def ensure_release_banner(text: str, release_tag: str, release_date: str) -> str:
    lines = text.lstrip().splitlines()
    cleaned: list[str] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        is_top_banner_line = (
            "Paper release:" in stripped
            or stripped.startswith("**Released:**")
            or stripped.startswith("Released:")
            or stripped in {'<div class="center">', "</div>"}
            or (stripped == "" and index < 8 and any("Paper release:" in candidate for candidate in lines[:8]))
        )
        if index < 12 and is_top_banner_line:
            continue
        cleaned.append(line)
    text = "\n".join(cleaned).lstrip()
    banner = f"**Paper release:** `{release_tag}`\n**Released:** {release_date}\n\n"
    return banner + text


def strip_trailing_whitespace(text: str) -> str:
    trailing_newline = "\n" if text.endswith("\n") else ""
    return "\n".join(line.rstrip() for line in text.splitlines()) + trailing_newline


def export_one(src: Path, dest: Path, pandoc_bin: str, release_tag: str, release_date: str) -> None:
    export_src = MARKDOWN_SOURCE_OVERRIDES.get(src, src)
    if export_src.is_relative_to(PAPER_DIR):
        pandoc_cwd = PAPER_DIR
        pandoc_input = str(export_src.relative_to(PAPER_DIR))
    else:
        pandoc_cwd = export_src.parent
        pandoc_input = export_src.name
    subprocess.run(
        [
            pandoc_bin,
            "-f",
            "latex",
            "-t",
            "gfm",
            "--wrap=none",
            pandoc_input,
            "-o",
            str(dest),
        ],
        check=True,
        cwd=pandoc_cwd,
    )
    text = postprocess_markdown(dest.read_text(encoding="utf-8"))
    text = ensure_release_banner(text, release_tag, release_date)
    dest.write_text(strip_trailing_whitespace(text), encoding="utf-8")
    if not dest.read_text(encoding="utf-8").strip():
        raise SystemExit(f"empty markdown export for {src}")


def current_release_metadata() -> tuple[str, str]:
    release_info = (PAPER_DIR / "release_info.tex").read_text(encoding="utf-8")
    id_match = re.search(r"\\newcommand\{\\OPHPaperReleaseID\}\{([^}]+)\}", release_info)
    date_match = re.search(r"\\newcommand\{\\OPHPaperReleaseDate\}\{([^}]+)\}", release_info)
    if not id_match:
        raise SystemExit("Could not read OPHPaperReleaseID from paper/release_info.tex")
    if not date_match:
        raise SystemExit("Could not read OPHPaperReleaseDate from paper/release_info.tex")
    return id_match.group(1), date_match.group(1)


def write_build_info(out_dir: Path, generated: list[str], release_tag: str) -> None:
    payload = {
        "release_tag": release_tag,
        "source_snapshot": "reverse-engineering-reality/paper, reverse-engineering-reality/extra, and released reverse-engineering-reality/cosmology papers",
        "generated_files": generated,
    }
    (out_dir / BUILD_INFO_NAME).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def resolve_source(name_or_path: str) -> Path:
    candidate = Path(name_or_path)
    if candidate.suffix == ".tex":
        if not candidate.is_absolute():
            candidate = (Path.cwd() / candidate).resolve()
        if candidate.is_file():
            return candidate

    basename = candidate.stem if candidate.suffix else name_or_path
    for directory in (PAPER_DIR, EXTRA_DIR, COSMOLOGY_DIR):
        source = directory / f"{basename}.tex"
        if source.is_file():
            return source
    raise SystemExit(f"missing paper source for {name_or_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export canonical TeX papers to Markdown.")
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT),
        help="Directory where markdown exports should be written.",
    )
    parser.add_argument(
        "--paper",
        action="append",
        default=[],
        help="Paper basename without .tex. Repeat to override the default set.",
    )
    parser.add_argument(
        "--pandoc",
        default=shutil.which("pandoc") or "/opt/homebrew/bin/pandoc",
        help="Pandoc binary to use.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    sources = [resolve_source(paper) for paper in args.paper] if args.paper else list(DEFAULT_SOURCES)
    pandoc_bin = args.pandoc

    if shutil.which(pandoc_bin) is None and not Path(pandoc_bin).exists():
        raise SystemExit(f"pandoc not found: {pandoc_bin}")

    release_tag, release_date = current_release_metadata()
    generated: list[str] = []
    for src in sources:
        dest = out_dir / f"{src.stem}.md"
        export_one(src, dest, pandoc_bin, release_tag, release_date)
        generated.append(dest.name)
        print(dest)

    write_build_info(out_dir, generated, release_tag)
    print(out_dir / BUILD_INFO_NAME)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
