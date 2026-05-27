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
DEFAULT_OUT = WORKSPACE_ROOT / "temp" / "markdown"
DEFAULT_CORE_PAPERS = [
    PAPER_DIR / "deriving_the_particle_zoo_from_observer_consistency.tex",
    PAPER_DIR / "observers_are_all_you_need.tex",
    PAPER_DIR / "paradise_as_fixed_point_consensus.tex",
    PAPER_DIR / "reality_as_consensus_protocol.tex",
    PAPER_DIR / "recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex",
    PAPER_DIR / "screen_microphysics_and_observer_synchronization.tex",
]
DEFAULT_SUPPLEMENTAL_PAPERS = [
]
DEFAULT_EXTRA_PAPERS = sorted(EXTRA_DIR.glob("*.tex"))
DEFAULT_SOURCES = [
    *DEFAULT_CORE_PAPERS,
    *DEFAULT_SUPPLEMENTAL_PAPERS,
    *DEFAULT_EXTRA_PAPERS,
]
BUILD_INFO_NAME = "_build_info.json"


def postprocess_markdown(text: str) -> str:
    text = re.sub(
        r"(\*\*Paper release:\*\* `[^`]+`)(?=\*\*Released:\*\*)",
        r"\1  \n",
        text,
        count=1,
    )
    return text


def export_one(src: Path, dest: Path, pandoc_bin: str) -> None:
    subprocess.run(
        [
            pandoc_bin,
            "-f",
            "latex",
            "-t",
            "gfm",
            "--wrap=preserve",
            src.name,
            "-o",
            str(dest),
        ],
        check=True,
        cwd=src.parent,
    )
    dest.write_text(postprocess_markdown(dest.read_text(encoding="utf-8")), encoding="utf-8")


def current_release_id() -> str:
    release_info = (PAPER_DIR / "release_info.tex").read_text(encoding="utf-8")
    match = re.search(r"\\newcommand\{\\OPHPaperReleaseID\}\{([^}]+)\}", release_info)
    if not match:
        raise SystemExit("Could not read OPHPaperReleaseID from paper/release_info.tex")
    return match.group(1)


def write_build_info(out_dir: Path, generated: list[str], release_tag: str) -> None:
    payload = {
        "release_tag": release_tag,
        "source_snapshot": "reverse-engineering-reality/paper and reverse-engineering-reality/extra",
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
    for directory in (PAPER_DIR, EXTRA_DIR):
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

    generated: list[str] = []
    for src in sources:
        dest = out_dir / f"{src.stem}.md"
        export_one(src, dest, pandoc_bin)
        generated.append(dest.name)
        print(dest)

    write_build_info(out_dir, generated, current_release_id())
    print(out_dir / BUILD_INFO_NAME)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
