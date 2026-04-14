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
DEFAULT_OUT = WORKSPACE_ROOT / "temp" / "markdown"
DEFAULT_PAPERS = [
    "deriving_the_particle_zoo_from_observer_consistency",
    "observers_are_all_you_need",
    "reality_as_consensus_protocol",
    "recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact",
    "screen_microphysics_and_observer_synchronization",
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
        "source_snapshot": "reverse-engineering-reality/paper",
        "generated_files": generated,
    }
    (out_dir / BUILD_INFO_NAME).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


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
    papers = args.paper or list(DEFAULT_PAPERS)
    pandoc_bin = args.pandoc

    if shutil.which(pandoc_bin) is None and not Path(pandoc_bin).exists():
        raise SystemExit(f"pandoc not found: {pandoc_bin}")

    generated: list[str] = []
    for paper_name in papers:
        src = PAPER_DIR / f"{paper_name}.tex"
        if not src.exists():
            raise SystemExit(f"missing paper source: {src}")
        dest = out_dir / f"{paper_name}.md"
        export_one(src, dest, pandoc_bin)
        generated.append(dest.name)
        print(dest)

    write_build_info(out_dir, generated, current_release_id())
    print(out_dir / BUILD_INFO_NAME)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
