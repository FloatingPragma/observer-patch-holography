#!/usr/bin/env python3
"""Rebuild the paper stack and its release manifest in one deterministic pass.

The recurring CI failure mode this tool removes: a paper PDF is rebuilt and
committed while ``paper/paper_release_manifest.json`` still records the old
hashes, so ``tools/validate_paper_release_manifest.py`` rejects the tree on
the next push. Building and manifest regeneration are one command here, in
the order the release pipeline requires:

  1. ``tools/build_tex_papers.py``            (tectonic, all registered papers)
  2. ``paper/tools/check_build_warnings.py``  (overfull/underfull gate)
  3. ``tools/generate_paper_release_manifest.py --allow-same-release``
  4. ``tools/validate_paper_release_manifest.py``

Release bumping is deliberately *not* part of this tool. Bumping the release
identifier (``tools/bump_paper_release.py``) is a separate editorial decision;
this tool refreshes artifacts for the identifier already recorded in
``paper/release_info.tex``.

Usage:
  python3 tools/refresh_paper_release.py            # full pass
  python3 tools/refresh_paper_release.py --no-gate  # skip the warnings gate
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def run(description: str, argv: list[str]) -> None:
    print(f"==> {description}")
    result = subprocess.run(argv, cwd=REPO_ROOT)
    if result.returncode != 0:
        print(f"refresh_paper_release: step failed: {description}", file=sys.stderr)
        raise SystemExit(result.returncode)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--no-gate",
        action="store_true",
        help="skip the build-warnings gate (use only while iterating on TeX)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    python = sys.executable

    run("build all registered papers", [python, "tools/build_tex_papers.py"])

    if not args.no_gate:
        logs = sorted(str(p) for p in REPO_ROOT.glob("paper/*.log"))
        logs += sorted(str(p) for p in REPO_ROOT.glob("extra/*.log"))
        run(
            "check build warnings against the allowlist",
            [python, "paper/tools/check_build_warnings.py", *logs],
        )

    run(
        "regenerate the release manifest",
        [python, "tools/generate_paper_release_manifest.py", "--allow-same-release"],
    )
    run(
        "validate the release manifest",
        [python, "tools/validate_paper_release_manifest.py"],
    )

    print("refresh_paper_release: papers, warnings gate, and manifest are consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
