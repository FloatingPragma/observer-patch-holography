#!/usr/bin/env python3
"""Rebuild the paper stack and its release manifest in one deterministic pass.

The recurring CI failure mode this tool removes: a paper PDF is rebuilt and
committed while ``paper/paper_release_manifest.json`` records different
hashes, so ``tools/validate_paper_release_manifest.py`` rejects the tree on
the next push. Building and manifest regeneration are one command here, in
the order the release pipeline requires:

  1. ``tools/build_tex_papers.py``            (tectonic, all registered papers)
  2. ``paper/tools/check_build_warnings.py``  (complete warning gate)
  3. ``tools/generate_paper_release_manifest.py --preview``
  4. ``tools/validate_paper_release_manifest.py``

Release bumping is deliberately *not* part of this tool. Bumping the release
identifier (``tools/bump_paper_release.py``) is a separate editorial decision;
the default mode refreshes review artifacts for the identifier already
recorded in ``paper/release_info.tex``. After review and a release bump,
``--publication`` applies the strict new-release checks.

Usage:
  python3 tools/refresh_paper_release.py               # local review pass
  python3 tools/refresh_paper_release.py --publication # strict release pass
  python3 tools/refresh_paper_release.py --no-gate     # skip warning gate
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))
import build_tex_papers as paper_sources  # noqa: E402


def run(description: str, argv: list[str]) -> None:
    print(f"==> {description}")
    result = subprocess.run(argv, cwd=REPO_ROOT)
    if result.returncode != 0:
        print(f"refresh_paper_release: step failed: {description}", file=sys.stderr)
        raise SystemExit(result.returncode)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--no-gate",
        action="store_true",
        help="skip the build-warnings gate (use only while iterating on TeX)",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--preview",
        action="store_true",
        help="prepare same-release PDFs and a manifest for local review (default)",
    )
    mode.add_argument(
        "--publication",
        action="store_true",
        help=(
            "prepare a publication candidate after a release bump; reject an "
            "existing local or remote tag with the selected release ID"
        ),
    )
    return parser.parse_args(argv)


def manifest_commands(
    python: str,
    *,
    publication: bool,
) -> tuple[list[str], list[str]]:
    """Return the generator and validator commands for one build mode."""

    generator_command = [python, "tools/generate_paper_release_manifest.py"]
    validator_command = [python, "tools/validate_paper_release_manifest.py"]
    if publication:
        validator_command.append("--publication")
    else:
        generator_command.append("--preview")
    return generator_command, validator_command


def main() -> int:
    args = parse_args()
    python = sys.executable

    run("build all registered papers", [python, "tools/build_tex_papers.py"])

    if not args.no_gate:
        logs = [
            str(tex_path.with_suffix(".log"))
            for tex_path in paper_sources.ALL_PAPERS.values()
        ]
        missing_logs = [path for path in logs if not Path(path).is_file()]
        if missing_logs:
            print(
                "refresh_paper_release: the current build did not emit every "
                "source-derived log:\n  " + "\n  ".join(missing_logs),
                file=sys.stderr,
            )
            return 1
        run(
            "check build warnings against the allowlist",
            [python, "paper/tools/check_build_warnings.py", *logs],
        )

    generator_command, validator_command = manifest_commands(
        python,
        publication=args.publication,
    )

    run("regenerate the release manifest", generator_command)
    run(
        "validate the release manifest",
        validator_command,
    )

    mode = "publication candidate" if args.publication else "local review preview"
    print(
        "refresh_paper_release: papers, warnings gate, and manifest are "
        f"consistent ({mode})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
