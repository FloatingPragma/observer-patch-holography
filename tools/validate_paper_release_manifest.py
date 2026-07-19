#!/usr/bin/env python3
"""Validate paper/paper_release_manifest.json against the source paper set.

Issue #514: release-manifest membership must be *derived from the source set*,
not a fixed entry count. This check imports the single source of truth for which
papers exist and which are release-tracked (``tools/build_tex_papers.py``) and
rejects the manifest when membership drifts:

  * papers            = RELEASE_TRACKED           (the release-tracked core set)
  * supplemental_papers = PAPERS - RELEASE_TRACKED (built but not release-tracked)
  * extra_papers      = EXTRA_PAPERS              (discovered from extra/*.tex)

For every section the manifest key set must equal the derived expected set (no
missing paper, no unexpected paper), and every listed PDF artifact must exist in
the checkout. Any mismatch exits non-zero. No fixed counts are hard-coded here;
add or remove a paper in build_tex_papers.py (or an extra/*.tex file) and the
expected set moves with it.

Usage:
  python3 tools/validate_paper_release_manifest.py
  python3 tools/validate_paper_release_manifest.py --manifest paper/paper_release_manifest.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = REPO_ROOT / "paper" / "paper_release_manifest.json"

# Import the source of truth. build_tex_papers only defines sets and globs
# extra/*.tex at import time; it has no import-time side effects beyond that.
sys.path.insert(0, str(REPO_ROOT / "tools"))
import build_tex_papers as source  # noqa: E402


def expected_sections() -> dict[str, set[str]]:
    release = set(source.RELEASE_TRACKED)
    supplemental = set(source.PAPERS) - release
    extra = set(source.EXTRA_PAPERS)
    return {
        "papers": release,
        "supplemental_papers": supplemental,
        "extra_papers": extra,
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_section(name: str, expected: set[str], section: dict, problems: list[str]) -> None:
    got = set(section)
    for missing in sorted(expected - got):
        problems.append(f"{name}: expected paper '{missing}' is missing from the manifest")
    for unexpected in sorted(got - expected):
        problems.append(f"{name}: manifest lists unexpected paper '{unexpected}' (not in the source set)")
    for paper_id, payload in section.items():
        pdf_rel = (payload or {}).get("pdf_path")
        if not pdf_rel:
            problems.append(f"{name}: '{paper_id}' has no pdf_path in the manifest")
            continue
        pdf_abs = REPO_ROOT / pdf_rel
        if not pdf_abs.exists():
            problems.append(f"{name}: artifact for '{paper_id}' is missing on disk: {pdf_rel}")
            continue
        # Content integrity: a listed artifact must match its declared sha256 + size_bytes,
        # so a silently-rebuilt / swapped / truncated PDF is rejected, not just an absent one.
        declared_sha = (payload or {}).get("sha256")
        if not declared_sha:
            problems.append(f"{name}: '{paper_id}' has no sha256 in the manifest")
        else:
            actual_sha = _sha256(pdf_abs)
            if actual_sha != declared_sha:
                problems.append(
                    f"{name}: '{paper_id}' sha256 mismatch for {pdf_rel}: "
                    f"manifest {declared_sha}, disk {actual_sha}"
                )
        declared_size = (payload or {}).get("size_bytes")
        if declared_size is None:
            problems.append(f"{name}: '{paper_id}' has no size_bytes in the manifest")
        else:
            actual_size = pdf_abs.stat().st_size
            try:
                declared_size_int = int(str(declared_size))
            except (TypeError, ValueError):
                problems.append(f"{name}: '{paper_id}' size_bytes is not an integer: {declared_size!r}")
            else:
                if declared_size_int != actual_size:
                    problems.append(
                        f"{name}: '{paper_id}' size_bytes mismatch for {pdf_rel}: "
                        f"manifest {declared_size_int}, disk {actual_size}"
                    )


def validate(manifest_path: Path) -> list[str]:
    problems: list[str] = []
    if not manifest_path.exists():
        return [f"manifest not found: {manifest_path}"]
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return [f"manifest is not valid JSON: {exc}"]
    for name, expected in expected_sections().items():
        check_section(name, expected, manifest.get(name, {}) or {}, problems)
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="path to the release manifest JSON")
    args = parser.parse_args()

    problems = validate(args.manifest)
    if problems:
        print("paper release manifest FAILED (derived from source set):")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    sections = expected_sections()
    counts = " + ".join(f"{len(v)} {k}" for k, v in sections.items())
    print(f"paper release manifest OK: {counts} (membership derived from build_tex_papers.py)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
