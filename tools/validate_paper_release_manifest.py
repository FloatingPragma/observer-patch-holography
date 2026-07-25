#!/usr/bin/env python3
"""Validate paper/paper_release_manifest.json against the source paper set.

Issue #514: release-manifest membership must be *derived from the source set*,
not a fixed entry count. This check imports the single source of truth for which
papers exist and which are release-tracked (``tools/build_tex_papers.py``) and
rejects the manifest when membership drifts:

  * papers            = RELEASE_TRACKED           (the release-tracked core set)
  * supplemental_papers = PAPERS - RELEASE_TRACKED (built but not release-tracked)
  * extra_papers      = EXTRA_PAPERS              (discovered from extra/*.tex)

For every section the manifest key set and paper-to-PDF mapping must equal the
derived source mapping (no missing paper, unexpected paper, or cross-paper
artifact substitution), and every listed PDF artifact must exist in the
checkout. Any mismatch exits non-zero. No fixed counts are hard-coded here; add
or remove a paper in build_tex_papers.py (or an extra/*.tex file) and the
expected mapping moves with it.

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


def _pdf_relative(tex_path: Path) -> str:
    return tex_path.with_suffix(".pdf").relative_to(REPO_ROOT).as_posix()


def expected_sections() -> dict[str, dict[str, str]]:
    release = set(source.RELEASE_TRACKED)
    return {
        "papers": {
            paper_id: _pdf_relative(source.PAPERS[paper_id])
            for paper_id in sorted(release)
        },
        "supplemental_papers": {
            paper_id: _pdf_relative(source.PAPERS[paper_id])
            for paper_id in sorted(set(source.PAPERS) - release)
        },
        "extra_papers": {
            paper_id: _pdf_relative(tex_path)
            for paper_id, tex_path in sorted(source.EXTRA_PAPERS.items())
        },
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_section(
    name: str,
    expected: dict[str, str],
    section: dict,
    problems: list[str],
) -> None:
    got = set(section)
    expected_ids = set(expected)
    for missing in sorted(expected_ids - got):
        problems.append(f"{name}: expected paper '{missing}' is missing from the manifest")
    for unexpected in sorted(got - expected_ids):
        problems.append(f"{name}: manifest lists unexpected paper '{unexpected}' (not in the source set)")
    for paper_id, payload in section.items():
        if not isinstance(payload, dict):
            problems.append(f"{name}: '{paper_id}' payload must be an object")
            continue
        pdf_rel = payload.get("pdf_path")
        if not pdf_rel:
            problems.append(f"{name}: '{paper_id}' has no pdf_path in the manifest")
            continue
        expected_pdf_rel = expected.get(paper_id)
        if expected_pdf_rel is not None and pdf_rel != expected_pdf_rel:
            problems.append(
                f"{name}: '{paper_id}' maps to {pdf_rel}, but its source "
                f"derives {expected_pdf_rel}"
            )
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


def check_release_surface(manifest: dict, problems: list[str]) -> None:
    """Require every release-surface PDF to have one registered source route."""
    expected = {
        str(Path(payload["pdf_path"]))
        for section_name in (
            "papers",
            "supplemental_papers",
            "extra_papers",
        )
        for payload in (
            (manifest.get(section_name, {}) or {}).values()
            if isinstance(manifest.get(section_name, {}) or {}, dict)
            else ()
        )
        if isinstance(payload, dict) and payload.get("pdf_path")
    }
    for pdf_rel, source_rel in source.NON_TEX_SOURCE_PDFS.items():
        if not (REPO_ROOT / source_rel).is_file():
            problems.append(
                f"registered source for {pdf_rel} is missing: {source_rel}"
            )
        if not (REPO_ROOT / pdf_rel).is_file():
            problems.append(f"registered non-TeX output is missing: {pdf_rel}")
        expected.add(str(pdf_rel))
    actual = {
        str(pdf.relative_to(REPO_ROOT))
        for directory in ("paper", "extra")
        for pdf in (REPO_ROOT / directory).glob("*.pdf")
    }
    for stray in sorted(actual - expected):
        problems.append(
            f"stray PDF is not implied by a registered source: {stray}"
        )


def validate(manifest_path: Path) -> list[str]:
    problems: list[str] = []
    if not manifest_path.exists():
        return [f"manifest not found: {manifest_path}"]
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return [f"manifest is not valid JSON: {exc}"]
    if not isinstance(manifest, dict):
        return ["manifest root must be an object"]
    for name, expected in expected_sections().items():
        section = manifest.get(name, {}) or {}
        if not isinstance(section, dict):
            problems.append(f"{name}: manifest section must be an object")
            continue
        check_section(name, expected, section, problems)
    check_release_surface(manifest, problems)
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
