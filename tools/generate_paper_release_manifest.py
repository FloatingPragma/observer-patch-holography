#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import build_tex_papers as paper_sources


RELEASE_INFO_RELATIVE = Path("paper/release_info.tex")
OUTPUT_RELATIVE = Path("paper/paper_release_manifest.json")
RELEASED_COSMOLOGY_TEX = ()
RELEASE_TRACKED_PDFS = {
    "deriving_the_particle_zoo_from_observer_consistency": Path(
        "paper/deriving_the_particle_zoo_from_observer_consistency.pdf"
    ),
    "observers_are_all_you_need": Path("paper/observers_are_all_you_need.pdf"),
    "paradise_as_fixed_point_consensus": Path("paper/paradise_as_fixed_point_consensus.pdf"),
    "reality_as_consensus_protocol": Path("paper/reality_as_consensus_protocol.pdf"),
    "recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency": Path(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf"
    ),
    "deriving_standard_model_gauge_structure_from_observer_overlap_consistency": Path(
        "paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf"
    ),
    "screen_microphysics_and_observer_synchronization": Path(
        "paper/screen_microphysics_and_observer_synchronization.pdf"
    ),
}
SUPPLEMENTAL_RELEASE_PDFS = {}

# Kept as a module-level alias for the focused regression tests. The canonical
# registry lives beside the source-derived TeX inventory.
NON_TEX_SOURCE_PDFS = paper_sources.NON_TEX_SOURCE_PDFS


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    verify_lean_theorem_count(repo_root)
    release_info = (repo_root / RELEASE_INFO_RELATIVE).read_text(encoding="utf-8")
    release_id = extract_macro(release_info, "OPHPaperReleaseID")
    release_date = extract_macro(release_info, "OPHPaperReleaseDate")

    manifest = {
        "release_id": release_id,
        "released_at": release_date,
        # Deterministic release metadata: a clean rebuild must not dirty the
        # manifest merely because it ran at a different wall-clock time.
        "generated_at": deterministic_generated_at(release_date),
        "papers": {},
        "supplemental_papers": {},
        "extra_papers": {},
    }
    fill_section(repo_root, manifest["papers"], RELEASE_TRACKED_PDFS)
    fill_section(repo_root, manifest["supplemental_papers"], SUPPLEMENTAL_RELEASE_PDFS)
    fill_section(repo_root, manifest["extra_papers"], discover_extra_pdfs(repo_root))

    output_path = repo_root / OUTPUT_RELATIVE
    previous_manifest = load_existing_manifest(output_path)
    tagged_manifest = load_tagged_manifest(repo_root, release_id)
    enforce_tag_immutability(tagged_manifest, manifest)
    enforce_release_bump(previous_manifest, manifest, args.allow_same_release)
    verify_no_stray_pdfs(repo_root, manifest)
    verify_pdf_release_lines(repo_root, manifest, args.skip_pdf_release_check)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output_path)
    return 0


def verify_lean_theorem_count(repo_root: Path) -> None:
    """Refuse to build the manifest while a stated Lean theorem-count floor drifts."""
    check = repo_root / "tools" / "check_lean_theorem_count.py"
    result = subprocess.run([sys.executable, str(check)], cwd=repo_root)
    if result.returncode != 0:
        raise SystemExit(
            "Lean theorem-count floor drifted; run "
            "`python3 tools/check_lean_theorem_count.py --fix` and re-run"
        )


def discover_extra_pdfs(repo_root: Path) -> dict[str, Path]:
    extra_dir = repo_root / "extra"
    discovered: dict[str, Path] = {}
    for tex_path in sorted(extra_dir.glob("*.tex")):
        pdf_path = tex_path.with_suffix(".pdf")
        if not pdf_path.is_file():
            raise SystemExit(
                f"missing PDF for extra paper {tex_path.name}: {pdf_path}. "
                "Run python3 tools/build_tex_papers.py --extra-only before regenerating the manifest."
            )
        discovered[tex_path.stem] = pdf_path.relative_to(repo_root)
    return discovered


def verify_no_stray_pdfs(repo_root: Path, manifest: dict) -> None:
    """Reject release-surface PDFs that no source implies.

    Membership is derived from the sources: the curated release set for
    ``paper/`` and the ``extra/*.tex`` glob for ``extra/``. Any other PDF in
    those directories would ship unhashed and unaudited, so its presence
    fails the manifest build (issue #514).
    """
    expected = {
        str(Path(payload["pdf_path"]))
        for payload in manifest_pdf_entries(manifest).values()
    }
    for pdf, source in NON_TEX_SOURCE_PDFS.items():
        if not (repo_root / source).is_file():
            raise SystemExit(
                f"registered source for {pdf} is missing: {source}. "
                "Fix the NON_TEX_SOURCE_PDFS registry."
            )
        if not (repo_root / pdf).is_file():
            raise SystemExit(
                f"registered non-TeX output is missing: {pdf}. "
                f"Rebuild it from {source} before regenerating the manifest."
            )
        expected.add(str(pdf))
    actual = {
        str(pdf.relative_to(repo_root))
        for directory in ("paper", "extra")
        for pdf in (repo_root / directory).glob("*.pdf")
    }
    strays = sorted(actual - expected)
    if strays:
        raise SystemExit(
            "stray PDFs on the release surface are not implied by any source: "
            f"{', '.join(strays)}. Remove them or register their source before "
            "regenerating the manifest."
        )


def fill_section(repo_root: Path, section: dict, pdfs: dict[str, Path]) -> None:
    for paper_id, relative_path in pdfs.items():
        pdf_path = repo_root / relative_path
        if not pdf_path.is_file():
            raise SystemExit(f"missing release PDF for {paper_id}: {pdf_path}")
        section[paper_id] = {
            "pdf_path": str(relative_path),
            "sha256": sha256(pdf_path),
            "size_bytes": pdf_path.stat().st_size,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write the current release manifest for the synced paper PDFs.",
    )
    parser.add_argument(
        "--allow-same-release",
        action="store_true",
        help=(
            "Allow PDF hash changes before the current release ID has been tagged. "
            "A tagged release is immutable and always requires a new release ID."
        ),
    )
    parser.add_argument(
        "--skip-pdf-release-check",
        action="store_true",
        help="Skip checking that each local PDF exposes the visible release line.",
    )
    return parser.parse_args()


def extract_macro(text: str, macro_name: str) -> str:
    pattern = re.compile(r"\\newcommand\{\\%s\}\{([^}]*)\}" % re.escape(macro_name))
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return match.group(1).strip()


def deterministic_generated_at(release_date: str) -> str:
    try:
        parsed = datetime.strptime(release_date, "%B %d, %Y")
    except ValueError as exc:
        raise SystemExit(
            f"paper release date must use 'Month D, YYYY': {release_date}"
        ) from exc
    return parsed.strftime("%Y-%m-%dT00:00:00Z")


def load_existing_manifest(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def load_tagged_manifest(repo_root: Path, release_id: str) -> dict | None:
    """Read the manifest committed at an existing local release tag.

    Same-release refreshes are useful while preparing a release. Once the
    release ID has a Git tag, that tag is the immutable local publication
    boundary and its PDF hashes cannot be replaced under the same ID.
    """

    tag_ref = f"refs/tags/{release_id}"
    exists = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", tag_ref],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    if exists.returncode != 0:
        return None

    tagged_path = f"{tag_ref}:{OUTPUT_RELATIVE.as_posix()}"
    result = subprocess.run(
        ["git", "show", tagged_path],
        cwd=repo_root,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or f"git show exited {result.returncode}"
        raise SystemExit(
            f"could not read the release manifest at {tag_ref}: {detail}"
        )
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"release manifest at {tag_ref} is not valid JSON: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"release manifest at {tag_ref} is not a JSON object")
    return payload


def enforce_tag_immutability(tagged_manifest: dict | None, manifest: dict) -> None:
    """Reject replacement of PDF bytes behind an existing release tag."""

    if tagged_manifest is None:
        return

    tagged_release_id = str(tagged_manifest.get("release_id", "")).strip()
    current_release_id = str(manifest.get("release_id", "")).strip()
    if tagged_release_id != current_release_id:
        return

    tagged_papers = manifest_pdf_entries(tagged_manifest)
    current_papers = manifest_pdf_entries(manifest)
    changed_papers = [
        paper_id
        for paper_id, payload in current_papers.items()
        if tagged_papers.get(paper_id, {}).get("sha256") != payload["sha256"]
    ]
    if not changed_papers:
        return

    changed_list = ", ".join(sorted(changed_papers))
    raise SystemExit(
        "PDF hashes differ from immutable Git tag "
        f"{current_release_id}: {changed_list}. "
        "Run python3 tools/bump_paper_release.py, rebuild every paper PDF, "
        "and regenerate the manifest. --allow-same-release cannot rewrite a "
        "tagged release."
    )


def enforce_release_bump(previous_manifest: dict | None, manifest: dict, allow_same_release: bool) -> None:
    if previous_manifest is None or allow_same_release:
        return

    previous_release_id = str(previous_manifest.get("release_id", "")).strip()
    current_release_id = str(manifest.get("release_id", "")).strip()
    if previous_release_id != current_release_id:
        return

    previous_papers = manifest_pdf_entries(previous_manifest)
    current_papers = manifest_pdf_entries(manifest)
    changed_papers = [
        paper_id
        for paper_id, payload in current_papers.items()
        if previous_papers.get(paper_id, {}).get("sha256") != payload["sha256"]
    ]
    if not changed_papers:
        return

    changed_list = ", ".join(sorted(changed_papers))
    raise SystemExit(
        "PDF hashes changed for the current release ID "
        f"{current_release_id}, but the release was not bumped first: {changed_list}. "
        "Run python3 tools/bump_paper_release.py, rebuild all current paper PDFs, and rerun this command. "
        "Use --allow-same-release only if the unchanged release ID is intentional."
    )


def verify_pdf_release_lines(repo_root: Path, manifest: dict, skip_pdf_release_check: bool) -> None:
    if skip_pdf_release_check:
        return

    release_id = str(manifest["release_id"]).strip()
    missing_release_line: list[str] = []
    tool_failures: list[str] = []
    for paper_id, payload in manifest_pdf_entries(manifest).items():
        pdf_path = repo_root / payload["pdf_path"]
        contains_release = pdf_contains_text(pdf_path, f"Paper release: {release_id}")
        if contains_release is True:
            continue
        if contains_release is False:
            missing_release_line.append(paper_id)
        else:
            tool_failures.append(f"{paper_id}: {contains_release}")

    if missing_release_line:
        raise SystemExit(
            "Local PDFs do not expose the current visible release line "
            f"{release_id}: {', '.join(sorted(missing_release_line))}. "
            "Rebuild all current paper PDFs after bumping paper/release_info.tex, then rerun this command. "
            "Every release bump must propagate to every paper PDF so the release IDs stay in sync."
        )

    if tool_failures:
        failures = "; ".join(tool_failures)
        raise SystemExit(
            "Could not verify the visible release line in the local PDFs: "
            f"{failures}. Install pdftotext or rerun with --skip-pdf-release-check."
        )


def pdf_contains_text(path: Path, needle: str) -> bool | str:
    if shutil.which("pdftotext") is None:
        return "pdftotext not installed"
    try:
        result = subprocess.run(
            ["pdftotext", str(path), "-"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
    except OSError as exc:
        return f"pdftotext failed: {exc}"

    if result.returncode != 0:
        stderr = result.stderr.strip() or f"exit code {result.returncode}"
        return f"pdftotext failed: {stderr}"
    return needle in result.stdout


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_pdf_entries(manifest: dict) -> dict:
    entries: dict = {}
    for section_name in ("papers", "supplemental_papers", "extra_papers"):
        entries.update(manifest.get(section_name, {}))
    return entries


if __name__ == "__main__":
    raise SystemExit(main())
