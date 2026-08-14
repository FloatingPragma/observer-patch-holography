#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path

from reproducible_build_env import build_environment


REPO_ROOT = Path(__file__).resolve().parent.parent
PAPER_DIR = REPO_ROOT / "paper"
EXTRA_DIR = REPO_ROOT / "extra"
FLAGSHIP_DIR = REPO_ROOT / "flagship"
COSMOLOGY_DIR = REPO_ROOT / "cosmology"

# Release-surface PDFs whose source is not a sibling TeX root. The manifest
# generator and validator share this registry so a source/output pair cannot
# disappear from one gate while remaining accepted by the other.
NON_TEX_SOURCE_PDFS = {}

PAPERS = {
    "from_observer_consensus_to_standard_physics": (
        FLAGSHIP_DIR / "from_observer_consensus_to_standard_physics.tex"
    ),
    "deriving_the_particle_zoo_from_observer_consistency": (
        PAPER_DIR / "deriving_the_particle_zoo_from_observer_consistency.tex"
    ),
    "observers_are_all_you_need": PAPER_DIR / "observers_are_all_you_need.tex",
    "paradise_as_fixed_point_consensus": PAPER_DIR / "paradise_as_fixed_point_consensus.tex",
    "reality_as_consensus_protocol": PAPER_DIR / "reality_as_consensus_protocol.tex",
    "recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency": (
        PAPER_DIR / "recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"
    ),
    "deriving_standard_model_gauge_structure_from_observer_overlap_consistency": (
        PAPER_DIR / "deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex"
    ),
    "screen_microphysics_and_observer_synchronization": PAPER_DIR / "screen_microphysics_and_observer_synchronization.tex",
}
EXTRA_PAPERS = {
    path.stem: path
    for path in sorted(EXTRA_DIR.glob("*.tex"))
}
COSMOLOGY_PAPERS = {
    path.stem: path
    for path in sorted(COSMOLOGY_DIR.glob("*.tex"))
}
# The compact proof remains manifest-bound for the challenge backend. No
# extra/ or cosmology/ PDF is attached to the public GitHub Release.
RELEASED_ADJUNCT_PAPERS = {
    "compact_proof_of_oph": EXTRA_PAPERS["compact_proof_of_oph"],
}
ALL_PAPERS = {**PAPERS, **EXTRA_PAPERS, **COSMOLOGY_PAPERS}

TECTONIC_FETCH_ATTEMPTS = 3
TECTONIC_TRANSIENT_FETCH_MARKERS = (
    "429 too many requests",
    "failed to retrieve",
    "connection reset by peer",
    "operation timed out",
)

RELEASE_TRACKED = (
    "from_observer_consensus_to_standard_physics",
    "recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency",
    "deriving_standard_model_gauge_structure_from_observer_overlap_consistency",
    "observers_are_all_you_need",
    "reality_as_consensus_protocol",
    "paradise_as_fixed_point_consensus",
    "screen_microphysics_and_observer_synchronization",
    "deriving_the_particle_zoo_from_observer_consistency",
)
RELEASE_TRACKED_SET = set(RELEASE_TRACKED)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the published OPH TeX papers with tectonic.",
    )
    parser.add_argument(
        "papers",
        nargs="*",
        help="Optional paper ids to build. Defaults to all known papers.",
    )
    parser.add_argument(
        "--release-only",
        action="store_true",
        help="Build only the release-tracked paper bundle used by paper_release_manifest.json.",
    )
    parser.add_argument(
        "--supplemental-only",
        action="store_true",
        help="Build only supplemental papers outside the release-tracked set.",
    )
    parser.add_argument(
        "--extra-only",
        action="store_true",
        help="Build only root-level extra/ papers.",
    )
    parser.add_argument(
        "--cosmology-only",
        action="store_true",
        help="Build only root-level cosmology/ papers.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print the known paper ids and exit.",
    )
    return parser.parse_args()


def resolve_targets(args: argparse.Namespace) -> list[str]:
    selected_modes = [
        args.release_only,
        args.supplemental_only,
        args.extra_only,
        args.cosmology_only,
    ]
    if sum(bool(mode) for mode in selected_modes) > 1:
        raise SystemExit(
            "choose at most one of --release-only, --supplemental-only, "
            "--extra-only, or --cosmology-only"
        )

    if args.list:
        for paper_id in sorted(ALL_PAPERS):
            if paper_id in RELEASE_TRACKED_SET:
                marker = "release"
            elif paper_id in EXTRA_PAPERS:
                marker = "extra"
            elif paper_id in COSMOLOGY_PAPERS:
                marker = "cosmology"
            else:
                marker = "supplemental"
            print(f"{paper_id}\t{marker}")
        raise SystemExit(0)

    if args.papers:
        unknown = [paper_id for paper_id in args.papers if paper_id not in ALL_PAPERS]
        if unknown:
            raise SystemExit(f"unknown paper ids: {', '.join(sorted(unknown))}")
        return args.papers

    if args.release_only:
        return list(RELEASE_TRACKED)
    if args.supplemental_only:
        return sorted(set(PAPERS) - RELEASE_TRACKED_SET)
    if args.extra_only:
        return sorted(EXTRA_PAPERS)
    if args.cosmology_only:
        return sorted(COSMOLOGY_PAPERS)
    return [
        *RELEASE_TRACKED,
        *sorted(set(ALL_PAPERS) - RELEASE_TRACKED_SET),
    ]


def build_one(paper_id: str) -> None:
    tex_path = ALL_PAPERS[paper_id]
    if not tex_path.exists():
        raise SystemExit(f"missing TeX source: {tex_path}")

    # Keep the log from this exact build.  The publication warning gate consumes
    # one log per discovered source; without --keep-logs an old ignored log can
    # make a clean-looking gate report describe a previous build instead.
    cmd = ["tectonic", "-X", "compile", tex_path.name, "--keep-logs"]
    for attempt in range(1, TECTONIC_FETCH_ATTEMPTS + 1):
        result = subprocess.run(
            cmd,
            cwd=tex_path.parent,
            text=True,
            capture_output=True,
            env=build_environment(REPO_ROOT),
        )
        if result.returncode == 0:
            break
        if attempt == TECTONIC_FETCH_ATTEMPTS or not transient_fetch_failure(result):
            break
        delay = 10 * attempt
        print(
            f"tectonic transient fetch failure for {paper_id}; "
            f"retrying in {delay}s ({attempt}/{TECTONIC_FETCH_ATTEMPTS})",
        )
        time.sleep(delay)
    if result.returncode != 0:
        if result.stdout.strip():
            print(result.stdout[-8000:])
        if result.stderr.strip():
            print(result.stderr[-8000:])
        raise SystemExit(f"tectonic failed for {paper_id}")

    print(tex_path.parent / f"{tex_path.stem}.pdf")


def transient_fetch_failure(result: subprocess.CompletedProcess[str]) -> bool:
    output = f"{result.stdout}\n{result.stderr}".lower()
    return any(marker in output for marker in TECTONIC_TRANSIENT_FETCH_MARKERS)


def main() -> int:
    if shutil.which("tectonic") is None:
        raise SystemExit("tectonic is required but was not found in PATH")

    args = parse_args()
    for paper_id in resolve_targets(args):
        build_one(paper_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
