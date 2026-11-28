"""Build or verify the four local reviewer drafts and their offline headline archive."""

import argparse
import hashlib
import io
from pathlib import Path
import re
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.relative_to(ROOT)
PUBLIC = BASE.parent
STUDIES = Path("model/candidate_laws/interval_echo/openwave")
DUST = Path("model/candidate_laws/interval_echo_published_records/spherical_dust")
COLLAR = PUBLIC / "row_19/selected_collar"
BAO = Path("code/cosmology/paid_trail_bao")
ARCHIVE = BASE / "gravity_review.zip"
MANIFEST = BASE / "MANIFEST.sha256"

NOTES = {
    PUBLIC / "row_19/GRAVITY_METHOD.md",
    PUBLIC / "row_26/README.md",
    PUBLIC / "row_27/README.md",
    PUBLIC / "row_28/README.md",
}
SHARED = {
    BASE / name for name in ("README.md", "PR_DRAFTS.md", "CHECKLIST.md", "REVIEW.md")
}
HEADLINE_INPUTS = {
    COLLAR / "parked_energy_charge.json",
    COLLAR / "parked_pair_flight.json",
    DUST / "results.json",
    DUST / "clock_pairs/regions.json",
    STUDIES / "row_26_newton/radial_release.json",
    STUDIES / "row_27_metric/dust_exterior/results.json",
    *(
        BAO / ("desi_gaussian_bao_ALL_GCcomb_" + name + ".txt")
        for name in ("mean", "cov")
    ),
}
NEW_EVIDENCE = {
    COLLAR / "matched_charge_transport.json",
    STUDIES / "row_26_newton/current_loop.json",
    STUDIES / "row_27_metric/parked_signals/common_delay.json",
    STUDIES / "row_28_acceleration/cut_messages.json",
}
EXTRA_SOURCES = {
    *(p.with_suffix(".py") for p in NEW_EVIDENCE),
    COLLAR / "parked_energy_charge.py",
    COLLAR / "parked_pair_flight.py",
    COLLAR / "parked_pair_flight.json",
    STUDIES / "row_26_newton/radial_release.py",
    STUDIES / "row_26_newton/written_source.py",
    STUDIES / "row_26_newton/written_source.json",
    STUDIES / "row_27_metric/dust_exterior/measure.py",
    STUDIES / "row_27_metric/dust_exterior/forecast.py",
    STUDIES / "row_27_metric/parked_signals/measure.py",
    STUDIES / "row_28_acceleration/capacity_distances.py",
    STUDIES / "row_28_acceleration/collar_abundance.py",
    DUST / "experiment.py",
    DUST / "clock_pairs/regions.py",
    DUST / "clock_pairs/REFINEMENT.md",
    Path("Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean"),
    Path("Lean/ObserverPatchHolography/EinsteinBranch/SphericalDust.lean"),
    Path("code/matched_einstein_limit/results.json"),
    Path(
        "model/candidate_laws/interval_echo_partition_repair/published_strings/archive.py"
    ),
    Path("model/experiments/interval_echo_clock/publication.json"),
    Path("cosmology/oph_dark_matter_paper.tex"),
    Path("cosmology/oph_cosmological_vacuum_and_structure_formation.tex"),
    Path(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"
    ),
}


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def manifest(contents):
    return "".join(
        f"{digest(blob)}  {name}\n" for name, blob in sorted(contents.items())
    ).encode()


def linked_sources():
    paths = set()
    notes = NOTES | SHARED
    for note in notes:
        markdown = (ROOT / note).read_text()
        links = re.findall(
            r"\]\(([^)]+)\)|^\[[^\]]+\]:\s*(\S+)", markdown, re.MULTILINE
        )
        for inline, reference in links:
            dest = (inline or reference).split("#")[0]
            if not dest or "://" in dest:
                continue
            path = (ROOT / note.parent / dest).resolve().relative_to(ROOT)
            if path in (ARCHIVE, MANIFEST):
                continue
            if not (ROOT / path).exists():
                raise FileNotFoundError(f"{note} links to missing {path}")
            if (ROOT / path).is_file():
                paths.add(path)
    return paths


def source_bytes(paths):
    return {p.as_posix(): (ROOT / p).read_bytes() for p in paths}


def make_archive():
    packet = source_bytes(
        NOTES
        | SHARED
        | HEADLINE_INPUTS
        | NEW_EVIDENCE
        | {BASE / "headline.py", BASE / "controls.py"}
    )
    packet["START_HERE.md"] = b"""# OPH gravity and dark matter reviewer drafts

Read openwave/interval_echo/gravity/README.md for the shared argument and four method notes.
From this extracted directory, run:

    python3 -I -B openwave/interval_echo/gravity/headline.py
    python3 -I -B openwave/interval_echo/gravity/controls.py
    LC_ALL=C shasum -a 256 -c REVIEW_MANIFEST.sha256

Python 3.11+ is sufficient. These reductions use the archived observer and DESI records.
The method notes also link to full generators, papers and larger experiments in the
research checkout. Those larger runs need that checkout and its scientific environment.
All four rows are local drafts on HOLD. No OpenWave task or cell award is claimed.
"""
    packet["REVIEW_MANIFEST.sha256"] = manifest(packet)
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as archive:
        for name, blob in sorted(packet.items()):
            info = zipfile.ZipInfo(name, (2026, 11, 28, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(
                info, blob, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9
            )
    return data.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="regenerate the local packet"
    )
    args = parser.parse_args()
    paths = (
        linked_sources()
        | NOTES
        | SHARED
        | HEADLINE_INPUTS
        | NEW_EVIDENCE
        | EXTRA_SOURCES
    )
    files = source_bytes(
        paths | {BASE / "bundle.py", BASE / "headline.py", BASE / "controls.py"}
    )
    if args.write:
        archive = make_archive()
        (ROOT / ARCHIVE).write_bytes(archive)
        files[ARCHIVE.as_posix()] = archive
        (ROOT / MANIFEST).write_bytes(manifest(files))
    else:
        archive = (ROOT / ARCHIVE).read_bytes()
        files[ARCHIVE.as_posix()] = archive
        if (ROOT / MANIFEST).read_bytes() != manifest(files):
            raise ValueError(
                "A source, receipt or reviewer file differs from the packet manifest"
            )
        # Compare entry contents, permitting different zlib versions across systems.
        with (
            zipfile.ZipFile(io.BytesIO(archive)) as actual,
            zipfile.ZipFile(io.BytesIO(make_archive())) as expected,
        ):
            if actual.namelist() != expected.namelist() or any(
                actual.read(p) != expected.read(p) for p in expected.namelist()
            ):
                raise ValueError(
                    "Archive contents differ from the current reviewer draft"
                )
    print(
        f"{len(files)} working-tree files bound; reviewer archive {len(archive):,} bytes"
    )


if __name__ == "__main__":
    main()
