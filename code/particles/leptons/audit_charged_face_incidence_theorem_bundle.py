#!/usr/bin/env python3
"""Audit the submitted conditional face-incidence theorem bundle.

The archive contains charged comparison data and is therefore compare-only.
This module never extracts or executes archive code.  It verifies member
safety, CRC, internal hashes, and agreement between the submitted conditional
receipt and the canonical reference-free-at-runtime conditional theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import stat
from pathlib import Path, PurePosixPath
from typing import Any
from zipfile import ZipFile, ZipInfo

import mpmath as mp


mp.mp.dps = 100

HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
REPO = HERE.parents[3]
WORKSPACE = HERE.parents[4]
DEFAULT_ARCHIVE = WORKSPACE / "correspondence" / "oph_charged_face_incidence_theorem.zip"
CANONICAL = (
    CODE_ROOT / "particles" / "runs" / "leptons"
    / "charged_face_incidence_conditional_theorem.json"
)
DEFAULT_OUT = (
    CODE_ROOT / "particles" / "runs" / "leptons"
    / "charged_face_incidence_conditional_theorem_review.json"
)

PREFIX = "oph_charged_face_incidence_theorem/"
EXPECTED_ARCHIVE_SHA256 = "cc550b311c5f2bc416d78f65b4ec9e1bf32f10f9b73fc4135fd05de626b48d4f"
EXPECTED_NARRATIVE_SHA256 = "69a2c494fdc06034d76915fa9b045e604a47aa2be77247f24a3cb2505e990bde"
EXPECTED_RECEIPT_SHA256 = "28e51c272ef322b0d4c2712e2ac651d5847a941f239fc3dd43544bd31d6c2316"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_member(info: ZipInfo) -> bool:
    path = PurePosixPath(info.filename)
    mode = info.external_attr >> 16
    return (
        not path.is_absolute()
        and ".." not in path.parts
        and not stat.S_ISLNK(mode)
    )


def parse_hashes(raw: bytes) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in raw.decode("utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        rows[name.strip()] = digest
    return rows


def maximum_abs_difference(left: list[str], right: list[str]) -> mp.mpf:
    return max(
        abs(mp.mpf(a) - mp.mpf(b))
        for a, b in zip(left, right, strict=True)
    )


def build_review(archive_path: Path, canonical_path: Path) -> dict[str, Any]:
    archive_raw = archive_path.read_bytes()
    canonical_raw = canonical_path.read_bytes()
    canonical = json.loads(canonical_raw)

    with ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        first_bad_crc = archive.testzip()
        hashes_raw = archive.read(PREFIX + "SHA256SUMS")
        declared_hashes = parse_hashes(hashes_raw)
        submitted_receipt_raw = archive.read(PREFIX + "conditional_theorem_receipt.json")
        submitted = json.loads(submitted_receipt_raw)
        hash_checks = {
            name: sha256(archive.read(PREFIX + name)) == digest
            for name, digest in declared_hashes.items()
        }

    expected_names = {
        PREFIX,
        PREFIX + "verify_conditional_theorem.py",
        PREFIX + "conditional_theorem_receipt.json",
        PREFIX + "PROOF.md",
        PREFIX + "README.md",
        PREFIX + "SHA256SUMS",
    }
    submitted_fixed = submitted["repair_map"]["unique_stable_fixed_point"]
    canonical_fixed = canonical["repair_map"]["unique_iteratively_stable_fixed_point"]
    submitted_masses = list(submitted["mass_readout"]["masses_mev"].values())
    canonical_masses = canonical["conditional_spectrum"]["masses_mev"]
    submitted_countermodels = submitted["non_entailment_countermodels"]
    canonical_countermodels = canonical["scoped_source_normalization_non_rigidity"]["witnesses"]

    numerical_agreement = {
        "contraction_absolute_difference": str(
            abs(
                mp.mpf(submitted["repair_map"]["contraction_constant"])
                - mp.mpf(canonical["repair_map"]["contraction_constant"])
            )
        ),
        "fixed_point_max_absolute_difference": str(
            maximum_abs_difference(submitted_fixed, canonical_fixed)
        ),
        "mass_max_absolute_difference_mev": str(
            maximum_abs_difference(submitted_masses, canonical_masses)
        ),
        "countermodel_mass_max_absolute_differences_mev": [
            str(maximum_abs_difference(left["masses_mev"], right["masses_mev"]))
            for left, right in zip(submitted_countermodels, canonical_countermodels, strict=True)
        ],
    }

    return {
        "artifact": "oph_charged_face_incidence_conditional_theorem_bundle_review",
        "status": "FROZEN_CONDITIONAL_SUPPLEMENT_VERIFIED_SOURCE_ENTAILMENT_OPEN",
        "compare_only": True,
        "runtime_archive_contains_charged_references": True,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "public_prediction_promotion_allowed": False,
        "runtime_dependency": {"mpmath": mp.__version__},
        "archive": {
            "path": (
                archive_path.relative_to(WORKSPACE).as_posix()
                if archive_path.is_relative_to(WORKSPACE)
                else archive_path.name
            ),
            "sha256": sha256(archive_raw),
            "expected_sha256": EXPECTED_ARCHIVE_SHA256,
            "archive_hash_matches": sha256(archive_raw) == EXPECTED_ARCHIVE_SHA256,
            "narrative_attachment_sha256": EXPECTED_NARRATIVE_SHA256,
            "safe_non_symlink_members": all(safe_member(info) for info in infos),
            "no_duplicate_member_names": len(names) == len(set(names)),
            "exact_member_set": set(names) == expected_names,
            "crc_check_pass": first_bad_crc is None,
            "crc_first_bad_member": first_bad_crc,
            "declared_hash_checks": hash_checks,
            "declared_hash_checks_pass": all(hash_checks.values()),
            "submitted_receipt_sha256": sha256(submitted_receipt_raw),
            "submitted_receipt_hash_matches": (
                sha256(submitted_receipt_raw) == EXPECTED_RECEIPT_SHA256
            ),
        },
        "canonical_conditional_receipt": {
            "path": canonical_path.relative_to(REPO).as_posix(),
            "sha256": sha256(canonical_raw),
            "status": canonical["status"],
            "checks_pass": canonical["checks_pass"],
        },
        "independent_numerical_agreement": numerical_agreement,
        "proved_conditionally": [
            "strict contraction and one globally iteration-stable fixed point for the declared affine map",
            "determinant-normalized circulant readout under the declared physical hypotheses",
            "three explicit source-normalization witnesses with the same Jacobian and contraction but different masses",
        ],
        "qualification": [
            "the submitted tuple remains branch-incoherent",
            "the repair sources and feedback are historically target-informed and not OPH-derived",
            "bare exact block balance and the nonzero endpoint amplitude correction require a staged physical bridge",
            "the countermodels prove scoped candidate-property non-rigidity, not non-entailment from every OPH axiom",
            "the submitted receipt combines source calculation and charged inverse comparison in one artifact",
            "physical family, phase, determinant, refinement, coherent-branch, and pole-scheme attachments remain open",
        ],
        "integration_decision": (
            "Retain as a reproducible conditional proof supplement. Mark fixed-point "
            "uniqueness and iterative stability closed after the map is declared; do not "
            "promote the map, mass values, or source entailment."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--canonical", type=Path, default=CANONICAL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    review = build_review(args.archive.resolve(), args.canonical.resolve())
    output = (json.dumps(review, indent=2, sort_keys=True) + "\n").encode()

    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == output
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(output)
    print(json.dumps({"status": review["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
