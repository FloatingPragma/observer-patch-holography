#!/usr/bin/env python3
"""Audit the submitted charged source-law rigidity archive as compare-only.

Archive code is never extracted or executed.  The audit verifies safe members,
CRC, internal hashes, the frozen receipt digest, and agreement with the
canonical conditional CFQ theorem and the prior declared-map receipt.  Because
the archive is historically target-informed and includes a downstream mass
readout, its output is forbidden as a candidate-builder ancestor.
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
DEFAULT_ARCHIVE = WORKSPACE / "correspondence" / "oph_charged_source_law_rigidity_v1.zip"
CANONICAL_CFQ = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_source_law_rigidity_conditional.json"
)
DECLARED_MAP = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_face_incidence_conditional_theorem.json"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_source_law_rigidity_v1_review.json"
)

PREFIX = "oph_charged_source_law_rigidity/"
EXPECTED_ARCHIVE_SHA256 = (
    "792dad5c1524984c52d7582e751e1123099978164def80fad34310a34e314dff"
)
EXPECTED_NARRATIVE_SHA256 = (
    "879f5afbdb704d8894fe826544e26dff90651040df45954485be422ea2c0779e"
)
EXPECTED_RECEIPT_SHA256 = (
    "d0acb1242714a38ea7ccc3b9f335bc9accf8e8bbdd9d6959d7e7e5dcf93be4d3"
)
EXPECTED_INTERNAL_HASHES = {
    "README.md": "3f81da1d7554cb6120447522ccbf5c41b0a3e29bcf76fb3ab9d2acf8317f2705",
    "SOURCE_LAW_RIGIDITY_THEOREM.md": "29e68d510033801d608c99d1b5ce92b2d36ca89ac65217b89182cac8e461e2f7",
    "source_law_rigidity_insertion.tex": "6cc5cf623e24882d693b6847f1c574f1663276f6c4b01cbe84f7fba5dfa171ab",
    "DEPENDENCY_AUDIT.md": "fd2ed35a651675d0662567e76796c1ec2bfa3efa2faa49dfd82e051d0f6be9cf",
    "cfq_carrier_certificate_template.json": "b89b6fed6fce08fa2392943f757f5c85ec9d4056bb35a329f046f044d7f5f6d6",
    "verify_source_law_rigidity.py": "56760d749feaf425b4a132337e402a54fd75fadc3a7a60d10fc938cd81633358",
    "source_law_rigidity_receipt.json": EXPECTED_RECEIPT_SHA256,
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_member(info: ZipInfo) -> bool:
    path = PurePosixPath(info.filename)
    mode = info.external_attr >> 16
    return (
        not path.is_absolute()
        and ".." not in path.parts
        and not stat.S_ISLNK(mode)
        and not (info.flag_bits & 0x1)
    )


def parse_hashes(raw: bytes) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for line in raw.decode("utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        hashes[name.strip()] = digest
    return hashes


def maximum_abs_difference(left: list[str], right: list[str]) -> mp.mpf:
    return max(
        abs(mp.mpf(left_value) - mp.mpf(right_value))
        for left_value, right_value in zip(left, right, strict=True)
    )


def build_review(
    archive_path: Path,
    canonical_cfq_path: Path,
    declared_map_path: Path,
) -> dict[str, Any]:
    archive_raw = archive_path.read_bytes()
    canonical_raw = canonical_cfq_path.read_bytes()
    declared_raw = declared_map_path.read_bytes()
    canonical = json.loads(canonical_raw)
    declared = json.loads(declared_raw)

    with ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        first_bad_crc = archive.testzip()
        submitted_hashes = parse_hashes(archive.read(PREFIX + "SHA256SUMS"))
        submitted_receipt_raw = archive.read(PREFIX + "source_law_rigidity_receipt.json")
        submitted = json.loads(submitted_receipt_raw)
        submitted_template = json.loads(
            archive.read(PREFIX + "cfq_carrier_certificate_template.json")
        )
        internal_hash_checks = {
            name: sha256(archive.read(PREFIX + name)) == digest
            for name, digest in submitted_hashes.items()
        }

    expected_names = {
        PREFIX,
        *(PREFIX + name for name in EXPECTED_INTERNAL_HASHES),
        PREFIX + "SHA256SUMS",
    }
    canonical_registers = canonical["conditional_cfq_packet"]["registers"]
    submitted_registers = submitted["normalized_trace_theorem"]["registers"]
    register_agreement = {
        name: (
            submitted_registers[name]["dimension"] == canonical_row["dimension"]
            and submitted_registers[name]["event_rank"] == canonical_row["event_rank"]
            and submitted_registers[name]["trace_weight"]
            == canonical_row["normalized_trace_weight"]
        )
        for name, canonical_row in canonical_registers.items()
    }

    canonical_paths = {
        row["name"]: row for row in canonical["conditional_cfq_packet"]["paths"]
    }
    submitted_paths = {
        row["name"]: row for row in submitted["primitive_path_catalogue"]["paths"]
    }
    path_agreement = {
        name: all(
            submitted_paths[name][field] == canonical_row[field]
            for field in (
                "block",
                "monomial",
                "sign",
                "registers",
                "multiplicity",
                "exact_trace_coefficient",
            )
        )
        for name, canonical_row in canonical_paths.items()
    }

    canonical_graphs = canonical["declared_finite_graph_models"]["graphs"]
    submitted_graphs = submitted["transition_graph_certificate"]["graphs"]
    graph_agreement = {
        name: all(
            submitted_graphs[name][field] == canonical_row[field]
            for field in ("nodes", "edges", "connected", "expected_nodes", "dimension_matches")
        )
        for name, canonical_row in canonical_graphs.items()
    }

    submitted_coefficients = {
        **submitted["repair_map"]["sources"],
        **submitted["repair_map"]["feedback"],
    }
    canonical_coefficients = canonical["declared_map_agreement"]["coefficient_values"]
    coefficient_differences = {
        name: str(mp.mpf(submitted_coefficients[name]) - mp.mpf(value))
        for name, value in canonical_coefficients.items()
    }

    submitted_masses = list(submitted["downstream_conditional_corollary"]["masses_mev"].values())
    declared_masses = declared["conditional_spectrum"]["masses_mev"]
    submitted_fixed = submitted["repair_map"]["unique_fixed_point"]
    declared_fixed = declared["repair_map"]["unique_iteratively_stable_fixed_point"]
    submitted_gates = submitted["physical_promotion_gate"]
    submitted_cfq_gates = {
        name: value for name, value in submitted_gates.items() if name.startswith("CFQ")
    }
    template_gates = submitted_template["gates"]

    archive_checks = {
        "archive_hash_matches": sha256(archive_raw) == EXPECTED_ARCHIVE_SHA256,
        "safe_non_symlink_unencrypted_members": all(safe_member(info) for info in infos),
        "no_duplicate_member_names": len(names) == len(set(names)),
        "exact_member_set": set(names) == expected_names,
        "crc_check_pass": first_bad_crc is None,
        "declared_hash_manifest_exact": submitted_hashes == EXPECTED_INTERNAL_HASHES,
        "all_declared_internal_hashes_match": all(internal_hash_checks.values()),
        "submitted_receipt_hash_matches": (
            sha256(submitted_receipt_raw) == EXPECTED_RECEIPT_SHA256
        ),
    }
    agreement_checks = {
        "all_registers_agree": all(register_agreement.values()),
        "all_paths_agree": all(path_agreement.values()),
        "all_graph_counts_agree": all(graph_agreement.values()),
        "all_scalar_coefficients_agree_to_emitted_precision": max(
            abs(mp.mpf(value)) for value in coefficient_differences.values()
        )
        < mp.mpf("1e-75"),
        "fixed_point_agrees_exactly_to_emitted_precision": maximum_abs_difference(
            submitted_fixed, declared_fixed
        )
        == 0,
        "conditional_mass_readout_agrees_exactly_to_emitted_precision": (
            maximum_abs_difference(submitted_masses, declared_masses) == 0
        ),
        "submitted_receipt_has_all_eight_physical_gates_false": (
            set(submitted_cfq_gates) == {f"CFQ{index}" for index in range(1, 9)}
            and not any(submitted_cfq_gates.values())
        ),
        "submitted_template_has_all_eight_physical_gates_false": (
            len(template_gates) == 8 and not any(template_gates.values())
        ),
        "canonical_promotion_boundary_retained": (
            canonical["historical_charged_target_informed"] is True
            and canonical["global_source_only"] is False
            and canonical["branch_tuple_coherent"] is False
            and canonical["public_prediction_promotion_allowed"] is False
        ),
    }

    return {
        "artifact": "oph_charged_source_law_rigidity_v1_review",
        "schema_version": 1,
        "status": (
            "FROZEN_CONDITIONAL_CFQ_SUPPLEMENT_VERIFIED_"
            "PHYSICAL_SOURCE_LAW_AND_PROMOTION_OPEN"
        ),
        "compare_only": True,
        "forbidden_as_candidate_ancestor": True,
        "runtime_archive_contains_downstream_charged_mass_values": True,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
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
            "narrative_attachment_sha256": EXPECTED_NARRATIVE_SHA256,
            "crc_first_bad_member": first_bad_crc,
            "submitted_receipt_sha256": sha256(submitted_receipt_raw),
            "internal_hash_checks": internal_hash_checks,
            "checks": archive_checks,
            "checks_pass": all(archive_checks.values()),
            "execution_boundary": (
                "Archive code was not extracted or executed by this canonical auditor. "
                "An independent one-time audit rebuilt the frozen receipt byte-for-byte "
                "with mpmath 1.3.0 and networkx 3.6.1."
            ),
        },
        "canonical_conditional_receipt": {
            "path": canonical_cfq_path.relative_to(REPO).as_posix(),
            "sha256": sha256(canonical_raw),
            "status": canonical["status"],
            "checks_pass": canonical["checks_pass"],
        },
        "independent_agreement": {
            "registers": register_agreement,
            "paths": path_agreement,
            "graph_counts": graph_agreement,
            "coefficient_differences": coefficient_differences,
            "fixed_point_max_absolute_difference": str(
                maximum_abs_difference(submitted_fixed, declared_fixed)
            ),
            "conditional_mass_max_absolute_difference_mev": str(
                maximum_abs_difference(submitted_masses, declared_masses)
            ),
            "checks": agreement_checks,
            "checks_pass": all(agreement_checks.values()),
        },
        "valid_conditional_result": (
            "CFQ1-CFQ8 imply one quotient-visible scalar affine response law "
            "inside the declared packet."
        ),
        "not_established": [
            "current OPH axioms imply CFQ1-CFQ8",
            "the declared connected graphs are physical charged response modules",
            "OPH central completed records descend to the noncentral rank-one CFQ projectors",
            "the eight paths, signs, clock, multiplicities, and exhaustion follow from dynamics",
            "the scalar law uniquely fixes one microscopic repair instrument",
            "the one-bit, 2/9 phase, or 6^-14 determinant baselines are physically attached",
            "the hybrid numerical branch is coherent or the endpoint eigenvalues are pole masses",
            "historical no-target ancestry",
        ],
        "integration_decision": (
            "Retain as a reproducible conditional finite-certificate theorem and a "
            "concrete promotion contract. Do not mark OPH source-law rigidity, charged "
            "mass prediction, or any physical CFQ gate closed."
        ),
        "checks_pass": all(archive_checks.values()) and all(agreement_checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--canonical-cfq", type=Path, default=CANONICAL_CFQ)
    parser.add_argument("--declared-map", type=Path, default=DECLARED_MAP)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    review = build_review(
        args.archive.resolve(),
        args.canonical_cfq.resolve(),
        args.declared_map.resolve(),
    )
    encoded = (json.dumps(review, indent=2, sort_keys=True) + "\n").encode()

    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == encoded
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": review["status"],
                "checks_pass": review["checks_pass"],
                "public_prediction_promotion_allowed": review[
                    "public_prediction_promotion_allowed"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
