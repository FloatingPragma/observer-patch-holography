#!/usr/bin/env python3
"""Adversarial compare-only audit of the submitted face-incidence candidate.

The original archive is preserved under the workspace ``correspondence``
surface.  This audit does not execute archive code.  It verifies the archive
and manifest, independently recomputes the candidate, inverts the comparison
triple into the three spectral coordinates, and evaluates coherent public and
source-audit branch tuples.

Because charged reference masses are deliberately consumed, this module and
its output are compare-only and may never become candidate ancestors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from functools import lru_cache
from pathlib import Path, PurePosixPath
from typing import Any
from zipfile import ZipFile


HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
REPO = HERE.parents[3]
WORKSPACE = HERE.parents[4]
DEFAULT_ARCHIVE = WORKSPACE / "correspondence" / "oph_charged_icosahedral_completion_v1.zip"
DEFAULT_OUT = (
    CODE_ROOT / "particles" / "runs" / "leptons"
    / "charged_icosahedral_completion_v1_review.json"
)
PUBLIC_PIXEL = (
    CODE_ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_public_pixel_certificate.json"
)
SOURCE_PIXEL = (
    CODE_ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
D10_GAP = (
    CODE_ROOT / "particles" / "runs" / "calibration"
    / "d10_ew_common_transport_gap_report.json"
)
KRAWCZYK = (
    CODE_ROOT / "particles" / "hierarchy" / "certificates"
    / "R_U_krawczyk_certificate.json"
)
LEGACY_D10_SOURCE = (
    CODE_ROOT / "particles" / "calibration" / "legacy_d10"
    / "particle_masses_paper_d10_d11.py"
)

PREFIX = "oph_charged_icosahedral_completion/"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def archive_member_json(archive: ZipFile, name: str) -> dict[str, Any]:
    return json.loads(archive.read(PREFIX + name))


def safe_archive_names(names: list[str]) -> bool:
    for name in names:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            return False
    return True


def verify_manifest(archive: ZipFile, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    for row in manifest["files"]:
        raw = archive.read(PREFIX + row["file"])
        checks.append(
            {
                "file": row["file"],
                "bytes_match": len(raw) == int(row["bytes"]),
                "sha256_match": sha256(raw) == row["sha256"],
            }
        )
    return checks


@lru_cache(maxsize=None)
def d10_branch(p_value: float) -> dict[str, float]:
    legacy_dir = str(LEGACY_D10_SOURCE.parent)
    if legacy_dir not in sys.path:
        sys.path.insert(0, legacy_dir)
    from particle_masses_paper_d10_d11 import (  # type: ignore[import-not-found]
        solve_alpha_u_from_p,
        unification_scale_gev,
    )

    alpha_u, report = solve_alpha_u_from_p(
        p_value,
        3,
        unification_scale_gev(p_value),
    )
    return {
        "P": p_value,
        "alpha_U": float(alpha_u),
        "v_gev": float(report["v"]),
    }


def face_candidate(p_value: float, v_gev: float, alpha_u: float) -> dict[str, Any]:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    alpha_p = (p_value - phi) / math.sqrt(math.pi)
    kappa = (alpha_u / 50.0) / (
        1.0 + alpha_u / 31.0 + alpha_u**2 / (31.0 * 10.0)
    )
    chi = (alpha_u**2 / 512.0) / (1.0 - alpha_u / 77.0)
    zeta = (alpha_p**2 / 21.0) / (
        1.0 - alpha_u / 27.0 - alpha_u**2 / (27.0 * 5.0)
    )
    rho = math.sqrt(2.0) * math.exp(-chi)
    delta = 2.0 / 9.0 + zeta
    roots = [
        1.0 + rho * math.cos(delta + 2.0 * math.pi * index / 3.0)
        for index in range(3)
    ]
    squares = sorted(root * root for root in roots)
    shape_geomean = math.prod(squares) ** (1.0 / 3.0)
    g0 = v_gev / (2.0 * 6.0**14) ** (1.0 / 3.0)
    geometric_mean = g0 * math.exp(kappa)
    masses = [geometric_mean * value / shape_geomean for value in squares]
    return {
        "P": p_value,
        "v_gev": v_gev,
        "alpha_U": alpha_u,
        "kappa": kappa,
        "chi": chi,
        "zeta": zeta,
        "delta": delta,
        "rho": rho,
        "masses_gev": masses,
    }


def inverse_coordinates(masses_gev: list[float], v_gev: float) -> dict[str, float]:
    electron, muon, tau = masses_gev
    geometric_mean = math.prod(masses_gev) ** (1.0 / 3.0)
    g0 = v_gev / (2.0 * 6.0**14) ** (1.0 / 3.0)
    kappa = math.log(geometric_mean / g0)

    invariant = sum(masses_gev) / sum(math.sqrt(value) for value in masses_gev) ** 2
    chi = -0.5 * math.log(3.0 * invariant - 1.0)
    rho = math.sqrt(2.0) * math.exp(-chi)
    root_sum = sum(math.sqrt(value) for value in masses_gev)
    root0 = 3.0 * math.sqrt(tau) / root_sum
    root1 = 3.0 * math.sqrt(electron) / root_sum
    root2 = 3.0 * math.sqrt(muon) / root_sum
    delta = math.atan2(
        (root2 - root1) / (math.sqrt(3.0) * rho),
        (root0 - 1.0) / rho,
    )
    return {
        "kappa": kappa,
        "chi": chi,
        "zeta": delta - 2.0 / 9.0,
        "delta": delta,
        "square_root_mass_invariant": invariant,
    }


def comparison(branch: dict[str, Any], references: list[float]) -> dict[str, Any]:
    masses = branch["masses_gev"]
    return {
        "masses_mev": [1000.0 * value for value in masses],
        "relative_residual_ppm": [
            1.0e6 * (value / target - 1.0)
            for value, target in zip(masses, references, strict=True)
        ],
        "max_absolute_residual_ppm": max(
            abs(1.0e6 * (value / target - 1.0))
            for value, target in zip(masses, references, strict=True)
        ),
    }


def build_review(archive_path: Path) -> dict[str, Any]:
    archive_raw = archive_path.read_bytes()
    with ZipFile(archive_path) as archive:
        names = archive.namelist()
        first_bad_crc_member = archive.testzip()
        manifest = archive_member_json(archive, "MANIFEST.json")
        candidate = archive_member_json(archive, "candidate.json")
        submitted_audit = archive_member_json(archive, "audit.json")
        reference = archive_member_json(archive, "compare_only_reference.json")
        manifest_checks = verify_manifest(archive, manifest)
        expected_names = {
            PREFIX,
            PREFIX + "MANIFEST.json",
            *(PREFIX + row["file"] for row in manifest["files"]),
        }

    names_ordered = ("electron", "muon", "tau")
    references = [float(reference["entries"][name]["value_gev"]) for name in names_ordered]
    submitted_inputs = candidate["source_inputs"]
    submitted = face_candidate(
        float(submitted_inputs["P"]),
        float(submitted_inputs["v_gev"]),
        float(submitted_inputs["alpha_U"]),
    )
    proposed_coordinates = {
        "kappa": float(candidate["finite_incidence_completion"]["kappa_affine"]),
        "chi": float(candidate["finite_incidence_completion"]["chi_amplitude"]),
        "zeta": float(candidate["finite_incidence_completion"]["zeta_phase"]),
    }
    target_coordinates = inverse_coordinates(references, float(submitted_inputs["v_gev"]))

    public = json.loads(PUBLIC_PIXEL.read_text(encoding="utf-8"))
    source = json.loads(SOURCE_PIXEL.read_text(encoding="utf-8"))
    d10_gap = json.loads(D10_GAP.read_text(encoding="utf-8"))
    krawczyk = json.loads(KRAWCZYK.read_text(encoding="utf-8"))
    public_p = float(public["public_endpoint_convention"]["P_C"])
    source_p = float(source["P_cand"])
    submitted_p = float(submitted_inputs["P"])
    submitted_v = float(submitted_inputs["v_gev"])
    truncated_d10 = d10_branch(submitted_p)
    public_d10 = d10_branch(public_p)
    source_d10 = d10_branch(source_p)
    public_alpha = public_d10["alpha_U"]
    source_alpha = source_d10["alpha_U"]
    public_v = public_d10["v_gev"]
    source_v = source_d10["v_gev"]
    krawczyk_alpha = float(krawczyk["center_c"])
    source_certificate_alpha = float(source["alpha_U_P_cand"])

    if not math.isclose(submitted_v, truncated_d10["v_gev"], rel_tol=0.0, abs_tol=2.0e-12):
        raise ValueError("submitted v is not the D10 value at submitted P")
    if not math.isclose(public_alpha, krawczyk_alpha, rel_tol=0.0, abs_tol=2.0e-15):
        raise ValueError("public D10 alpha_U does not match the Krawczyk center")
    if not math.isclose(source_alpha, source_certificate_alpha, rel_tol=0.0, abs_tol=2.0e-15):
        raise ValueError("source-audit D10 alpha_U does not match its certificate")

    branch_rows = {
        "submitted_hybrid_truncated": submitted,
        "coherent_truncated_P": face_candidate(
            submitted_p, submitted_v, truncated_d10["alpha_U"]
        ),
        "coherent_canonical_public": face_candidate(public_p, public_v, public_alpha),
        "coherent_source_audit": face_candidate(source_p, source_v, source_alpha),
    }
    branch_audit = {
        name: {
            "inputs": {key: row[key] for key in ("P", "v_gev", "alpha_U")},
            **comparison(row, references),
        }
        for name, row in branch_rows.items()
    }

    candidate_masses = [float(candidate["readout"]["masses_gev"][name]) for name in names_ordered]
    independent_rebuild_max_abs = max(
        abs(left - right)
        for left, right in zip(submitted["masses_gev"], candidate_masses, strict=True)
    )
    coordinate_differences = {
        key: proposed_coordinates[key] - target_coordinates[key]
        for key in ("kappa", "chi", "zeta")
    }
    coordinate_relative_differences = {
        key: proposed_coordinates[key] / target_coordinates[key] - 1.0
        for key in ("kappa", "chi", "zeta")
    }

    tau_api_mev = float(reference["entries"]["tau"]["api_value"])
    tau_used_mev = 1000.0 * references[2]
    tau_candidate_mev = 1000.0 * candidate_masses[2]
    try:
        archive_display_path = str(archive_path.relative_to(WORKSPACE))
    except ValueError:
        archive_display_path = archive_path.name

    return {
        "artifact": "oph_charged_icosahedral_completion_v1_adversarial_review",
        "status": "FROZEN_RETROSPECTIVE_TARGET_INFORMED_FACE_CANDIDATE_NOT_COMPLETION",
        "compare_only": True,
        "charged_reference_data_consumed": True,
        "runtime_builder_reference_separated": True,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "archive": {
            "path": archive_display_path,
            "sha256": sha256(archive_raw),
            "safe_member_paths": safe_archive_names(names),
            "no_duplicate_member_names": len(names) == len(set(names)),
            "exact_manifest_member_set": set(names) == expected_names,
            "crc_check_pass": first_bad_crc_member is None,
            "crc_first_bad_member": first_bad_crc_member,
            "manifest_checks": manifest_checks,
            "manifest_checks_pass": all(
                row["bytes_match"] and row["sha256_match"] for row in manifest_checks
            ),
            "submitted_candidate_status": candidate["status"],
            "submitted_audit_status": submitted_audit["status"],
            "independent_rebuild_max_abs_mass_difference_gev": independent_rebuild_max_abs,
        },
        "valid_geometric_result": {
            "receipt": "charged_icosahedral_face_carrier_frontier.json",
            "statement": (
                "The OPH icosahedron supplies the A5/C3 face orbit and local three-corner "
                "C3 fibers. It does not canonically supply one global physical generation "
                "space or the spectral value laws."
            ),
        },
        "upstream_receipts": {
            str(PUBLIC_PIXEL.relative_to(REPO)): sha256(PUBLIC_PIXEL.read_bytes()),
            str(SOURCE_PIXEL.relative_to(REPO)): sha256(SOURCE_PIXEL.read_bytes()),
            str(D10_GAP.relative_to(REPO)): sha256(D10_GAP.read_bytes()),
            str(KRAWCZYK.relative_to(REPO)): sha256(KRAWCZYK.read_bytes()),
            str(LEGACY_D10_SOURCE.relative_to(REPO)): sha256(LEGACY_D10_SOURCE.read_bytes()),
        },
        "three_coordinate_inverse_audit": {
            "coordinate_count": 3,
            "target_inferred": target_coordinates,
            "proposed_equalizers": proposed_coordinates,
            "proposed_minus_target": coordinate_differences,
            "relative_differences": coordinate_relative_differences,
            "conclusion": (
                "The three post-hoc equalizers match the three independently target-inferred "
                "spectral coordinates to roughly one part per million or better. Runtime "
                "separation does not remove this historical discrete formula selection."
            ),
        },
        "branch_provenance": {
            "submitted_P": float(submitted_inputs["P"]),
            "submitted_alpha_U": float(submitted_inputs["alpha_U"]),
            "alpha_U_at_submitted_P_from_D10_recompute": truncated_d10["alpha_U"],
            "submitted_v_matches_D10_recompute": True,
            "submitted_P_matches_compare_only_D10_probe": (
                float(submitted_inputs["P"]) == float(d10_gap["p"])
            ),
            "d10_probe_artifact_status": d10_gap["status"],
            "canonical_public_P": public_p,
            "canonical_public_P_status": public["status"],
            "public_alpha_U_matches_Krawczyk_center": True,
            "public_P_uses_measured_Thomson_endpoint": True,
            "source_audit_P": source_p,
            "source_audit_P_status": source["status"],
            "source_alpha_U_matches_source_certificate": True,
            "branch_rows": branch_audit,
        },
        "tau_precision_audit": {
            "comparison_field_used_mev": tau_used_mev,
            "packet_api_value_mev": tau_api_mev,
            "candidate_mev": tau_candidate_mev,
            "residual_ppm_against_used_rounded_field": 1.0e6 * (tau_candidate_mev / tau_used_mev - 1.0),
            "residual_ppm_against_packet_api_value": 1.0e6 * (tau_candidate_mev / tau_api_mev - 1.0),
        },
        "proved_and_open": {
            "proved": [
                "archive safety, hashes, and numerical reproducibility",
                "icosahedral face incidence and local C3 circulant algebra",
                "algebraic mass readout conditional on all proposed laws",
            ],
            "open": [
                "quotient-visible attachment to the physical charged generation space",
                "ln(2) MaxEnt action gap and ensemble-to-Yukawa bridge",
                "charged connection and phase-transport law",
                "normalized Z6 charged determinant character",
                "derivation of all three equalizers from a declared repair operator with uniqueness, stability, and refinement",
                "one coherent source branch and a running-to-pole mass-scheme map",
            ],
        },
        "integration_decision": (
            "Preserve as the leading frozen icosahedral face-incidence completion conjecture "
            "and retrospective numerical candidate. Do not replace public theorem mass rows "
            "or the current non-identifiability result until the open source theorems close."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    review = build_review(args.archive.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(review, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": review["status"], "out": str(args.out)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
