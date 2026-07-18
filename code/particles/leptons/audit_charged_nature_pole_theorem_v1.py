#!/usr/bin/env python3
"""Audit the submitted charged nature/pole bridge theorem.

The archive is treated as a compare-only theorem packet.  This auditor reads
its hashed payloads without executing submitted programs, independently
reconstructs the finite face-response spectrum, and separates three claims:

1. the finite algebraic response calculation;
2. the conditional positive-square-root and pole-zero transport lemmas; and
3. the missing physical OPH nature-identification and interacting-kernel
   certificates.

The load-bearing NI6 and RP4 premises already identify the physical Yukawa
response and singularity readout with the face operator.  Consequently the
downstream equality is valid but is not a derivation of those premises from
OPH dynamics.
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


HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
REPO = HERE.parents[3]
WORKSPACE = HERE.parents[4]
DEFAULT_ARCHIVE = (
    WORKSPACE / "correspondence" / "oph_charged_nature_pole_theorem_v1.zip"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_nature_pole_theorem_v1_review.json"
)
PARENT_REVIEW = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_icosahedral_completion_v1_review.json"
)

PREFIX = "oph_charged_nature_pole_theorem_v1/"
EXPECTED_ARCHIVE_SHA256 = (
    "73e8f5a8d0ca2a5a0c306b0f4f5f319c3000c98a999a240021bc13879dcbf897"
)
EXPECTED_NARRATIVE_SHA256 = (
    "19054f289bf363d850bda4ff299d9334894dbfbda40fdaac0feea4a5e441a7ee"
)
EXPECTED_THEOREM_SHA256 = (
    "a6439d0080c5b97f02d8909eed20f336ec49d8d78718d1c03a53e609f318667e"
)
EXPECTED_BUILDER_SHA256 = (
    "c7c2c25b7d337d3181cac7082f39bf4eac58987e3ca9b63330f0d09523acadba"
)
EXPECTED_VERIFIER_SHA256 = (
    "eafcc13c17db2d4065cb9da138fd8d426eef9e84f60303fb1248c73094e3b646"
)
EXPECTED_THEOREM_RECEIPT_SHA256 = (
    "63e65606f6fef573c849002c94816eb0f29b53bad23bd6c9b17d29b655a1b96f"
)
EXPECTED_VERIFICATION_RECEIPT_SHA256 = (
    "fe7333565819c4faa480ed5107adef59188a6da9a01b1929687953d10f2e7ce7"
)
EXPECTED_NEGATIVE_RECEIPT_SHA256 = (
    "1c8e6595738629f17a8429d01a0cf95fbd220094185ad7657d74552146cee206"
)
EXPECTED_FINAL_STATUS_SHA256 = (
    "57842db63817351c7500aa97489f0ce2447bcce22418e9bd2c6e34cc36ac80d7"
)

mp.mp.dps = 90


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
    result: dict[str, str] = {}
    for line in raw.decode("utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        result[name.strip()] = digest
    return result


def close(left: mp.mpf, right: mp.mpf, tolerance: mp.mpf = mp.mpf("1e-65")) -> bool:
    return abs(left - right) <= tolerance * max(mp.mpf(1), abs(left), abs(right))


def reconstruct_face_endpoint(nature: dict[str, Any]) -> dict[str, Any]:
    source = nature["source_inputs"]
    P = mp.mpf(source["P"])
    v = mp.mpf(source["v_gev"])
    alpha_u = mp.mpf(source["alpha_U"])
    phi = (1 + mp.sqrt(5)) / 2
    alpha_p = (P - phi) / mp.sqrt(mp.pi)

    q_kappa = alpha_u / 31 + alpha_u**2 / 310
    q_chi = alpha_u / 77
    q_zeta = alpha_u / 27 + alpha_u**2 / 135
    kappa = (alpha_u / 50) / (1 + q_kappa)
    chi = (alpha_u**2 / 512) / (1 - q_chi)
    zeta = (alpha_p**2 / 21) / (1 - q_zeta)
    delta = mp.mpf(2) / 9 + zeta
    rho = mp.sqrt(2) * mp.e ** (-chi)
    roots = [1 + rho * mp.cos(delta + 2 * mp.pi * k / 3) for k in range(3)]
    squares = sorted(root * root for root in roots)
    shape_mean = mp.power(mp.fprod(squares), mp.mpf(1) / 3)
    g_end = v / mp.power(2 * mp.power(6, 14), mp.mpf(1) / 3) * mp.e**kappa
    masses = [g_end * square / shape_mean for square in squares]

    fixed = nature["fixed_point"]
    submitted = nature["physical_mass_and_yukawa_operators"]["masses_gev"]
    submitted_masses = [
        mp.mpf(submitted[label]) for label in ("electron", "muon", "tau")
    ]
    fixed_checks = {
        "alpha_P": close(mp.mpf(fixed["alpha_P"]), alpha_p),
        "kappa": close(mp.mpf(fixed["kappa"]), kappa),
        "chi": close(mp.mpf(fixed["chi_rho"]), chi),
        "zeta": close(mp.mpf(fixed["zeta_delta"]), zeta),
        "strict_contraction": mp.mpf(fixed["contraction_constant"]) < 1,
    }
    mass_checks = [close(left, right) for left, right in zip(masses, submitted_masses)]
    return {
        "fixed_point_checks": fixed_checks,
        "mass_checks": mass_checks,
        "reconstructed_masses_mev": [mp.nstr(1000 * mass, 30) for mass in masses],
        "submitted_masses_mev": [mp.nstr(1000 * mass, 30) for mass in submitted_masses],
        "g_end_gev": mp.nstr(g_end, 30),
        "checks_pass": all(fixed_checks.values()) and all(mass_checks),
    }


def build_review(archive_path: Path, parent_review_path: Path = PARENT_REVIEW) -> dict[str, Any]:
    archive_raw = archive_path.read_bytes()
    parent = json.loads(parent_review_path.read_bytes())

    with ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        first_bad_crc = archive.testzip()
        hashes = parse_hashes(archive.read(PREFIX + "SHA256SUMS"))

        def raw(name: str) -> bytes:
            return archive.read(PREFIX + name)

        def packet(name: str) -> dict[str, Any]:
            return json.loads(raw(name))

        internal_hash_checks = {
            name: sha256(raw(name)) == digest for name, digest in hashes.items()
        }
        theorem_raw = raw("NATURE_IDENTIFICATION_AND_POLE_THEOREM.md")
        builder_raw = raw("build_nature_pole_packet.py")
        verifier_raw = raw("verify_nature_pole_theorem.py")
        theorem_receipt_raw = raw("nature_pole_theorem_receipt.json")
        verification_raw = raw("nature_pole_verification_receipt.json")
        negative_raw = raw("negative_control_receipt.json")
        final_raw = raw("FINAL_STATUS.json")
        nature = packet("nature_identification_packet.json")
        pole = packet("pole_operator_packet.json")
        certificate = packet("nature_pole_certificate_template.json")
        theorem_receipt = json.loads(theorem_receipt_raw)
        verification = json.loads(verification_raw)
        negative = json.loads(negative_raw)
        final_status = json.loads(final_raw)

    expected_names = {PREFIX + name for name in hashes} | {PREFIX + "SHA256SUMS"}
    archive_checks = {
        "archive_hash_matches": sha256(archive_raw) == EXPECTED_ARCHIVE_SHA256,
        "safe_non_symlink_unencrypted_members": all(safe_member(info) for info in infos),
        "no_duplicate_member_names": len(names) == len(set(names)),
        "exact_member_set": set(names) == expected_names,
        "crc_check_pass": first_bad_crc is None,
        "nineteen_internal_hashes_declared": len(hashes) == 19,
        "all_internal_hashes_match": all(internal_hash_checks.values()),
        "theorem_hash_matches": sha256(theorem_raw) == EXPECTED_THEOREM_SHA256,
        "builder_hash_matches": sha256(builder_raw) == EXPECTED_BUILDER_SHA256,
        "verifier_hash_matches": sha256(verifier_raw) == EXPECTED_VERIFIER_SHA256,
        "theorem_receipt_hash_matches": (
            sha256(theorem_receipt_raw) == EXPECTED_THEOREM_RECEIPT_SHA256
        ),
        "verification_receipt_hash_matches": (
            sha256(verification_raw) == EXPECTED_VERIFICATION_RECEIPT_SHA256
        ),
        "negative_receipt_hash_matches": (
            sha256(negative_raw) == EXPECTED_NEGATIVE_RECEIPT_SHA256
        ),
        "final_status_hash_matches": sha256(final_raw) == EXPECTED_FINAL_STATUS_SHA256,
    }

    reconstruction = reconstruct_face_endpoint(nature)
    physical_status_checks = {
        "submitted_conditional_theorem_pass": verification["conditional_theorem_pass"] is True,
        "submitted_unconditional_claim_false": (
            verification["unconditional_nature_pole_claim"] is False
        ),
        "final_nature_identification_open": (
            final_status["physical_nature_identification_pass"] is False
        ),
        "final_interacting_kernel_open": (
            final_status["interacting_renormalized_pole_operator_pass"] is False
        ),
        "promotion_forbidden": final_status["promotion_allowed"] is False,
        "all_NI_physical_gates_open": all(
            str(status).startswith("OPEN")
            for key, status in theorem_receipt["physical_gate_status"].items()
            if key.startswith("NI")
        ),
        "all_RP_physical_gates_open": all(
            str(status).startswith("OPEN")
            for key, status in theorem_receipt["physical_gate_status"].items()
            if key.startswith("RP")
        ),
        "minimal_kernel_is_free_identity_witness": (
            pole["minimal_exact_realization"]["self_energy"] == "Sigma_0(s)=0"
            and pole["minimal_exact_realization"]["nature_status"]
            == "EXISTENCE_WITNESS_ONLY"
        ),
        "interacting_kernel_contract_open": (
            pole["interacting_kernel_contract"]["status"]
            == "OPEN_UNTIL_SOURCE_1PI_KERNEL_IS_EMITTED"
        ),
        "infrared_contract_open": pole["infrared_contract"]["status"].startswith("OPEN"),
    }

    theorem_text = theorem_raw.decode("utf-8")
    logical_boundary = {
        "regular_C3_character_lemma_valid": True,
        "positive_square_root_transport_valid": True,
        "regular_left_right_field_factors_preserve_zero_set": True,
        "Nielsen_zero_set_lemma_valid_conditional_on_regular_factorization": True,
        "NI6_assumes_physical_Yukawa_response_equals_face_response": (
            "Charged response isometry" in theorem_text
            and "widehat Y_e\\widehat Y_e^\\dagger" in theorem_text
        ),
        "RP4_assumes_singularity_readout_equals_face_response": (
            "CFQ–Dyson fixed-point intertwiner" in theorem_text
            and "mathcal P(\\mathfrak B(x))=J_LM_F(x)J_L^\\dagger" in theorem_text
        ),
        "substantive_necessary_and_sufficient_OPH_derivation_proved": False,
    }

    submitted_masses = [
        mp.mpf(nature["physical_mass_and_yukawa_operators"]["masses_gev"][label])
        for label in ("electron", "muon", "tau")
    ]
    alternate_masses = [2 * submitted_masses[0], submitted_masses[1] / 2, submitted_masses[2]]
    alternate_same_product = close(
        mp.fprod(alternate_masses), mp.fprod(submitted_masses), mp.mpf("1e-60")
    )
    adversarial_verifier_boundary = {
        "submitted_verifier_reconstructs_fixed_point": True,
        "submitted_verifier_reconstructs_response_shape": True,
        "submitted_verifier_checks_mass_matrix_internal_spectrum": True,
        "submitted_verifier_checks_mass_determinant": True,
        "submitted_verifier_binds_M_to_g_end_times_response_shape": False,
        "submitted_verifier_binds_parent_artifact_hashes": False,
        "same_product_alternate_spectrum_exists": (
            alternate_same_product
            and 0 < alternate_masses[0] < alternate_masses[1] < alternate_masses[2]
        ),
        "alternate_spectrum_factors": ["2", "1/2", "1"],
        "alternate_spectrum_mev": [mp.nstr(1000 * mass, 24) for mass in alternate_masses],
        "consequence": (
            "A coherent replacement of M, Y, their three mass fields, and the pole roots by "
            "this same-product ordered spectrum satisfies the submitted verifier's tested "
            "mass/projector/determinant/pole relations because it never checks M=g_end*S."
        ),
    }

    provenance_checks = {
        "runtime_external_files_empty": (
            nature["provenance"]["runtime_source_separation"] is True
        ),
        "historical_target_exposure_declared": (
            nature["provenance"]["historical_target_exposure"] is True
        ),
        "prospective_prediction_status_rejected": (
            nature["provenance"]["prospective_prediction_status"] == "NOT_PROSPECTIVE"
        ),
        "parent_artifacts_hash_bound_in_bundle": False,
        "source_coordinates_are_builder_defaults_not_loaded_parent_receipts": True,
        "parent_candidate_branch_coherent": parent["branch_tuple_coherent"] is True,
    }

    checks_pass = (
        all(archive_checks.values())
        and reconstruction["checks_pass"]
        and all(physical_status_checks.values())
        and all(value for key, value in logical_boundary.items() if key != "substantive_necessary_and_sufficient_OPH_derivation_proved")
        and adversarial_verifier_boundary["same_product_alternate_spectrum_exists"]
        and negative["all_negative_controls_rejected"] is True
    )

    return {
        "artifact": "oph_charged_nature_pole_theorem_v1_review",
        "schema_version": 1,
        "status": (
            "CONDITIONAL_NATURE_AND_POLE_TRANSPORT_LEMMAS_VALID_"
            "PHYSICAL_NI_RP_CERTIFICATES_SOURCE_BINDING_AND_PROMOTION_OPEN"
        ),
        "compare_only": True,
        "forbidden_as_candidate_ancestor": True,
        "runtime_charged_reference_consumed": False,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "physical_nature_identification_proved": False,
        "interacting_pole_kernel_proved": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "provenance": {
            "archive_path": (
                archive_path.relative_to(WORKSPACE).as_posix()
                if archive_path.is_relative_to(WORKSPACE)
                else archive_path.name
            ),
            "archive_sha256": sha256(archive_raw),
            "narrative_attachment_sha256": EXPECTED_NARRATIVE_SHA256,
            "parent_candidate_review": parent_review_path.relative_to(REPO).as_posix(),
        },
        "archive": {
            "member_count": len(infos),
            "crc_first_bad_member": first_bad_crc,
            "internal_hash_checks": internal_hash_checks,
            "checks": archive_checks,
            "checks_pass": all(archive_checks.values()),
            "execution_boundary": "Submitted programs are not executed by this canonical auditor.",
        },
        "finite_face_endpoint_reconstruction": reconstruction,
        "conditional_theorem_boundary": {
            "checks": logical_boundary,
            "valid_consequence": (
                "Given NI6, the unique positive square root transports M_F to the physical "
                "left charged response. Given RP4 and the exact regular renormalized-kernel "
                "receipts, the declared singularity readout has the same operator."
            ),
            "nonconsequence": (
                "The theorem does not derive NI6, RP4, the chiral carrier, the physical "
                "family attachment, or the interacting charged 1PI kernel from OPH."
            ),
            "promotion_equivalence_boundary": (
                "The advertised necessity direction permits any equivalent packet that "
                "restates the desired physical equality. It is a complete gate checklist, "
                "not an independent source-selection or uniqueness theorem."
            ),
        },
        "physical_gate_status": {
            "checks": physical_status_checks,
            "submitted_NI_RP_status": theorem_receipt["physical_gate_status"],
        },
        "submitted_verifier_boundary": adversarial_verifier_boundary,
        "provenance_and_branch_boundary": {
            "checks": provenance_checks,
            "reason": (
                "The builder hard-codes P, v, and alpha_U defaults and supplies no hashed "
                "parent receipts. The same target-informed hybrid candidate remains the "
                "numerical ancestor, so a clean runtime invocation does not repair source "
                "or historical ancestry."
            ),
        },
        "independent_reproduction_observation": {
            "submitted_verifier_returned_conditional_pass": True,
            "submitted_negative_controls_returned_pass": True,
            "builder_semantically_reproduced_source_coordinates_and_masses": True,
            "builder_byte_reproduced_nature_packet": False,
            "builder_byte_reproduced_pole_and_certificate_packets": True,
            "audit_environment": "Python 3.13.7, NumPy 2.3.5, mpmath 1.3.0, SymPy 1.14.0",
            "byte_drift_cause": (
                "Floating eigensystem and matrix serialization changed while the analytic "
                "source coordinates and masses remained equivalent. No dependency lock is supplied."
            ),
        },
        "external_reference_audit": {
            "mixed_fermion_pole_algebra": "consistent with arXiv:1308.3140",
            "Nielsen_complex_pole_gauge_independence": "consistent with arXiv:hep-ph/9907254",
            "charged_infraparticle_boundary": (
                "consistent with Buchholz, Phys. Lett. B 174 (1986) 331, and arXiv:2307.06114"
            ),
            "OPH_specific_bridge_supplied_by_references": False,
        },
        "still_open": [
            "derive the chiral charged-lepton carrier from the OPH reconstruction",
            "derive a quotient-visible refinement-natural regular C3 family action and intertwiners",
            "derive NI6 rather than assume the charged response-isometry identity",
            "derive the physical determinant line and remove the hybrid source branch",
            "emit the BRST/Ward-complete interacting charged 1PI kernel",
            "derive RP4 rather than assume the CFQ-to-Dyson singularity readout",
            "supply QED infrared dressing, thresholds, running, and pole/threshold conventions",
            "bind every upstream source receipt before comparison",
        ],
        "integration_decision": (
            "Integrate the regular-C3, positive-square-root, field-redefinition, Nielsen, and "
            "fixed-point transport implications as conditional bridge lemmas. Retain the "
            "minimal zero-self-energy kernel as an existence witness only. Do not promote "
            "physical nature identification, an interacting pole theorem, or charged masses."
        ),
        "checks_pass": checks_pass,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--parent-review", type=Path, default=PARENT_REVIEW)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    review = build_review(args.archive.resolve(), args.parent_review.resolve())
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
                "physical_nature_identification_proved": review[
                    "physical_nature_identification_proved"
                ],
                "interacting_pole_kernel_proved": review[
                    "interacting_pole_kernel_proved"
                ],
                "public_prediction_promotion_allowed": review[
                    "public_prediction_promotion_allowed"
                ],
            },
            indent=2,
        )
    )
    return 0 if review["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
