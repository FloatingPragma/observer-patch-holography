#!/usr/bin/env python3
"""Independent verifier for the FZ-12 formal observation-map packet.

This verifier intentionally does not import the producer.  It uses only the
Python standard library, fixed custody digests, and a separate Fraction
derivation of every reported coefficient and truncation identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE / "runtime" / "seam_current_edge_prediction_receipt.json"
DEFAULT_CUSTODY = HERE / "runtime" / "fz12_custody_projection.json"
DEFAULT_RECEIPT = HERE / "runtime" / "fz12_observation_map_receipt.json"

SOURCE_SHA256 = "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
SOURCE_BYTES = 9296
CUSTODY_SHA256 = "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
CUSTODY_BYTES = 3624
CUSTODY_SELF_SHA256 = (
    "sha256:334d9877023d4eaa2ac64bca51d9e2193f0eaad81afd664b22cada2b890a8527"
)
RECEIPT_SCHEMA = "oph.fz12.formal_observation_map.v1"
RECEIPT_STATUS = (
    "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
    "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
)


class VerificationError(ValueError):
    """The independent verifier rejected an input."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def check_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    """Reject missing and unknown fields on a typed receipt surface."""

    check(set(value) == expected, f"{label} keys drift")


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def self_digest(value: dict[str, Any], field: str) -> str:
    body = {key: item for key, item in value.items() if key != field}
    return "sha256:" + sha256(canonical_json_bytes(body))


def load_json(path: Path, *, canonical: bool = True) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError(f"invalid JSON: {path}") from error
    check(isinstance(value, dict), f"JSON root is not an object: {path}")
    if canonical:
        check(raw == canonical_json_bytes(value), f"noncanonical JSON: {path}")
    return raw, value


def verify_source(path: Path) -> dict[str, Any]:
    raw, source = load_json(path)
    check(sha256(raw) == SOURCE_SHA256, "source fixed byte digest drift")
    check(len(raw) == SOURCE_BYTES, "source byte-count drift")
    check(
        source.get("schema") == "oph.seam_current_edge_prediction_candidate.v1",
        "source schema drift",
    )
    check(
        source.get("status")
        == (
            "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
            "PHYSICAL_PRODUCER_OPEN"
        ),
        "source status drift",
    )
    check(
        source.get("receipt_sha256") == self_digest(source, "receipt_sha256"),
        "source self-digest drift",
    )
    candidate = source.get("conditional_physical_candidate", {})
    check(
        candidate.get("coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "source coefficient drift",
    )
    check(
        candidate.get("scale_free_relations")
        == {
            "B0_over_C4_squared": "10/21",
            "B6_over_C4_squared": "-2/63",
            "B6_over_B0": "-1/15",
        },
        "source scale-free relation drift",
    )
    boundary = source.get("exposure_and_custody_boundary", {})
    for key in (
        "target_values_read",
        "comparison_data_read",
        "public_measurement_read",
        "comparison_permitted",
    ):
        check(boundary.get(key) is False, f"source exposure drift: {key}")
    return source


def verify_custody(path: Path) -> dict[str, Any]:
    raw, custody = load_json(path)
    check(sha256(raw) == CUSTODY_SHA256, "custody fixed byte digest drift")
    check(len(raw) == CUSTODY_BYTES, "custody byte-count drift")
    check(
        custody.get("schema") == "oph.fz12.custody_projection.v1",
        "custody schema drift",
    )
    check(
        custody.get("status")
        == "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
        "custody status drift",
    )
    check(
        custody.get("projection_sha256") == self_digest(custody, "projection_sha256"),
        "custody self-digest drift",
    )
    check(custody.get("projection_sha256") == CUSTODY_SELF_SHA256, "custody pin drift")
    scope = custody.get("projection_scope", {})
    check(scope.get("register_row") == "FZ-12", "custody row drift")
    for key in (
        "includes_measurement_values",
        "includes_comparison_values",
        "includes_likelihood_values",
        "includes_other_campaign_rows",
    ):
        check(scope.get(key) is False, f"custody data-scope drift: {key}")
    source = custody.get("source_receipt", {})
    check(source.get("bytes") == SOURCE_BYTES, "custody source byte-count drift")
    check(source.get("sha256") == SOURCE_SHA256, "custody source digest drift")
    row = custody.get("frozen_row_projection", {})
    check(row.get("id") == "FZ-12", "custody projected row drift")
    check(row.get("content_sha256") == SOURCE_SHA256, "custody content pin drift")
    check(
        row.get("canonical_full_row_sha256")
        == "4b83abb7c65c1a5f87d0839952a1e62bd88c8cd49f4694013cebe683637fc9c7",
        "custody full-row pin drift",
    )
    contract = custody.get("external_custody_contract", {})
    check(contract.get("rows") == ["FZ-12"], "custody contract row drift")
    check(
        contract.get("source_commit") == "bc5595f8dbb2d2886e2a64ddf447f69fbb00eb3f",
        "source commit drift",
    )
    check(
        contract.get("custody_commit") == "54b450af0bb5bd0fee4842f5c5f654d08d6baa2d",
        "custody commit drift",
    )
    check(
        contract.get("decision_rule_custody_commit")
        == "25da61a800226e0232336ccc86de8dec7d6b51c6",
        "decision-rule custody commit drift",
    )
    check(
        contract.get("artifact_sha256", {}).get(
            "seam_current_edge_prediction_frozen_2026-08-02.json"
        )
        == "d6d8f7e299dc8b38efd88c3e27135c5e4ca8eddc7bc11425f0b6844908ef76df",
        "frozen prediction artifact pin drift",
    )
    check(
        contract.get("in_repo_artifact_sha256", {}).get(
            "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json"
        )
        == SOURCE_SHA256,
        "custody source ancestry drift",
    )
    return custody


def independent_formulas() -> tuple[dict[str, Any], dict[str, str]]:
    c4 = Fraction(-1, 20)
    b0 = Fraction(1, 840)
    b6 = Fraction(-1, 12600)
    p = c4 / 2
    q0 = b0 / 2 - c4 * c4 / 8
    q6 = b6 / 2
    check(2 * p == c4, "independent k4 square-root identity failed")
    check(p * p + 2 * q0 == b0, "independent isotropic k6 identity failed")
    check(2 * q6 == b6, "independent anisotropic k6 identity failed")

    r2 = 3 * p
    r40 = 5 * q0
    r46 = 5 * q6
    s2 = -r2
    s40 = r2 * r2 - r40
    s46 = -r46
    check(r2 + s2 == 0, "independent inverse k2 identity failed")
    check(r40 + s40 + r2 * s2 == 0, "independent inverse isotropic k4 failed")
    check(r46 + s46 == 0, "independent inverse anisotropic k4 failed")
    check(q6 == b6 / 2, "independent transverse coefficient failed")

    coefficients = {
        "omega": {
            "k": "1",
            "a2_k3": str(p),
            "a4_k5_isotropic": str(q0),
            "a4_k5_I6": str(q6),
        },
        "radial_group_velocity": {
            "constant": "1",
            "a2_k2": str(r2),
            "a4_k4_isotropic": str(r40),
            "a4_k4_I6": str(r46),
        },
        "inverse_radial_speed": {
            "constant": "1",
            "a2_k2": str(s2),
            "a4_k4_isotropic": str(s40),
            "a4_k4_I6": str(s46),
        },
        "formal_inverse_radial_speed_excess": {
            "a2_k2": str(s2),
            "a4_k4_isotropic": str(s40),
            "a4_k4_I6": str(s46),
        },
        "transverse_group_velocity": {"a4_k4_grad_S_I6": str(q6)},
    }
    c4_squared = c4 * c4
    scale_free = {
        "omega_k5_isotropic_over_C4_squared": str(q0 / c4_squared),
        "omega_k5_I6_over_C4_squared": str(q6 / c4_squared),
        "radial_k4_isotropic_over_C4_squared": str(r40 / c4_squared),
        "radial_k4_I6_over_C4_squared": str(r46 / c4_squared),
        "inverse_speed_k4_isotropic_over_C4_squared": str(s40 / c4_squared),
        "inverse_speed_k4_I6_over_C4_squared": str(s46 / c4_squared),
        "transverse_grad_S_I6_over_C4_squared": str(q6 / c4_squared),
    }
    return coefficients, scale_free


def verify_receipt(path: Path) -> dict[str, Any]:
    raw, receipt = load_json(path)
    check_keys(
        receipt,
        {
            "ancestry",
            "certificate_scope",
            "exact_formal_result",
            "exposure_boundary",
            "formal_input",
            "interpretation_boundary",
            "issue",
            "receipt_sha256",
            "schema",
            "status",
        },
        "receipt",
    )
    check(receipt.get("schema") == RECEIPT_SCHEMA, "receipt schema drift")
    check(receipt.get("status") == RECEIPT_STATUS, "receipt status drift")
    check(receipt.get("issue") == 666, "receipt issue drift")
    check(
        receipt.get("receipt_sha256") == self_digest(receipt, "receipt_sha256"),
        "receipt self-digest drift",
    )
    ancestry = receipt.get("ancestry", {})
    check_keys(
        ancestry,
        {"source_receipt", "custody_projection", "frozen_artifacts_modified"},
        "ancestry",
    )
    source = ancestry.get("source_receipt", {})
    check(
        source
        == {
            "path": "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
            "bytes": SOURCE_BYTES,
            "sha256": f"sha256:{SOURCE_SHA256}",
            "schema": "oph.seam_current_edge_prediction_candidate.v1",
            "status": (
                "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
                "PHYSICAL_PRODUCER_OPEN"
            ),
            "receipt_sha256": (
                "sha256:3b1344dbce713545deba54b48e027a88faedac9b3a8776795b9abe33cafde7c3"
            ),
        },
        "receipt source ancestry drift",
    )
    custody = ancestry.get("custody_projection", {})
    check(
        custody
        == {
            "path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
            "bytes": CUSTODY_BYTES,
            "sha256": f"sha256:{CUSTODY_SHA256}",
            "schema": "oph.fz12.custody_projection.v1",
            "status": "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
            "projection_sha256": CUSTODY_SELF_SHA256,
            "row": "FZ-12",
            "canonical_full_row_sha256": (
                "4b83abb7c65c1a5f87d0839952a1e62bd88c8cd49f4694013cebe683637fc9c7"
            ),
            "source_commit": "bc5595f8dbb2d2886e2a64ddf447f69fbb00eb3f",
            "custody_commit": "54b450af0bb5bd0fee4842f5c5f654d08d6baa2d",
            "decision_rule_custody_commit": (
                "25da61a800226e0232336ccc86de8dec7d6b51c6"
            ),
        },
        "receipt custody ancestry drift",
    )
    check(
        ancestry.get("frozen_artifacts_modified") is False, "frozen mutation flag drift"
    )
    serialized = json.dumps(receipt, sort_keys=True)
    check(
        "claims/frozen_prediction_register.json" not in serialized,
        "broad register dependency",
    )
    check(
        "inverse_speed_delay_per_unit_coordinate_length" not in serialized,
        "retired delay field survived",
    )

    formal_input = receipt.get("formal_input", {})
    check_keys(
        formal_input,
        {
            "dispersion_symbol",
            "dimensionless_normalization",
            "frozen_coefficients",
            "frozen_scale_free_relations",
            "extra_map_premise",
            "extra_map_premise_discharged_here",
        },
        "formal input",
    )
    check(
        formal_input.get("dispersion_symbol")
        == "Lambda = k^2 + C4 k^4 + (B0 + B6 I6(n)) k^6",
        "formal input symbol drift",
    )
    check(
        formal_input.get("dimensionless_normalization")
        == {
            "q": "q = a k",
            "lambda_hat": "lambda_hat = a^2 Lambda",
            "omega_hat": "omega_hat = a omega",
            "normalized_symbol": (
                "lambda_hat = q^2 - q^4/20 + (1/840 - I6/12600) q^6 + O(q^8)"
            ),
        },
        "dimensionless normalization drift",
    )
    check(
        formal_input.get("frozen_coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "formal input coefficient drift",
    )
    check(
        formal_input.get("frozen_scale_free_relations")
        == {
            "B0_over_C4_squared": "10/21",
            "B6_over_C4_squared": "-2/63",
            "B6_over_B0": "-1/15",
        },
        "formal input ray drift",
    )
    check(
        formal_input.get("extra_map_premise")
        == (
            "lambda_hat is identified with omega_hat squared and omega_hat is "
            "the positive formal branch with leading term q"
        ),
        "frequency-map premise statement drift",
    )

    coefficients, scale_free = independent_formulas()
    result = receipt.get("exact_formal_result", {})
    check_keys(
        result,
        {
            "coefficients_in_a_units",
            "formal_inverse_radial_speed_excess",
            "frequency",
            "inverse_radial_speed",
            "machine_checks",
            "radial_group_velocity",
            "ring_for_square_root",
            "ring_for_velocity_and_inverse",
            "scale_free_coefficients",
            "transverse_group_velocity",
        },
        "formal result",
    )
    check(
        result.get("coefficients_in_a_units") == coefficients,
        "result coefficient drift",
    )
    check(
        result.get("scale_free_coefficients") == scale_free, "scale-free result drift"
    )
    check(
        result.get("ring_for_square_root") == "Q[I6][q] modulo q^8", "square ring drift"
    )
    check(
        result.get("ring_for_velocity_and_inverse") == "Q[I6][q] modulo q^6",
        "velocity ring drift",
    )
    expected_display = {
        "frequency": (
            "omega = k - (a^2/40) k^3 + a^4 (19/67200 - I6/25200) k^5 + O(a^6 k^7)"
        ),
        "radial_group_velocity": (
            "d omega/dk = 1 - (3a^2/40) k^2 + a^4 (19/13440 - I6/5040) k^4 + O(a^6 k^6)"
        ),
        "inverse_radial_speed": (
            "1/(d omega/dk) = 1 + (3a^2/40) k^2 + a^4 "
            "(283/67200 + I6/5040) k^4 + O(a^6 k^6)"
        ),
        "formal_inverse_radial_speed_excess": (
            "1/(d omega/dk) - 1 = (3a^2/40) k^2 + a^4 "
            "(283/67200 + I6/5040) k^4 + O(a^6 k^6)"
        ),
        "transverse_group_velocity": (
            "(1/k) grad_S omega = -(a^4/25200) k^4 grad_S I6 + O(a^6 k^6)"
        ),
    }
    for key, expected in expected_display.items():
        check(result.get(key) == expected, f"displayed formal series drift: {key}")
    check("O(a^6 k^7)" in result.get("frequency", ""), "frequency order term missing")
    for key in (
        "radial_group_velocity",
        "inverse_radial_speed",
        "formal_inverse_radial_speed_excess",
        "transverse_group_velocity",
    ):
        check(
            "O(a^6 k^6)" in result.get(key, ""), f"velocity order term missing: {key}"
        )
    check(
        result.get("machine_checks")
        == {
            "omega_hat_squared_equals_lambda_hat_mod_q8": True,
            "radial_times_inverse_radial_equals_one_mod_q6": True,
            "formal_inverse_radial_speed_excess_is_inverse_minus_one": True,
            "transverse_coefficient_is_one_over_k_times_spherical_gradient": True,
        },
        "machine-check declaration drift",
    )

    scope = receipt.get("certificate_scope", {})
    check(
        scope
        == {
            "type": "exact derived formal-series map conditional on the frozen FZ-12 ray",
            "derived_only": True,
            "new_prediction": False,
            "physical_prediction": False,
            "comparison_result": False,
            "statement": (
                "the certificate changes coordinates from the frozen spatial "
                "symbol coefficients to formal frequency and velocity coefficients"
            ),
        },
        "certificate scope drift",
    )
    check(
        receipt.get("formal_input", {}).get("extra_map_premise_discharged_here")
        is False,
        "frequency-map premise promotion drift",
    )
    boundary = receipt.get("interpretation_boundary", {})
    check(
        boundary
        == {
            "formal_truncation_only": True,
            "analytic_remainder_bound": False,
            "spatial_symbol_to_frequency_squared_proved": False,
            "physical_sector_identified": False,
            "photon_attachment": False,
            "polarization_attachment": False,
            "physical_clock_or_units": False,
            "preferred_frame_or_boost_map": False,
            "cosmological_propagation_map": False,
            "SME_or_other_EFT_normalization": False,
            "line_of_sight_delay_integral": False,
            "wave_packet_dynamics_proved": False,
            "time_evolution_proved": False,
            "trajectory_or_time_of_flight_map": False,
            "likelihood_or_nuisance_model": False,
            "statement": (
                "the inverse-speed expression is the reciprocal of the radial "
                "formal group velocity; the nonzero transverse drift means this "
                "scalar excess is not a trajectory or time-of-flight map"
            ),
        },
        "interpretation boundary drift",
    )
    exposure = receipt.get("exposure_boundary", {})
    check(
        exposure
        == {
            "target_values_read": False,
            "comparison_data_read": False,
            "public_measurement_read": False,
            "comparison_permitted": False,
            "comparison_inputs": [],
            "inputs_read": [
                "target-free frozen FZ-12 source receipt",
                "FZ-12-only custody projection without measurement values",
            ],
            "no_data_inputs_verified": True,
            "comparison_state": "INELIGIBLE_DERIVED_FORMAL_MAP_ONLY",
        },
        "exposure boundary drift",
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--custody-projection", type=Path, default=DEFAULT_CUSTODY)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        verify_source(args.source)
        verify_custody(args.custody_projection)
        verify_receipt(args.receipt)
    except (OSError, VerificationError, UnicodeError) as error:
        print(f"FZ12_OBSERVATION_MAP_INDEPENDENT_VERIFY_FAIL: {error}", file=sys.stderr)
        return 1
    print("FZ12_OBSERVATION_MAP_INDEPENDENT_VERIFY_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
