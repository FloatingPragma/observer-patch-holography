#!/usr/bin/env python3
"""Derive the exact no-data observation map of the frozen FZ-12 ray.

FZ-12 freezes a conditional spatial symbol

    Lambda = k^2 + C4 k^4 + (B0 + B6 I6(n)) k^6.

This certificate performs one algebraic operation only.  It introduces the
dimensionless variables ``q = a k``, ``lambda_hat = a^2 Lambda``, and
``omega_hat = a omega``.  Conditional on the additional identification
``lambda_hat = omega_hat^2`` and the positive formal branch, it computes the
truncated frequency, radial group velocity, inverse radial speed, and
transverse spherical-gradient coefficient in exact rational formal-series
rings.  It supplies no convergence or remainder bound and makes no
physical-sector, photon, clock, frame, cosmological, SME, or observation-data
identification.

The producer reads the immutable FZ-12 source receipt and a canonical FZ-12-only
custody projection.  It refuses schema, status, fixed byte digest, coefficient,
custody, or exposure drift.  It does not parse the multi-campaign register and
never writes either input.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
SOURCE_RECEIPT_PATH = HERE / "runtime" / "seam_current_edge_prediction_receipt.json"
CUSTODY_PROJECTION_PATH = HERE / "runtime" / "fz12_custody_projection.json"
RECEIPT_PATH = HERE / "runtime" / "fz12_observation_map_receipt.json"

SCHEMA = "oph.fz12.formal_observation_map.v1"
STATUS = (
    "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
    "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
)
SOURCE_SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
SOURCE_STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)
CUSTODY_SCHEMA = "oph.fz12.custody_projection.v1"
CUSTODY_STATUS = "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA"
EXPECTED_SOURCE_BYTES = 9296
EXPECTED_SOURCE_SHA256 = (
    "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
)
EXPECTED_CUSTODY_BYTES = 3624
EXPECTED_CUSTODY_SHA256 = (
    "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
)
EXPECTED_CUSTODY_SELF_SHA256 = (
    "sha256:334d9877023d4eaa2ac64bca51d9e2193f0eaad81afd664b22cada2b890a8527"
)

# A polynomial is sparse in (power of k, power of the formal scalar I6).
Monomial = tuple[int, int]
Polynomial = dict[Monomial, Fraction]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def add(*polynomials: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return {monomial: value for monomial, value in out.items() if value}


def scale(polynomial: Polynomial, factor: Fraction) -> Polynomial:
    return {
        monomial: factor * coefficient
        for monomial, coefficient in polynomial.items()
        if factor * coefficient
    }


def multiply(left: Polynomial, right: Polynomial, *, max_k: int) -> Polynomial:
    out: Polynomial = {}
    for (left_k, left_i6), left_coefficient in left.items():
        for (right_k, right_i6), right_coefficient in right.items():
            k_power = left_k + right_k
            if k_power > max_k:
                continue
            monomial = (k_power, left_i6 + right_i6)
            out[monomial] = (
                out.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
            )
    return {monomial: value for monomial, value in out.items() if value}


def truncate(polynomial: Polynomial, *, max_k: int) -> Polynomial:
    return {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if monomial[0] <= max_k and coefficient
    }


def radial_derivative(polynomial: Polynomial) -> Polynomial:
    return {
        (k_power - 1, i6_power): Fraction(k_power) * coefficient
        for (k_power, i6_power), coefficient in polynomial.items()
        if k_power > 0 and coefficient
    }


def transverse_gradient_prefactor(polynomial: Polynomial) -> Polynomial:
    """Return (1/k) d/dI6, the coefficient multiplying grad_S I6."""

    return {
        (k_power - 1, i6_power - 1): Fraction(i6_power) * coefficient
        for (k_power, i6_power), coefficient in polynomial.items()
        if k_power > 0 and i6_power > 0 and coefficient
    }


def coefficient(polynomial: Polynomial, k_power: int, i6_power: int = 0) -> Fraction:
    return polynomial.get((k_power, i6_power), Fraction(0))


def formal_map(c4: Fraction, b0: Fraction, b6: Fraction) -> dict[str, Any]:
    """Build and internally verify the formal map through the frozen orders."""

    dispersion: Polynomial = {
        (2, 0): Fraction(1),
        (4, 0): c4,
        (6, 0): b0,
        (6, 1): b6,
    }
    omega: Polynomial = {
        (1, 0): Fraction(1),
        (3, 0): c4 / 2,
        (5, 0): b0 / 2 - c4 * c4 / 8,
        (5, 1): b6 / 2,
    }
    require(
        multiply(omega, omega, max_k=7) == dispersion,
        "formal positive-branch square-root identity modulo q^8 failed",
    )

    radial = truncate(radial_derivative(omega), max_k=4)
    delta = add(radial, {(0, 0): Fraction(-1)})
    inverse_radial = truncate(
        add(
            {(0, 0): Fraction(1)},
            scale(delta, Fraction(-1)),
            multiply(delta, delta, max_k=4),
        ),
        max_k=4,
    )
    require(
        multiply(radial, inverse_radial, max_k=5) == {(0, 0): Fraction(1)},
        "formal inverse radial-speed identity modulo q^6 failed",
    )
    inverse_speed_excess = add(inverse_radial, {(0, 0): Fraction(-1)})

    transverse = truncate(transverse_gradient_prefactor(omega), max_k=4)
    require(
        transverse == {(4, 0): b6 / 2},
        "formal transverse spherical-gradient identity failed",
    )

    expected = {
        "omega": {
            "k": "1",
            "a2_k3": "-1/40",
            "a4_k5_isotropic": "19/67200",
            "a4_k5_I6": "-1/25200",
        },
        "radial_group_velocity": {
            "constant": "1",
            "a2_k2": "-3/40",
            "a4_k4_isotropic": "19/13440",
            "a4_k4_I6": "-1/5040",
        },
        "inverse_radial_speed": {
            "constant": "1",
            "a2_k2": "3/40",
            "a4_k4_isotropic": "283/67200",
            "a4_k4_I6": "1/5040",
        },
        "formal_inverse_radial_speed_excess": {
            "a2_k2": "3/40",
            "a4_k4_isotropic": "283/67200",
            "a4_k4_I6": "1/5040",
        },
        "transverse_group_velocity": {
            "a4_k4_grad_S_I6": "-1/25200",
        },
    }
    actual = {
        "omega": {
            "k": str(coefficient(omega, 1)),
            "a2_k3": str(coefficient(omega, 3)),
            "a4_k5_isotropic": str(coefficient(omega, 5)),
            "a4_k5_I6": str(coefficient(omega, 5, 1)),
        },
        "radial_group_velocity": {
            "constant": str(coefficient(radial, 0)),
            "a2_k2": str(coefficient(radial, 2)),
            "a4_k4_isotropic": str(coefficient(radial, 4)),
            "a4_k4_I6": str(coefficient(radial, 4, 1)),
        },
        "inverse_radial_speed": {
            "constant": str(coefficient(inverse_radial, 0)),
            "a2_k2": str(coefficient(inverse_radial, 2)),
            "a4_k4_isotropic": str(coefficient(inverse_radial, 4)),
            "a4_k4_I6": str(coefficient(inverse_radial, 4, 1)),
        },
        "formal_inverse_radial_speed_excess": {
            "a2_k2": str(coefficient(inverse_speed_excess, 2)),
            "a4_k4_isotropic": str(coefficient(inverse_speed_excess, 4)),
            "a4_k4_I6": str(coefficient(inverse_speed_excess, 4, 1)),
        },
        "transverse_group_velocity": {
            "a4_k4_grad_S_I6": str(coefficient(transverse, 4)),
        },
    }
    require(actual == expected, "FZ-12 formal observation-map coefficient drift")

    c4_squared = c4 * c4
    scale_free = {
        "omega_k5_isotropic_over_C4_squared": str(coefficient(omega, 5) / c4_squared),
        "omega_k5_I6_over_C4_squared": str(coefficient(omega, 5, 1) / c4_squared),
        "radial_k4_isotropic_over_C4_squared": str(coefficient(radial, 4) / c4_squared),
        "radial_k4_I6_over_C4_squared": str(coefficient(radial, 4, 1) / c4_squared),
        "inverse_speed_k4_isotropic_over_C4_squared": str(
            coefficient(inverse_radial, 4) / c4_squared
        ),
        "inverse_speed_k4_I6_over_C4_squared": str(
            coefficient(inverse_radial, 4, 1) / c4_squared
        ),
        "transverse_grad_S_I6_over_C4_squared": str(
            coefficient(transverse, 4) / c4_squared
        ),
    }
    require(
        scale_free
        == {
            "omega_k5_isotropic_over_C4_squared": "19/168",
            "omega_k5_I6_over_C4_squared": "-1/63",
            "radial_k4_isotropic_over_C4_squared": "95/168",
            "radial_k4_I6_over_C4_squared": "-5/63",
            "inverse_speed_k4_isotropic_over_C4_squared": "283/168",
            "inverse_speed_k4_I6_over_C4_squared": "5/63",
            "transverse_grad_S_I6_over_C4_squared": "-1/63",
        },
        "FZ-12 scale-free formal observation-map drift",
    )
    return {
        "coefficients_in_a_units": actual,
        "scale_free_coefficients": scale_free,
        "machine_checks": {
            "omega_hat_squared_equals_lambda_hat_mod_q8": True,
            "radial_times_inverse_radial_equals_one_mod_q6": True,
            "formal_inverse_radial_speed_excess_is_inverse_minus_one": True,
            "transverse_coefficient_is_one_over_k_times_spherical_gradient": True,
        },
    }


def load_source_receipt(path: Path | None = None) -> tuple[bytes, dict[str, Any]]:
    if path is None:
        path = SOURCE_RECEIPT_PATH
    raw = path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid FZ-12 source receipt JSON") from error
    require(isinstance(receipt, dict), "FZ-12 source receipt root is not an object")
    require(
        raw == base.canonical_json_bytes(receipt), "noncanonical FZ-12 source receipt"
    )
    require(
        base.tagged_sha256(raw) == f"sha256:{EXPECTED_SOURCE_SHA256}",
        "FZ-12 source fixed byte digest drift",
    )
    require(len(raw) == EXPECTED_SOURCE_BYTES, "FZ-12 source byte-count drift")
    require(receipt.get("schema") == SOURCE_SCHEMA, "FZ-12 source schema drift")
    require(receipt.get("status") == SOURCE_STATUS, "FZ-12 source status drift")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    require(
        receipt.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "FZ-12 source self-digest drift",
    )
    candidate = receipt.get("conditional_physical_candidate", {})
    require(
        candidate.get("coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "FZ-12 source coefficient drift",
    )
    require(
        candidate.get("scale_free_relations")
        == {
            "B0_over_C4_squared": "10/21",
            "B6_over_C4_squared": "-2/63",
            "B6_over_B0": "-1/15",
        },
        "FZ-12 source scale-free ray drift",
    )
    require(
        "not automatically a frequency squared" in candidate.get("symbol_name", ""),
        "FZ-12 source frequency boundary drift",
    )
    boundary = receipt.get("exposure_and_custody_boundary", {})
    for key in (
        "target_values_read",
        "comparison_data_read",
        "public_measurement_read",
        "comparison_permitted",
    ):
        require(boundary.get(key) is False, f"FZ-12 source exposure drift: {key}")
    require(
        receipt.get("promotion_gates", {}).get("all_discharged") is False,
        "FZ-12 source physical gates unexpectedly discharged",
    )
    return raw, receipt


def load_custody_projection(
    source_raw: bytes, path: Path | None = None
) -> tuple[bytes, dict[str, Any]]:
    if path is None:
        path = CUSTODY_PROJECTION_PATH
    raw = path.read_bytes()
    try:
        projection = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid FZ-12 custody projection JSON") from error
    require(isinstance(projection, dict), "FZ-12 custody projection is not an object")
    require(
        raw == base.canonical_json_bytes(projection),
        "noncanonical FZ-12 custody projection",
    )
    require(
        base.tagged_sha256(raw) == f"sha256:{EXPECTED_CUSTODY_SHA256}",
        "FZ-12 custody fixed byte digest drift",
    )
    require(len(raw) == EXPECTED_CUSTODY_BYTES, "FZ-12 custody byte-count drift")
    require(projection.get("schema") == CUSTODY_SCHEMA, "FZ-12 custody schema drift")
    require(projection.get("status") == CUSTODY_STATUS, "FZ-12 custody status drift")
    body = {
        key: value for key, value in projection.items() if key != "projection_sha256"
    }
    require(
        projection.get("projection_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "FZ-12 custody self-digest drift",
    )
    require(
        projection.get("projection_sha256") == EXPECTED_CUSTODY_SELF_SHA256,
        "FZ-12 custody fixed self-digest drift",
    )
    scope = projection.get("projection_scope", {})
    require(scope.get("register_row") == "FZ-12", "FZ-12 custody row drift")
    for key in (
        "includes_measurement_values",
        "includes_comparison_values",
        "includes_likelihood_values",
        "includes_other_campaign_rows",
    ):
        require(scope.get(key) is False, f"FZ-12 custody scope drift: {key}")
    source_pin = projection.get("source_receipt", {})
    require(
        source_pin.get("bytes") == EXPECTED_SOURCE_BYTES, "custody source size drift"
    )
    require(
        source_pin.get("sha256") == EXPECTED_SOURCE_SHA256,
        "custody source digest drift",
    )
    require(
        source_pin.get("sha256")
        == base.tagged_sha256(source_raw).removeprefix("sha256:"),
        "custody projection no longer binds source bytes",
    )
    frozen_row = projection.get("frozen_row_projection", {})
    require(frozen_row.get("id") == "FZ-12", "custody frozen row id drift")
    require(
        frozen_row.get("content_sha256") == EXPECTED_SOURCE_SHA256,
        "custody frozen row content pin drift",
    )
    contract = projection.get("external_custody_contract", {})
    require(contract.get("rows") == ["FZ-12"], "custody contract row drift")
    require(
        contract.get("in_repo_artifact_sha256", {}).get(
            "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json"
        )
        == EXPECTED_SOURCE_SHA256,
        "custody contract source pin drift",
    )
    return raw, projection


def build_receipt() -> dict[str, Any]:
    source_raw, source = load_source_receipt()
    custody_raw, custody = load_custody_projection(source_raw)
    coefficients = source["conditional_physical_candidate"]["coefficients"]
    result = formal_map(
        Fraction(coefficients["C4_over_a2"]),
        Fraction(coefficients["B0_over_a4"]),
        Fraction(coefficients["B6_over_a4"]),
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 666,
        "certificate_scope": {
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
        "ancestry": {
            "source_receipt": {
                "path": "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
                "bytes": len(source_raw),
                "sha256": base.tagged_sha256(source_raw),
                "schema": SOURCE_SCHEMA,
                "status": SOURCE_STATUS,
                "receipt_sha256": source["receipt_sha256"],
            },
            "custody_projection": {
                "path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
                "bytes": len(custody_raw),
                "sha256": base.tagged_sha256(custody_raw),
                "schema": CUSTODY_SCHEMA,
                "status": CUSTODY_STATUS,
                "projection_sha256": custody["projection_sha256"],
                "row": custody["frozen_row_projection"]["id"],
                "canonical_full_row_sha256": custody["frozen_row_projection"][
                    "canonical_full_row_sha256"
                ],
                "source_commit": custody["external_custody_contract"]["source_commit"],
                "custody_commit": custody["external_custody_contract"][
                    "custody_commit"
                ],
                "decision_rule_custody_commit": custody["external_custody_contract"][
                    "decision_rule_custody_commit"
                ],
            },
            "frozen_artifacts_modified": False,
        },
        "formal_input": {
            "dispersion_symbol": "Lambda = k^2 + C4 k^4 + (B0 + B6 I6(n)) k^6",
            "dimensionless_normalization": {
                "q": "q = a k",
                "lambda_hat": "lambda_hat = a^2 Lambda",
                "omega_hat": "omega_hat = a omega",
                "normalized_symbol": (
                    "lambda_hat = q^2 - q^4/20 + (1/840 - I6/12600) q^6 + O(q^8)"
                ),
            },
            "frozen_coefficients": coefficients,
            "frozen_scale_free_relations": source["conditional_physical_candidate"][
                "scale_free_relations"
            ],
            "extra_map_premise": (
                "lambda_hat is identified with omega_hat squared and omega_hat is "
                "the positive formal branch with leading term q"
            ),
            "extra_map_premise_discharged_here": False,
        },
        "exact_formal_result": {
            "ring_for_square_root": "Q[I6][q] modulo q^8",
            "ring_for_velocity_and_inverse": "Q[I6][q] modulo q^6",
            "frequency": (
                "omega = k - (a^2/40) k^3 + a^4 (19/67200 - I6/25200) k^5 + O(a^6 k^7)"
            ),
            "radial_group_velocity": (
                "d omega/dk = 1 - (3a^2/40) k^2 + a^4 "
                "(19/13440 - I6/5040) k^4 + O(a^6 k^6)"
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
            **result,
        },
        "interpretation_boundary": {
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
        "exposure_boundary": {
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
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(
            "invalid committed FZ-12 observation-map receipt"
        ) from error
    require(isinstance(committed, dict), "FZ-12 observation-map root is not an object")
    require(
        raw == base.canonical_json_bytes(committed),
        "noncanonical FZ-12 observation-map receipt",
    )
    require(committed.get("schema") == SCHEMA, "FZ-12 observation-map schema drift")
    require(committed.get("status") == STATUS, "FZ-12 observation-map status drift")
    body = {key: value for key, value in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "FZ-12 observation-map self-digest drift",
    )
    rebuilt = build_receipt()
    require(
        raw == base.canonical_json_bytes(rebuilt),
        "FZ-12 observation-map ancestry or result drift",
    )
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write and args.verify:
        parser.error("choose --write or --verify")
    receipt = verify_committed_receipt() if args.verify else build_receipt()
    if args.write:
        RECEIPT_PATH.parent.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["exact_formal_result"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
