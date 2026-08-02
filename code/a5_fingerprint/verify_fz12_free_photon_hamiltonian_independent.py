#!/usr/bin/env python3
"""Separately replay the conditional FZ-12 transverse-oscillator packet.

The verifier does not import the producer.  It separately reconstructs the thirty seam
directions over ``Q(sqrt(5))``, checks the basis-free transverse projector,
repeats the dimensionless exact-symbol wave-packet calculation, and verifies
the fixed Lean source and parent receipts byte for byte.

The replay is separate code, not an independent scientific implementation.
The verified object remains a mathematical transverse-oscillator candidate.
This file does not attach physical position, time, fields, electrons,
interactions, frames, sources, detectors, or comparison data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
DEFAULT_SOURCE = RUNTIME / "seam_current_edge_prediction_receipt.json"
DEFAULT_REMAINDER = RUNTIME / "fz12_full_symbol_remainder_receipt.json"
DEFAULT_FREE_PHOTON_LEAN = (
    REPO_ROOT / "Lean" / "Screen" / "SeamCurrentFreePhotonLift.lean"
)
DEFAULT_RECEIPT = RUNTIME / "fz12_free_photon_hamiltonian_receipt.json"

SCHEMA = "oph.fz12.transverse_oscillator_candidate.v1"
STATUS = (
    "EXACT_BASIS_FREE_RANK_TWO_TRANSVERSE_OSCILLATOR_CANDIDATE_AND_"
    "AUXILIARY_WAVEPACKET__PHYSICAL_FIELD_ATTACHMENTS_OPEN"
)

PARENT_SPECS = (
    {
        "path": "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
        "bytes": 9296,
        "sha256": (
            "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
        ),
        "schema": "oph.seam_current_edge_prediction_candidate.v1",
        "status": (
            "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
            "PHYSICAL_PRODUCER_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "frozen target-free FZ-12 source ray",
    },
    {
        "path": "code/a5_fingerprint/runtime/fz12_full_symbol_remainder_receipt.json",
        "bytes": 5367,
        "sha256": (
            "ce4052f633c891515ebc319e9d7ff2bc0044bb8bb9e76a890e8e22c8e882dcfa"
        ),
        "schema": "oph.fz12.full_symbol_remainder.v1",
        "status": (
            "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
            "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "exact full spatial symbol and q<=1 remainder",
    },
)

LEAN_BYTES = 13119
LEAN_SHA256 = "5c816286127815938c7c391d1764a67c899ce544b0b2bbb6768388de4fb5699c"


class VerificationError(ValueError):
    """A transverse-oscillator packet contract failed closed."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


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


def load_parent(path: Path, spec: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError(f"invalid parent JSON: {path.name}") from error
    check(isinstance(value, dict), f"parent root is not an object: {path.name}")
    check(raw == canonical_json_bytes(value), f"noncanonical parent: {path.name}")
    check(len(raw) == spec["bytes"], f"parent byte-count drift: {path.name}")
    check(sha256(raw) == spec["sha256"], f"parent raw hash drift: {path.name}")
    check(value.get("schema") == spec["schema"], f"parent schema drift: {path.name}")
    check(value.get("status") == spec["status"], f"parent status drift: {path.name}")
    check(
        value.get(spec["self_field"]) == self_digest(value, spec["self_field"]),
        f"parent self-digest drift: {path.name}",
    )
    return raw, value


def load_lean(path: Path) -> bytes:
    raw = path.read_bytes()
    check(len(raw) == LEAN_BYTES, "transverse-oscillator Lean byte-count drift")
    check(sha256(raw) == LEAN_SHA256, "transverse-oscillator Lean raw hash drift")
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise VerificationError(
            "transverse-oscillator Lean source is not UTF-8"
        ) from error
    required = (
        "theorem transversePolarization_finrank",
        "theorem transverseProjector_idempotent",
        "theorem photonSpatialAction_commutes_projector",
        "theorem photonModeFrequency_sq",
        "theorem photonModeFrequency_zero",
        "theorem photonMode_second_order",
        "theorem photonModeEnergy_firstVariation_generator_zero",
        "does not identify\nthe carrier completion with physical position",
        "#print axioms photonModeEnergy_firstVariation_generator_zero",
    )
    check(
        all(fragment in source for fragment in required),
        "transverse-oscillator Lean claim drift",
    )
    check("sorry" not in source and "admit" not in source, "Lean placeholder found")
    return raw


# Exact carrier arithmetic.  Scalars are represented as a + b sqrt(5).
Q5 = tuple[Fraction, Fraction]
Vec3 = tuple[Q5, Q5, Q5]
ZERO: Q5 = (Fraction(0), Fraction(0))
ONE: Q5 = (Fraction(1), Fraction(0))
PHI: Q5 = (Fraction(1, 2), Fraction(1, 2))


def q5(value: int) -> Q5:
    return (Fraction(value), Fraction(0))


def q5_add(left: Q5, right: Q5) -> Q5:
    return (left[0] + right[0], left[1] + right[1])


def q5_sub(left: Q5, right: Q5) -> Q5:
    return (left[0] - right[0], left[1] - right[1])


def q5_mul(left: Q5, right: Q5) -> Q5:
    return (
        left[0] * right[0] + 5 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def q5_scale(value: Q5, scalar: Fraction) -> Q5:
    return (scalar * value[0], scalar * value[1])


def q5_div(left: Q5, right: Q5) -> Q5:
    denominator = right[0] * right[0] - 5 * right[1] * right[1]
    check(denominator != 0, "division by zero in Q(sqrt(5))")
    inverse = (right[0] / denominator, -right[1] / denominator)
    return q5_mul(left, inverse)


def q5_float(value: Q5) -> float:
    return float(value[0]) + float(value[1]) * math.sqrt(5.0)


def vec_sub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(  # type: ignore[return-value]
        q5_sub(x, y) for x, y in zip(left, right, strict=True)
    )


def vec_scale(value: Vec3, scalar: Fraction) -> Vec3:
    return tuple(q5_scale(x, scalar) for x in value)  # type: ignore[return-value]


def dot(left: Vec3, right: Vec3) -> Q5:
    result = ZERO
    for x, y in zip(left, right, strict=True):
        result = q5_add(result, q5_mul(x, y))
    return result


def port_vectors() -> tuple[Vec3, ...]:
    result = []
    for sign_one in (1, -1):
        for sign_two in (1, -1):
            a = q5_scale(ONE, Fraction(sign_one))
            b = q5_scale(PHI, Fraction(sign_two))
            result.extend(((ZERO, a, b), (a, b, ZERO), (b, ZERO, a)))
    check(len(result) == 12, "port census drift")
    return tuple(result)


def unit_edge_directions() -> np.ndarray:
    ports = port_vectors()
    norm_squared = q5_add(q5(2), PHI)
    inverse_norm_squared = q5_div(ONE, norm_squared)
    inverse_sqrt_five: Q5 = (Fraction(0), Fraction(1, 5))
    exact = []
    for left in range(12):
        for right in range(left + 1, 12):
            normalized_dot = q5_mul(
                dot(ports[left], ports[right]), inverse_norm_squared
            )
            if normalized_dot == inverse_sqrt_five:
                exact.append(
                    vec_scale(vec_sub(ports[right], ports[left]), Fraction(1, 2))
                )
    exact = tuple(exact)
    check(len(exact) == 30, "seam census drift")
    check(all(dot(row, row) == ONE for row in exact), "nonunit exact seam direction")
    result = np.asarray(
        [[q5_float(component) for component in row] for row in exact],
        dtype=np.float64,
    )
    check(
        float(np.max(np.abs(np.sum(result**2, axis=1) - 1.0))) <= 1.0e-12,
        "floating seam normalization drift",
    )
    return result


def projector_contract() -> dict[str, Any]:
    directions = (
        ("axis_x", np.asarray([1.0, 0.0, 0.0])),
        ("diagonal_yz", np.asarray([0.0, 1.0, 1.0]) / math.sqrt(2.0)),
        ("normalize_123", np.asarray([1.0, 2.0, 3.0]) / math.sqrt(14.0)),
    )
    rows = []
    for label, direction in directions:
        projector = np.eye(3) - np.outer(direction, direction)
        eigenvalues = np.linalg.eigvalsh(projector)
        idempotence_residual = float(
            np.linalg.norm(projector @ projector - projector)
        )
        transversality_residual = float(np.linalg.norm(direction @ projector))
        trace_residual = abs(float(np.trace(projector)) - 2.0)
        spectrum_residual = float(
            np.max(np.abs(eigenvalues - np.asarray([0.0, 1.0, 1.0])))
        )
        numerical_rank = int(np.linalg.matrix_rank(projector, tol=1.0e-12))
        check(idempotence_residual <= 1.0e-10, "projector idempotence drift")
        check(transversality_residual <= 1.0e-10, "projector transverse drift")
        check(trace_residual <= 1.0e-10, "projector trace drift")
        check(spectrum_residual <= 1.0e-10, "projector spectrum drift")
        check(numerical_rank == 2, "projector rank drift")
        rows.append(
            {
                "momentum_direction_fixture": label,
                "expected_spectrum": ["0", "1", "1"],
                "expected_trace": "2",
                "numerical_rank": numerical_rank,
                "idempotence_check_passed": True,
                "transversality_check_passed": True,
                "trace_check_passed": True,
                "spectrum_check_passed": True,
            }
        )
    return {
        "formula": "P_T(k) = I - khat tensor khat",
        "global_basis_selected": False,
        "numerical_check_tolerance": "1e-10",
        "cross_platform_serialization": (
            "theorem-exact spectrum and tolerance-pass booleans; raw ULPs omitted"
        ),
        "lean_rank_theorem": "transversePolarization_finrank",
        "lean_projector_theorems": [
            "transverseProjector_dot_zero",
            "transverseProjector_idempotent",
            "photonSpatialAction_commutes_projector",
        ],
        "numerical_rows": rows,
    }


def lambda_hat(q_values: np.ndarray | float, projections: np.ndarray) -> np.ndarray:
    q = np.asarray(q_values, dtype=np.float64)
    return 6.0 * np.mean(1.0 - np.cos(q[..., np.newaxis] * projections), axis=-1)


def lambda_hat_derivative(
    q_values: np.ndarray | float, projections: np.ndarray
) -> np.ndarray:
    q = np.asarray(q_values, dtype=np.float64)
    return 6.0 * np.mean(
        projections * np.sin(q[..., np.newaxis] * projections), axis=-1
    )


def wave_packet_contract(directions: np.ndarray) -> dict[str, Any]:
    direction = np.asarray([1.0, 2.0, 3.0], dtype=np.float64)
    direction /= np.linalg.norm(direction)
    projections = directions @ direction

    sample_count = 4096
    domain_length = 400.0
    spacing = domain_length / sample_count
    q_modes = 2.0 * math.pi * np.fft.fftfreq(sample_count, d=spacing)
    q_center = 0.5
    q_width = 0.05
    x_center = 100.0
    auxiliary_time = 30.0

    spectral = np.exp(-((q_modes - q_center) ** 2) / (4.0 * q_width**2))
    spectral = spectral.astype(np.complex128) * np.exp(-1j * q_modes * x_center)
    spectral[(q_modes <= 0.0) | (q_modes > 1.0)] = 0.0

    symbol = lambda_hat(q_modes, projections)
    frequency = np.sqrt(np.maximum(symbol, 0.0))
    derivative = lambda_hat_derivative(q_modes, projections)
    group_velocity = np.zeros_like(q_modes)
    support = symbol > 1.0e-14
    group_velocity[support] = derivative[support] / (2.0 * frequency[support])

    initial = np.fft.ifft(spectral)
    evolved = np.fft.ifft(spectral * np.exp(-1j * frequency * auxiliary_time))
    x_grid = np.arange(sample_count, dtype=np.float64) * spacing

    def centroid(field: np.ndarray) -> float:
        density = np.abs(field) ** 2
        return float(np.dot(x_grid, density) / np.sum(density))

    initial_center = centroid(initial)
    evolved_center = centroid(evolved)
    weights = np.abs(spectral) ** 2
    expected_velocity = float(np.dot(weights, group_velocity) / np.sum(weights))
    expected_shift = expected_velocity * auxiliary_time
    observed_shift = evolved_center - initial_center
    residual = observed_shift - expected_shift

    seed = np.asarray([1.0, 0.0, 0.0])
    polarization_one = seed - direction * float(np.dot(direction, seed))
    polarization_one /= np.linalg.norm(polarization_one)
    polarization_two = np.cross(direction, polarization_one)
    polarization_two /= np.linalg.norm(polarization_two)
    initial_one = initial[:, np.newaxis] * polarization_one
    initial_two = initial[:, np.newaxis] * polarization_two
    evolved_one = evolved[:, np.newaxis] * polarization_one
    evolved_two = evolved[:, np.newaxis] * polarization_two
    norm_residual = abs(
        float(np.vdot(evolved_one, evolved_one).real)
        / float(np.vdot(initial_one, initial_one).real)
        - float(np.vdot(evolved_two, evolved_two).real)
        / float(np.vdot(initial_two, initial_two).real)
    )
    orthogonality_residual = abs(float(np.dot(polarization_one, polarization_two)))
    transversality_residual = max(
        abs(float(np.dot(direction, polarization_one))),
        abs(float(np.dot(direction, polarization_two))),
    )

    check(abs(initial_center - x_center) <= 1.0e-10, "packet center drift")
    check(abs(residual) <= 1.0e-8, "packet transport drift")
    check(norm_residual <= 1.0e-12, "transverse-degeneracy drift")
    check(orthogonality_residual <= 1.0e-12, "transverse orthogonality drift")
    check(transversality_residual <= 1.0e-12, "transverse constraint drift")
    check(
        np.all(
            (q_modes[np.abs(spectral) > 1.0e-12] > 0.0)
            & (q_modes[np.abs(spectral) > 1.0e-12] <= 1.0)
        ),
        "packet support left q<=1 contract",
    )

    center_symbol = float(lambda_hat(q_center, projections))
    center_derivative = float(lambda_hat_derivative(q_center, projections))
    center_velocity = center_derivative / (2.0 * math.sqrt(center_symbol))
    return {
        "fixture_role": (
            "dimensionless auxiliary transverse-oscillator wave-packet check"
        ),
        "physical_wave_packet_claimed": False,
        "direction": "normalize(1,2,3)",
        "sample_count": sample_count,
        "periodic_domain_length": "400",
        "q_center": "0.5",
        "q_width": "0.05",
        "q_support": "0 < q <= 1",
        "initial_center": "100",
        "auxiliary_time": "30",
        "center_symbol_6_digits": format(center_symbol, ".6g"),
        "center_auxiliary_group_velocity_6_digits": format(center_velocity, ".6g"),
        "spectral_auxiliary_group_velocity_6_digits": format(
            expected_velocity, ".6g"
        ),
        "expected_center_shift_6_digits": format(expected_shift, ".6g"),
        "observed_center_shift_6_digits": format(observed_shift, ".6g"),
        "center_shift_tolerance": "1e-8",
        "center_shift_check_passed": True,
        "center_shift_residual_clamped": "0",
        "two_local_transverse_norm_ratio_tolerance": "1e-12",
        "two_local_transverse_norm_ratio_check_passed": True,
        "local_transverse_orthogonality_check_passed": True,
        "local_transverse_constraint_check_passed": True,
        "ulp_residuals_clamped": True,
        "cross_platform_serialization": (
            "theorem-zero and tolerance-passing ULP residuals serialize as zero; "
            "nonzero diagnostics are rounded to six significant digits"
        ),
        "local_transverse_witness_only": True,
        "global_polarization_basis_claimed": False,
        "chosen_diagnostic_evolution": (
            "spectral phase exp(-i sqrt(lambda_hat(q,n)) tau)"
        ),
        "time_reversed_positive_phase_is_also_solution": True,
        "auxiliary_group_velocity_formula": (
            "d omega_hat/dq = [sum_e mu_e sin(q mu_e)] / "
            "[5*2*sqrt(lambda_hat)]"
        ),
    }


def expected_receipt(
    parents: list[tuple[bytes, dict[str, Any]]], lean_raw: bytes
) -> dict[str, Any]:
    parent_pins = []
    for (raw, value), spec in zip(parents, PARENT_SPECS, strict=True):
        parent_pins.append(
            {
                "path": spec["path"],
                "role": spec["role"],
                "bytes": len(raw),
                "sha256": "sha256:" + sha256(raw),
                "schema": spec["schema"],
                "status": spec["status"],
                "self_digest": value[spec["self_field"]],
            }
        )
    parent_pins.append(
        {
            "path": "Lean/Screen/SeamCurrentFreePhotonLift.lean",
            "role": "kernel-checked basis-free transverse oscillator lift",
            "bytes": len(lean_raw),
            "sha256": "sha256:" + sha256(lean_raw),
            "theorems": [
                "transversePolarization_finrank",
                "transverseProjector_idempotent",
                "photonModeFrequency_sq",
                "photonModeFrequency_zero",
                "photonMode_second_order",
                "photonModeEnergy_firstVariation_generator_zero",
            ],
            "sorry_free": True,
        }
    )
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 667,
        "scope": {
            "type": (
                "conditional basis-free transverse oscillator/free-vector candidate"
            ),
            "declared_oscillator_relation_proved": True,
            "physical_prediction": False,
            "comparison_permitted": False,
            "statement": (
                "the packet proves the rank-two transverse and declared "
                "oscillator consequences of choosing the exact FZ-12 spatial "
                "symbol as a scalar stiffness"
            ),
        },
        "parent_pins": parent_pins,
        "transverse_sector": projector_contract(),
        "hamiltonian_contract": {
            "state": "(A_T, Pi_T) in ker(k dot -) x ker(k dot -)",
            "transverse_dimension_for_nonzero_k": 2,
            "spatial_action": "A_T maps to Lambda_a(k) A_T",
            "generator": "(A_T,Pi_T) maps to (Pi_T,-Lambda_a(k) A_T)",
            "canonical_energy": "H_k=(|Pi_T|^2+Lambda_a(k)|A_T|^2)/2",
            "exact_mode_equation": "A_T''+Lambda_a(k)A_T=0",
            "frequency_branch": "omega_a(k)=sqrt(Lambda_a(k))>=0",
            "exact_frequency_relation": "omega_a(k)^2=Lambda_a(k)",
            "zero_mode": "Lambda_a(0)=omega_a(0)=0",
            "polarization_degeneracy": (
                "Lambda_a(k) is scalar on the complete rank-two transverse fiber"
            ),
            "transversality_preserved": True,
            "declared_energy_first_variation_on_generator_zero": True,
            "trajectory_energy_conservation_proved": False,
            "auxiliary_parameter_only": True,
        },
        "exact_symbol_contract": {
            "lambda_hat": "(1/5) sum_e [1-cos(q (w_e dot n))]",
            "omega_hat": "sqrt(lambda_hat)",
            "direction_count": 30,
            "q_domain_consumed": "0 <= q <= 1",
            "positive_bound_consumed": "lambda_hat >= (19/20)q^2",
            "target_free": True,
        },
        "wave_packet_check": wave_packet_contract(unit_edge_directions()),
        "physical_boundary": {
            "carrier_completion_identified_with_physical_position": False,
            "auxiliary_parameter_identified_with_physical_clock": False,
            "transverse_fiber_identified_with_laboratory_photon": False,
            "gauge_or_gauss_constraint_quotient_proved": False,
            "u1_or_maxwell_field_identification_proved": False,
            "physical_dilation_or_SI_scale_identified": False,
            "carrier_rest_frame_selected": False,
            "observer_frame_or_boost_transport_proved": False,
            "lorentz_covariance_proved": False,
            "continuum_locality_proved": False,
            "physical_causality_proved": False,
            "cofinal_refinement_or_global_gluing_proved": False,
            "reality_pairing_k_and_minus_k_proved": False,
            "positive_residue_quantization_or_hilbert_space_proved": False,
            "canonical_symplectic_form_or_kinetic_normalization_selected": False,
            "trajectory_existence_or_flow_proved": False,
            "trajectory_energy_conservation_proved": False,
            "physical_masslessness_proved": False,
            "physical_wave_packet_or_detector_readout_proved": False,
            "source_selected_hamiltonian": False,
            "spatial_symbol_identified_with_physical_frequency_squared": False,
            "electron_dispersion_selected": False,
            "positron_dispersion_selected": False,
            "electromagnetic_interaction_or_pair_kinematics_proved": False,
            "comparison_nuisance_model_proved": False,
            "statement": (
                "the exact omega-squared relation belongs to the declared "
                "conditional oscillator generator; identifying it with a "
                "physical field requires every open attachment"
            ),
        },
        "exposure_boundary": {
            "public_measurement_read": False,
            "comparison_data_read": False,
            "target_values_read": False,
            "comparison_inputs": [],
            "comparison_permitted": False,
            "comparison_budget_consumed": False,
            "score_emitted": False,
            "evidence_claimed": False,
            "verdict_emitted": False,
        },
        "semantic_receipt_mutation_guards": {
            "guard_role": (
                "resigned semantic receipt mutations exercised by the test suite"
            ),
            "producer_executes_counterfactuals": False,
            "test_suite_exercises_guards": True,
            "guards": [
                "longitudinal rank-three fiber",
                "non-idempotent projector",
                "transverse-direction-split stiffness",
                "negative square-root branch",
                "mutated FZ-12 symbol sign",
                "q support above one",
                "physical clock promotion",
                "electron or interaction promotion",
            ],
        },
    }
    receipt["receipt_sha256"] = self_digest(receipt, "receipt_sha256")
    return receipt


def verify(
    source_path: Path,
    remainder_path: Path,
    lean_path: Path,
    receipt_path: Path,
) -> dict[str, Any]:
    parents = [
        load_parent(source_path, PARENT_SPECS[0]),
        load_parent(remainder_path, PARENT_SPECS[1]),
    ]
    lean_raw = load_lean(lean_path)
    raw = receipt_path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid transverse-oscillator receipt JSON") from error
    check(
        isinstance(receipt, dict),
        "transverse-oscillator receipt root is not an object",
    )
    check(
        raw == canonical_json_bytes(receipt),
        "transverse-oscillator receipt is noncanonical",
    )
    check(receipt.get("schema") == SCHEMA, "transverse-oscillator schema drift")
    check(receipt.get("status") == STATUS, "transverse-oscillator status drift")
    check(
        receipt.get("receipt_sha256") == self_digest(receipt, "receipt_sha256"),
        "transverse-oscillator self-digest drift",
    )
    check(
        receipt == expected_receipt(parents, lean_raw),
        "transverse-oscillator semantic drift",
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--remainder", type=Path, default=DEFAULT_REMAINDER)
    parser.add_argument(
        "--free-photon-lean", type=Path, default=DEFAULT_FREE_PHOTON_LEAN
    )
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        verify(args.source, args.remainder, args.free_photon_lean, args.receipt)
    except (OSError, VerificationError) as error:
        print(f"FZ12_TRANSVERSE_OSCILLATOR_REPLAY_FAIL: {error}", file=sys.stderr)
        return 1
    print("FZ12_TRANSVERSE_OSCILLATOR_SEPARATE_REPLAY_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
