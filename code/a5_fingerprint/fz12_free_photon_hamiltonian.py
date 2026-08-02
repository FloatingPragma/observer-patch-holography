#!/usr/bin/env python3
"""Build the conditional basis-free transverse-oscillator lift of FZ-12.

The frozen FZ-12 source supplies a spatial symbol.  This packet makes one
standard mathematical completion explicit: on the momentum-orthogonal
rank-two fiber, the scalar symbol is used as the stiffness of a declared
oscillator generator.  Its nonnegative frequency has ``omega^2 = Lambda``
and acts equally on both transverse directions.  A deterministic narrow-band
packet checks exact-symbol group transport in dimensionless coordinates.

The construction is not a physical field derivation.  It does not identify
the carrier with position, the auxiliary parameter with a clock, the
rank-two fiber with a photon or Maxwell field, or the scale, frame, source,
interaction, quantization, and readout with physical objects.  No public
comparison input is read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
SOURCE_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
REMAINDER_PATH = RUNTIME / "fz12_full_symbol_remainder_receipt.json"
LEAN_PATH = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentFreePhotonLift.lean"
RECEIPT_PATH = RUNTIME / "fz12_free_photon_hamiltonian_receipt.json"

SCHEMA = "oph.fz12.transverse_oscillator_candidate.v1"
STATUS = (
    "EXACT_BASIS_FREE_RANK_TWO_TRANSVERSE_OSCILLATOR_CANDIDATE_AND_"
    "AUXILIARY_WAVEPACKET__PHYSICAL_FIELD_ATTACHMENTS_OPEN"
)

SOURCE_CONTRACT = {
    "path": SOURCE_PATH,
    "schema": "oph.seam_current_edge_prediction_candidate.v1",
    "status": (
        "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
        "PHYSICAL_PRODUCER_OPEN"
    ),
    "bytes": 9296,
    "sha256": "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915",
    "self_field": "receipt_sha256",
}
REMAINDER_CONTRACT = {
    "path": REMAINDER_PATH,
    "schema": "oph.fz12.full_symbol_remainder.v1",
    "status": (
        "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
        "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
    ),
    "bytes": 5367,
    "sha256": "ce4052f633c891515ebc319e9d7ff2bc0044bb8bb9e76a890e8e22c8e882dcfa",
    "self_field": "receipt_sha256",
}


class FreePhotonError(ValueError):
    """The conditional transverse-oscillator packet failed closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreePhotonError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return base.canonical_json_bytes(value)


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    path = contract["path"]
    raw = path.read_bytes()
    require(len(raw) == contract["bytes"], f"parent byte-count drift: {path.name}")
    require(sha256(raw) == contract["sha256"], f"parent digest drift: {path.name}")
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise FreePhotonError(f"invalid parent JSON: {path.name}") from error
    require(isinstance(value, dict), f"parent root is not an object: {path.name}")
    require(raw == canonical_json_bytes(value), f"noncanonical parent: {path.name}")
    require(value.get("schema") == contract["schema"], f"schema drift: {path.name}")
    require(value.get("status") == contract["status"], f"status drift: {path.name}")
    self_field = contract["self_field"]
    body = {key: item for key, item in value.items() if key != self_field}
    require(
        value.get(self_field) == base.tagged_sha256(canonical_json_bytes(body)),
        f"parent self-digest drift: {path.name}",
    )
    return raw, value


def load_lean_proof() -> bytes:
    raw = LEAN_PATH.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise FreePhotonError("non-UTF-8 transverse-oscillator Lean proof") from error
    normalized = re.sub(r"\s+", " ", text)
    required = (
        (
            "theorem transversePolarization_finrank {k : FreePhotonVec3} "
            "(hk : k \u2260 0) : "
            "Module.finrank \u211d (TransversePolarization k) = 2 := by"
        ),
        (
            "noncomputable def transverseProjector (k v : FreePhotonVec3) : "
            "FreePhotonVec3 :="
        ),
        "theorem photonModeFrequency_sq (a : \u211d) (k : FreePhotonVec3) :",
        "theorem photonModeEnergy_firstVariation_generator_zero",
        "does not identify the carrier completion with physical position",
    )
    for fragment in required:
        require(
            fragment in normalized,
            "claim-bearing transverse-oscillator Lean proof drift",
        )
    require("sorry" not in text and "admit" not in text, "placeholder in Lean proof")
    return raw


def q5_dot(left: tuple, right: tuple) -> tuple:
    total = base.ZERO
    for x, y in zip(left, right, strict=True):
        total = base.q5_add(total, base.q5_mul(x, y))
    return total


def q5_sub(left: tuple, right: tuple) -> tuple:
    return tuple(base.q5_sub(x, y) for x, y in zip(left, right, strict=True))


def q5_float(value: tuple) -> float:
    return float(value[0]) + float(value[1]) * math.sqrt(5.0)


def unit_edge_directions() -> np.ndarray:
    """Reconstruct the thirty exact edge directions, converted to floats."""

    vertices = base.cartesian_vertices()
    inverse_norm = base.q5_div(base.ONE, base.NORM_SQ)
    differences = []
    for left in range(12):
        for right in range(left + 1, 12):
            normalized_dot = base.q5_mul(
                q5_dot(vertices[left], vertices[right]), inverse_norm
            )
            if normalized_dot == base.INV_SQRT5:
                differences.append(q5_sub(vertices[right], vertices[left]))
    require(len(differences) == 30, "edge-direction census drift")
    directions = np.asarray(
        [[q5_float(component) / 2.0 for component in row] for row in differences],
        dtype=np.float64,
    )
    norm_residual = float(np.max(np.abs(np.sum(directions**2, axis=1) - 1.0)))
    require(norm_residual <= 1.0e-12, "edge directions are not unit normalized")
    return directions


def lambda_hat(q: np.ndarray | float, projections: np.ndarray) -> np.ndarray:
    values = np.asarray(q, dtype=np.float64)
    return 6.0 * np.mean(
        1.0 - np.cos(values[..., np.newaxis] * projections), axis=-1
    )


def lambda_hat_derivative(
    q: np.ndarray | float, projections: np.ndarray
) -> np.ndarray:
    values = np.asarray(q, dtype=np.float64)
    return 6.0 * np.mean(
        projections * np.sin(values[..., np.newaxis] * projections), axis=-1
    )


def projector_checks() -> dict[str, Any]:
    """Tolerance checks for the exact basis-free projector theorem."""

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
        require(idempotence_residual <= 1.0e-10, "projector idempotence drift")
        require(transversality_residual <= 1.0e-10, "projector transverse drift")
        require(trace_residual <= 1.0e-10, "projector trace drift")
        require(spectrum_residual <= 1.0e-10, "projector spectrum drift")
        require(numerical_rank == 2, "projector rank drift")
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


def wave_packet_check(directions: np.ndarray) -> dict[str, Any]:
    """Propagate an auxiliary packet by the declared positive symbol root."""

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
    supported = symbol > 1.0e-14
    group_velocity[supported] = derivative[supported] / (2.0 * frequency[supported])

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
    shift_residual = observed_shift - expected_shift

    seed = np.asarray([1.0, 0.0, 0.0])
    polarization_one_vector = seed - direction * float(np.dot(direction, seed))
    polarization_one_vector /= np.linalg.norm(polarization_one_vector)
    polarization_two_vector = np.cross(direction, polarization_one_vector)
    polarization_two_vector /= np.linalg.norm(polarization_two_vector)
    polarization_one = initial[:, np.newaxis] * polarization_one_vector
    polarization_two = initial[:, np.newaxis] * polarization_two_vector
    evolved_one = evolved[:, np.newaxis] * polarization_one_vector
    evolved_two = evolved[:, np.newaxis] * polarization_two_vector
    norm_one = float(np.vdot(evolved_one, evolved_one).real)
    norm_two = float(np.vdot(evolved_two, evolved_two).real)
    input_norm_one = float(np.vdot(polarization_one, polarization_one).real)
    input_norm_two = float(np.vdot(polarization_two, polarization_two).real)
    polarization_norm_ratio_residual = abs(
        norm_one / input_norm_one - norm_two / input_norm_two
    )
    polarization_orthogonality_residual = abs(
        float(np.dot(polarization_one_vector, polarization_two_vector))
    )
    polarization_transversality_residual = max(
        abs(float(np.dot(direction, polarization_one_vector))),
        abs(float(np.dot(direction, polarization_two_vector))),
    )

    require(abs(initial_center - x_center) <= 1.0e-10, "initial packet center drift")
    require(abs(shift_residual) <= 1.0e-8, "wave-packet group transport drift")
    require(
        polarization_norm_ratio_residual <= 1.0e-12,
        "polarization-degenerate evolution drift",
    )
    require(
        polarization_orthogonality_residual <= 1.0e-12
        and polarization_transversality_residual <= 1.0e-12,
        "local transverse-polarization witness drift",
    )
    require(
        np.all((q_modes[np.abs(spectral) > 1.0e-12] > 0.0)
               & (q_modes[np.abs(spectral) > 1.0e-12] <= 1.0)),
        "wave packet left the certified q domain",
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


def parent_pin(
    raw: bytes, payload: dict[str, Any], contract: dict[str, Any], role: str
) -> dict[str, Any]:
    return {
        "path": contract["path"].relative_to(REPO_ROOT).as_posix(),
        "role": role,
        "bytes": len(raw),
        "sha256": f"sha256:{sha256(raw)}",
        "schema": contract["schema"],
        "status": contract["status"],
        "self_digest": payload[contract["self_field"]],
    }


def build_receipt() -> dict[str, Any]:
    source_raw, source = load_parent(SOURCE_CONTRACT)
    remainder_raw, remainder = load_parent(REMAINDER_CONTRACT)
    lean_raw = load_lean_proof()
    require(
        source["exposure_and_custody_boundary"]["public_measurement_read"] is False,
        "source exposure boundary drift",
    )
    require(
        remainder["physical_boundary"]["spatial_symbol_to_frequency_squared_proved"]
        is False,
        "parent spatial symbol was physically promoted",
    )
    require(
        remainder["geometry_contract"]["q_domain"]
        == {"minimum": "0", "maximum": "1", "inclusive": True},
        "parent q domain drift",
    )
    directions = unit_edge_directions()
    projector = projector_checks()
    wave_packet = wave_packet_check(directions)

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
        "parent_pins": [
            parent_pin(
                source_raw,
                source,
                SOURCE_CONTRACT,
                "frozen target-free FZ-12 source ray",
            ),
            parent_pin(
                remainder_raw,
                remainder,
                REMAINDER_CONTRACT,
                "exact full spatial symbol and q<=1 remainder",
            ),
            {
                "path": LEAN_PATH.relative_to(REPO_ROOT).as_posix(),
                "role": "kernel-checked basis-free transverse oscillator lift",
                "bytes": len(lean_raw),
                "sha256": f"sha256:{sha256(lean_raw)}",
                "theorems": [
                    "transversePolarization_finrank",
                    "transverseProjector_idempotent",
                    "photonModeFrequency_sq",
                    "photonModeFrequency_zero",
                    "photonMode_second_order",
                    "photonModeEnergy_firstVariation_generator_zero",
                ],
                "sorry_free": True,
            },
        ],
        "transverse_sector": projector,
        "hamiltonian_contract": {
            "state": "(A_T, Pi_T) in ker(k dot -) x ker(k dot -)",
            "transverse_dimension_for_nonzero_k": 2,
            "spatial_action": "A_T maps to Lambda_a(k) A_T",
            "generator": "(A_T,Pi_T) maps to (Pi_T,-Lambda_a(k) A_T)",
            "canonical_energy": (
                "H_k=(|Pi_T|^2+Lambda_a(k)|A_T|^2)/2"
            ),
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
            "lambda_hat": (
                "(1/5) sum_e [1-cos(q (w_e dot n))]"
            ),
            "omega_hat": "sqrt(lambda_hat)",
            "direction_count": 30,
            "q_domain_consumed": "0 <= q <= 1",
            "positive_bound_consumed": "lambda_hat >= (19/20)q^2",
            "target_free": True,
        },
        "wave_packet_check": wave_packet,
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
    receipt["receipt_sha256"] = base.tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise FreePhotonError("invalid committed transverse-oscillator receipt") from error
    require(
        isinstance(committed, dict),
        "transverse-oscillator receipt root is not an object",
    )
    require(
        raw == canonical_json_bytes(committed),
        "noncanonical transverse-oscillator receipt",
    )
    require(committed.get("schema") == SCHEMA, "transverse-oscillator schema drift")
    require(committed.get("status") == STATUS, "transverse-oscillator status drift")
    body = {key: item for key, item in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(canonical_json_bytes(body)),
        "transverse-oscillator self-digest drift",
    )
    require(
        raw == canonical_json_bytes(build_receipt()),
        "transverse-oscillator ancestry or result drift",
    )
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if not args.write and not args.verify:
        parser.error("choose --write or --verify")
    if args.write:
        RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    if args.verify:
        verify_committed_receipt()
        print("FZ12_TRANSVERSE_OSCILLATOR_CANDIDATE_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
