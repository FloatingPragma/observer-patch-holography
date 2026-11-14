"""Target-blind local dynamics for a twelve-port echosahedral carrier.

The regular icosahedral incidence is used only as a hidden local coupling
presentation.  No global S2/H3 coordinates or clock candidate enter this
module.  The finite unitary channel is an exact A5 intertwiner, so independently
reorienting a carrier commutes with its local propagation.

This is a carrier-level construction, not a BW/KMS clock.  Its dimensionless
evolution parameter has no physical angular normalization, and all physical
clock/emergence receipts remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from typing import Any

import numpy as np
from scipy.linalg import expm

from oph_fpe.core.icosahedral import (
    build_geodesic_icosahedral_tower,
    icosahedral_a5_port_permutations,
)


@dataclass(frozen=True)
class LocalRecurrentCarrierState:
    """Complex twelve-channel state for a batch of local carriers."""

    amplitudes: np.ndarray
    intrinsic_phase: np.ndarray

    def __post_init__(self) -> None:
        amplitudes = np.array(self.amplitudes, dtype=np.complex128, copy=True)
        amplitudes = _validated_amplitudes(amplitudes)
        phase = np.array(self.intrinsic_phase, dtype=float, copy=True)
        if phase.shape != (amplitudes.shape[0],) or not np.all(np.isfinite(phase)):
            raise ValueError("intrinsic_phase must contain one finite value per carrier")
        if np.any((phase < 0.0) | (phase >= 1.0)):
            raise ValueError("intrinsic_phase must be represented in [0, 1)")
        amplitudes.setflags(write=False)
        phase.setflags(write=False)
        object.__setattr__(self, "amplitudes", amplitudes)
        object.__setattr__(self, "intrinsic_phase", phase)

    @property
    def carrier_count(self) -> int:
        return int(self.amplitudes.shape[0])


def reference_icosahedral_coupling() -> np.ndarray:
    """Return the 12x12 combinatorial Laplacian of the local carrier."""

    base = build_geodesic_icosahedral_tower(0).levels[0]
    adjacency = np.zeros((12, 12), dtype=float)
    edges = np.asarray(base.edges, dtype=np.int64)
    adjacency[edges[:, 0], edges[:, 1]] = 1.0
    adjacency[edges[:, 1], edges[:, 0]] = 1.0
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def _reference_adjacency() -> np.ndarray:
    base = build_geodesic_icosahedral_tower(0).levels[0]
    adjacency = np.zeros((12, 12), dtype=np.int64)
    for left, right in np.asarray(base.edges, dtype=np.int64):
        adjacency[int(left), int(right)] = 1
        adjacency[int(right), int(left)] = 1
    return adjacency


def _graph_distances(adjacency: np.ndarray) -> np.ndarray:
    distances = np.full((12, 12), -1, dtype=np.int64)
    for start in range(12):
        distances[start, start] = 0
        queue = [start]
        while queue:
            node = queue.pop(0)
            for neighbor in np.flatnonzero(adjacency[node]):
                neighbor = int(neighbor)
                if distances[start, neighbor] < 0:
                    distances[start, neighbor] = distances[start, node] + 1
                    queue.append(neighbor)
    return distances


def _solve_unique_fraction_system(
    rows: list[list[int]], rhs: list[int]
) -> tuple[Fraction, ...]:
    """Solve an overdetermined exact system with four unknowns."""

    matrix = [
        [Fraction(value) for value in row] + [Fraction(target)]
        for row, target in zip(rows, rhs, strict=True)
    ]
    pivot_rows: list[int] = []
    active = 0
    for column in range(4):
        pivot = next(
            (index for index in range(active, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[active], matrix[pivot] = matrix[pivot], matrix[active]
        scale = matrix[active][column]
        matrix[active] = [entry / scale for entry in matrix[active]]
        for index in range(len(matrix)):
            if index == active or not matrix[index][column]:
                continue
            factor = matrix[index][column]
            matrix[index] = [
                matrix[index][j] - factor * matrix[active][j]
                for j in range(5)
            ]
        pivot_rows.append(column)
        active += 1
    if pivot_rows != [0, 1, 2, 3]:
        raise ValueError("impulse/readback polynomial is not uniquely determined")
    if any(all(not row[column] for column in range(4)) and row[4] for row in matrix):
        raise ValueError("impulse/readback polynomial constraints are inconsistent")
    return tuple(matrix[index][4] for index in range(4))


def reference_impulse_readback_report() -> dict[str, Any]:
    """Recover the maximal-distance echo from target-blind impulse readback.

    A delta is injected at every port. The local adjacency recurrence is read
    for k=0..diameter. One homogeneous polynomial filter is solved from those
    observations by requiring cancellation on every nearer distance shell and
    unit response on the maximal shell. No antipode map or gauge label is an
    input to the solve.
    """

    adjacency = _reference_adjacency()
    distances = _graph_distances(adjacency)
    diameter = int(distances.max())
    if diameter != 3:
        raise ValueError(f"reference carrier diameter drifted to {diameter}")
    farthest_counts = np.sum(distances == diameter, axis=1)
    if not np.all(farthest_counts == 1):
        raise ValueError("every impulse source must have one maximal-distance port")
    powers = [np.eye(12, dtype=np.int64)]
    for _ in range(diameter):
        powers.append(powers[-1] @ adjacency)
    target = (distances == diameter).astype(np.int64)
    rows: list[list[int]] = []
    rhs: list[int] = []
    for source in range(12):
        for readback in range(12):
            rows.append(
                [int(powers[k][source, readback]) for k in range(diameter + 1)]
            )
            rhs.append(int(target[source, readback]))
    coefficients = _solve_unique_fraction_system(rows, rhs)
    reconstructed = np.zeros((12, 12), dtype=object)
    for k, coefficient in enumerate(coefficients):
        reconstructed += coefficient * powers[k]
    if any(
        reconstructed[i, j] != Fraction(int(target[i, j]))
        for i in range(12)
        for j in range(12)
    ):
        raise ValueError("the solved impulse filter does not isolate the farthest shell")
    farthest_map = tuple(int(np.flatnonzero(target[source])[0]) for source in range(12))
    if len(set(farthest_map)) != 12 or any(
        farthest_map[farthest_map[index]] != index for index in range(12)
    ):
        raise ValueError("the maximal-distance echo is not a permutation involution")
    return {
        "diameter": diameter,
        "source_count": 12,
        "recurrence_steps": [0, 1, 2, 3],
        "homogeneous_filter_coefficients": [str(value) for value in coefficients],
        "polynomial_identity": "J = (A^3 - 4*A^2 - 5*A + 10*I)/10",
        "unique_solution_rank": 4,
        "unique_farthest_port_per_source": True,
        "nearer_shells_cancelled": True,
        "farthest_port_map": list(farthest_map),
        "target_labels_used": False,
        "downstream_labels_used": False,
    }


def reference_negative_antipode_response() -> np.ndarray:
    """Return ``R=-J`` from the operational maximal-distance echo protocol."""

    report = reference_impulse_readback_report()
    response = np.zeros((12, 12), dtype=np.complex128)
    response[
        np.arange(12),
        np.asarray(report["farthest_port_map"], dtype=np.int64),
    ] = -1.0
    return response


def apply_local_charged_response(
    state: LocalRecurrentCarrierState,
) -> LocalRecurrentCarrierState:
    """Apply the impulse/readback-derived reversible response operator once."""

    amplitudes = _validated_amplitudes(state.amplitudes)
    response = reference_negative_antipode_response()
    transformed = amplitudes @ response.T
    return LocalRecurrentCarrierState(
        amplitudes=np.asarray(transformed, dtype=np.complex128),
        intrinsic_phase=np.asarray(state.intrinsic_phase, dtype=float),
    )


def initialize_local_recurrent_carriers(
    carrier_count: int,
    *,
    seed: int,
) -> LocalRecurrentCarrierState:
    """Initialize normalized source states without a candidate clock scale."""

    count = int(carrier_count)
    if count <= 0:
        raise ValueError("carrier_count must be positive")
    rng = np.random.default_rng(int(seed))
    amplitudes = rng.normal(size=(count, 12)) + 1j * rng.normal(size=(count, 12))
    norms = np.linalg.norm(amplitudes, axis=1, keepdims=True)
    amplitudes = np.asarray(amplitudes / norms, dtype=np.complex128)
    phase = rng.random(count, dtype=np.float64)
    return LocalRecurrentCarrierState(amplitudes=amplitudes, intrinsic_phase=phase)


def propagate_local_recurrent_carriers(
    state: LocalRecurrentCarrierState,
    *,
    intrinsic_step: float,
    coupling_strength: float = 1.0,
) -> LocalRecurrentCarrierState:
    """Apply one A5-equivariant local unitary and an R/Z phase increment."""

    amplitudes = _validated_amplitudes(state.amplitudes)
    step = float(intrinsic_step)
    strength = float(coupling_strength)
    if not math.isfinite(step) or not math.isfinite(strength):
        raise ValueError("intrinsic_step and coupling_strength must be finite")
    coupling = reference_icosahedral_coupling()
    unitary = expm(-1j * step * strength * coupling)
    propagated = amplitudes @ unitary.T
    phase = np.mod(np.asarray(state.intrinsic_phase, dtype=float) + step, 1.0)
    if phase.shape != (amplitudes.shape[0],) or not np.all(np.isfinite(phase)):
        raise ValueError("intrinsic_phase must contain one finite value per carrier")
    return LocalRecurrentCarrierState(
        amplitudes=np.asarray(propagated, dtype=np.complex128),
        intrinsic_phase=phase,
    )


def local_port_statistics(state: LocalRecurrentCarrierState) -> np.ndarray:
    """Return quotient-visible per-port intensities; hidden XYZ is not consulted."""

    amplitudes = _validated_amplitudes(state.amplitudes)
    return np.asarray(np.abs(amplitudes) ** 2, dtype=float)


def local_a5_dynamics_report(
    *,
    intrinsic_step: float = 0.137,
    coupling_strength: float = 1.0,
    tolerance: float = 5.0e-12,
) -> dict[str, Any]:
    """Recompute A5 covariance, unitarity, spectrum, and clock firewalls."""

    step = float(intrinsic_step)
    strength = float(coupling_strength)
    if not all(math.isfinite(value) for value in (step, strength, tolerance)):
        raise ValueError("audit parameters must be finite")
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    coupling = reference_icosahedral_coupling()
    response = reference_negative_antipode_response()
    unitary = expm(-1j * step * strength * coupling)
    unitary_one_period_later = expm(-1j * (step + 1.0) * strength * coupling)
    identity = np.eye(12, dtype=np.complex128)
    unitarity_residual = float(
        np.linalg.norm(unitary.conj().T @ unitary - identity, ord="fro")
    )
    commutator_residuals: list[float] = []
    equivariance_residuals: list[float] = []
    probe = initialize_local_recurrent_carriers(3, seed=0xA5).amplitudes
    for permutation in icosahedral_a5_port_permutations():
        matrix = _permutation_matrix(permutation)
        commutator_residuals.append(
            float(np.linalg.norm(matrix @ coupling - coupling @ matrix, ord="fro"))
        )
        left = (probe @ matrix.T) @ unitary.T
        right = (probe @ unitary.T) @ matrix.T
        equivariance_residuals.append(float(np.linalg.norm(left - right, ord="fro")))
    eigenvalues = np.linalg.eigvalsh(coupling)
    multiplicities = _spectral_multiplicities(eigenvalues)
    maximum_commutator = max(commutator_residuals, default=math.inf)
    maximum_equivariance = max(equivariance_residuals, default=math.inf)
    nontriviality_norm = float(np.linalg.norm(unitary - identity, ord="fro"))
    response_involution_residual = float(
        np.linalg.norm(response @ response - identity, ord="fro")
    )
    response_self_adjoint_residual = float(
        np.linalg.norm(response.conj().T - response, ord="fro")
    )
    response_coupling_commutator_residual = float(
        np.linalg.norm(response @ coupling - coupling @ response, ord="fro")
    )
    response_probe = LocalRecurrentCarrierState(
        amplitudes=probe,
        intrinsic_phase=np.zeros(probe.shape[0], dtype=float),
    )
    response_after = apply_local_charged_response(response_probe)
    response_norm_residual = float(
        np.linalg.norm(
            np.linalg.norm(response_after.amplitudes, axis=1)
            - np.linalg.norm(probe, axis=1)
        )
    )
    period_descent_residual = float(
        np.linalg.norm(unitary_one_period_later - unitary, ord="fro")
    )
    phase_probe = LocalRecurrentCarrierState(
        amplitudes=np.eye(1, 12, dtype=np.complex128),
        intrinsic_phase=np.asarray([0.95], dtype=float),
    )
    phase_after = propagate_local_recurrent_carriers(
        phase_probe,
        intrinsic_step=0.1,
        coupling_strength=strength,
    ).intrinsic_phase
    phase_register_wrap_receipt = bool(
        phase_after.shape == (1,) and abs(float(phase_after[0]) - 0.05) <= tolerance
    )
    dynamics_descends_to_r_mod_z = period_descent_residual <= tolerance
    nontrivial = nontriviality_norm > tolerance
    receipt = bool(
        unitarity_residual <= tolerance
        and maximum_commutator <= tolerance
        and maximum_equivariance <= tolerance
        and sorted(row["multiplicity"] for row in multiplicities) == [1, 3, 3, 5]
        and response_involution_residual <= tolerance
        and response_self_adjoint_residual <= tolerance
        and response_coupling_commutator_residual <= tolerance
        and response_norm_residual <= tolerance
    )
    payload = {
        "schema": "oph.echosahedral_local_recurrent_dynamics.v1",
        "carrier_port_count": 12,
        "local_edge_count": 30,
        "local_face_count": 20,
        "a5_action_count": len(icosahedral_a5_port_permutations()),
        "coupling": "combinatorial_icosahedral_graph_laplacian",
        "evolution": "finite_unitary_exp_minus_i_u_gL",
        "intrinsic_phase_space": "R_mod_Z",
        "intrinsic_step": step,
        "coupling_strength": strength,
        "coupling_spectrum": multiplicities,
        "unitarity_residual": unitarity_residual,
        "maximum_a5_coupling_commutator_residual": maximum_commutator,
        "maximum_a5_propagation_equivariance_residual": maximum_equivariance,
        "propagation_nontriviality_norm": nontriviality_norm,
        "charged_response_operator": "negative_graph_antipode_involution",
        "charged_response_source": "target_blind_maximal_distance_impulse_readback",
        "charged_response_impulse_readback": reference_impulse_readback_report(),
        "charged_response_involution_residual": response_involution_residual,
        "charged_response_self_adjoint_residual": response_self_adjoint_residual,
        "charged_response_coupling_commutator_residual": (
            response_coupling_commutator_residual
        ),
        "charged_response_norm_residual": response_norm_residual,
        "CHARGED_RESPONSE_OPERATOR_RECEIPT": bool(
            response_involution_residual <= tolerance
            and response_self_adjoint_residual <= tolerance
            and response_coupling_commutator_residual <= tolerance
            and response_norm_residual <= tolerance
        ),
        "one_phase_period_propagator_residual": period_descent_residual,
        "hidden_xyz_coordinates_used_by_dynamics": False,
        "global_support_chart_used_by_dynamics": False,
        "candidate_clock_scale_used_by_dynamics": False,
        "LOCAL_A5_EQUIVARIANT_PROPAGATION_RECEIPT": receipt,
        "LOCAL_NONTRIVIAL_REVERSIBLE_PROPAGATION_RECEIPT": bool(
            receipt and nontrivial
        ),
        "LOCAL_PHASE_REGISTER_MOD_ONE_RECEIPT": phase_register_wrap_receipt,
        "LOCAL_DYNAMICS_DESCENDS_TO_R_MOD_Z_RECEIPT": bool(
            receipt and nontrivial and dynamics_descends_to_r_mod_z
        ),
        "LOCAL_INTRINSIC_R_MOD_Z_RECURRENCE_RECEIPT": False,
        "PHYSICAL_2PI_CLOCK_SELECTION_RECEIPT": False,
        "BW_KMS_CLOCK_RECEIPT": False,
        "PHYSICAL_H3_KMS_EMERGENCE_RECEIPT": False,
        "claim_boundary": (
            "This certifies a target-blind finite local A5-equivariant unitary "
            "carrier channel and separately audits a mod-one phase register. The "
            "generic propagator does not descend to that circle, and neither object "
            "chooses radians, 2pi, a KMS temperature, an H3 frame, or physical time."
        ),
    }
    payload["dynamics_sha256"] = "sha256:" + hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return payload


def _validated_amplitudes(value: np.ndarray) -> np.ndarray:
    amplitudes = np.asarray(value, dtype=np.complex128)
    if amplitudes.ndim != 2 or amplitudes.shape[1] != 12:
        raise ValueError("carrier amplitudes must have shape (carrier_count, 12)")
    if not np.all(np.isfinite(amplitudes.real)) or not np.all(
        np.isfinite(amplitudes.imag)
    ):
        raise ValueError("carrier amplitudes must be finite")
    return amplitudes


def _permutation_matrix(permutation: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((12, 12), dtype=float)
    matrix[np.asarray(permutation, dtype=np.int64), np.arange(12)] = 1.0
    return matrix


def _spectral_multiplicities(eigenvalues: np.ndarray) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for value in np.asarray(eigenvalues, dtype=float):
        if rows and abs(value - float(rows[-1]["eigenvalue"])) <= 1.0e-10:
            rows[-1]["multiplicity"] += 1
        else:
            rows.append({"eigenvalue": float(value), "multiplicity": 1})
    return rows


__all__ = [
    "LocalRecurrentCarrierState",
    "initialize_local_recurrent_carriers",
    "local_a5_dynamics_report",
    "local_port_statistics",
    "propagate_local_recurrent_carriers",
    "reference_impulse_readback_report",
    "reference_icosahedral_coupling",
    "reference_negative_antipode_response",
]
