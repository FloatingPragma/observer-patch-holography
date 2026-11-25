"""Native rail rotations and a first-derivative-matched right-head analyzer."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import root

from .. import signed_transfers as rails
from ..full_algebra import station


def nominal_experiment():
    record = json.loads(Path(station.__file__).with_name("results.json").read_text())
    chi = np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    effects = station.radar_menus()[6]

    def left(blocks):
        return np.trace(blocks.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)

    right = rails.right_coin
    paulis = rails.phase.AXES[1:]
    eye = np.eye(2)

    directions = np.zeros((4, 2, 4, 4), dtype=complex)
    for b, (offset, scale) in enumerate(((0.2, 0.1), (0.05, 0.03))):
        h = (chi[b] - offset * np.eye(4)) / scale
        for j, p in enumerate((h, h @ h)):
            directions[2 * b + j, b] = p - np.trace(p) * np.eye(4) / 4

    def log_slopes(reference, directions):
        return -np.einsum("bij,kbji->k", station.log_blocks(reference), directions).real

    constraint_rows = np.array(
        [
            station.range_gradient(station.pair(effects, chi))[1]
            @ np.array(
                [station.pair(effects, direction) for direction in directions]
            ).T,
            log_slopes(left(chi), [left(d) for d in directions]),
            log_slopes(right(chi)[None], [right(d)[None] for d in directions]),
        ]
    )
    dependent = [0, 2, 3]
    coeff = np.ones(4)
    coeff[dependent] = np.linalg.solve(
        constraint_rows[:, dependent], -constraint_rows[:, 1]
    )
    tangent = np.einsum("k,kbij->bij", coeff, directions)

    def entropy(blocks):
        p = np.linalg.eigvalsh(blocks)
        return -np.sum(p * np.log(p))

    def constraints(blocks):
        return np.array(
            [
                station.geometry(station.pair(effects, blocks))[1],
                entropy(left(blocks)),
                entropy(right(blocks)),
            ]
        )

    base = constraints(chi)
    step = 1e-4

    def displaced(dependent_coeff):
        x = np.copy(coeff)
        x[dependent] = dependent_coeff
        return chi + step * np.einsum("i,ibjk->bjk", x, directions)

    fit = root(lambda x: (constraints(displaced(x)) - base) / step, coeff[dependent])
    active = displaced(fit.x)
    residual = np.max(np.abs(constraints(active) - base))
    assert residual < 2e-13
    prep = station.preparation(active)
    assert (
        prep["min_kept"] > 0 and prep["min_dephased"] > 0 and prep["max_error"] < 5e-14
    )

    def bloch(matrix):
        return np.array([np.trace(matrix @ p).real for p in paulis])

    def pauli(v):
        return sum(c * p for c, p in zip(v, paulis))

    axis = bloch(right(chi))
    gap = np.linalg.norm(axis)
    axis /= gap
    m = bloch(right(tangent))
    m /= np.linalg.norm(m)
    assert abs(axis @ m) < 1e-11
    omega = np.log((1 + gap) / (1 - gap))
    p_axis = (omega * m + 2 * np.pi * np.cross(axis, m)) / np.hypot(omega, 2 * np.pi)
    effect = (eye + pauli(p_axis)) / 2
    assert np.linalg.norm(effect @ effect - effect) < 1e-13
    diff_r = right(tangent)
    impulse = np.trace(diff_r @ effect).real

    # Boost the observable by U_s. An equivalent forward laboratory
    # experiment rotates the state by U_s^* before applying the effect.
    boost = np.log(2) / (2 * np.pi)
    endpoint = []
    for s in (-boost, 0.0, boost, 2 * boost):
        u = expm(-0.5j * s * omega * pauli(axis))
        signal = np.trace(effect @ u @ diff_r @ u.conj().T).real
        gate_error = rails.check_rail(pauli(axis), s * omega)
        endpoint.append(
            dict(
                modular_time=float(s),
                measured_ratio=float(signal / impulse),
                affine_ratio=float(np.exp(2 * np.pi * s)),
                rail_operator_error=gate_error,
            )
        )

    du = -0.5j * omega * pauli(axis)
    first = np.trace(effect @ (du @ diff_r - diff_r @ du)).real
    return dict(
        tangent_coefficients=coeff.tolist(),
        entropy_slope=float(log_slopes(chi, [tangent])[0]),
        area_slope=float(
            station.range_gradient(station.pair(effects, chi))[0]
            @ station.pair(effects, tangent)
        ),
        right_Z_slope=float(bloch(right(tangent))[2]),
        finite_step=step,
        constraint_max_abs=float(residual),
        finite_step_left_state_difference=float(
            np.linalg.norm(left(active) - left(chi))
        ),
        finite_step_right_state_difference=float(
            np.linalg.norm(right(active) - right(chi))
        ),
        finite_step_commutator_norm=float(np.linalg.norm(chi @ active - active @ chi)),
        minimum_preparation_weight=min(prep["min_kept"], prep["min_dephased"]),
        preparation_max_error=prep["max_error"],
        analyzer_reference_slope=float(impulse),
        normalized_analyzer_modular_derivative=float(first / impulse),
        analyzer_endpoints=endpoint,
    )
