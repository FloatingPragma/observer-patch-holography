"""Compare one-pole calibration with all native seam and cap effects on held stars."""

import argparse
import json
from copy import copy
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm

from model.experiments.interval_echo_angular import modular
from model.experiments.interval_echo_einstein_station import cap
from model.experiments.interval_echo_null_balance.experiment import check_native_gate
from model.experiments.interval_echo_space.shell_source import prepare
from model.executable_model import federation
from .certificate import verify
from .stationary import waiting_protocol


def acquire(rho, operator, ports, rng, shots):
    """Read the star using the port and balanced two-port terminal settings."""
    ids = list(ports)
    p = rho.diagonal()[ids].real
    rates = np.r_[p, max(0.0, 1 - p.sum())]
    counts = rng.multinomial(shots, rates / rates.sum())
    records = [dict(ports=ids, counts=counts.tolist())]
    result = operator.diagonal()[ids] @ counts[:-1] / shots
    variance = operator.diagonal()[ids] ** 2 @ p - (operator.diagonal()[ids] @ p) ** 2
    for a, b in combinations(ports, 2):
        base = (rho[a, a] + rho[b, b]).real / 2
        overlap = rho[a, b].real
        rates = np.maximum([base + overlap, base - overlap, 1 - 2 * base], 0.0)
        draws = rng.multinomial(shots, rates / rates.sum())
        records.append(dict(pair=[a, b], counts=draws.tolist()))
        result += operator[a, b] * (draws[0] - draws[1]) / shots
        variance += operator[a, b] ** 2 * (2 * base - 4 * overlap**2)
    return float(result), float(variance / shots), records


def energy_read(rho, edges, rng, shots):
    value, variance = 0.0, 0.0
    records = []
    for a, b in edges:
        # The complete 0,8,9 read consists of a miss, symmetric click or
        # antisymmetric click of the native two-port half-exchange analyzer.
        total = (rho[a, a] + rho[b, b]).real
        overlap = rho[a, b].real
        rates = np.maximum([1 - total, total / 2 + overlap, total / 2 - overlap], 0.0)
        draws = rng.multinomial(shots, rates / rates.sum())
        levels = np.array([0.0, 8.0, 9.0])
        records.append(dict(pair=[a, b], counts=draws.tolist()))
        value += levels @ draws / shots
        variance += (levels**2 @ rates - (levels @ rates) ** 2) / shots
    return float(value), float(variance), records


def exchange_only():
    # On the one-load filter ancilla the contrast is scalar. At the
    # read-basis exposure its phase is unity, and at each multiplexed
    # angle it contributes the same phase for every control word.
    full = cap.NativeFilter()
    changed = copy(full)
    exchange, _ = cap.law.local_term_parts()
    changed.binary = full.injection.T @ exchange @ full.injection
    phase = expm(-0.5j * np.pi * changed.binary) @ expm(0.5j * np.pi * full.binary)
    changed.read_basis = phase @ full.read_basis
    odds = 2.0 ** (-np.array([8.0, 9.0]) / 2)
    before, after = [f.success(odds, 1 / 16) for f in (full, changed)]
    # The uniform one-load column is annihilated by every exchange.
    uniform = np.array([0.0, 1.0, 1.0, 0.0]) / np.sqrt(12)
    return dict(
        filter_effect_error=float(np.linalg.norm(before - after, ord=2)),
        uniform_column_norm=float(np.linalg.norm(changed.binary @ uniform)),
    )


def run():
    bounds = verify()
    mid = {
        key: np.mean(json.loads(value))
        for key, value in bounds.items()
        if isinstance(value, str)
    }
    n = cap.chord_coordinates()
    basis = np.column_stack((np.ones(12) / np.sqrt(12), n / 2))
    parent, untwirled, prob = modular.common_reference()
    stars, logs, _ = modular.regional_operators(parent)
    _, changed, _ = modular.regional_operators(untwirled)
    sector = basis.T @ logs @ basis
    edges = cap.cell().seams
    all_effects = []
    for a, b in edges:
        h = np.zeros((12, 12))
        h[a, a] = h[b, b] = 8.5
        h[a, b] = h[b, a] = -0.5
        all_effects.append(h)
    all_effects = np.array(all_effects)
    total = all_effects.sum(axis=0)
    total_low = basis.T @ total @ basis
    spokes = []
    rims = []
    for i, star in enumerate(stars):
        spokes.append([j for j, e in enumerate(edges) if i in e])
        rims.append(
            [j for j, e in enumerate(edges) if i not in e and set(e) <= set(star)]
        )
    assert all(len(g) == 5 for g in spokes + rims)
    positive_native = mid["spoke_weight"] * all_effects[spokes].sum(axis=1) + mid[
        "rim_weight"
    ] * all_effects[rims].sum(axis=1)
    pos = basis.T @ positive_native @ basis
    predictions = (
        pos + mid["identity"] * np.eye(4) + mid["total_energy_weight"] * total_low
    )

    def norm(a):
        return float(np.max(np.linalg.norm(a, ord=2, axis=(-2, -1))))

    residual = norm(predictions - sector)
    angle_residual = norm((pos - pos.mean(axis=0)) - (sector - sector.mean(axis=0)))
    direction = n[0]
    axial = n @ direction / 2
    triplet = basis[:, 1:]
    transverse = (triplet @ triplet.T - np.outer(axial, axial)) / 2
    rng = np.random.default_rng(9150347)
    # Difference of two independently prepared shell states. The axial and
    # transverse mixtures have equal total energy throughout this sector.
    board = federation(0)
    carrier = SimpleNamespace(car=board, active=board.seams)
    _, axes = np.linalg.eigh(transverse)
    pilots = [axial, *axes[:, -2:].T]
    achieved = []
    pulse_counts = []
    errors = []
    for pilot in pilots:
        wanted = np.zeros(board.n, dtype=complex)
        wanted[:12] = pilot
        actual, events = prepare(carrier, wanted)
        errors.append(float(1 - abs(np.vdot(wanted, actual)) ** 2))
        pulse_counts.append(events)
        achieved.append(np.outer(actual[:12], actual[:12].conj()))
    states = [achieved[0], (achieved[1] + achieved[2]) / 2]
    shots = 2**20
    readings = {}
    for label, reader in (
        ("modular", lambda r: acquire(r, logs[0], stars[0], rng, shots)),
        ("spokes", lambda r: energy_read(r, [edges[j] for j in spokes[0]], rng, shots)),
        ("rim", lambda r: energy_read(r, [edges[j] for j in rims[0]], rng, shots)),
    ):
        (ax, va, axial_counts), (tr, vt, transverse_counts) = map(reader, states)
        readings[label] = {
            "difference": ax - tr,
            "standard_error": np.sqrt(va + vt),
            "axial": axial_counts,
            "transverse": transverse_counts,
        }
    # The operator witness probes every individual seam even with asymmetric weights.
    pole = direction
    delta = np.zeros((4, 4))
    delta[0, 0] = 1
    delta[1:, 1:] = (
        -np.eye(3) / 3
        - mid["separating_quadrupole_weight"]
        * (3 * np.outer(pole, pole) - np.eye(3))
        / 2
    )
    delta[0, 1:] = delta[1:, 0] = -mid["separating_dipole_weight"] * pole / 2
    tangent = basis @ delta @ basis.T
    witness = np.einsum("ij,nji->n", tangent, all_effects)
    hamiltonian = cap.native_energy(
        [cap.seam_indices(e) for e in edges], np.ones(len(edges))
    )[0]
    stationary, waits = waiting_protocol(parent, hamiltonian)
    _, stationary_logs, _ = modular.regional_operators(stationary)
    waiting_rng = np.random.default_rng(9150348)
    ax, tr = [
        acquire(state, stationary_logs[0], stars[0], waiting_rng, shots)
        for state in states
    ]
    return dict(
        stationary_waits=waits,
        stationary_triplet=dict(
            expectation=float(
                np.trace((states[0] - states[1]) @ stationary_logs[0]).real
            ),
            difference=ax[0] - tr[0],
            standard_error=float(np.sqrt(ax[1] + tr[1])),
            axial=ax[2],
            transverse=tr[2],
        ),
        rounded_stationary_source_trace_change=float(
            np.linalg.norm(stationary - parent, ord="nuc")
        ),
        rounded_stationary_source_commutator=norm(
            hamiltonian @ stationary - stationary @ hamiltonian
        ),
        rounded_stationary_charge_change=norm(stationary_logs - logs),
        true_bounds=bounds,
        stationary_bounds=verify(stationary=True),
        reference_acceptance=prob,
        unchanged_contact_error=check_native_gate(),
        exchange_only_control=exchange_only(),
        held_stars_operator_max=residual,
        centered_stars_operator_max=angle_residual,
        unchanged_weights_untwirled_error=norm(predictions - basis.T @ changed @ basis),
        unconstrained_witness_seam_range=[witness.min(), witness.max()],
        unconstrained_witness_modular=float(np.trace(tangent @ logs[0])),
        equal_energy_triplet=dict(
            preparation_gates=pulse_counts,
            max_infidelity=max(errors),
            shots_per_setting=shots,
            ordinary_energy_difference=float(
                np.trace((states[0] - states[1]) @ total).real
            ),
            expectation={
                label: float(np.trace((states[0] - states[1]) @ eff).real)
                for label, eff in (
                    ("modular", logs[0]),
                    ("spokes", all_effects[spokes[0]].sum(axis=0)),
                    ("rim", all_effects[rims[0]].sum(axis=0)),
                )
            },
            sampled=readings,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(output)
    else:
        print(output, end="")
