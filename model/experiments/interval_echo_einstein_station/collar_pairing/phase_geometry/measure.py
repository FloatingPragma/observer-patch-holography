"""Read a freely formed two-head state with material position detectors."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..measure import log_density
from ..ranging.measure import parity_blocks
from ...shared_geometry_probe.flight import prepared_states, formation
from .remote import AXES, PAULI, port_records


SHOTS = 2**18
REFERENCE_SHOTS = 2**74
SEED = 526018
IDENTITY, X, Y, Z = AXES


def read_dimers(board, route):
    dimers = []
    for a, b in (route[:2], route[2:]):
        dimers.append(
            min(
                (x, y)
                for x in board.neighbours[a]
                for y in board.neighbours[b]
                if x not in route
                and y not in route
                and y in board.neighbours[x]
                and board.degree[x] == 5
            )
        )
    assert len(set(route).union(*dimers)) == 8
    return dimers


def marked_gate(q, k):
    blocks = parity_blocks(q, k, 1)
    h = blocks[1:, 1:].transpose(0, 2, 1, 3).reshape(4, 4)
    bell = np.array([0, 1, 1, 0] if q == 1 else [1, 0, 0, 1])
    np.testing.assert_array_equal(2 * h, k * (2 * np.eye(4) - np.outer(bell, bell)))
    return expm(-1j * np.pi * h / 8).reshape(2, 2, 2, 2)


def transfer(q, k, rho):
    gate = marked_gate(q, k)[:, :, 0, :]
    # Each initial rail is loaded at its first port. Retain both head output
    # indices until tracing, so loss or residual entanglement would be visible.
    embed = np.einsum("iax,jby->ijabxy", gate, gate).reshape(4, 4, 4)
    return np.einsum("phx,xy,qhy->pq", embed, rho, embed.conj())


def rail_analyzers(q):
    # Restrict the native blank seam to its two possible material launch sites.
    native = formation.emission.four_h(8)
    law = formation.emission.LAW
    cols = [law.SRC.index((1 + q, 1)), law.SRC.index((1, 1 + q))]
    h = native[np.ix_(cols, cols)] / 4
    np.testing.assert_array_equal(h, 8 * IDENTITY + (IDENTITY - X) / 2)
    np.testing.assert_allclose(expm(-1j * np.pi * h), X, atol=1e-13)
    half = expm(-1j * np.pi * h / 2)
    reference_gap = 8 * ((1 + q) ** 2 - 1)
    turn = np.diag([np.exp(-1j * np.pi * np.sign(reference_gap) / 2), 1])
    outputs = []
    # The native phase and blank half-exchange realize X,Y,Z up to outcome labels.
    for unitary in (half @ turn, half, IDENTITY):
        outputs.append(unitary.conj().T @ Z @ unitary)
    w = Z if q == 1 else Z @ X
    for i, (detector, axis) in enumerate(zip(outputs, (X, Y, Z))):
        mapped = w.conj().T @ detector @ w
        orientation = np.trace(mapped @ axis).real / 2
        np.testing.assert_allclose(mapped, np.sign(orientation) * axis, atol=1e-13)
        outputs[i] = np.sign(orientation) * detector
    return np.array(
        [np.kron(a, b) for a in (IDENTITY, *outputs) for b in (IDENTITY, *outputs)]
    )


def readings(q, k, states):
    detectors = rail_analyzers(q)
    return {
        name: np.einsum("kij,ji->k", detectors, transfer(q, k, rho)).real
        for name, rho in states.items()
    }


def counts(rng, shots, means):
    # NumPy's binomial takes a signed machine count. Aggregate exact integer counts.
    chunk = min(shots, 2**60)
    assert shots % chunk == 0
    out = []
    for mean in means:
        probability = np.clip((1 + mean) / 2, 0, 1)
        draws = rng.binomial(chunk, probability, size=shots // chunk)
        out.append(sum(map(int, draws)))
    return out


def signs(successes, shots):
    return np.array([(2 * c - shots) / shots for c in successes])


def active_record(states):
    rng = np.random.default_rng(SEED)
    all_means = {q: readings(q, 8, states) for q in (1, -1)}
    alpha_reference, alpha_trial = 0.004, 0.003
    ref = counts(rng, REFERENCE_SHOTS, all_means[1]["formed"][1:])
    estimated = (
        PAULI[0] + np.einsum("i,ijk->jk", signs(ref, REFERENCE_SHOTS), PAULI[1:])
    ) / 4
    smallest = float(np.linalg.eigvalsh(estimated)[0])
    epsilon = np.sqrt(2 * np.log(30 / alpha_reference) / REFERENCE_SHOTS)
    operator_error = np.sqrt(15) * epsilon / 2
    assert smallest > operator_error
    log_error = -np.log1p(-operator_error / smallest)
    coefficients = np.einsum("kij,ji->k", PAULI[1:], -log_density(estimated)).real / 4
    sampling = 2 * np.sqrt(
        np.log(2 / alpha_trial) * (coefficients @ coefficients) / SHOTS
    )
    exact = np.trace(
        (states["formed"] - states["independent"]) @ (-log_density(states["formed"]))
    ).real
    runs = []
    for q, means in all_means.items():
        for name, m in means.items():
            np.testing.assert_allclose(
                m, np.einsum("kij,ji->k", PAULI, states[name]).real, atol=1e-13
            )
        changed = readings(q, 0, states)
        np.testing.assert_allclose(
            changed["formed"], changed["independent"], atol=1e-13
        )
        yes = {
            name: counts(rng, SHOTS, m[1:])
            for name, m in means.items()
            if name != "tracial"
        }
        contrast = signs(yes["formed"], SHOTS) - signs(yes["independent"], SHOTS)
        center = float(coefficients @ contrast)
        radius = float(sampling + 2 * log_error)
        runs.append(
            dict(
                probe_charge=q,
                successes=yes,
                estimate=center,
                confidence_interval=[center - radius, center + radius],
                changed_law_max_population_contrast=float(
                    np.max(np.abs(changed["formed"] - changed["independent"]))
                ),
            )
        )
    return dict(
        seed=SEED,
        word_order=[a + b for a in "IXYZ" for b in "IXYZ" if a + b != "II"],
        reference_shots_per_word=REFERENCE_SHOTS,
        reference_successes=ref,
        reference_lambda_min=smallest,
        reference_operator_error_bound=operator_error,
        reference_log_observable_error_bound=log_error,
        held_shots_per_word_and_preparation=SHOTS,
        held_sampling_bound=sampling,
        family_confidence=1 - alpha_reference - 2 * alpha_trial,
        formed_minus_independent_ideal_nats=float(exact),
        held=runs,
    )


def run():
    board = formation.federation(1)
    route = formation.region(board)
    states = dict(
        zip(("formed", "independent", "tracial"), prepared_states(board, route, 1))
    )
    # Remove scalar roundoff in the propagated trace before forming Born probabilities.
    states = {name: rho / np.trace(rho).real for name, rho in states.items()}
    return dict(
        cells=board.n // 12,
        route=route,
        material_read_dimers=read_dimers(board, route),
        remote=port_records(board, route, states),
        gated=active_record(states),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(result)
    else:
        print(result, end="")
