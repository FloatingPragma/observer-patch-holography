"""Check the phase predictions from raw carrier data and draw the six registered trials."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import fisher_exact

from tools.strict_json import load
from ..submission.row_04_check import screen
from . import pair_phase as protocol

BOARD = Path(__file__).resolve().parents[4] / "executable_model/carriers/FEDERATION_LEVEL0.npz"
COUNTS = Path(__file__).with_name("pair_phase_counts.json")
SHOTS = 2**26


def raw_flights(generator, initial):
    state = initial.astype(complex)
    previous = 0
    for t in (0, 16, 64, 256):
        if t != previous:
            state = expm_multiply(-1j * (t - previous) * generator, state)
        yield state
        previous = t


def bosons(one, x, z, pair_ports):
    # Expand two normalized, symmetrized product states, including their diagonal.
    modes = []
    for i, j in (x, z):
        outer = np.outer(one[:, i], one[:, j])
        modes.append((outer + outer.T) / np.sqrt(2))
    waves = np.stack(modes, axis=-1)
    waves = waves @ (np.array([[1, 1], [1, -1]]) / 2)
    return np.sqrt(2) * waves[tuple(pair_ports.T)], waves[np.arange(len(one)), np.arange(len(one))]


def detector(port_pairs, lookup, receiver, states, doubles):
    u, v, w, x = np.array([0, 1, 5, 9]) + 12 * receiver
    # Give each input word one position in the ordered first-refusal instrument.
    first = ((port_pairs == u) | (port_pairs == v)).sum(axis=1)
    second = ((port_pairs == w) | (port_pairs == x)).sum(axis=1)
    status = np.where(first != 1, 1, np.where(second != 1, 2, 3))
    X, Z = lookup[u, w], lookup[v, x]
    status[[X, Z]] = 4
    weight = abs(states) ** 2
    bins = [
        np.full(2, 0.5),
        weight[status == 1].sum(axis=0) + (abs(doubles) ** 2).sum(axis=0),
        weight[status == 2].sum(axis=0),
        weight[status == 3].sum(axis=0),
        abs(states[X] - states[Z]) ** 2 / 2,
        abs(states[X] + states[Z]) ** 2 / 2,
    ]
    return np.array(bins).T


def draw(source, dest, bins):
    if bins.shape != (2, 6) or not np.isfinite(bins).all() or (bins < 0).any():
        raise ValueError("Invalid multinomial")
    np.testing.assert_allclose(bins.sum(axis=1), 1.0, atol=protocol.TOL, rtol=0)
    return np.array(
        [
            np.random.Generator(np.random.PCG64([41029, source, dest, i])).multinomial(
                SHOTS, row / row.sum()
            )
            for i, row in enumerate(bins)
        ]
    )


def admission_floor(bins):
    """Lower bound on the weakest squared singular value, including source acceptance."""
    admitted = np.asarray(bins)[:, -2:]
    q = admitted.sum(axis=1)
    if min(q) == 0:
        return 0.0
    f = admitted[:, 0] / q
    return float(np.prod(q) / sum(q) * (f[0] - f[1]) ** 2)


def decision(counts, forecast):
    if (
        counts.shape != (2, 6)
        or counts.dtype.kind not in "iu"
        or (counts < 0).any()
        or not np.all(counts.sum(axis=1) == SHOTS)
    ):
        raise ValueError("Expected two enrolled cohorts")
    finals = counts[:, 4:]
    admission = finals.sum(axis=1)
    expected = np.asarray(forecast)[:, 4:]
    if min(admission) == 0 or min(expected.sum(axis=1)) == 0:
        return 1.0, False
    fraction = finals[:, 0] / admission
    delta = fraction[0] - fraction[1]
    f = expected[:, 0] / expected.sum(axis=1)
    p = float(fisher_exact(finals, alternative="two-sided").pvalue)
    return p, bool(p < 0.01 / 6 and delta * (f[0] - f[1]) > 0)


def reproduce(prediction):
    edges, lookup, exchange, adjacent = screen(BOARD)
    pairs = np.transpose(np.triu_indices(len(lookup), 1))
    np.testing.assert_array_equal(lookup[tuple(pairs.T)], np.arange(len(pairs)))
    adjacency = sparse.csr_matrix((np.ones(len(edges)), edges.T), shape=lookup.shape)
    adjacency += adjacency.T
    single = (sparse.diags(np.ravel(adjacency.sum(axis=1))) - adjacency) / 2
    actual_receiver, actual_association, held = [], [], []
    floors = {}
    for source in (0, 2, 3):
        ports = 12 * source + np.array([0, 1, 5, 9])
        X, Z = lookup[ports[0], ports[2]], lookup[ports[1], ports[3]]
        initial = np.zeros((len(pairs), 2))
        initial[X], initial[Z] = [0.5, 0.5], [0.5, -0.5]
        floors[str(source)] = (1 - 2 * np.sum(initial * (exchange @ initial), axis=0) / 16).tolist()
        local = np.any(edges // 12 == source, axis=1)
        neighbors = sorted(set((edges[local] // 12).ravel().tolist()) - {source})
        if len(neighbors) != 3:
            raise ValueError("Wrong held-out routes")
        for law, strength in (("native", 8), ("exchange", 0), ("independent", None)):
            if strength is None:
                columns = np.eye(len(lookup))[:, ports]
                flights = (
                    bosons(one, (0, 2), (1, 3), pairs) for one in raw_flights(single, columns)
                )
            else:
                h = exchange - sparse.diags(2 * strength * adjacent)
                flights = ((state, np.zeros((0, 2))) for state in raw_flights(h, initial))
            for t, (states, doubles) in zip((0, 16, 64, 256), flights, strict=True):
                weight = abs(states) ** 2
                global_bins = np.array(
                    [
                        np.full(2, 0.5),
                        weight[adjacent == 1].sum(axis=0),
                        weight[adjacent == 0].sum(axis=0),
                        (abs(doubles) ** 2).sum(axis=0),
                    ]
                ).T
                actual_association.append([source, law, t, *global_bins.tolist()])
                # Undo the sign basis before applying the receiver to each sharp word.
                basis_change = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
                sharp_states, sharp_doubles = states @ basis_change, doubles @ basis_change
                for dest in [source, *neighbors]:
                    bins = detector(pairs, lookup, dest, states, doubles)
                    actual_receiver.append([source, dest, law, t, *bins.tolist()])
                    sharp = detector(pairs, lookup, dest, sharp_states, sharp_doubles)
                    source_control, receiver_control = protocol.dephased_controls(bins)
                    np.testing.assert_allclose(
                        sharp.mean(axis=0), source_control[0], atol=1e-12, rtol=0
                    )
                    a, b = lookup[12 * dest, 12 * dest + 5], lookup[12 * dest + 1, 12 * dest + 9]
                    diagonal = (abs(states[a]) ** 2 + abs(states[b]) ** 2) / 2
                    np.testing.assert_allclose(receiver_control[:, 4], diagonal, atol=1e-12, rtol=0)
                    minimum = np.linalg.svd(states[[a, b]], compute_uv=False)[-1] ** 2
                    if admission_floor(bins) > minimum + 1e-14:
                        raise ArithmeticError(
                            "Receiver bound exceeds direct channel singular value"
                        )
                    if source != 0 and dest != source and law == "native" and t == 256:
                        held.append((source, dest, bins))
    # Verify all bins and appointment identities before drawing from the second calculation.
    expected_rows = prediction["receiver"]
    calculated = dict(
        board=sha256(json.dumps(edges.tolist()).encode()).hexdigest(),
        ports=[0, 1, 5, 9],
        signs=["plus", "minus"],
        outcomes=[
            "source refused",
            "first rail refused",
            "second rail refused",
            "unequal source-side ports",
            "equal",
            "unequal",
        ],
        association_outcomes=["source refused", "adjacent", "separated", "double"],
        all_time_adjacent_given_launch=floors,
        receiver=actual_receiver,
        association=actual_association,
    )
    worst = protocol.verify_forecast(calculated, prediction)
    predicted = {
        (s, d): np.array(rows)
        for s, d, law, t, *rows in expected_rows
        if law == "native" and t == 256
    }
    trials = []
    for source, dest, bins in held:
        counts = draw(source, dest, bins)
        p, passed = decision(counts, predicted[source, dest])
        trials.append(
            dict(source=source, receiver=dest, counts=counts.tolist(), p_value=p, passed=passed)
        )
    return trials, worst


def check_trials(trials, archived):
    for row, saved in zip(trials, archived, strict=True):
        if set(saved) != set(row):
            raise ValueError("Wrong trial schema")
        for field in row.keys() - {"p_value"}:
            if json.dumps(row[field]) != json.dumps(saved[field]):
                raise ValueError(f"Different trial {field}")
        if type(saved["p_value"]) is not float or not np.isfinite(saved["p_value"]):
            raise ValueError("Invalid Fisher probability")
        np.testing.assert_allclose(row["p_value"], saved["p_value"], rtol=1e-9, atol=1e-300)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new count archive")
    args = parser.parse_args()
    trials, worst = reproduce(load(protocol.PREDICTIONS))
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(
                "[\n" + ",\n".join(json.dumps(r, allow_nan=False) for r in trials) + "\n]\n"
            )
    else:
        check_trials(trials, load(COUNTS))
    print("Held-out detections:", [(r["source"], r["receiver"], r["passed"]) for r in trials])
    print("Maximum absolute probability discrepancy:", worst)
