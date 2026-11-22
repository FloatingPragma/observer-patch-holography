"""Distinguish the stored energy of a prepared collar from its action on material probes."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

law = candidate.native_law
KAPPA = candidate.KAPPA
RECEIPT = Path(__file__).with_name("backreaction.json")


def coin_basis(magnitude, parity):
    """Two unnormalized integer eigenvectors on the occupied memory, one for each coin."""
    vectors = np.zeros((law.N_LAB, 2), dtype=int)
    for coin in (0, 1):
        vectors[law.LABELS.index((coin, magnitude)), coin] = 1
        if magnitude:
            vectors[law.LABELS.index((coin, -magnitude)), coin] = parity
    return vectors


def collar_vector(magnitude, parity, coin_phase):
    lab = coin_basis(magnitude, parity) @ np.array([1, coin_phase])
    return np.r_[0, lab] / np.linalg.norm(lab)


def classify():
    """Compute the complete scalar-action subspaces on occupied collars at the pin."""
    import sympy as sp

    W = law._write_pattern()  # sqrt(2) times the native write, exactly integral
    reflection, contrast = law._label_ops()
    twice_record_h = (
        np.eye(law.N_LAB, dtype=int) - reflection + int(2 * KAPPA) * contrast
    )
    rows = []
    exchange, agreement = law.local_term_parts()
    native = exchange + KAPPA * agreement
    for magnitude, parity in [(0, 1), (1, 1), (1, -1), (2, 1), (2, -1)]:
        basis = coin_basis(magnitude, parity)
        energy = int(KAPPA * magnitude**2 + (1 - parity) / 2)
        if np.any(twice_record_h @ basis != 2 * energy * basis):
            raise ValueError("Record eigenbasis failed")
        # For an eigenvector with positive E, native scalar action on every
        # load is equivalent to W* (record tensor load) = 0 for every load.
        restrictions = np.column_stack(
            [
                (W.T @ np.kron(column[:, None], np.eye(law.N_SRC, dtype=int))).ravel()
                for column in basis.T
            ]
        )
        gram = restrictions.T @ restrictions
        constraints = sp.Matrix(energy * gram)
        nullspace = [list(map(int, v)) for v in constraints.nullspace()]
        for coins in nullspace:
            record = np.r_[0, basis @ coins]
            columns = np.kron(record[:, None], np.eye(law.N_SRC))
            np.testing.assert_allclose(
                native @ columns, energy * columns, atol=2e-13, rtol=0
            )
        rows.append(
            dict(
                energy=energy,
                write_overlap_gram=gram.tolist(),
                scalar_action_coin_vectors=nullspace,
            )
        )
    if [len(r["scalar_action_coin_vectors"]) for r in rows] != [2, 0, 0, 1, 1]:
        raise ValueError("Unexpected scalar-action spectrum at the pinned coupling")
    return rows


def source_reads():
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    vacant = np.eye(law.N_SRC)[:, law.SRC.index((0, 0))]
    zero = np.kron(collar_vector(0, 1, -1), vacant)
    records = []
    for magnitude, parity in [(1, 1), (1, -1), (2, 1), (2, -1)]:
        record = np.kron(collar_vector(magnitude, parity, -1), vacant)
        preparation = (zero + record) / np.sqrt(2)
        arms = {}
        for name, generator, kappa in (
            ("native", h, KAPPA),
            ("exchange_only", exchange, 0),
        ):
            readings = []
            for duration in (1 / 32, 3 / 64):
                state = expm_multiply(-1j * duration * generator, preparation)
                x = abs(np.vdot(preparation, state)) ** 2
                y = abs(np.vdot((zero - 1j * record) / np.sqrt(2), state)) ** 2
                gap = np.arctan2(2 * y - 1, 2 * x - 1) / duration
                np.testing.assert_allclose(
                    gap, kappa * magnitude**2 + (1 - parity) / 2, atol=2e-12
                )
                readings.append(
                    dict(exposure=duration, x=float(x), y=float(y), gap=float(gap))
                )
            arms[name] = readings
        records.append(dict(magnitude=magnitude, sign_parity=parity, arms=arms))
    return records


def lift_one_load(board, edge):
    """Keep all eleven states on the marked seam during a full-sphere unit-load flight."""
    a, b = edge
    memories = 1 + law.N_LAB
    free = candidate.generator(board).tocsr()
    unit = [law.SRC.index(pair) for pair in [(1, 0), (0, 1)]]
    edge_h = law.seam_term(KAPPA)[np.ix_(unit, unit)]
    marked_ports = np.array(edge)
    seam = sparse.coo_matrix(
        (edge_h.ravel(), (np.repeat(marked_ports, 2), np.tile(marked_ports, 2))),
        shape=free.shape,
    )
    external = sparse.kron(free - seam, sparse.eye(memories), format="csr")
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    offsets = np.arange(memories) * law.N_SRC
    zero = h[np.ix_(offsets, offsets)]
    # The local vacuum block acts also when the probe is elsewhere.
    total = external + sparse.kron(sparse.eye(board.n), zero, format="csr")
    local = np.concatenate([offsets + i for i in unit])
    correction = h[np.ix_(local, local)] - np.kron(np.eye(2), zero)
    indices = np.r_[
        a * memories + np.arange(memories), b * memories + np.arange(memories)
    ]
    total += sparse.coo_matrix(
        (
            correction.ravel(),
            (np.repeat(indices, len(indices)), np.tile(indices, len(indices))),
        ),
        shape=total.shape,
    ).tocsr()
    return total


def flight(level):
    board = federation(level)
    edge = min(board.external)
    receiver = list(range((edge[1] // 12) * 12, (edge[1] // 12 + 1) * 12))
    memories = 1 + law.N_LAB
    sources = {
        "blank": np.eye(memories)[:, 0],
        "zero": collar_vector(0, 1, -1),
        "even_double": collar_vector(2, 1, -1),
        "odd_double": collar_vector(2, -1, -1),
        "ordinary_even_double": collar_vector(2, 1, 1),
    }
    states = np.zeros((board.n, memories, len(sources)), dtype=complex)
    states[edge[0]] = np.column_stack(list(sources.values()))
    states = states.reshape(board.n * memories, len(sources))
    total = lift_one_load(board, edge)
    # Remove a common scalar to reduce the cost of exponentiation.
    total = total - sparse.eye(total.shape[0]) * (6 * KAPPA)
    appointments = []
    for duration in (3 / 4, 3 / 2):
        evolved = expm_multiply(-1j * duration * total, states).reshape(
            board.n, memories, len(sources)
        )
        loads = (abs(evolved) ** 2).sum(axis=1)
        np.testing.assert_allclose(loads.sum(axis=0), 1, atol=3e-12, rtol=0)
        patch = loads[receiver].T
        outcomes = np.column_stack((patch, 1 - patch.sum(axis=1)))
        np.testing.assert_allclose(
            outcomes[1:],
            np.broadcast_to(outcomes[1], outcomes[1:].shape),
            rtol=0,
            atol=3e-12,
        )
        distance = 0.5 * np.abs(outcomes[0] - outcomes[1]).sum()
        if distance < 0.005:
            raise ValueError("The blank-memory sensitivity control lost its signal")
        appointments.append(
            dict(
                exposure=duration,
                blank_vs_marked_tv=float(distance),
                patch_plus_quiet={
                    key: list(row) for key, row in zip(sources, outcomes)
                },
            )
        )
    return dict(
        cells=board.cells,
        launch=edge[0],
        seam=list(edge),
        receiver_ports=receiver,
        appointments=appointments,
    )


def run():
    return dict(
        scalar_sectors=classify(),
        collar_ramsey=source_reads(),
        flights=[flight(level) for level in (0, 1, 2)],
    )


def compare(actual, saved):
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError("Receipt fields differ")
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, list):
        for x, y in zip(actual, saved, strict=True):
            compare(x, y)
    elif isinstance(actual, float):
        np.testing.assert_allclose(actual, saved, atol=1e-13, rtol=2e-10)
    elif actual != saved:
        raise ValueError((actual, saved))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    observed = run()
    if args.write:
        RECEIPT.write_text(json.dumps(observed, indent=2, allow_nan=False) + "\n")
    else:
        compare(observed, json.loads(RECEIPT.read_text()))
    print(
        "Scalar-action dimensions:",
        [len(r["scalar_action_coin_vectors"]) for r in observed["scalar_sectors"]],
    )
    print(
        "Patch sensitivity:",
        [
            [p["blank_vs_marked_tv"] for p in f["appointments"]]
            for f in observed["flights"]
        ],
    )
