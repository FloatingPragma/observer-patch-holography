"""Exact graph witnesses and a phase-sensitive observer of a public repair."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from ..__main__ import local_ops
from ..asynchronous import priority
from ..flight_confluence import native_double
from ..strict_pair import verify_columns


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def components(n, edges):
    left, right = zip(*edges)
    matrix = sparse.csr_matrix(
        (np.ones(len(left), dtype=int), (left, right)), shape=(n, n)
    )
    return connected_components(matrix, directed=False)


def fibers_and_certificate(board):
    owners = [[] for _ in range(board.n)]
    for edge in board.seams:
        for p in edge:
            owners[p].append(edge)
    winner = [min(incident, key=priority) for incident in owners]
    fibers = defaultdict(list)
    for p, edge in enumerate(winner):
        fibers[edge].append(p)

    witnesses = {}
    witness_counts = []
    for edge, block in sorted(fibers.items()):
        require(set(block) <= set(edge) and len(block) <= 2, "invalid winning fiber")
        if len(block) == 1:
            continue
        candidates = []
        for x, y in (edge, edge[::-1]):
            for z in board.neighbours[x]:
                if (
                    z != y
                    and len(fibers[winner[z]]) == 1
                    and z not in board.neighbours[y]
                ):
                    candidates.append((x, y, z))
        require(candidates, "missing single-neighbor certificate")
        witnesses[edge] = min(candidates)
        witness_counts.append(len(candidates))
    cross = [(p, q) for p, q in board.seams if winner[p] != winner[q]]
    count, label = components(board.n, cross)
    require(
        count == board.cells,
        "first derivative separates a different number of components",
    )
    require(
        all(
            len(set(label[start : start + 12])) == 1 for start in range(0, board.n, 12)
        ),
        "first derivative components differ from observer cells",
    )
    require(components(board.n, board.seams)[0] == 1, "native board is disconnected")
    pair_count = len(witnesses)
    n = board.n
    summary = dict(
        cells=board.cells,
        ports=n,
        binary_words=n * (n - 1) // 2,
        singles=len(fibers) - pair_count,
        pairs=pair_count,
        min_singleton_witnesses=min(witness_counts),
        max_singleton_witnesses=max(witness_counts),
        double_effect_dimensions=[n + 2 * pair_count, count, 1],
    )
    return summary, witnesses, winner


def observer(board, witnesses, winner):
    echo, _, write, _ = local_ops()
    exchange, contrast = native_law.local_term_parts()
    seam_h = exchange + KAPPA * contrast
    half = expm(-0.5j * np.pi * seam_h)
    full = expm(-1j * np.pi * seam_h)
    origins = [native_law.SRC.index(word) for word in ((2, 0), (0, 2))]
    mix = (half @ echo.T @ write)[np.ix_(origins, origins)]
    np.testing.assert_allclose(mix.conj().T @ mix, np.eye(2), atol=2e-13, rtol=0)
    np.testing.assert_allclose(
        full[np.ix_(origins, origins)], [[0, 1], [1, 0]], atol=1e-13, rtol=0
    )
    x, y, z = next(
        triple
        for edge, triple in sorted(witnesses.items())
        if edge in board.external
        and board.degree[triple[2]] == 5
        and edge == triple[:2]
    )
    require(
        winner[x] == winner[y] != winner[z],
        "public transcripts did not separate the phase",
    )

    launch = np.zeros((board.n, 2), dtype=complex)
    # Half exchange from y on (x,y), full swap x to z, then the addressed vacancy reference of z.
    launch[y] = half[origins[1], origins[1]]
    for col, duration in enumerate((np.pi / 64, 3 * np.pi / 64)):
        launch[z, col] = half[origins[0], origins[1]] * full[origins[1], origins[0]]
        launch[z, col] *= np.exp(-1j * duration * 4 * KAPPA)
    np.testing.assert_allclose(launch[z] / launch[y], [1, -1], rtol=0, atol=2e-13)

    def terminals(wave):
        xy = abs(mix @ wave[[x, y]]) ** 2
        complement = abs(np.delete(wave, [x, y], axis=0)) ** 2
        probs = np.vstack((xy, complement.sum(axis=0))).T
        np.testing.assert_allclose(probs.sum(axis=1), 1, atol=4e-13, rtol=0)
        return probs

    # A scalar contrast/reference term on all blank doubles provides only a phase.
    lap = board.laplacian()
    native = terminals(expm_multiply(-0.125j * lap, launch))
    control = terminals(launch)
    np.testing.assert_allclose(control, [[0.25, 0.25, 0.5]] * 2, rtol=0, atol=2e-13)
    verify_exchange = np.zeros(2)
    for col in range(2):
        # The conditional states after an immediate complete read are identical for the
        # two inputs. Each branch supports just one of their initial occupied positions.
        verify_exchange[col] = abs(launch[y, col]) ** 2 + abs(launch[z, col]) ** 2
    np.testing.assert_allclose(verify_exchange, 1, rtol=0, atol=2e-13)
    for p in (y, z):
        np.testing.assert_allclose(
            abs(launch[p, 0]) ** 2, abs(launch[p, 1]) ** 2, atol=2e-13, rtol=0
        )
    rng = np.random.default_rng(270926)
    shots = 2**18

    def bank(prob):
        return dict(
            probabilities=prob.round(12).tolist(),
            counts=[rng.multinomial(shots, row).tolist() for row in prob],
        )

    return dict(
        seam=[x, y],
        singleton=z,
        shots_per_phase=shots,
        flight_exposure="1/4",
        outcomes=["origin_x_analyzer", "origin_y_analyzer", "other_winner"],
        input_phase_order=["+", "-"],
        flight_first=bank(native),
        exchange_off=bank(control),
        sweep_first_phase_difference_exact=0,
    )


def run():
    verify_columns()
    native_double()
    boards = [federation(level) for level in range(3)]
    tables = [fibers_and_certificate(board) for board in boards]
    return dict(
        boards=[row[0] for row in tables], observer=observer(boards[0], *tables[0][1:])
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    results = run()
    if args.record is not None:
        with args.record.open("x") as output:
            output.write(json.dumps(results, indent=2) + "\n")
    elif results != json.loads(Path(__file__).with_name("results.json").read_text()):
        raise ValueError("record or certificate differs")
    print(json.dumps(results, indent=2))
