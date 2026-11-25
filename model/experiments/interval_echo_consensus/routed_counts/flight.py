"""One vacancy and one double propagating on the fully open registered screen."""

from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.candidate import KAPPA
from .material import Routing


def open_generator(board):
    n = board.n
    i, j = np.where(~np.eye(n, dtype=bool))
    size = n * n
    rows = []
    cols = []
    adjacent = np.zeros((n, n), int)
    neighbours = np.full((n, 6), -1, int)
    for a, nb in enumerate(board.neighbours):
        neighbours[a, : len(nb)] = nb
        adjacent[a, nb] = 1
    for k in range(6):
        dest = neighbours[i, k]
        mask = dest >= 0
        rows.append((i * n + j)[mask])
        cols.append((dest * n + np.where(dest == j, i, j))[mask])
        dest = neighbours[j, k]
        mask = (dest >= 0) & (dest != i)
        rows.append((i * n + j)[mask])
        cols.append((i * n + dest)[mask])
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    hopping = coo_matrix(
        (-np.ones(len(rows)) / 2, (rows, cols)), shape=(size, size)
    ).tocsr()
    diag = np.zeros(size)
    # Relative to the all-one material, including the degree-five private references.
    diag[i * n + j] = KAPPA * (
        12 + 2 * (board.degree[i] - board.degree[j]) + 2 * adjacent[i, j]
    )
    classical = diag.copy()
    diag += -np.asarray(hopping.sum(axis=1)).ravel()
    h = hopping + diags(diag)
    if (h - h.T).nnz:
        raise ArithmeticError("Asymmetric native exchange")
    return h, classical


def read_effects(board, reverse, pairs):
    """Pull the strict chord predicate back through every actual routed permutation."""
    routing = Routing(board)
    chord = routing.chords[-1] if reverse else routing.chords[0]
    orientation = np.empty((len(pairs), 2), dtype=int)
    available = {tuple(p) for p in pairs}
    for row in reversed(range(len(pairs))) if reverse else range(len(pairs)):
        route = routing.legs(pairs[row], chord)
        a, b = chord
        for i, j in reversed(route):
            if a == i:
                a = j
            elif a == j:
                a = i
            if b == i:
                b = j
            elif b == j:
                b = i
        # Earlier quiet probes have removed all configurations that would have
        # written at an earlier request. Check that nothing is counted twice.
        key = tuple(sorted((a, b)))
        if key not in available:
            raise ValueError("Route repeated a previously tested material pair")
        available.remove(key)
        orientation[row] = (a, b)
    if available:
        raise ValueError("The no-mark branch contains an untested discrepancy")
    return chord, orientation.T


def cohort(board, rng):
    h, classical = open_generator(board)
    n = board.n
    v = 0
    pair = next(
        (i, j)
        for i, j in board.seams
        if i // 12 != 0
        and j // 12 == i // 12
        and board.degree[i] == board.degree[j]
        and i not in board.neighbours[v]
        and j not in board.neighbours[v]
    )
    start = np.zeros((n * n, 2), complex)
    for c, d in enumerate(pair):
        start[n * v + d, c] = 1
    t = 0.25
    native = expm_multiply(-1j * t * h, start)
    disabled = np.exp(-1j * t * classical[:, None]) * start
    offers = list(combinations(range(n), 2))
    u, w = np.array(offers).T
    effects = {r: read_effects(board, r, offers) for r in (False, True)}
    target = int(np.flatnonzero((u == v) & (w == pair[0]))[0])
    shots = 2**20
    arrays = {}
    summaries = []
    for mode, amplitudes in [
        ("native", native),
        ("which_source", native),
        ("exchange_parked_during_flight", disabled),
    ]:
        for reverse in (False, True):
            for basis in ("Z", "X"):
                # Third component has all ports at unit load and reference |+>.
                # It never writes; its projection is disjoint from every successful record.
                amplitudes0, amplitudes1 = amplitudes.T
                if basis == "Z" or mode == "which_source":
                    p0, p1 = abs(amplitudes0) ** 2 / 3, abs(amplitudes1) ** 2 / 3
                    if basis == "X":
                        p0 = p1 = (p0 + p1) / 2
                else:
                    p0 = abs(amplitudes0 + 1j * amplitudes1) ** 2 / 6
                    p1 = abs(amplitudes0 - 1j * amplitudes1) ** 2 / 6
                columns = []
                for a in (p0, p1):
                    mat = a.reshape(n, n)
                    a0, a1 = effects[reverse][1]
                    columns.append(mat[a0, a1] + mat[a1, a0])
                probs = np.column_stack(columns)
                no_mark = [1 / 6, 1 / 6] if basis == "Z" else [1 / 3, 0.0]
                probs = np.vstack((probs, no_mark))
                np.testing.assert_allclose(probs.sum(), 1.0, atol=1e-12)
                name = f"{mode}_{int(reverse)}_{basis}"
                counts = rng.multinomial(shots, probs.ravel() / probs.sum()).reshape(
                    probs.shape
                )
                arrays[name] = counts
                summaries.append(
                    dict(
                        receiving_chord=list(effects[reverse][0]),
                        mode=mode,
                        reverse=reverse,
                        reference=basis,
                        marker=int(counts[:-1].sum()),
                        quiet=counts[-1].tolist(),
                        selected_pair_reference=counts[target].tolist(),
                        complete_marker_reference=counts[:-1].sum(axis=0).tolist(),
                    )
                )
    return (
        dict(
            cells=board.cells,
            ports=n,
            starts_per_bank=shots,
            wait=t,
            vacancy=v,
            double_origins=list(map(int, pair)),
            target_row=target,
            output_order="lexicographic unordered port pairs, then no mark; reference outcomes 0,1",
            cohorts=summaries,
        ),
        arrays,
    )
