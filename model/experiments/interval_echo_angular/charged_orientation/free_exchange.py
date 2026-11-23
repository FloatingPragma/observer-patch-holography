"""Unrestricted free return of two identically signed orbital packets on the 20-cell host."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from model.executable_model import federation
from .measure import (
    HERE,
    EPS,
    SEED,
    SHOTS,
    candidate,
    code_columns,
    counts,
    design,
    step,
)

TIMES = (0, 1, 2, 4, 8, 16)


def pair_basis(n):
    left, right = np.triu_indices(n, 1)
    lookup = np.full((n, n), -1, dtype=int)
    labels = np.arange(len(left))
    lookup[left, right] = labels
    lookup[right, left] = labels
    return left, right, lookup


def generator(board, sign, odd, basis, kappa=candidate.KAPPA):
    """Project blank-seam H onto two equal defects relative to the unit background.

    Each moving defect has the degree-dependent load correction from the one-defect
    generator. A seam shared by the two equal defects has neither a contrast nor a hop.
    Relative to adding their two individual contrasts, subtract twice κ there.
    The uniform background energy and 12κ fixed-number offsets are common to both arms.
    """
    left, right, lookup = basis
    indices = np.arange(len(left))
    is_seam = np.any(board.table[left] == right[:, None], axis=1)
    degree = board.degree[left] + board.degree[right] - 2 * is_seam
    diagonal = (
        degree / 2
        + 2 * sign * kappa * (12 - board.degree[left] - board.degree[right])
        - 2 * kappa * is_seam
    )
    rows, cols, values = [indices], [indices], [diagonal]
    for mover, parked in ((left, right), (right, left)):
        for neighbor in board.table[mover].T:
            valid = (neighbor >= 0) & (neighbor != parked)
            offset = (mover[valid] < parked[valid]) != (neighbor[valid] < parked[valid])
            phase = 1 - 2 * offset if odd else np.ones(valid.sum())
            rows.append(lookup[neighbor[valid], parked[valid]])
            cols.append(indices[valid])
            values.append(-0.5 * phase)
    h = sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(left), len(left)),
    ).tocsr()
    if (h - h.T).nnz:
        raise ArithmeticError("The free generator must be self-adjoint")
    return h


def prepare(src, dst, basis):
    left, right, lookup = basis
    modes = code_columns([0, 1, 2], 3)
    initial = np.zeros((len(left), 4), complex)
    for i, a in enumerate(src):
        for j, b in enumerate(dst):
            initial[lookup[a, b]] = np.kron(modes[i], modes[j])
    return initial


def read(k, rng):
    swap = np.eye(4)[[0, 2, 1, 3]]
    accepted = float(np.linalg.norm(k) ** 2 / 4)
    own, cross = [
        float(abs(np.vdot(pattern, k) / 4) ** 2) for pattern in (np.eye(4), swap)
    ]
    chars = abs(k) ** 2
    return dict(
        accepted=accepted,
        own_reference=own,
        crossed_reference=cross,
        reference_counts=[
            counts([p, accepted - p, 1 - accepted], rng) for p in (own, cross)
        ],
        character_counts=[counts([*row, 1 - sum(row)], rng) for row in chars.T],
    )


def acquire(board, basis, initial, sign, odd):
    h = generator(board, sign, odd, basis)
    states = initial.copy()
    previous = 0
    kernels = []
    for t in TIMES:
        states = step(h, states, t - previous)
        previous = t
        np.testing.assert_allclose(
            states.conj().T @ states, np.eye(4), atol=2e-11, rtol=0
        )
        kernels.append(initial.conj().T @ states)
    return kernels


def run():
    rng = np.random.default_rng(SEED + 1)
    board = federation(0)
    owner, destinations, eligible = design(board)
    src, dst = eligible[owner][1], eligible[destinations[0]][1]
    basis = pair_basis(board.n)
    initial = prepare(src, dst, basis)
    signs = []
    for sign in (-1, 1):
        native = acquire(board, basis, initial, sign, False)
        fermion = acquire(board, basis, initial, sign, True)
        results = []
        for t, n, f in zip(TIMES, native, fermion):
            registered, comparator = read(n, rng), read(f, rng)
            deltas = []
            for panel in ("reference_counts", "character_counts"):
                delta = (
                    np.asarray(registered[panel]) - np.asarray(comparator[panel])
                ) / SHOTS
                for (setting, outcome), value in np.ndenumerate(delta):
                    deltas.append((abs(value), panel, int(setting), int(outcome)))
            largest, panel, setting, outcome = max(deltas)
            results.append(
                dict(
                    exposure=t,
                    registered=registered,
                    odd_worldlines=comparator,
                    even_return_weight=float(np.linalg.norm((n + f) / 2) ** 2 / 4),
                    odd_return_weight=float(np.linalg.norm((n - f) / 2) ** 2 / 4),
                    largest_count_difference=float(largest),
                    separated=bool(largest > 2 * EPS),
                    setting=[panel, setting, outcome],
                )
            )
        signs.append(dict(sign=sign, appointments=results))
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "FREE_EXCHANGE.md").read_bytes()
        ).hexdigest(),
        cells=board.cells,
        ports=board.n,
        basis_size=len(basis[0]),
        source=src,
        receiver=dst,
        shots=SHOTS,
        seed=SEED + 1,
        difference_radius=2 * EPS,
        arms=signs,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    result = run()
    with path.open("x") as out:
        json.dump(result, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
