"""Unconditioned physical sign observations and reported modular energy under erasure."""

import argparse
from itertools import product
import json

import numpy as np
from scipy.linalg import expm
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import measure


def entropy(state):
    eig = np.linalg.eigvalsh(state)
    positive = eig[eig > 1e-13]
    return -float(positive @ np.log(positive))


def kmean(state, ref):
    vals, vecs = np.linalg.eigh(ref)
    diag = np.diag(vecs.conj().T @ state @ vecs).real
    assert sum(diag[vals <= 1e-13]) < 1e-10
    ok = vals > 1e-13
    return -float(diag[ok] @ np.log(vals[ok]))


def run():
    law = measure.trail.measure.LAW
    vac = measure.edge_operators(8)
    ref = expm(-measure.BETA * vac)
    z = np.trace(ref)
    ref /= z
    blank = np.diag(np.eye(11)[0])
    ids = [1 + law.LABELS.index((0, d)) for d in (1, -1)]
    source = np.diag(np.eye(11)[ids[0]])
    shots = 2**18
    rng = np.random.default_rng(191223)
    thinning = []
    for p in (0.0, 0.25, 0.5, 0.75, 1.0):
        new_ref = p * ref + (1 - p) * blank
        alt = p * source + (1 - p) * blank
        unflagged = kmean(alt - blank, new_ref)
        flagged = p * measure.BETA * 8.5
        expected = p * (measure.BETA * 8.5 + np.log1p(z * (1 - p) / p)) if p else 0
        np.testing.assert_allclose(unflagged, expected, atol=2e-12)
        # Mean difference in relative entropy needs each state's entropy as well.
        relative = unflagged - entropy(alt) + entropy(blank)
        thinning.append(
            dict(
                retained=p,
                flagged_increment=flagged,
                missing_as_blank_increment=unflagged,
                relative_entropy_contrast_missing_as_blank=relative,
            )
        )
    g = vac[np.ix_(ids, ids)]
    timing = []
    for n in (1, 2, 4, 16, 64):
        U = expm(-1j * np.pi * g / n)
        kernel = abs(U) ** 2
        state = np.linalg.matrix_power(kernel, n) @ np.array([1.0, 0.0])
        np.testing.assert_allclose(
            state[0], (1 + np.cos(np.pi / n) ** n) / 2, atol=4e-13
        )
        receipts = rng.multinomial(shots, state / state.sum())
        spectral = rng.binomial(shots, 0.5)
        timing.append(
            dict(
                read_intervals=n,
                native_sign_same=float(state[0]),
                sampled_sign_same=int(receipts[0]),
                sampled_energy_9=int(spectral),
                energy=float(np.diag(g) @ state),
            )
        )
    # Reflection on the signed difference is a conserved local observable
    # at the full native contact, including all nine possible ordinary words.
    full = measure.trail.measure.kernels()[2]
    reflect = np.eye(11, dtype=int)
    reflect[1:, 1:] = law._label_ops()[0]
    R = np.kron(reflect, np.eye(9, dtype=int))
    np.testing.assert_array_equal(R @ full, full @ R)
    # The full three-head solver can change its signed preparation without
    # recalibrating propagation. Equivalent reduced ordinary states for all
    # eight sign outcomes follow by independent applications of reflection.
    from . import flight

    board, path, h, _, _, tags, _, _ = flight.cell_operators()
    port = path[-1]
    reads = []
    for signs in product((1, -1), repeat=len(path) - 1):
        mem = np.ones(1)
        for u, v, s in zip(path, path[1:], signs):
            coin = 0 if u < v else 1
            mem = np.kron(mem, np.eye(4)[tags.index((coin, s))])
        init = np.kron(np.eye(board.n)[port], mem)
        evolved = expm_multiply(-4j * (h - 72 * sparse.eye(h.shape[0])), init).reshape(
            board.n, -1
        )
        reads.append(evolved @ evolved.conj().T)
    deviation = max(np.linalg.norm(r - reads[0], ord=2) for r in reads)
    assert deviation < 1e-10
    return dict(
        trials=shots,
        reference_partition=z,
        thinning=thinning,
        physical_observations=timing,
        max_ordinary_density_difference=deviation,
        sign_branches_verified=len(reads),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    path = measure.HERE / "observation.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        measure.check(json.loads(path.read_text()), result)
