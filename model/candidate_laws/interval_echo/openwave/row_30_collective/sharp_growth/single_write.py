"""All-start response of a locally repaired strand, with every contact read once."""

import argparse
from itertools import combinations, product

import numpy as np

from model.executable_model import federation
from .. import experiment
from .growth import SEED, MATERIAL_SEED
from .measure import HERE, write_data

SHOTS = 2**34
CLOCK = 2**30
TIMES = np.array([0.0, 0.5, 1.0, 2.0, 4.0])


def enroll():
    for level in range(3):
        board = federation(level)
        for closed, route in experiment.grow(board, MATERIAL_SEED):
            if not closed and len(route) == 6:
                return level, route
    raise ValueError("No open six-dimer component")


def columns(n):
    """Physical matrices of i gamma_a gamma_b in the two fixed end conventions."""
    states = np.arange(2**n)
    spin = 1 - 2 * ((states[:, None] >> np.arange(n)) & 1)
    ext = np.column_stack([np.ones(len(spin), int), spin, -np.ones(len(spin), int)])
    wall = -ext[:, :-1] * ext[:, 1:]
    operators = []
    for a, b in combinations(range(2 * (n + 1)), 2):
        w = wall.copy()
        coefficient = np.full(len(spin), 1j)
        for j, kind in (divmod(b, 2), divmod(a, 2)):
            coefficient *= np.prod(w[:, :j], axis=1)
            if kind:
                coefficient *= 1j * w[:, j]
            w[:, j] *= -1
        target = (1 - np.cumprod(-w, axis=1)[:, :n]) // 2 @ (1 << np.arange(n))
        op = np.zeros((len(spin), len(spin)), complex)
        op[target, states] = coefficient
        operators.append(op)
    ops = np.stack(operators)
    np.testing.assert_allclose(
        ops.reshape(len(ops), -1).conj() @ ops.reshape(len(ops), -1).T / len(spin),
        np.eye(len(ops)),
        atol=1e-13,
    )
    return spin, ops


def traces():
    n = 6
    spin, ops = columns(n)
    source = np.eye(1 << n)[np.arange(1 << n) ^ (1 << (n // 2))]
    # Check the physical dimer flip is contained in this basis with unit length.
    np.testing.assert_allclose(
        np.sum(abs(ops.reshape(len(ops), -1).conj() @ source.ravel() / len(spin)) ** 2),
        1,
    )
    patterns = list(product((0, 1), repeat=n - 1))
    transfers = []
    state = np.arange(len(spin))
    for bits in patterns:
        gaps = 4 * experiment.K + np.asarray(bits)
        # Facing loads coincide at 11 only for z_left=-1,z_right=+1.
        facing = (spin[:, :-1] == -1) & (spin[:, 1:] == 1)
        energy = (1 - facing) @ gaps
        ham = np.diag(energy)
        for i in range(n):
            ham[state, state ^ (1 << i)] -= 0.5
        e, v = np.linalg.eigh(ham)
        source_e = v.T @ source @ v
        rows = []
        for t in TIMES:
            op = v @ (source_e * np.exp(-1j * t * (e[:, None] - e[None, :]))) @ v.T
            cs = ops.reshape(len(ops), -1).conj() @ op.ravel() / len(spin)
            np.testing.assert_allclose(cs.imag, 0, atol=2e-14)
            rows.append(cs.real)
        transfers.append(rows)
    return patterns, np.asarray(transfers)


def summarize(counts):
    totals = counts.sum(axis=-1)
    diff = counts[..., 0] + counts[..., 3] - counts[..., 1] - counts[..., 2]
    signals = diff / totals
    pooled = diff.sum(axis=0) / totals.sum(axis=0)
    pooled_norm = np.sum(pooled**2, axis=-1)
    norm = np.sum(signals**2, axis=-1)
    # Joint Hoeffding + union bound. All underlying signal vectors are
    # orthogonal projections of a norm-one observable (or convex mixtures).
    delta = np.sqrt(
        2 * np.log(2 * (signals.size + pooled.size) / 0.01) / np.min(totals)
    )
    radius = np.sqrt(signals.shape[-1]) * delta
    return dict(
        pattern_norm=norm.tolist(),
        mean_conditional_norm=norm.mean(axis=0).tolist(),
        pooled_norm=pooled_norm.tolist(),
        norm_error_bound=2 * radius + radius**2,
    )


def run():
    level, route = enroll()
    patterns, signals = traces()
    rng = np.random.default_rng([SEED, 3, 19])
    # Each measurement has its own complete first-coin census. All 32 records
    # are equiprobable on independent sharp writes.
    totals = rng.multinomial(
        SHOTS, np.ones(len(patterns)) / len(patterns), size=signals.shape[1:]
    ).transpose(2, 0, 1)
    first = rng.binomial(totals, 0.5)
    agree = np.clip((1 + signals) / 2, 0, 1)
    pp = rng.binomial(first, agree)
    mm = rng.binomial(totals - first, agree)
    counts = np.stack([pp, first - pp, totals - first - mm, mm], axis=-1)
    hands = rng.binomial(
        CLOCK, (1 + np.stack([np.cos(TIMES), np.sin(TIMES)], axis=1)) / 2
    )
    angles = np.unwrap(np.arctan2(*(2 * hands.T[::-1] / CLOCK - 1)))
    res = summarize(counts)
    return (counts, hands), dict(
        level=level,
        route=route,
        patterns=patterns,
        times_from_hand=angles.tolist(),
        shots=SHOTS,
        clock_shots=CLOCK,
        **res,
        exact_pooled_norm=np.sum(signals.mean(axis=0) ** 2, axis=-1).tolist(),
        exact_mean_conditional_norm=np.mean(
            np.sum(signals**2, axis=-1), axis=0
        ).tolist(),
    )


if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--record", action="store_true")
    opt = args.parse_args()
    (counts, hands), result = run()
    path = HERE / "single_write.npz"
    rec = HERE / "single_write.json"
    if opt.record:
        with path.open("xb") as f:
            np.savez_compressed(f, counts=counts, hands=hands)
        with rec.open("xb") as f:
            f.write(write_data(result) + b"\n")
    else:
        with np.load(path, allow_pickle=False) as data:
            np.testing.assert_array_equal(counts, data["counts"])
            np.testing.assert_array_equal(hands, data["hands"])
        assert rec.read_bytes() == write_data(result) + b"\n"
    print(
        result["pooled_norm"],
        result["mean_conditional_norm"],
        result["norm_error_bound"],
    )
