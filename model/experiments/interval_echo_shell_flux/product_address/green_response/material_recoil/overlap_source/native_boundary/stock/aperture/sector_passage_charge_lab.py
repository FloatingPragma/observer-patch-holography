"""Count the connected two-home source on complete, independently prepared registers."""

import math

import numpy as np

from . import peer_register_lab, sector_sampling
from .sector_passage_source import mul


def read(words, subset, full, mats, log_odds):
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(words))
    tables = peer_register_lab.acquisition.returns.ticket_tables(mats)
    x, y = (words[:, [0, 3]] > 0).astype(int).T
    cells = 2 * x + y
    flags = np.column_stack((x, y, (x - x.mean()) * (y - y.mean())))
    ranges = np.ptp(flags, axis=0)
    rows = []
    counts = np.zeros((2, 4), dtype=np.int64)
    rng = np.random.default_rng(390415)

    def observe(members, original):
        b = len(members)
        selected = members[
            np.arange(b)[:, None], rng.integers(members.shape[1], size=(b, 128))
        ].ravel()
        for k, start in enumerate((selected, original[: len(selected)])):
            end, _ = sector_sampling.report_from(tables, subset, native, rng, 8, start)
            counts[k] += np.bincount(cells[end], minlength=4)
            if k == 0:
                rows.append(flags[end].reshape(b, 128, 3).mean(1))

    for pulse, reg in (
        (601331, 955201),
        (771125, 944109),
        (210927, 917451),
        (534971, 451143),
    ):
        peer_register_lab.run(observe=observe, pulse_seed=pulse, register_seed=reg)
    scores = np.concatenate(rows)
    n = len(scores) * 128
    assert np.all(counts.sum(1) == n)
    by_cell = np.array([flags[np.flatnonzero(cells == i)[0]] for i in range(4)])
    ordinary = counts[1] @ by_cell / n
    var = counts[1] @ ((by_cell - ordinary) ** 2) / (n - 1)
    diff = scores.mean(0) - ordinary
    ell = math.log(4 * 6 / 0.006)

    def bound(v, k):
        return np.sqrt(2 * ell * v / k) + 7 * ranges * ell / (3 * (k - 1))

    pad = bound(scores.var(0, ddof=1), len(scores)) + bound(var, n)
    ends = np.stack((diff - pad, diff + pad), axis=1)
    cov = ends[2] - mul(ends[0], ends[1])[::-1]
    q = mul(-np.array(log_odds)[::-1], cov)
    return {
        "registers": len(scores),
        "trials_per_arm": n,
        "joint_counts": counts.tolist(),
        "local_contrasts": diff.tolist(),
        "covariance": float(diff[2] - diff[0] * diff[1]),
        "covariance_interval": cov.tolist(),
        "connected_modular_nats": q.tolist(),
        "combined_coverage": 0.993,
    }
