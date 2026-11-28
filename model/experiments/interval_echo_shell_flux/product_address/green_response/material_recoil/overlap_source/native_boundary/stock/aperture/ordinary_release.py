"""Compare the same cold probe on a written seam and on one mobile ordinary load."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import cold_release as cold
from .cold_release import probability_band, difference

HERE = cold.HERE


def two_loads(board, k, marked=()):
    # i carries two units while j carries one. The blank native exchange
    # exactly preserves those two different occupations, including swaps.
    n = board.n
    i = np.repeat(np.arange(n), n - 1)
    j = np.tile(np.arange(n - 1), n)
    j += j >= i

    def index(a, b):
        return a * (n - 1) + b - (b > a)

    np.testing.assert_array_equal(index(i, j), np.arange(len(i)))
    left, right = [], []
    all_ports = np.arange(n)
    adjacent = np.zeros(len(i), dtype=int)
    incidence = np.zeros(n, dtype=int)
    for a, b in board.seams:
        if (a, b) in marked:
            incidence[[a, b]] += 1
            continue
        spectator = all_ports[(all_ports != a) & (all_ports != b)]
        left.extend((index(a, spectator), index(spectator, a), np.array([index(a, b)])))
        right.extend(
            (index(b, spectator), index(spectator, b), np.array([index(b, a)]))
        )
        adjacent[[index(a, b), index(b, a)]] = 1
    u, v = np.concatenate(left), np.concatenate(right)
    adj = sparse.csr_matrix(
        (np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(len(i), len(i))
    )
    degree = np.asarray(adj.sum(axis=1)).ravel()
    # Subtract the common isolated cost 6*k*(1²+2²). The contrast
    # cross term on a seam occupied by both loads is -2*k*(1*2).
    saved = k * (4 * incidence[i] + incidence[j] + 4 * adjacent)
    h = sparse.diags(degree / 2 - saved, format="csr") - adj / 2
    return h, i, j


def release_acceleration(h, mask, indices, shell):
    twice = (2 * h).astype(np.int64)
    assert (twice != 2 * h).nnz == 0
    a = mask.astype(np.int64)
    b = twice @ a
    numerator = b * b - a * (twice @ b)
    per_shell = np.zeros(shell.max() + 1, dtype=np.int64)
    np.add.at(per_shell, shell[indices], numerator)
    return per_shell / (2 * a.sum())


def check_joint_columns(k):
    law = cold.emission.LAW
    four = cold.emission.four_h(k)
    heads = 9 * np.array([i + 1 for i, (_, d) in enumerate(law.LABELS) if abs(d) == 2])
    native = four[np.ix_(heads, heads)]
    for a, b in ((0, 0), (2, 0), (0, 2), (1, 0), (0, 1), (2, 1), (1, 2)):
        word = law.SRC.index((a, b))
        cols = heads + word
        target = np.zeros_like(four[:, cols])
        target[cols] = native
        np.testing.assert_array_equal(four[:, cols], target)
        # Blank contact: unit exchange and the full squared contrast.
        expected = np.zeros(99, dtype=int)
        expected[word] += 2 + 4 * k * (a - b) ** 2
        expected[law.SRC.index((b, a))] -= 2
        np.testing.assert_array_equal(four[:, word], expected)


def experiment():
    board = cold.federation(0)
    edge = tuple(sorted(min(board.external)))
    h0, _, _ = cold.motion(board, [], 2, 8)
    hd, _, _ = cold.motion(board, [edge], 2, 8)
    ho, i, j = two_loads(board, 8)
    hj, ji, jj = two_loads(board, 8, [edge])
    np.testing.assert_array_equal(i, ji)
    np.testing.assert_array_equal(j, jj)
    check_joint_columns(8)
    shell = cold.distances(board, edge)
    n = board.n
    rng = np.random.default_rng(501928)
    t = np.array(cold.TIMES)
    rows = {}
    for name, h in (("blank", h0), ("ordinary", ho), ("written", hd), ("joint", hj)):
        mobile = name in ("ordinary", "joint")
        positions = i if mobile else np.arange(n)
        laws = np.zeros((len(t), shell.max() + 1))
        source_stays = np.zeros(len(t))
        accel = np.zeros(shell.max() + 1)
        starts = []
        for z in edge:
            start = (j == z) if mobile else (np.arange(n) != z)
            # Both lots have the same conditional probe marginal and no
            # imaginary coherences. The ordinary source subsequently moves.
            assert start.sum() == n - 1
            starts.append(start / np.sqrt(start.sum()))
            accel += release_acceleration(h, start, positions, shell) / 2
        state = np.array(starts).T.astype(complex)
        last = 0
        for r, time in enumerate(t):
            state = expm_multiply((-1j * (time - last)) * h, state)
            last = time
            prob = (abs(state) ** 2).mean(axis=1)
            laws[r] = np.bincount(shell[positions], weights=prob)
            if mobile:
                source_stays[r] = sum(prob[np.isin(j, edge)])
        np.testing.assert_allclose(laws.sum(axis=1), 1, atol=1e-11)
        rows[name] = dict(
            counts=[rng.multinomial(cold.STARTS, p / p.sum()).tolist() for p in laws],
            probabilities=laws.tolist(),
            acceleration=accel.tolist(),
        )
        if mobile:
            rows[name]["ordinary_on_initial_seam"] = source_stays.tolist()
    return dict(
        board_ports=n,
        seam=edge,
        seed=501928,
        times=t.tolist(),
        shots=cold.STARTS,
        arms=rows,
    )


def decode(record):
    arms = record["arms"]
    bands = {}
    for k, v in arms.items():
        assert all(sum(row) == record["shots"] and min(row) >= 0 for row in v["counts"])
        bands[k] = probability_band(v["counts"], record["shots"])
    ordinary = difference(bands["ordinary"], bands["blank"])
    held = difference(bands["written"], bands["blank"])
    a, b = ordinary[0, 0]
    x, y = held[0, 0]
    gain = None
    misses = []
    if a > 0 or b < 0:
        ratios = [v / u for v in (x, y) for u in (a, b)]
        gain = [min(ratios), max(ratios)]
        for i in range(len(held)):
            for r in range(held.shape[1]):
                forecasts = [f * v for f in gain for v in ordinary[i, r]]
                if max(forecasts) < held[i, r, 0] or min(forecasts) > held[i, r, 1]:
                    misses.append([i, r])
    opposite = []
    for i, j in np.ndindex(held.shape[:2]):
        if (ordinary[i, j, 0] > 0 and held[i, j, 1] < 0) or (
            ordinary[i, j, 1] < 0 and held[i, j, 0] > 0
        ):
            opposite.append([i, j])
    added = difference(bands["joint"], bands["ordinary"])
    residual = difference(added, held)
    nonadditive = np.argwhere((residual[..., 0] > 0) | (residual[..., 1] < 0))
    return dict(
        joint_minus_ordinary=added.tolist(),
        joint_additivity_residual=residual.tolist(),
        resolved_nonadditivity=nonadditive.tolist(),
        ordinary_contrasts=ordinary.tolist(),
        written_contrasts=held.tolist(),
        first_support_ratio=gain,
        held_disagreements=misses,
        opposite_signs=opposite,
    )


if __name__ == "__main__":
    cmd = argparse.ArgumentParser(description=__doc__)
    modes = cmd.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cmd.parse_args()
    path = HERE / "ordinary_release.json"
    if args.counts_only:
        data = json.loads(path.read_bytes())
        cold.apparatus.check(data["decoded"], decode(data))
    else:
        record = experiment()
        data = dict(
            **record,
            protocol_sha256=sha256((HERE / "COLD_RELEASE.md").read_bytes()).hexdigest(),
            mixed_protocol_sha256=sha256(
                (HERE / "MIXED_RELEASE.md").read_bytes()
            ).hexdigest(),
            decoded=decode(record),
        )
        if args.write:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            cold.apparatus.check(
                json.loads(path.read_bytes()), json.loads(json.dumps(data))
            )
    print(data["decoded"]["first_support_ratio"], data["decoded"]["opposite_signs"])
