"""Compare clock and release responses of an ordinary medium to a written seam."""

import argparse
from hashlib import sha256
import json
from math import exp, factorial

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from . import cold_release as cold
from . import ordinary_release as ordinary

HERE = cold.HERE
K = 8
SHOTS = 2**44
CLOCK_STEP = 1 / 512
RELEASE_TIMES = (0, 1 / 32, 1 / 16)
# Simultaneous Hoeffding allowance (over all enrolled binary events).
RISK = np.log(2 * 10000 / 0.01)


def thermal(h):
    values, modes = eigh(h.toarray())
    weights = np.exp(-cold.apparatus.BETA * values)
    diagonal = modes**2 @ weights
    return diagonal / sum(diagonal), diagonal


def release_bank(h, i, j, shells):
    n = len(shells)
    rows = np.arange(len(i))
    start = sparse.csr_matrix((np.ones(len(i)), (rows, j)), shape=(len(i), n))
    twice = (2 * h).astype(np.int64)
    assert (twice != 2 * h).nnz == 0
    wave = start.astype(np.int64)
    hw = twice @ wave
    per_address = hw.multiply(hw) - wave.multiply(twice @ hw)
    project = sparse.csr_matrix((np.ones(len(i), dtype=np.int64), (shells[i], rows)))
    numerators = (project @ per_address).toarray().cumsum(axis=0)
    # Independent finite-time propagation with all ordinary positions kept.
    banks = np.empty((len(RELEASE_TIMES), int(shells.max()) + 1, n))
    for lo in range(0, n, 12):
        hi = min(n, lo + 12)
        psi = start[:, lo:hi].toarray().astype(complex) / np.sqrt(n - 1)
        last = 0
        for t, time in enumerate(RELEASE_TIMES):
            if time != last:
                psi = expm_multiply((-1j * (time - last)) * h, psi)
            banks[t, :, lo:hi] = project @ (abs(psi) ** 2)
            last = time
    np.testing.assert_allclose(banks.sum(axis=1), 1, atol=2e-13, rtol=0)
    return numerators, banks


def ramsey_moments(h, hb, i, j, site, order=8):
    # V maps each successful bath address to the newly inserted double at x.
    n = hb.shape[0]
    rows = np.flatnonzero(i == site)
    v = sparse.csr_matrix((np.ones(len(rows)), (rows, j[rows])), shape=(len(i), n))
    powers = [v]
    for _ in range((order + 1) // 2):
        powers.append(h @ powers[-1])
    pinned = [
        (powers[k // 2].T @ powers[(k + 1) // 2]).toarray() for k in range(order + 1)
    ]
    free = [np.eye(n)]
    for _ in range(order):
        free.append(hb @ free[-1])
    # Only diagonal bath elements remain after the declared sharp address read.
    out = np.zeros((order + 1, n), complex)
    for k in range(order + 1):
        for a in range(k + 1):
            b = k - a
            diagonal = np.einsum("ij,ji->i", free[a], pinned[b])
            out[k] += (1j**a * (-1j) ** b) / (factorial(a) * factorial(b)) * diagonal
    return out


def degree_clock(board, mask, source, x):
    near = [v for v in board.neighbours[x] if tuple(sorted((x, v))) not in mask]
    adjacent = source[near].sum() / (1 - source[x])
    # The reference holds the same ordinary unit; one shared exchange costs
    # half a unit less than two independent incident exchanges.
    m = len(board.neighbours[x]) - len(near)
    return len(near) / 2 - 4 * K * m - (4 * K + 0.5) * adjacent


def clock_samples(moments, weights, x, rng):
    # Attempted insertions up to the fixed number of successful Ramsey trials.
    failure = int(rng.negative_binomial(SHOTS, 1 - weights[x])) if weights[x] else 0
    p = weights.copy()
    p[x] = 0
    p /= p.sum()
    c = moments @ p
    result = []
    for dt in (CLOCK_STEP, CLOCK_STEP / 2):
        fringe = -np.polynomial.polynomial.polyval(dt, c).imag
        result.append(int(rng.binomial(SHOTS, (1 + fringe) / 2)))
    return dict(plus=result, failed=failure, success=SHOTS)


def rate_interval(record):
    # Richardson extrapolation uses |f^(5)| <= M^5 for the two-arm overlap.
    # Native Gershgorin gives joint norm <=64 and single norm <=8 on this
    # one-seam board. Check the bound against both actual sparse generators.
    M = 80
    h = CLOCK_STEP
    signs = 2 * np.asarray(record["plus"]) / SHOTS - 1
    estimate = (8 * signs[1] - signs[0]) / (3 * h)
    # In bounding the two remainders separately the coefficient is 1/288.
    trunc = M**5 * h**4 / 288
    poly = 2 * exp(h * M) * (h * M) ** 9 / factorial(9)
    error = 3 * (np.sqrt(2 * RISK / SHOTS) + poly) / h + trunc
    return np.array([estimate - error, estimate + error])


def band(value, error):
    return np.array([max(0, value - error), min(1, value + error)])


def conditional_interval(count, x, near):
    n = sum(count)
    a = sum(count[i] for i in near)
    # Conditional on the number of successful (j != x) draws; denominator
    # is observed exactly. Bound the neighbor indicator in that subexperiment.
    successes = n - count[x]
    return band(a / successes, np.sqrt(RISK / (2 * successes)))


def delta(a, b):
    return np.array([a[0] - b[1], a[1] - b[0]])


def signed_product(a, b):
    p = np.outer(a, b)
    return np.array([p.min(), p.max()])


def decode(record):
    site = record["stations"]
    scans = record["clock"]
    gain = delta(
        rate_interval(scans["far_standard"]), rate_interval(scans["near_standard"])
    )
    board = cold.federation(0)
    marks, _ = cold.read_map(board, record["spectral_tape"])
    shells = cold.distances(board, [v for edge in marks for v in edge])
    counts = record["source_counts"]
    for bank in counts.values():
        assert len(bank) == board.n and min(bank) >= 0 and sum(bank) == 2**40
    predicted = []
    for x, near in enumerate(board.neighbours):
        live = [v for v in near if tuple(sorted((x, v))) not in marks]
        incidence = len(near) - len(live)
        ordinary = delta(
            conditional_interval(counts["written"], x, live),
            conditional_interval(counts["blank"], x, near),
        )
        predicted.append(-signed_product(gain, incidence + ordinary)[::-1])
    out = []
    for row in site:
        x = row["address"]
        observed = delta(
            rate_interval(scans["written"][str(x)]),
            rate_interval(scans["blank"][str(x)]),
        )
        out.append(
            dict(
                hops=row["hops"],
                frequency_change=observed.tolist(),
                from_ordinary=predicted[x].tolist(),
            )
        )
    # Use point estimates for the astronomical continuation. The fitted
    # astronomical amplitude absorbs the overall frequency unit.
    shell_read = np.zeros(shells.max() + 1)
    for x, near in enumerate(board.neighbours):
        live = [v for v in near if tuple(sorted((x, v))) not in marks]

        def fraction(k, vertices):
            return sum(counts[k][v] for v in vertices) / (sum(counts[k]) - counts[k][x])

        shell_read[shells[x]] += gain.mean() * (
            len(near) - len(live) + fraction("written", live) - fraction("blank", near)
        )
    kernel = np.asarray(record["ordinary_integer_kernel"]) / (2 * (board.n - 1))
    d = np.array(counts["written"]) / sum(counts["written"]) - np.array(
        counts["blank"]
    ) / sum(counts["blank"])
    center = kernel @ d
    width = 2 * np.ptp(kernel, axis=1) * np.sqrt(RISK / (2 * 2**40))
    exterior = np.stack([center - width, center + width], axis=1)[1:]
    inside = {}
    for name, tape in record["release_counts"].items():
        assert all(sum(row) == SHOTS and min(row) >= 0 for row in tape)
        values = np.cumsum(np.asarray(tape), axis=1) / SHOTS
        inside[name] = values - values[0]
    radius = np.sqrt(RISK / (2 * SHOTS))
    shifts = inside["thermal"] - inside["blank"]
    imitation = inside["imitated"] - inside["blank"]
    return dict(
        ordinary_neighbor_gain=gain.tolist(),
        exterior_acceleration_from_counts=exterior.tolist(),
        shell_clock_slowdown_from_counts=shell_read.tolist(),
        clocks=out,
        inside_changes=shifts.tolist(),
        imitated_changes=imitation.tolist(),
        two_arm_change_error=4 * radius,
        direct_changes=(inside["direct"] - inside["blank"]).tolist(),
    )


def run():
    rng = np.random.default_rng(914683)
    board = cold.federation(0)
    seam = tuple(sorted(min(board.external)))
    shells = cold.distances(board, seam)
    selections = [int(np.flatnonzero(shells == r)[0]) for r in range(shells.max() + 1)]
    source_counts = {}
    hamiltonians = {}
    diag = {}
    weights = {}
    for name, mask in (("blank", []), ("written", [seam])):
        hb, _, _ = cold.motion(board, mask, 1, K)
        h, i, j = ordinary.two_loads(board, K, mask)
        hamiltonians[name] = (hb, h)
        weights[name], diag[name] = thermal(hb)
        source_counts[name] = rng.multinomial(2**40, weights[name]).tolist()
    ordinary.check_joint_columns(K)
    tape = cold.source_tape(board, list(seam), 2, K, rng)
    detected, reading = cold.read_map(board, tape)
    assert detected == [seam]
    clock = {}
    numerator = {}
    release_counts = {}
    accelerations = {}
    exact_clock = {}
    for name, (hb, h) in hamiltonians.items():
        mask = [] if name == "blank" else [seam]
        assert (
            max(np.asarray(abs(h).sum(axis=1)).ravel())
            + max(np.asarray(abs(hb).sum(axis=1)).ravel())
            <= 80
        )
        numer, banks = release_bank(h, i, j, shells)
        numerator[name] = numer
        pairs = (
            (
                ("blank", weights["blank"]),
                ("imitated", np.array(source_counts["written"]) / 2**40),
            )
            if name == "blank"
            else (("direct", weights["blank"]), ("thermal", weights["written"]))
        )
        for tag, p in pairs:
            probs = banks @ p
            release_counts[tag] = [
                rng.multinomial(SHOTS, v / v.sum()).tolist() for v in probs
            ]
            accelerations[tag] = (numer @ p / (2 * (board.n - 1))).tolist()
        exact_clock[name] = [
            degree_clock(board, mask, weights[name], x) for x in range(board.n)
        ]
        clock[name] = {}
        for x in selections:
            data = ramsey_moments(h, hb, i, j, x)
            p = weights[name]
            measured = clock_samples(data, p, x, rng)
            conditioned = p.copy()
            conditioned[x] = 0
            conditioned /= conditioned.sum()
            np.testing.assert_allclose(
                -(data[1] @ conditioned).imag, exact_clock[name][x], atol=1e-10
            )
            clock[name][str(x)] = measured
            if name == "blank" and x == selections[0]:
                for tag, u in (
                    ("near_standard", min(board.neighbours[x])),
                    ("far_standard", selections[-1]),
                ):
                    point = np.zeros(board.n)
                    point[u] = 1
                    clock[tag] = clock_samples(data, point, x, rng)
    # An exact common kernel for every exterior cut, all possible source
    # address distributions and every normalized mixture of recorded releases.
    np.testing.assert_array_equal(numerator["blank"][1:], numerator["written"][1:])
    distance = np.array(shells)
    clock_change = np.asarray(exact_clock["written"]) - exact_clock["blank"]
    native_shell = np.bincount(distance, weights=-clock_change)
    beta = cold.apparatus.BETA
    ratios = diag["written"] / diag["blank"] - 1
    locality = []
    # Shift both bath generators by 3I. Their absolute row sums after the
    # shift are bounded by 11 (blank: 3), uniformly in federation size.
    for r in sorted(set(shells)):
        order = 2 * int(r) + 1
        tail = 2 * exp(11 * beta) * (11 * beta) ** order / factorial(order)
        locality.append(
            dict(
                hops=int(r),
                sites=int(sum(shells == r)),
                max_weight_ratio_change=float(np.max(abs(ratios[shells == r]))),
                uniform_bound=tail,
            )
        )
    return dict(
        protocol_sha256=sha256((HERE / "MEDIATED_RELEASE.md").read_bytes()).hexdigest(),
        seed=914683,
        ports=board.n,
        seam=seam,
        beta=beta,
        shots_per_release=SHOTS,
        times=list(RELEASE_TIMES),
        spectral_tape=tape,
        modular_reading=reading,
        source_counts=source_counts,
        stations=[
            dict(
                address=x,
                hops=int(shells[x]),
                blank_neighbors=list(map(int, board.neighbours[x])),
            )
            for x in selections
        ],
        clock=clock,
        clock_step=CLOCK_STEP,
        release_counts=release_counts,
        exact_cumulative_accelerations=accelerations,
        ordinary_integer_kernel=numerator["blank"].tolist(),
        max_exterior_integer_kernel_difference=int(
            np.abs(numerator["blank"][1:] - numerator["written"][1:]).max()
        ),
        source_probabilities={k: v.tolist() for k, v in weights.items()},
        gibbs_weight_range=locality,
        shell_clock_slowdown=native_shell.tolist(),
        spatial_clock_change=clock_change.tolist(),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    g = cli.add_mutually_exclusive_group()
    g.add_argument("--write", action="store_true")
    g.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "mediated_release.json"
    if args.counts_only:
        receipt = json.loads(path.read_bytes())
        cold.apparatus.check(receipt["decoded"], decode(receipt))
    else:
        receipt = run()
        receipt["decoded"] = decode(receipt)
        if args.write:
            path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        else:
            cold.apparatus.check(
                json.loads(path.read_bytes()), json.loads(json.dumps(receipt))
            )
    print(
        receipt["decoded"]["ordinary_neighbor_gain"],
        receipt["max_exterior_integer_kernel_difference"],
    )
