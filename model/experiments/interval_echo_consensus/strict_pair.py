"""Free release of a two-load discrepancy after the Interval Echo strict repair.

Run with --full to propagate both instruments at every comparison appointment.
Add --cross-check for direct sign-basis propagation of the full two-load sector.
"""

import argparse
from fractions import Fraction
from itertools import combinations
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law as law
from model.executable_model import federation
from model.executable_model.carrier import Carrier
from .__main__ import local_ops

START_PORT = 156
TIMES = (0, 1, 7, 37)
STRENGTHS = (0, 0.5, KAPPA)


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def verify_columns():
    """Reconstruct the echo and both free generators exactly on the native basis."""
    w = sp.Matrix(law._write_pattern()) / sp.sqrt(2)
    u = sp.zeros(9).row_join(w.T).col_join(w.row_join(sp.eye(90) - w * w.T))
    blank = sp.diag(*([1] * 9 + [0] * 90))
    d = [None] * 9 + [d for _, d in law.LABELS for _ in law.SRC]
    neutral = sp.diag(*[int(x == 0) for x in d])
    strict_output = sp.diag(
        *[int(abs(x or 0) == 2 and law.SRC[i % 9] == (1, 1)) for i, x in enumerate(d)]
    )
    pair = sp.diag(*[int(law.SRC[i % 9] == (1, 1)) for i in range(99)])
    contrast_one = sp.diag(*[int(abs(x or 0) == 1) for x in d])
    hit = strict_output * u * blank
    quiet = u.T * (u * blank - hit)
    occupied = sp.eye(99) - blank
    ready = u.T * hit
    r, dd = map(sp.Matrix, law._label_ops())
    native_echo, _, commit, probes = local_ops()
    check(np.allclose(np.array(u, float), native_echo, atol=1e-14, rtol=0), "echo mismatch")
    check(np.allclose(np.array(hit, float), commit, atol=1e-14, rtol=0), "commit mismatch")
    for exact, (_, actual) in zip([occupied, quiet, ready], probes, strict=True):
        check(np.allclose(np.array(exact, float), actual, atol=1e-14, rtol=0), "probe mismatch")
    generators = []
    for coef, actual in zip([(sp.eye(10) - r) / 2, dd], law.local_term_parts(), strict=True):
        h = u * sp.diag(sp.zeros(9), sp.kronecker_product(coef, sp.eye(9))) * u.T
        generators.append(h)
        check(np.allclose(np.array(h, float), actual, atol=1e-14, rtol=0), "free-law mismatch")
        check(neutral * h == sp.zeros(99), "free term can feed neutral records")
    for op in [hit, occupied, quiet, ready]:
        check(
            neutral * op * (sp.eye(99) - neutral) == sp.zeros(99),
            "strict operation creates neutral",
        )
    check(
        hit.T * hit + quiet.T * quiet + occupied.T * occupied == sp.eye(99),
        "incomplete instrument",
    )
    check(
        (neutral * u * sp.eye(99)[:, 0]).norm() ** 2 == 1,
        "raw empty echo must create a neutral record",
    )
    for word in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        x = sp.eye(99)[:, law.SRC.index(word)]
        check(hit * x == sp.zeros(99, 1) and quiet * x == x, "strict action on binary blank input")
        check(u.T * contrast_one * u * x == int(word[0] != word[1]) * x, "binary XOR reader")
        check(u.T * pair * u * x == int(word == (1, 1)) * x, "blank pair flag")
        swapped = sp.eye(99)[:, law.SRC.index(word[::-1])]
        check(generators[0] * x == (x - swapped) / 2, "blank exchange column")
        check(generators[1] * x == (word[0] - word[1]) ** 2 * x, "blank contrast column")
    captured = {}
    for word in [(2, 0), (0, 2)]:
        y = hit * sp.eye(99)[:, law.SRC.index(word)]
        check(
            y.norm() ** 2 == 1 and (pair * y).norm() ** 2 == 1,
            "strict capture did not produce two units",
        )
        captured[word] = y
        check((pair * u * y).norm() == 0, "written-seam flag at release")
    parities = []
    for odd in (0, 1):
        memory = sp.zeros(11, 1)
        for coin in (0, 1):
            for sign in (-1, 1):
                memory[1 + law.LABELS.index((coin, 2 * sign))] = sp.Rational(sign if odd else 1, 2)
        parities.append(sp.kronecker_product(memory, sp.eye(9)[:, law.SRC.index((1, 1))]))
        for word in ((0, 0), (0, 1), (1, 0), (1, 1)):
            col = sp.kronecker_product(memory, sp.eye(9)[:, law.SRC.index(word)])
            for h, coefficient in zip(generators, (odd, 4), strict=True):
                check(h * col == int(word != (1, 1)) * coefficient * col, "exact signed block")
            check(pair * u * col == sp.zeros(99, 1), "native pair flag on strict memory")
    for word, y in captured.items():
        sign = 1 if word == (2, 0) else -1
        check(y == (parities[0] + sign * parities[1]) / sp.sqrt(2), "signed launch mismatch")


def strict_scheduler(board):
    """Resolve every sharp double without passing a destination to the scheduler."""
    from model.experiments.interval_echo_consensus.asynchronous import SharpService, settle

    service = SharpService(range(board.n), board.seams)
    selected = []
    for port in range(board.n):
        initial = [0] * board.n, [0] * len(board.seams)
        initial[0][port] = 2
        final = settle(service, initial, "reference", 0)
        for order in ("forward", "reverse", "shuffled"):
            check(settle(service, initial, order, 31) == final, "repair address depends on visits")
        loads, collars = final
        marks = [i for i, x in enumerate(collars) if x]
        check(len(marks) == 1, "repair wrote the wrong number of memories")
        edge = board.seams[marks[0]]
        check(
            collars[marks[0]] == (2 if port == edge[0] else -2)
            and port in edge
            and {i for i, x in enumerate(loads) if x} == set(edge)
            and sum(loads) == 2
            and all(loads[i] == 1 for i in edge),
            "wrong repaired pair",
        )
        selected.append([*edge, collars[marks[0]]])
    internal = sum(a // 12 == b // 12 for a, b, _ in selected)
    check(
        len(selected) == 240 and internal == 180 and tuple(selected[156][:2]) == (156, 157),
        "unexpected scheduler output panel",
    )
    return selected


def return_floor(moves, strength, odd=0):
    """Exact bound for a positive separation; zero when the estimate is vacuous."""
    if moves <= 0 or odd not in (0, 1) or strength < 0:
        raise ValueError("expected a positive exit count, sign parity and nonnegative coupling")
    v2 = Fraction(moves, 4)
    gap = 3 * Fraction(strength) + odd - Fraction(moves, 2)
    if gap <= 0 or gap**2 <= v2:
        return Fraction(0)
    return ((gap**2 - v2) / (gap**2 + v2)) ** 2


def reader_census(board):
    """Enumerate the six contrast flags on each possible binary pair, for every seam."""
    left, right = np.triu_indices(board.n, 1)
    occupied = [(left == p) | (right == p) for p in range(board.n)]
    counts = {"internal": {}, "intercell": {}}
    for a, b in board.seams:
        hit = np.ones(left.size, dtype=bool)
        for port, partner in ((a, b), (b, a)):
            probes = [q for q in board.neighbours[port] if q != partner and q // 12 == port // 12][
                :3
            ]
            check(len(set(probes)) == 3, "endpoint needs three distinct blank probe seams")
            for q in probes:
                hit &= occupied[port] ^ occupied[q]
        check(np.array_equal(hit, (left == a) & (right == b)), "aliased six-probe reading")
        moves = int(board.degree[a] + board.degree[b] - 2)
        table = counts["internal" if a // 12 == b // 12 else "intercell"]
        table[moves] = table.get(moves, 0) + 1
    return counts


def flight_parts(board, tag):
    live = Carrier("blank seams", board.n, [e for e in board.seams if e != tag])
    sector = PairSector(live)
    exchange = sector.exchange()
    # Keep the degree references of the intact board when separating the written seam.
    weight = live.degree + 6 - board.degree
    cost = weight[sector.aa] + weight[sector.bb] - 2 * sector.adjacent
    start = sector.index(*tag)
    initial = np.zeros(sector.P)
    initial[start] = 1
    moves = int(sum(board.degree[p] for p in tag) - 2)
    off_diagonal = exchange[:, start].toarray().ravel()[initial == 0]
    check(
        exchange[start, start] == moves / 2
        and cost[start] == 10
        and cost[initial == 0].min() >= 9
        and off_diagonal @ off_diagonal == moves / 4,
        "launch energy, complementary cost or off-diagonal norm",
    )
    return sector, exchange, cost, initial


def complete_bins(visible):
    bins = np.r_[visible, 1 - visible.sum()]
    check(
        np.isfinite(bins).all() and bins.min() >= -1e-11 and bins.max() <= 1 + 1e-11,
        "invalid detector distribution",
    )
    return bins


def strict_flight(board, tag, home, strengths=STRENGTHS, times=TIMES):
    sector, exchange, cost, initial = flight_parts(board, tag)
    tag_slot = home.index(tag)
    detector = [sector.index(*e) for e in home]
    rows = []
    for strength in strengths:
        occupied = np.zeros((len(times), len(home)))
        for odd in (0, 1):
            potential = strength * cost + (4 * strength + odd) * (1 - initial)
            h = exchange + sparse.diags(potential - potential[initial == 1].item(), format="csr")
            state = initial
            previous = 0
            for i, t in enumerate(times):
                if t != previous:
                    state = expm_multiply(-1j * (t - previous) * h, state)
                previous = t
                check(abs(np.vdot(state, state).real - 1) < 1e-10, "lost flight norm")
                occupied[i] += abs(state[detector]) ** 2 / 2
        floor = return_floor(int(sum(board.degree[p] for p in tag) - 2), strength)
        for t, visible in zip(times, occupied, strict=True):
            check(visible[tag_slot] >= float(floor) - 1e-10, "violated all-time bound")
            original = visible.copy()
            original[tag_slot] = 0
            rows.append(
                {
                    "kappa": strength,
                    "t": t,
                    "pair_flag": complete_bins(original).tolist(),
                    "six_probe": complete_bins(visible).tolist(),
                }
            )
    return rows


def cross_check(board, tag, home, rows):
    """Construct the joint signed-memory Hamiltonian from individual allowed load hops."""
    from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission.neutral_pair import (
        chebyshev_flight,
    )

    words = list(combinations(range(board.n), 2))
    indices = {w: i for i, w in enumerate(words)}
    n = len(words)
    edges = set(board.seams)
    src, dst = [], []
    cost = []
    for i, word in enumerate(words):
        cost.append(12 - 2 * (word in edges) - int(len(set(word) & set(tag)) == 1))
        a, b = word
        for p, other in ((a, b), (b, a)):
            for q in board.neighbours[p]:
                if q != other and tuple(sorted((p, q))) != tag:
                    src.append(i)
                    dst.append(indices[tuple(sorted((q, other)))])
    hops = sparse.csr_matrix((np.full(len(src), -0.5), (dst, src)), shape=(n, n))
    exchange = hops + sparse.diags(np.bincount(src, minlength=n) / 2)
    out = np.ones(n)
    out[indices[tag]] = 0
    flip = sparse.diags(-out / 2)
    initial = np.zeros(2 * n)
    initial[indices[tag]] = 1
    side = [indices[e] for e in home]
    worst = 0.0
    for k in STRENGTHS:
        diagonal = exchange + sparse.diags(k * (np.asarray(cost) + 4 * out) + out / 2)
        h = sparse.bmat([[diagonal, flip], [flip, diagonal]], format="csr")
        check((h != h.T).nnz == 0, "joint generator is not self-adjoint")
        off_radius = np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal())
        for t in TIMES:
            state = initial
            if t:
                state, _, tail = chebyshev_flight(h, off_radius, initial, t)
                check(tail < 1e-14, "Chebyshev series did not stabilize")
            check(abs(np.vdot(state, state).real - 1) < 1e-10, "joint flight norm")
            visible = (abs(state.reshape(2, n)[:, side]) ** 2).sum(axis=0)
            row = next(r for r in rows if r["kappa"] == k and r["t"] == t)
            for reader in ("six_probe", "pair_flag"):
                if reader == "pair_flag":
                    visible[home.index(tag)] = 0
                worst = max(worst, float(np.max(abs(complete_bins(visible) - row[reader]))))
    check(worst < 2e-10, "joint propagation disagrees with parity sectors")
    return {
        "dimension": 2 * n,
        "bins_checked": 2 * (len(home) + 1) * len(TIMES) * len(STRENGTHS),
        "maximum_absolute_difference": worst,
    }


def run(full=False, compare=False):
    board = federation(0)
    verify_columns()
    launches = strict_scheduler(board)
    a, b, sign = launches[START_PORT]
    tag = (a, b)
    home = [e for e in board.seams if e[0] // 12 == e[1] // 12 == a // 12]
    check(len(home) == 30 and tag in home and sign == 2, "unexpected launch and home")
    census = reader_census(board)
    check(
        census == {"internal": {8: 308, 9: 284, 10: 8}, "intercell": {10: 30}},
        "unexpected seam census",
    )
    report = {
        "source_port": START_PORT,
        "launches_by_source_port": launches,
        "seams_by_exit_count": census,
        "home_order": [*home, "quiet"],
        "uniform_floor_at_eight": str(return_floor(10, KAPPA)),
        "uniform_floor_at_six": str(return_floor(10, 6)),
    }
    if full or compare:
        report["outcomes"] = strict_flight(board, tag, home)
    if compare:
        report["independent_replay"] = cross_check(board, tag, home, report["outcomes"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--cross-check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.full, args.cross_check), allow_nan=False))
