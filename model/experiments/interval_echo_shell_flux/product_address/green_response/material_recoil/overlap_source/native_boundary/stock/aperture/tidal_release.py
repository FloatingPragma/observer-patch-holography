"""Two-probe coincidence and chord readings of the closed monitored material."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
import json

import numpy as np

from . import material_release as flight
from . import relay_range as relay

lab = flight.lab
HERE = flight.HERE


def pair_paths(A, B, layers):
    paths = {(A, B): 1}
    for p in layers:
        new = defaultdict(int)
        for (i, j), count in paths.items():
            u, v = int(p[i]), int(p[j])
            if u == j:
                new[i, j] += 2 * count
                new[u, v] += 2 * count
            else:
                for x in (i, u):
                    for y in (j, v):
                        new[x, y] += count
        paths = new
    assert sum(paths.values()) == 4 ** len(layers) and all(i != j for i, j in paths)
    return paths


def category_weights(A, B, layers):
    n = len(layers[0])
    weights = defaultdict(lambda: np.zeros((3, n), dtype=int))
    for (i, j), count in pair_paths(A, B, layers).items():

        def cat(k):
            return 0 if k == i else 1 if k == j else 2

        row = weights[cat(A), cat(B)]
        row[2] += count
        for idx, k in enumerate((i, j)):
            row[2, k] -= count
            row[idx, k] += count
    groups = []
    for (u, v), mat in weights.items():
        total = int(mat[:, 0].sum())
        assert np.all(mat.sum(0) == total)
        mat[:, [A, B]] = 0
        groups.append((u, v, total, mat))
    return groups


@lru_cache(maxsize=None)
def draws(pop, count):
    if count == 0:
        return [((), Q(1))]
    ans = []
    for t, n in enumerate(pop):
        if n:
            rest = list(pop)
            rest[t] -= 1
            for tail, p in draws(tuple(rest), count - 1):
                ans.append(((t, *tail), Q(n, sum(pop)) * p))
    return ans


@lru_cache(maxsize=None)
def current(m, a, b, t):
    return sum(
        p * (int(out % 3 == t) - int(b == t))
        for out, p in flight.native(9 * m + 3 * a + b)
    )


@lru_cache(maxsize=None)
def coincidence(pop, row, t, cats):
    m, ab = divmod(row, 9)
    a, b = divmod(ab, 3)
    val = Q(0)
    for loads, p in draws(pop, cats.count(2)):
        tail = iter(loads)
        x = [next(tail) if c == 2 else (a, b)[c] for c in cats]
        if len(x) == 2 or x[2] == t:
            val += p * current(m, x[0], x[1], t)
    return val


def currents(dist, t, ends, layers, groups, launch):
    """Return J and E[J Y_k] for each service, with Y_k a second sharp flag."""
    A, B = ends
    n = len(layers[0])
    j1, c1, j2 = Q(0), Q(0), Q(0)
    c2 = np.full(n, Q(0), object)
    resident = defaultdict(Q)
    for (pop, labels), weight in dist.items():
        for m in labels[-1]:
            prob = weight / len(labels[-1])
            j = prob * current(m, *launch, t)
            j1 += j
            c1 += j * Q(pop[t], n)
            for row, w in flight.native(9 * m + 3 * launch[0] + launch[1]):
                resident[pop, row] += w * prob
    den = 4 ** len(layers)
    for u, v, total, mat in groups:
        j2 += Q(total, den) * sum(
            p * coincidence(pop, row, t, (u, v)) for (pop, row), p in resident.items()
        )
        for k in range(3):
            if np.any(mat[k]):
                val = sum(
                    p * coincidence(pop, row, t, (u, v, k))
                    for (pop, row), p in resident.items()
                )
                c2 += (val / den) * mat[k]
    cs = np.tile(c1, n)
    cs[[A, B]] = 0
    return np.array([j1, j2], object), np.array([cs, c2], object)


def pair_back(w, layers):
    """Backward expected distance for two walkers that share each swap coin."""
    for p in reversed(layers):
        together = p != np.arange(len(p))
        i = np.flatnonzero(together)
        j = p[i]
        new = w + w[p, :] + w[:, p] + w[p[:, None], p[None, :]]
        new[i, j] = 2 * (w[i, j] + w[j, i])
        np.fill_diagonal(new, 0)
        w = new
    return w


@lru_cache(maxsize=1)
def fastest_relay():
    """Optimal unambiguous masked bit among all two-slot load codes."""
    from itertools import product, combinations

    words = list(product(range(3), repeat=2))
    full = lab.gates()[0]
    table = np.zeros((2, 99), int)
    for col in range(99):
        m, ab = divmod(col, 9)
        a, b = divmod(ab, 3)
        for direction in range(2):
            inputs = (
                (3 * np.arange(3) + b) if direction == 0 else (3 * a + np.arange(3))
            )
            marg = (
                full[:, 9 * m + inputs].reshape(11, 3, 3, 3).sum(0).sum(axis=direction)
            )
            outcomes = [np.outer(marg[:, i], marg[:, j]) for i, j in words]
            table[direction, col] = max(
                x[y == 0].sum() + y[x == 0].sum() for x, y in combinations(outcomes, 2)
            )
    limited = relay.success_tickets()
    assert np.all(table >= 2 * limited[0]) and np.all(table >= 2 * limited[1])
    return table


def chart(h, ends, cost):
    A, B = ends
    via = np.minimum(h[:, A, None] + h[B, None, :], h[:, B, None] + h[A, None, :])
    use = h - via > cost.numerator // cost.denominator
    return np.where(use, via, h), use.astype(int)


def pair_read(w, ends, layers, C):
    A, B = ends
    den = 4 ** len(layers)
    first = pair_back(w, layers)
    second = pair_back(w, layers[::-1])
    both = pair_back(second, layers)

    def delta(kernel, c):
        return (kernel[B] - kernel[A]) @ c

    initial = int(w[A, B])
    mid = (int(first[A, B]) + delta(first, C[0])) / den
    late = (int(both[A, B]) + delta(both, C[0])) / den**2 + delta(second, C[1]) / den
    return np.array([Q(initial), mid, late], object)


def observation(h, ends, layers, parts, cost):
    A, B = ends
    ja, ca = parts
    C = ca[3] - ca[1] - ca[2] + ca[0]
    # Launch at B minus launch at A measures the signed projected chord.
    J = ja[2] - ja[1]
    basis = np.full(len(h), Q(0), object)
    basis[A], basis[B] = -1, 1
    mid = flight.forward(basis * (1 + J[0]), layers)
    late = flight.forward(mid + basis * J[1], layers[::-1])
    base, use = chart(h, ends, cost)
    path = pair_read(base, ends, layers, C) + cost * pair_read(use, ends, layers, C)
    dA, dB = [base[i] + cost * use[i] for i in ends]
    axis = (dA**2 - dB**2) / (2 * path[0] ** 2)
    projected = np.array([basis @ axis, mid @ axis, late @ axis], object)
    assert projected[0] == 1 and sum(mid) == sum(late) == 0
    assert np.all(C[:, [A, B]] == 0)
    n = len(h)
    equilibrium = (int(base.sum()) + int(use.sum()) * cost) / (n * (n - 1))
    return dict(
        separation=list(map(str, path)),
        normalized_pair_second=str((path[2] - 2 * path[1] + path[0]) / path[0]),
        normalized_projected=list(map(str, projected)),
        projected_second=str(projected[2] - 2 * projected[1] + projected[0]),
        equilibrium_separation=str(equilibrium),
    )


def native_source(dist, t, ends, layers, groups):
    arms = [
        currents(dist, t, ends, layers, groups, x)
        for x in ((0, 0), (t, 0), (0, t), (t, t))
    ]
    return tuple(np.array(x, object) for x in zip(*arms, strict=True))


def conditional(dist, t):
    ok = {k: p for k, p in dist.items() if k[0][t] == 0}
    prob = sum(ok.values(), Q(0))
    return prob, {k: p / prob for k, p in ok.items()}


def minimax(matrix):
    """Certify the least relative discrepancy for a unit inward chord tide."""
    from itertools import combinations
    from scipy.optimize import linprog
    from sympy import Matrix, eye, zeros

    n = matrix.cols
    units = Matrix.ones(matrix.rows, 1)
    system = (
        (matrix.row_join(-units))
        .col_join((-matrix).row_join(-units))
        .col_join(-eye(n + 1))
    )
    rhs = (-units).col_join(units).col_join(zeros(n + 1, 1))
    objective = zeros(n + 1, 1)
    objective[-1] = 1
    approx = np.asarray(system, float)
    cap = np.asarray(rhs, float).ravel()
    fit = linprog(
        np.asarray(objective, float).ravel(),
        A_ub=approx,
        b_ub=cap,
        bounds=[(None, None)] * (n + 1),
        method="highs",
    )
    assert fit.success
    tight = np.flatnonzero(abs(approx @ fit.x - cap) < 1e-6).tolist()
    for indices in combinations(tight, n + 1):
        basis = system[list(indices), :]
        if basis.det() == 0:
            continue
        inv = basis.inv()
        primal = inv * rhs[list(indices), :]
        dual = -inv.T * objective
        if all(v <= 0 for v in system * primal - rhs) and all(v >= 0 for v in dual):
            assert primal[-1] == -(dual.T * rhs[list(indices), :])[0]
            return dict(
                error=str(primal[-1]),
                coefficients=list(map(str, primal[:-1])),
                active=list(indices),
                dual=list(map(str, dual)),
            )
    raise ValueError("No exact active primal and dual certificate")


def clock_test(panel, treatment, reference):
    """Allow two arbitrary positive time intervals in each preparation."""
    from sympy import Matrix, Rational

    rows = []
    for t in ("1", "2"):
        for ruler in ("single_relay", "double_relay"):
            increments = []
            for name in (treatment, reference):
                series = list(
                    map(Q, panel[t][name]["charts"][ruler]["normalized_projected"])
                )
                increments += [series[2] - series[1], series[1] - series[0]]
            a, b, c, d = increments
            rows.append([a, -b, -c, d, Q(1)])
    matrix = Matrix(
        [[Rational(v.numerator, v.denominator) for v in row] for row in rows]
    )
    rays = matrix.nullspace()
    assert len(rays) == 1 and rays[0][-1] != 0
    ray = rays[0] / rays[0][-1]
    assert matrix * ray == Matrix.zeros(4, 1)
    coeff = matrix[:, :4]
    stationary = (coeff[:, 0] + coeff[:, 1]).row_join(coeff[:, 2] + coeff[:, 3])
    inv = coeff.inv()
    bounds = [max(0, -ray[i]) / sum(abs(v) for v in inv.row(i)) for i in range(4)]
    return dict(
        minimum_fractional_inward_error_lower=str(max(bounds)),
        arbitrary_intervals=minimax(coeff),
        equal_intervals=minimax(stationary),
        interval_ratios=[str(ray[i] / ray[i + 1]) for i in (0, 2)],
        stationary_additive_clock_solution=bool(
            ray[0] == ray[1] and ray[2] == ray[3] and all(v > 0 for v in ray)
        ),
        rank=matrix.rank(),
        unit_inward_solution=list(map(str, ray)),
        positive_clocks_possible=all(v > 0 for v in ray),
        unit_outward_positive_clocks_possible=all(v < 0 for v in list(ray)[:-1]),
    )


def schedule(board, marks, h, preparations, reverse):
    ends = marks[-1]
    A, B = ends
    layers = flight.matching_layers(board.n, set(board.seams) - set(marks))
    if reverse:
        layers.reverse()
    groups = category_weights(A, B, layers)
    tickets = relay.success_tickets()
    data = {}
    # An autonomous clock on a separately prepared stationary copy counts
    # each changed shared word, whichever edge produced it.
    half = lab.gates()[1]
    numbers = np.array([len(board.seams) - len(marks), 1, 1, 1])
    clock_rng = np.random.default_rng(62760)
    source_clocks = {}
    for name, dist in preparations.items():
        p = lab.before(dist) @ (8 - half.diagonal()) / 8
        exact = numbers @ p
        bank = clock_rng.binomial(2**32, np.array(p, float))
        measured = numbers @ bank / 2**32
        width = float(
            np.linalg.norm(numbers)
            * np.sqrt(np.log(2 * len(preparations) / 0.01) / (2 * 2**32))
        )
        assert abs(measured - float(exact)) <= width
        source_clocks[name] = dict(
            expected_ticks=str(exact),
            regional_samples=bank.tolist(),
            shots_per_type=2**32,
            measured_ticks=float(measured),
            simultaneous_99_halfwidth=width,
        )
    data["stationary_clocks"] = source_clocks
    for clean in (False, True):
        panel = {}
        for t in (1, 2):
            native = {}
            source_parts = {}
            for name, dist in preparations.items():
                selection, arm = conditional(dist, t) if clean else (Q(1), dist)
                if not selection:
                    continue
                pre = lab.before(arm)
                parts = native_source(arm, t, ends, layers, groups)
                source_parts[name] = parts
                clocks = pre[3] @ tickets.reshape(4, 99).T / 16
                assert clocks[0] == clocks[1] and clocks[2] == clocks[3]
                fast = pre @ fastest_relay().T / 32
                assert (
                    np.all(fast[0] == 1)
                    and np.all(fast[1:3] == 0)
                    and fast[3, 0] == fast[3, 1] > 0
                )
                costs = {
                    "fastest_relay": 1 / fast[3, 0],
                    "blank": Q(int(h[A, B])),
                    "single_relay": 1 / clocks[0],
                    "double_relay": 1 / clocks[2],
                }
                readings = {
                    k: observation(h, ends, layers, parts, cost)
                    for k, cost in costs.items()
                }
                C = parts[1][3] - parts[1][1] - parts[1][2] + parts[1][0]
                native[name] = dict(
                    acceptance=str(selection),
                    chart_cost={k: str(v) for k, v in costs.items()},
                    coincident_current=[[str(x) for x in c] for c in C],
                    relative_currents=list(map(str, parts[0][2] - parts[0][1])),
                    charts=readings,
                )
                if name in ("source", "original") and clean:
                    native[name]["direct_clean"] = clean_replay(
                        h, ends, layers, arm, t, costs, readings, reverse
                    )
            for a, b in [("source", "original"), ("positive", "reference")]:
                native[a]["fixed_reference_chart"] = {
                    c: observation(h, ends, layers, source_parts[a], Q(cost))
                    for c, cost in native[b]["chart_cost"].items()
                }
            panel[str(t)] = native
        data["clean" if clean else "all"] = panel
    # Every pulse of the changed transport law keeps each material reading
    # on the third contact in place, regardless of its stochastic head.
    # This erases J and the coincident currents for every arm and launch.
    for col in range(99):
        assert all(row % 9 == col % 9 for row, _ in flight.native(col, True))
    empty_parts = (
        np.full((4, 2), Q(0), object),
        np.full((4, 2, board.n), Q(0), object),
    )
    data["hold_material_common_blank_chart"] = observation(
        h, ends, layers, empty_parts, Q(int(h[A, B]))
    )
    data["clock_freedom"] = {
        tag: clock_test(data["all"], a, b)
        for tag, a, b in [
            ("repair", "source", "original"),
            ("matched", "positive", "reference"),
        ]
    }
    common = {}
    for a, b in [
        ("source", "original"),
        ("positive", "reference"),
        ("ordinary", "original"),
        ("single", "original"),
        ("reference_plus_double", "reference"),
        ("reference_plus_single", "reference"),
    ]:
        rows = [
            [
                Q(data["all"][t][name]["charts"]["fastest_relay"]["projected_second"])
                for t in ("1", "2")
            ]
            for name in (a, b)
        ]
        q = (rows[1][0] - rows[1][1]) / (rows[0][0] - rows[0][1])
        gain = rows[1][0] - q * rows[0][0]
        observed = (
            Q(source_clocks[b]["expected_ticks"])
            / Q(source_clocks[a]["expected_ticks"])
        ) ** 2
        common[a] = dict(
            target_time_ratio_squared=str(q),
            inward_response_at_equal_projection=str(gain),
            measured_time_ratio_squared=str(observed),
            equal_sweep=list(map(str, [x - y for x, y in zip(*rows, strict=True)])),
            word_clock=list(
                map(str, [observed * x - y for x, y in zip(*rows, strict=True)])
            ),
        )
    data["common_ruler"] = common
    return data


def clean_replay(h, ends, layers, dist, t, costs, readings, reverse):
    # Independent oracle. Both loads are t, so the first marked appointment
    # leaves their positions fixed. Transport the two physical probes forward
    # as a joint distribution. At the second appointment read the partner
    # against the sharp background, retaining the head and composition.
    A, B = ends
    paths = pair_paths(A, B, layers)
    den = 4 ** len(layers)
    passed = {key: Q(count, den) for key, count in paths.items()}
    outgoing = defaultdict(Q)
    jump_counts = np.zeros(len(h), dtype=np.int64)
    rng = np.random.default_rng(729034 + 10 * t + int(reverse))
    local_bank = rng.multinomial(2**38, np.array(list(passed.values()), float))
    rates = [Q(0), Q(0)]
    for (pop, labels), p in dist.items():
        assert pop[t] == 0
        for m in labels[-1]:
            for row, w in flight.native(9 * m + 4 * t):
                for other, count in enumerate(pop):
                    if not count:
                        continue
                    for d, (a, b) in enumerate(((t, other), (other, t))):
                        prob = p * w * Q(count, len(h) * len(labels[-1]))
                        for code, s in flight.native(9 * (row // 9) + 3 * a + b):
                            if code % 9 != 3 * a + b:
                                rates[d] += prob * s
    for ((i, j), p), starts in zip(passed.items(), local_bank, strict=True):
        occupied = {i, j}
        if (A in occupied) == (B in occupied):
            outgoing[i, j] += p
            continue
        origin, dest, rate = (A, B, rates[0]) if A in occupied else (B, A, rates[1])
        partner = j if i == origin else i
        sign = 1 if origin == A else -1
        jump_counts[partner] += sign * rng.binomial(starts, float(rate))
        outgoing[dest if i == origin else i, dest if j == origin else j] += p * rate
        outgoing[i, j] += p * (1 - rate)
    # Exact forward joint aggregation. Integer transport separately from
    # the rational scale avoids multiplying full Fraction matrices.
    import math

    common = math.lcm(*(p.denominator for p in outgoing.values()))
    assert common * den < 2**63
    arr = np.zeros(h.shape, np.int64)
    for k, p in outgoing.items():
        arr[k] = int(p * common)
    # Forwards through reversed layers equals the backward action of
    # the forward order since every matching is self-adjoint.
    arr = pair_back(arr, layers)
    assert int(arr.sum()) == common * den
    index = np.flatnonzero(arr)
    probs = arr.ravel()[index].astype(float) / (common * den)
    counts = rng.multinomial(2**38, probs / probs.sum())
    result = {}
    for name, cost in costs.items():
        base, use = chart(h, ends, cost)
        d = base.ravel()[index] + cost * use.ravel()[index]
        expected = sum(
            Q(int(weight), common * den) * v
            for weight, v in zip(arr.ravel()[index], d, strict=True)
        )
        assert expected == Q(readings[name]["separation"][2])
        measured = counts @ np.asarray(d, float) / counts.sum()
        # A fresh bank counted signed jumps with the other trace at k.
        # Source-weighted prediction and independent arrivals share only
        # the native kernel and the frozen schedule and chart.
        q0, q1 = [pair_back(m, layers[::-1]) for m in (base, use)]
        delta = ((q0[B] - q0[A]) + float(cost) * (q1[B] - q1[A])) / den
        delta[[A, B]] = 0
        base_free = sum(
            p * (int(q0[i, j]) + cost * int(q1[i, j])) / den
            for (i, j), p in passed.items()
        )
        local = float(base_free) + delta @ (jump_counts / 2**38)
        local_band = (
            2 * np.max(abs(delta)) * np.sqrt(np.log(2 * 64 / 0.01) / (2 * 2**38))
        )
        assert abs(local - float(expected)) <= local_band
        band = float(max(d) - min(d)) * np.sqrt(
            np.log(2 * 64 / 0.01) / (2 * counts.sum())
        )
        assert abs(measured - float(expected)) <= band
        result[name] = dict(
            mean=float(measured),
            simultaneous_99_halfwidth=band,
            local_prediction=float(local),
            local_99_halfwidth=float(local_band),
        )
    return dict(
        signed_jump_counts_by_partner=jump_counts.tolist(),
        supported_ordered_landings=len(index),
        draws=int(counts.sum()),
        charts=result,
    )


def run():
    _, words, cells, *_ = lab.energy.pair.setup()
    board, _, marks, *_ = lab.energy.returned(words)
    _, h, _ = flight.records.thermal_radar.chart(board, marks)
    u = np.full(len(words), Q(1, len(words)), object)
    arms = dict(
        original=flight.records.ensemble(words, u),
        source=flight.records.ensemble(
            words, lab.energy.least_change(cells) / len(words)
        ),
        ordinary=flight.records.ensemble(words, u, extra=1),
        single=flight.records.ensemble(words, u, single=1),
    )
    w = json.loads((HERE / "material_release.json").read_bytes())["matched_mixtures"]
    arms.update(positive=flight.mix(arms, w[0]), reference=flight.mix(arms, w[1]))
    # Repeat the ordinary intervention on the actual reference used in
    # the fixed-history/volume contrast, including its pre-existing donor.
    assert set(w[1]) == {"source", "ordinary"}
    for name, extra, single in [
        ("reference_plus_double", 1, 0),
        ("reference_plus_single", 0, 1),
    ]:
        augmented = dict(
            source=flight.records.ensemble(
                words,
                lab.energy.least_change(cells) / len(words),
                extra=extra,
                single=single,
            ),
            ordinary=flight.records.ensemble(words, u, extra=extra + 1, single=single),
        )
        arms[name] = flight.mix(augmented, w[1])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "TIDAL_RELEASE.md",
                HERE / "material_release.json",
                HERE / "relay_range.json",
            )
        },
        ends=marks[-1],
        center=219,
        endpoint_radii=h[219, list(marks[-1])].tolist(),
        schedules={
            tag: schedule(board, marks, h, arms, rev)
            for tag, rev in [("forward", False), ("reverse", True)]
        },
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    raw = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = HERE / "tidal_release.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(
        "Replayed joint arrivals and opposite-launch projections on both service orders."
    )
