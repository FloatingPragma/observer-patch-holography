"""Native arrival before an ordinary return as a receiver of a source fluctuation."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

import networkx as nx
import numpy as np
from scipy.linalg import solve_banded

from . import closed_energy as energy
from .measure import check

HERE = energy.HERE
stock = energy.pair.stock
SHOTS = 2**16
LIMIT = 4096


def ring(level):
    board = stock.observations.federation(level)
    graph = nx.Graph()
    graph.add_edges_from(board.seams)
    cycles = sorted(
        (tuple(c) for c in nx.cycle_basis(graph, root=0) if len(c) % 2 == 0),
        key=lambda c: (-len(c), c),
    )
    for cycle in cycles:
        for route in (cycle, cycle[1:] + cycle[:1]):
            even = set(route[::2])
            sides = nx.Graph()
            sides.add_nodes_from(sorted(even))
            sides.add_edges_from(
                (x, y) for x in sorted(even) for y in sorted(set(graph[x]) - set(route))
            )
            matching = nx.bipartite.maximum_matching(sides, top_nodes=sorted(even))
            if all(x in matching for x in even):
                leaves = [int(matching[x]) for x in route[::2]]
                assert len(set(leaves)) == len(even) and set(leaves).isdisjoint(route)
                return list(map(int, route)), leaves
    raise ValueError("No even cycle with disjoint material sidings")


def response(half_length, hazard):
    """Exact homogeneous answer from finite binomial polynomials."""
    p = F(hazard)
    even = sum(
        F(comb(half_length, k)) * p ** (k // 2) for k in range(0, half_length + 1, 2)
    )
    odd = sum(
        F(comb(half_length, k)) * p ** ((k + 1) // 2)
        for k in range(1, half_length + 1, 2)
    )
    return odd / even


def exact_first_step(n, h):
    """Rational Dirichlet solve for surviving to the first return on the even sites."""
    q = 1 - h
    # U_0=q on arrival. Off the starting site U_i=q*(U_{i-1}+2U_i+U_{i+1})/4.
    off = -q / 4
    diag = [1 - q / 2] * (n - 1)
    rhs = [F(0)] * (n - 1)
    rhs[0] = rhs[-1] = q * q / 4
    for i in range(1, n - 1):
        ratio = off / diag[i - 1]
        diag[i] -= ratio * off
        rhs[i] -= ratio * rhs[i - 1]
    sol = [F(0)] * (n - 1)
    for i in reversed(range(n - 1)):
        sol[i] = (rhs[i] - (off * sol[i + 1] if i + 1 < n - 1 else 0)) / diag[i]
    cycle = [q, *sol]
    assert all(
        4 * cycle[i] == q * (cycle[i - 1] + 2 * cycle[i] + cycle[(i + 1) % n])
        for i in range(1, n)
    )
    return 1 - (sol[0] + 2 * q + sol[-1]) / 4


def local_response(h):
    n = len(h)
    probs = []
    for start in range(n):
        a = np.roll(h, -start)
        miss = (1 - a[1:]) / 4
        matrix = np.zeros((3, n - 1))
        matrix[1] = 1 - 2 * miss
        matrix[0, 1:] = -miss[:-1]
        matrix[2, :-1] = -miss[1:]
        rhs = a[1:].copy()
        rhs[[0, -1]] += miss[[0, -1]] * a[0]
        rest = solve_banded((1, 1), matrix, rhs)
        probs.append((rest[0] + 2 * a[0] + rest[-1]) / 4)
    return np.mean(probs)


def local_cap(n, k):
    if not k:
        return F(0)
    gap, rem = divmod(n, k)
    harm = sum((F(1, i) for i in range(1, gap)), start=F())
    return F(k, n) + (k * harm + F(rem, gap)) / (2 * n)


def native():
    W, V, _, h0 = stock.kernels()
    full = (2 * np.eye(99, dtype=int) - h0) ** 2
    tickets = np.array(
        [np.repeat(np.arange(99), col) for col in full.T], dtype=np.int16
    )
    write = np.array(
        [np.repeat(np.arange(99), col) for col in W[:, :9].T], dtype=np.int16
    )
    d = np.array([0, *[x[1] for x in stock.LAW.LABELS]])
    for j, (a, b) in enumerate(stock.LAW.SRC):
        np.testing.assert_array_equal(d[write[j] // 9], a - b)
    # One blank half exchange on (1,0) gives the fair direction; a blank
    # full contact moves the walker exactly one chosen physical edge.
    np.testing.assert_array_equal(np.flatnonzero(V[:, 3]), [1, 3])
    np.testing.assert_array_equal(V[[1, 3], 3], 32)
    np.testing.assert_array_equal(tickets[3], 1)
    # Probability of crossing conditional on every *current* sharp head,
    # with a fair direction. Magnitude and written status survive either port read.
    chances = {}
    for load in (1, 2):
        moving = []
        for head in range(11):
            hits = 0
            for old in (load, 3 * load):
                out = tickets[9 * head + old]
                np.testing.assert_array_equal(abs(d[out // 9]), abs(d[head]))
                np.testing.assert_array_equal(out // 9 == 0, head == 0)
                assert set(out % 9) <= {old, 4 * load - old}
                hits += int(sum(out % 9 != old))
            moving.append(F(hits, 8))
        chances[load] = moving
    assert chances[1] == [
        F(1) if i == 0 else F(int(abs(v) == 1), 4) for i, v in enumerate(d)
    ]
    assert chances[2] == [F(1)] + [F(0)] * 10
    return tickets, write, d, chances


def moving_response(sites):
    """First local transit of a forward courier versus an independent return."""
    n = len(sites)
    if not np.any(sites):
        return 0.0
    survival = np.ones(n)
    first = []
    for t in range(1, n + 1):
        h = np.roll(sites, -t) / 4
        first.append(float(np.mean(survival * h)))
        survival *= 1 - h
    cycle = (3 / 4) ** int(sites.sum())
    np.testing.assert_allclose(survival, cycle)
    eigen = np.cos(np.pi * np.arange(1, n, 2) / (2 * n)) ** 2
    # Pr(T>=j)=(2/n)sum eigen**(j-1) for j>=2, and Pr(T>=1)=1.
    polys = np.polynomial.polynomial.polyval(eigen, first)
    return first[0] * (1 - 2 * len(eigen) / n) + (2 / n) * sum(
        polys / (1 - cycle * eigen**n)
    )


def trial_bank(n, active, mode, seed, load=1, no_write=False, alignment="independent"):
    tickets, write, contrast, _ = native()
    rng = np.random.default_rng(seed)
    # Independent source preparations use actual write outcomes; all subsequent
    # sharp pulses on a resident contact update that same trial's stored head.
    size = len(active)
    assert size == n
    side = rng.integers(2, size=(SHOTS, n), dtype=np.int8)
    if alignment == "common":
        side[:] = side[:, 0, None]
    if alignment == "opposite":
        ids = np.flatnonzero(active)
        assert len(ids) % 2 == 0
        side[:, ids[1::2]] = 1 - side[:, ids[::2]]
    word = active[None, :] * (1 + 2 * side)
    heads = (write[word, rng.integers(64, size=word.shape, dtype=np.int8)] // 9).astype(
        np.int8
    )
    if no_write:
        heads[:] = 0
    sums = contrast[heads].sum(1)
    # The full blank swap removes either sharp material output to a
    # private apparatus rail and restores it upon reversal.
    swap = tickets[:9, 0]
    assert np.array_equal(swap[3 * np.arange(3)], np.arange(3))
    assert np.array_equal(swap[np.arange(3)], 3 * np.arange(3))
    start = rng.integers(n, size=SHOTS)
    at = start.copy()
    ids = np.arange(SHOTS)
    passed = returned = total = 0
    for tick in range(1, LIMIT + 1):
        # In units of two steps the same-site bin has two of four native paths.
        hop = rng.integers(4, size=len(ids), dtype=np.int8)
        at = (at + (hop == 0).astype(np.int8) - (hop == 3)) % n
        if mode == "pooled":
            pos = rng.integers(n, size=len(ids))
        elif mode == "moving":
            pos = (start[ids] + tick) % n
        else:
            assert mode == "local"
            pos = at
        old = load * (1 + 2 * rng.integers(2, size=len(ids)))
        out = tickets[
            9 * heads[ids, pos].astype(int) + old, rng.integers(4, size=len(ids))
        ]
        heads[ids, pos] = out // 9
        yes = out % 9 != old
        home = at == start[ids]
        passed += int(yes.sum())
        returned += int((home & ~yes).sum())
        keep = ~(yes | home)
        total += len(ids)
        ids, at = ids[keep], at[keep]
        if not len(ids):
            break
    assert passed + returned + len(ids) == SHOTS
    delta = np.sqrt(np.log(2 * 500 / 0.01) / (2 * SHOTS))
    return dict(
        trials=SHOTS,
        transits=passed,
        returns=returned,
        pending=len(ids),
        appointments=total,
        probability_band=[
            max(0, passed / SHOTS - delta),
            min(1, (passed + len(ids)) / SHOTS + delta),
        ],
        source_mean=float(sums.mean()),
        source_second_moment=float(np.mean(sums**2)),
    )


def unrecovered(w, weights, keys):
    groups = {}
    y = w[:, 1] != 1
    for word, v, z in zip(w, y, weights, strict=True):
        key = tuple(word[keys])
        mass, hit = groups.get(key, (F(0), F(0)))
        groups[key] = (mass + z, hit + z * bool(v))
    return sum(a * (1 - a / m) for m, a in groups.values()) / sum(weights)


def middle_meter(w, weights, cells, n, seed, source_interval):
    tickets, write, contrast, _ = native()
    # Port 20 to the smaller numbered blank neighbor, on the actual board.
    chart = json.loads((HERE / "closed_energy.json").read_bytes())
    home = chart["ordinary_home_ports"]
    contact = min(chart["ordinary_edges"])
    assert home[1] in contact
    left = contact[0] == home[1]
    # Orient the measured head so that a positive reading means n>1.
    load = w[:, 1]
    inputs = 3 * load + 1 if left else 3 + load
    possible = write[inputs] // 9
    np.testing.assert_array_equal(
        contrast[possible] * (1 if left else -1),
        np.broadcast_to((load - 1)[:, None], possible.shape),
    )
    d = w[:, :4] - 1
    np.testing.assert_array_equal(d.sum(1), 0)
    assert all(v == 0 for v in weights @ d)
    covariance = d.T @ (weights[:, None] * d)
    assert all(v == 0 for v in covariance.sum(1))
    f = sum(weights[cells == 0])
    assert covariance[1, 1] == f
    out = {}
    for attempts in (1, 2):
        alpha = 1 - F(3, 4) ** attempts
        predicted = response(n, alpha * f)
        # A real source word is sampled at each appointment. The second full
        # pulse keeps the actual first head after its failed sharp port read.
        rng = np.random.default_rng(seed + attempts)
        loc = np.zeros(2**18, int)
        trials = len(loc)
        hits = returned = 0
        for _ in range(LIMIT):
            hop = rng.integers(4, size=len(loc))
            loc = (loc + (hop == 0).astype(int) - (hop == 3)) % n
            drawn = rng.choice(len(w), size=len(loc), p=np.asarray(weights, float))
            mark = write[inputs[drawn], rng.integers(64, size=len(loc))] // 9
            yes = np.zeros(len(loc), bool)
            for _ in range(attempts):
                live = np.flatnonzero(~yes)
                arrival = 1 + 2 * rng.integers(2, size=len(live))
                dst = tickets[9 * mark[live] + arrival, rng.integers(4, size=len(live))]
                yes[live] = dst % 9 != arrival
                mark[live] = dst // 9
            hits += int(yes.sum())
            returned += int(((loc == 0) & ~yes).sum())
            loc = loc[(loc != 0) & ~yes]
            if not len(loc):
                break
        assert hits + returned + len(loc) == trials
        width = np.sqrt(np.log(2 * 500 / 0.01) / (2 * trials))
        band = [hits / trials - width, (hits + len(loc)) / trials + width]
        forecast = [float(response(n, alpha * F(float(x)))) for x in source_interval]
        assert band[0] <= float(predicted) <= band[1] and max(
            band[0], forecast[0]
        ) <= min(band[1], forecast[1])
        out[str(attempts)] = dict(
            ordinary_gain=str(alpha),
            pooled=float(predicted),
            one_reused=float(f * response(n, alpha)),
            trials=trials,
            transits=hits,
            returns=returned,
            pending=len(loc),
            interval=band,
            from_source=forecast,
        )
    return dict(
        contact=contact,
        all_source_words=len(w),
        signed_middle_variance=str(f),
        covariance=[[str(x) for x in row] for row in covariance],
        clocks=out,
    )


def source_experiment(n, seed):
    _, w, cells, *_ = energy.pair.setup()
    _, write, _, chances = native()
    z = energy.least_change(cells)
    rows = {}
    for arm, weights in (("original", np.full(len(w), F(1), object)), ("repaired", z)):
        weights = weights / sum(weights)
        third = w[:, 6]
        hazard = sum(
            weight * chances[1][x] for weight, x in zip(weights, third, strict=True)
        )
        persistent = sum(
            sum(weights[third == i]) * response(n, p) for i, p in enumerate(chances[1])
        )
        row = dict(
            original_third=dict(
                pooled_hazard=str(hazard),
                pooled=float(response(n, hazard)),
                one_reused=float(persistent),
            )
        )
        for name, keys in (("unresolved", []), ("third", [6]), ("middle", [1])):
            var = unrecovered(w, weights, keys)
            # Two independent records conditional on the same witness disagree
            # with probability 2 E Var(Y|witness). Feed both observed bits as
            # ordinary (zero or one) load rails of the new native write.
            pred = response(n, var / 2)
            row[name] = dict(
                variance=str(var), hazard=str(var / 2), transit=float(pred)
            )
        # Independent simulation with fresh sampled source pairs each appointment.
        # Actual binary rails, native write labels and probe pulse determine hits.
        rng = np.random.default_rng(seed + len(rows))
        tickets, _, _, _ = native()
        loc = np.zeros(SHOTS, int)
        passed = returned = 0
        for tick in range(LIMIT):
            step = rng.integers(4, size=len(loc))
            loc = (loc + (step == 0).astype(int) - (step == 3)) % n
            left, right = rng.choice(
                len(w), size=(2, len(loc)), p=np.asarray(weights, float)
            )
            word = 3 * (w[left, 1] != 1) + (w[right, 1] != 1)
            mark = write[word, rng.integers(64, size=len(word))] // 9
            inp = 1 + 2 * rng.integers(2, size=len(word))
            out = tickets[9 * mark + inp, rng.integers(4, size=len(word))]
            yes = out % 9 != inp
            passed += int(yes.sum())
            returned += int(((loc == 0) & ~yes).sum())
            loc = loc[(loc != 0) & ~yes]
            if not len(loc):
                break
        assert passed + returned + len(loc) == SHOTS
        delta = np.sqrt(np.log(2 * 500 / 0.01) / (2 * SHOTS))
        band = [passed / SHOTS - delta, (passed + len(loc)) / SHOTS + delta]
        assert band[0] < row["unresolved"]["transit"] < band[1]
        reader = np.random.default_rng(seed + 100 + len(rows))
        nsource = 2**20
        s = reader.choice(len(w), size=nsource, p=np.asarray(weights, float))
        count = int((w[s, 1] != 1).sum())
        eps = np.sqrt(np.log(800) / (2 * nsource))
        source_band = [count / nsource - eps, count / nsource + eps]
        assert 0 < source_band[0] < source_band[1] < 0.5
        transformed = [
            float(response(n, F(float(f * (1 - f) / 2)))) for f in source_band
        ]
        assert max(transformed[0], band[0]) <= min(transformed[1], band[1])
        row["source_meter"] = dict(
            trials=nsource,
            nonunit=count,
            rate_band=source_band,
            transit_band=transformed,
        )
        row["original_middle_write"] = middle_meter(
            w, weights, cells, n, seed + 300 + 2 * len(rows), source_band
        )
        row["simulated_pair"] = dict(
            trials=SHOTS, transits=passed, returns=returned, pending=len(loc), band=band
        )
        rows[arm] = row
    return rows


def run():
    boards = {}
    seed = 3490872
    for level in range(3):
        route, sidings = ring(level)
        n = len(sidings)
        rows = []
        for k in sorted({0, 1, 4, n // 4, n}):
            for layout in ("even", "cluster") if 1 < k < n else ("even",):
                ids = np.arange(k) * n // k if k and layout == "even" else np.arange(k)
                sites = np.zeros(n, dtype=np.int8)
                sites[ids] = 1
                fresh = response(n, F(k, 4 * n))
                assert exact_first_step(n, F(k, 4 * n)) == fresh
                local = float(local_response(sites / 4))
                cap = local_cap(n, k)
                assert local <= float(cap) + 1e-12
                assert (
                    abs(local_response(np.full(n, k / (4 * n))) - float(fresh)) < 1e-10
                )
                moving = float(moving_response(sites))
                gap = int(max(np.diff(np.r_[ids, ids[0] + n]))) if k else 0
                moving_floor = (
                    (F(1, 4) * F(comb(2 * (gap - 1), gap - 1), 4 ** (gap - 1)))
                    if k
                    else F(0)
                )
                assert float(moving_floor) <= moving + 1e-13
                row = dict(
                    moving=moving,
                    longest_gap=gap,
                    moving_floor=float(moving_floor),
                    k=k,
                    layout=layout,
                    mixed=float(fresh),
                    resident=local,
                    ideal_amplitude=np.sqrt(k / n) / 2,
                    all_placement_upper=str(cap),
                    mixed_mean_appointments=float(fresh / F(k, 4 * n)) if k else n,
                )
                if level <= 1 and (layout == "even" or k == 4) and k > 0:
                    for name, pred in (
                        ("pooled", float(fresh)),
                        ("local", local),
                        ("moving", moving),
                    ):
                        seed += 1
                        records = trial_bank(n, sites, name, seed)
                        lo, hi = records["probability_band"]
                        assert lo <= pred <= hi
                        row[name + "_records"] = records
                rows.append(row)
        boards[str(level)] = dict(
            route=route, sidings=sidings, half_length=n, rows=rows
        )
    n = boards["0"]["half_length"]
    active = np.ones(n, np.int8)
    control = trial_bank(n, active, "pooled", seed + 20, no_write=True)
    double = trial_bank(n, active, "pooled", seed + 21, load=2)
    assert control["transits"] == SHOTS and double["transits"] == 0
    sites = np.zeros(n, np.int8)
    sites[np.arange(4) * n // 4] = 1
    alternatives = {}
    for j, alignment in enumerate(("common", "opposite")):
        r = trial_bank(n, sites, "pooled", seed + 23 + j, alignment=alignment)
        assert r["source_second_moment"] == (16 if j == 0 else 0)
        low, hi = r["probability_band"]
        assert low < float(response(n, F(1, n))) < hi
        alternatives[alignment] = r
    return dict(
        orientation_controls=alternatives,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "FLUCTUATION_RETURN.md",
                HERE / "closed_pairing.json",
                HERE / "closed_passage.json",
            )
        },
        boards=boards,
        closed_source=source_experiment(n, seed + 22),
        exchange_only_preparation=control,
        double_probe=double,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dst = HERE / "fluctuation_return.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(json.loads(dst.read_bytes()), result)
    print({k: r["half_length"] for k, r in result["boards"].items()})
