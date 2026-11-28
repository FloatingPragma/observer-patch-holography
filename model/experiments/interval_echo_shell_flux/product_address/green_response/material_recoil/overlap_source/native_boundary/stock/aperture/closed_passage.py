"""Predict ordinary first-passage clocks from the local records of the closed source."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog
from scipy.stats import beta
from scipy.sparse.csgraph import laplacian, shortest_path
from sympy import Matrix, linsolve, symbols

from . import closed_energy as energy

HERE = energy.HERE
pair = energy.pair
SLOTS = 6
CLOCKS = 32768
SOURCE_SHOTS = 2**24


def blank_range(w, cells, codes, h, ids, groups, returned, z):
    rows = np.vstack(
        ([cells == i for i in range(4)], w[:, :4].T, h.diagonal()[codes], returned)
    )
    A = pair.group_matrix(rows, ids, groups)
    b = rows @ z
    blank = pair.group_matrix((w[:, 6] == 0)[None, :], ids, groups)[0]
    records = {}
    for key, sign in (("min", 1), ("max", -1)):
        opt = linprog(sign * blank, A_eq=A, b_eq=b.astype(float), bounds=(0, None))
        assert opt.success
        active = np.flatnonzero(opt.x > 1e-9)
        (exact,) = linsolve((Matrix(A[:, active]), Matrix(b.tolist())))
        weights = np.full(groups, Q(0), dtype=object)
        weights[active] = [Q(int(v.p), int(v.q)) for v in exact]
        dual = np.array(pair.rational(opt.eqlin.marginals), dtype=object)
        assert all(weights >= 0) and np.all(A @ weights == b)
        assert np.all(dual @ A <= sign * blank) and dual @ b == sign * (blank @ weights)
        records[key] = dict(
            blank_fraction=str(blank @ weights / len(w)),
            component_weights=list(map(str, weights)),
            dual=list(map(str, dual)),
        )
    return records


def full_pulses(heads):
    _, half, h, h0 = pair.stock.kernels()
    identity = np.eye(99, dtype=int)
    contrast = h - h0
    # In the native four-H units the commuting contrast spectrum is
    # {0,32,128}, so all its phases are one at the full exchange pulse.
    np.testing.assert_array_equal(h0 @ h0, 4 * h0)
    np.testing.assert_array_equal(h0 @ contrast, contrast @ h0)
    np.testing.assert_array_equal(
        contrast @ (contrast - 32 * identity) @ (contrast - 128 * identity), 0
    )
    tickets = (2 * identity - h0) ** 2
    np.testing.assert_array_equal(tickets, half // 8 - 4 * identity)
    assert np.all(half % 8 == 0) and tickets.min() == 0
    np.testing.assert_array_equal(tickets.sum(0), 4)
    np.testing.assert_allclose(
        tickets / 4, abs(expm(-1j * np.pi * h / 4)) ** 2, atol=1e-13
    )
    # Each appointment ends in a sharp material and memory read. Written
    # status and magnitude persist, including pulses on empty routing rails.
    old, new = np.nonzero(tickets)
    np.testing.assert_array_equal(old // 9 == 0, new // 9 == 0)
    np.testing.assert_array_equal(pair.stock.MAG[old // 9], pair.stock.MAG[new // 9])
    transfers = {}
    for load in (1, 2):
        rows = []
        for m in range(11):
            sides = []
            for s in (3 * load, load):
                t = tickets[:, 9 * m + s].reshape(11, 9).sum(0)
                assert t[s] + t[3 * load + load - s] == 4
                sides.append(int(4 - t[s]))
            rows.append(sides)
        transfers[str(load)] = rows
    assert all(
        row == ([4, 4] if m == 0 else [0, 0]) for m, row in enumerate(transfers["2"])
    )
    assert set(
        np.flatnonzero(np.any(np.array(transfers["1"]) != transfers["2"], axis=1))
    ) == set(np.flatnonzero(pair.stock.MAG == 1))
    # For each enrolled label the parked survey sees an empty or one-probe sector.
    empty = tickets[:, 9 * np.unique(heads)].reshape(11, 9, -1).sum(0)
    np.testing.assert_array_equal(empty[0], 4)
    return transfers


def contacts(board, seams):
    adj = board.adjacency().toarray().astype(int)
    for u, v in seams:
        adj[u, v] = adj[v, u] = 0
    neighbors = np.tile(np.arange(board.n)[:, None], (1, SLOTS))
    for x in range(board.n):
        live = np.flatnonzero(adj[x])
        neighbors[x, : len(live)] = live
    return adj, neighbors


def dipole(adj, u, v):
    # A grounded integer Poisson problem. Sparse rational elimination and
    # substitution certify the potential and the shell averages exactly.
    L = laplacian(adj)
    n = len(adj)
    rhs = np.zeros(n, int)
    rhs[u], rhs[v] = 1, -1
    x = symbols(f"v:{n - 1}")
    eq = [
        int(L[i, i]) * x[i]
        - sum(x[j] for j in np.flatnonzero(adj[i, :-1]))
        - int(rhs[i])
        for i in range(n - 1)
    ]
    (sol,) = linsolve(eq, x)
    phi = np.array([Q(int(f.p), int(f.q)) for f in sol] + [Q(0)], dtype=object)
    assert np.array_equal(L @ phi, rhs)
    return phi - sum(phi) / n


def passage_moments(adj, start, stop):
    # In actual local appointment ticks P=I-L/6. Dirichlet first-step
    # equations give E[T] and E[T²]; two freshly randomized legs add.
    live = np.arange(len(adj)) != stop
    L = laplacian(adj)[np.ix_(live, live)]
    mean = np.linalg.solve(L, np.full(len(L), SLOTS))
    second = np.linalg.solve(L, SLOTS * (2 * mean - 1))
    k = start - int(start > stop)
    return mean[k], second[k] - mean[k] ** 2


def clock_trips(table, a, b, rng):
    where = np.full(CLOCKS, a)
    returning = np.zeros(CLOCKS, bool)
    count, ticks, total = 0, 0, 0
    while len(where):
        ticks += 1
        where = table[where, rng.integers(SLOTS, size=len(where))]
        returning |= where == b
        done = returning & (where == a)
        count += int(done.sum())
        total += ticks * int(done.sum())
        where, returning = where[~done], returning[~done]
    assert count == CLOCKS
    return dict(completions=count, last_tick=ticks, sum_ticks=total)


def short_gate(table, a, b, rng):
    # No source label or transition probabilities are supplied to the terminal reader.
    n = 2**19
    state = np.full(n, a)
    arrived = np.zeros(n, bool)
    for _ in range(3):
        state = table[state, rng.integers(SLOTS, size=n)]
        arrived |= state == b
    x = np.eye(len(table), dtype=np.int64)[b]
    for k in range(1, 4):
        x = x[table].sum(1)
        x[b] = SLOTS**k
    p = Q(int(x[a]), SLOTS**3)
    assert abs(arrived.mean() - float(p)) < np.sqrt(np.log(1600) / (2 * n))
    return dict(trials=n, arrived=int(arrived.sum()), probability=str(p))


def pulse_table(board, seams, load, tickets, blind):
    """24 equally likely tickets per port and head: six slots and four pulse outcomes."""
    n = board.n
    table = np.empty((11 * n, 4 * SLOTS), int)
    slots = contacts(board, [])[1]
    for m in range(11):
        for x, options in enumerate(slots):
            for slot, y in enumerate(options):
                edge = tuple(sorted((x, y)))
                dest = np.full(4, m * n + y)
                if not blind and edge in seams[:2]:
                    dest[:] = m * n + x
                elif not blind and edge == seams[2]:
                    left, right = edge
                    word = 3 * load if x == left else load
                    pulse = np.repeat(np.arange(99), tickets[:, 9 * m + word])
                    head, material = np.divmod(pulse, 9)
                    assert set(material) <= {3 * load, load}
                    dest = head * n + np.where(material == 3 * load, left, right)
                table[m * n + x, 4 * slot : 4 * (slot + 1)] = dest
    return table


def held_gate(board, seams, outside, w, z, transfers):
    # Evolve integer counts through three actual addressed appointments. Every
    # load is detected from the outside. The terminal observer uses no head read.
    _, half, _, _ = pair.stock.kernels()
    tickets = half // 8 - 4 * np.eye(99, dtype=int)
    n = board.n
    a, b = outside
    reports = {}
    for k, (label, load, blind) in enumerate(
        (("single", 1, False), ("double", 2, False), ("memory_blind", 2, True))
    ):
        table = pulse_table(board, seams, load, tickets, blind)
        # Absorb at the first detector arrival.
        for h in range(11):
            table[h * n + b] = h * n + b
        exact = (np.arange(11 * n) % n == b).astype(np.int64)
        for _ in range(3):
            exact = exact[table].sum(1)
        numerators = exact[a::n]
        if not blind:
            # Only one length-three itinerary crosses the third seam.
            np.testing.assert_array_equal(
                numerators, 16 * np.array(transfers[str(load)])[:, 0]
            )
        else:
            assert len(set(numerators)) == 1
        rows = {}
        for j, (tag, weights) in enumerate(
            (("original", np.full(len(w), Q(1), object)), ("repaired", z))
        ):
            shots = 2**30
            rng = np.random.default_rng(738061 + 2 * k + j)
            drawn = rng.multinomial(shots, np.asarray(weights, float) / len(w))
            counts = np.zeros(11 * n, np.int64)
            np.add.at(counts, w[:, 6] * n + a, drawn)
            for _ in range(3):
                after = np.zeros_like(counts)
                for x in np.flatnonzero(counts):
                    split = rng.multinomial(int(counts[x]), np.full(24, 1 / 24))
                    np.add.at(after, table[x], split)
                counts = after
            assert counts.sum() == shots
            arrived = int(counts[b::n].sum())
            p = weights @ numerators[w[:, 6]] / (len(w) * 24**3)
            lo, hi = (
                beta.ppf(0.0005, arrived, shots - arrived + 1),
                beta.isf(0.0005, arrived + 1, shots - arrived),
            )
            assert lo <= float(p) <= hi
            rows[tag] = dict(
                trials=shots,
                arrived=arrived,
                expected=str(p),
                interval=[float(lo), float(hi)],
            )
        bands = [r["interval"] for r in rows.values()]
        rows["difference_interval"] = [
            bands[1][0] - bands[0][1],
            bands[1][1] - bands[0][0],
        ]
        reports[label] = rows
    return reports


def run():
    _, w, cells, _, codes, h, ids, groups = pair.setup()
    board, home, seams, _, _, prices = energy.returned(w)
    z = energy.least_change(cells)
    blank = w[:, 6] == 0
    transfers = full_pulses(w[:, 4:])
    source = {}
    for seed, (label, weights) in enumerate(
        (("original", np.ones(len(w), int)), ("repaired", z))
    ):
        rng = np.random.default_rng(733061 + seed)
        p = sum(Q(weight) for weight in weights[blank]) / len(w)
        draw = rng.multinomial(SOURCE_SHOTS, np.asarray(weights, float) / len(w))
        count = int(draw[blank].sum())
        d = np.sqrt(np.log(800) / (2 * SOURCE_SHOTS))
        assert abs(count / SOURCE_SHOTS - float(p)) < d
        active = np.isin(w[:, 6], np.flatnonzero(np.array(transfers["1"])[:, 0] == 2))
        source[label] = dict(
            blank_fraction=str(p),
            trials=SOURCE_SHOTS,
            blanks=count,
            single_half_crossings=int(draw[active].sum()),
            halfwidth=float(d),
        )
    # Use the same blank apparatus contacts as the protected material echoes.
    # Store and restore all four loads on every accepted input word, outside
    # the timed survey. These four side rails belong to the apparatus.
    blank_swap = (2 * np.eye(9, dtype=int) - pair.stock.kernels()[3][:9, :9]).argmax(0)
    parked_words = blank_swap[3 * w[:, :4]]
    assert np.array_equal(parked_words, w[:, :4])
    restored = blank_swap[parked_words]
    np.testing.assert_array_equal(restored, 3 * w[:, :4])
    bounds = blank_range(w, cells, codes, h, ids, groups, prices, z)
    closed, closed_slots = contacts(board, seams)
    opened, open_slots = contacts(board, seams[:2])
    u, v = seams[2]
    for table in (closed_slots, open_slots):
        # A uniform launch remains uniform after each actual appointment.
        np.testing.assert_array_equal(np.bincount(table.ravel()), SLOTS)
    phi = dipole(closed, u, v)
    R = phi[u] - phi[v]
    factor = SLOTS * board.n / (1 + R)
    hops = shortest_path(closed, indices=[u, v], unweighted=True).min(0).astype(int)
    shell = [
        sum(phi[hops == i] ** 2) / int(sum(hops == i)) * factor
        for i in range(max(hops) + 1)
    ]
    forbidden = set(home) | {u, v}
    outside = [min(set(np.flatnonzero(closed[x])) - forbidden) for x in (u, v)]
    gate = held_gate(board, seams, outside, w, z, transfers)
    for name, bank in source.items():
        for receiver, extra in (("double", 0), ("single", 1)):
            estimated = (bank["blanks"] + extra * bank["single_half_crossings"] / 2) / (
                216 * bank["trials"]
            )
            width = bank["halfwidth"] * (1 + extra / 2) / 216
            forecast = [estimated - width, estimated + width]
            measured = gate[receiver][name]["interval"]
            assert max(forecast[0], measured[0]) <= min(forecast[1], measured[1])
            gate[receiver][name]["from_independent_source_bank"] = forecast
    # All exterior vertices are harmonic. A squared harmonic voltage has
    # nonpositive graph Laplacian there. Check the effective gravitational
    # flux on every enclosing graph ball with exact fractions.
    L = laplacian(closed)
    mass = L @ (phi**2)
    assert sum(mass) == 0 and all(mass[hops > 0] <= 0)
    enclosed = [factor * sum(mass[hops <= i]) for i in range(max(hops) + 1)]
    assert all(
        m1 >= m2 >= 0 for m1, m2 in zip(enclosed[:-1], enclosed[1:], strict=True)
    )
    reports = {}
    for i, (a, b) in enumerate(((u, v), outside)):
        difference = factor * (phi[a] - phi[b]) ** 2
        rows = {}
        for j, (tag, adj, table) in enumerate(
            (("closed", closed, closed_slots), ("open", opened, open_slots))
        ):
            m, var = np.sum(
                [passage_moments(adj, a, b), passage_moments(adj, b, a)], axis=0
            )
            rng = np.random.default_rng(734061 + 2 * i + j)
            record = clock_trips(table, a, b, rng)
            # Chebyshev + union for four means, family coverage at least 99%.
            width = 20 * np.sqrt(var / CLOCKS)
            assert abs(record["sum_ticks"] / CLOCKS - m) < width
            rows[tag] = dict(
                mean=float(m),
                variance=float(var),
                halfwidth=float(width),
                **record,
                short_gate=short_gate(table, a, b, rng),
            )
        assert (
            abs(rows["closed"]["mean"] - rows["open"]["mean"] - float(difference))
            < 1e-8
        )
        reports[f"{a},{b}"] = dict(closed_minus_open=str(difference), arms=rows)
    dp = Q(source["repaired"]["blank_fraction"]) - Q(
        source["original"]["blank_fraction"]
    )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOSED_PASSAGE.md",
                HERE / "closed_steering.json",
                HERE / "closed_energy.json",
            )
        },
        accepted=len(w),
        parked_and_restored=len(parked_words),
        ports=board.n,
        participating_heads=[list(e) for e in seams],
        full_exchange_quarters=transfers,
        single_exception_words=int(np.sum(pair.stock.MAG[w[:, 6]] == 1)),
        local_source=source,
        repaired_minus_native_blank=str(dp),
        fixed_price_blank_extremes=bounds,
        unit_current_closed_resistance=str(R),
        shell_populations=np.bincount(hops).tolist(),
        closed_minus_open_mean_partner_time_up_to_constant=list(map(str, shell)),
        clocks=reports,
        direct_exterior_gate=gate,
        enclosed_flux=list(map(str, enclosed)),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    result = run()
    path = HERE / "closed_passage.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        old = json.loads(path.read_bytes())
        from . import measure

        measure.check(old, result)
    print(
        {
            k: v["blank_fraction"]
            for k, v in result["fixed_price_blank_extremes"].items()
        }
    )
