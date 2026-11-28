"""An exterior release triggered only by several locally measured repairs."""

import argparse
from collections import deque
from dataclasses import replace
from fractions import Fraction as Q
from hashlib import sha256
from itertools import pairwise
import json
import math

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import dijkstra

from . import collective_circuit as circuit
from . import collective_geometry as geom
from . import collective_source as source
from . import collective_tide as tide

HERE = source.tidal.HERE


def enroll(board, s, ports, start, L):
    circ = ports + ports
    X, V, Y = circ[start], circ[start + L - 1], circ[start + L]
    center = circ[start + (L - 1) // 2]
    active = {e: j for j, e in enumerate(s.contacts)}
    plus = lil_matrix((board.n, board.n), dtype=float)
    chosen = {}
    for u, v in list(pairwise(ports)) + [(ports[-1], ports[0])]:
        seam = circuit.edge(u, v)
        if seam in active:
            chosen[seam] = active[seam]
        else:
            assert s.blank[u, v]
            plus[u, v] = plus[v, u] = 1
    visited = set(ports)
    for u in (X, V):
        for v in board.neighbours[u]:
            assert s.blank[u, v]
            plus[u, v] = plus[v, u] = 1
            visited.add(v)
    queue = deque(sorted(visited))
    while queue:
        u = queue.popleft()
        for v in s.blank[u].indices:
            v = int(v)
            if v not in visited:
                plus[u, v] = plus[v, u] = 1
                visited.add(v)
                queue.append(v)
    assert len(visited) == board.n and len(chosen) == 4
    upper = plus.copy()
    for a, b in chosen:
        upper[a, b] = upper[b, a] = 4
    radii = dijkstra(upper.tocsr(), directed=False, indices=center).astype(int)
    assert len(ports) + len(visited - set(ports)) == board.n
    active_chart = replace(
        s, blank=plus.tocsr(), contacts=list(chosen), radii=radii, center=center
    )

    def witness(chart):
        stations = np.array(sorted({center, X, V}), int)
        index = np.full(board.n, -1, int)
        index[stations] = np.arange(len(stations))
        return replace(
            chart, pairs=np.array([[X, V]]), starts=stations, index=index, center=center
        )

    return witness(active_chart), witness(s), list(chosen.values()), (X, V, Y, center)


def measured_sum(volume, tags, singles, chosen, rng, comparisons, trials=2**26):
    p = source.chances(volume, tags, singles)
    stats = source.records()[-1]
    N = trials
    eps = np.sqrt(np.log(2 * comparisons / 0.01) / (2 * N))
    banks = []
    for j in chosen:
        a, b, _, _ = stats[tags[j]]
        c = 4 * (p[j] - a) - b
        assert 0 <= c <= b
        counts = rng.multinomial(N, np.array([a, b - c, c, 1 - a - b], float))
        score = (counts[0] + counts[1] / 4 + counts[2] / 2) / N
        banks.append(
            dict(
                contact=j,
                counts=counts.tolist(),
                p_lower=score - eps,
                p_upper=score + eps,
            )
        )
    return dict(
        trials_per_contact=N,
        banks=banks,
        combined_cost=[
            sum(1 / v["p_upper"] for v in banks),
            sum(1 / v["p_lower"] for v in banks),
        ],
    )


def run():
    board, layout, s, mode, tried, (ports, i, L) = circuit.find()
    selected, whole, chosen, (X, V, Y, C) = enroll(board, s, ports, i, L)
    m = len(layout)
    original = ["original"] * m
    repaired = ["source"] * m
    p = source.chances(board.n, original)[0]
    q = source.chances(board.n, repaired)[0]
    b = source.records()[2]["original"][1]
    step = b / (4 * board.n)
    exact_n = (q - p) / step
    threshold = math.floor((Q(4, 13) - p) / step) + 1
    assert 0 < threshold < exact_n and math.ceil(exact_n) <= board.n - 4 * m
    settings = {
        "original": (original, 0),
        "all": (repaired, 0),
        "ordinary_one_per_patch": (original, m),
        "ordinary_double_per_patch": (original, 0),
        "ordinary_below": (original, threshold - 1),
        "ordinary_first": (original, threshold),
        "ordinary_equivalent": (original, exact_n),
    }
    for k in range(m):
        one = original.copy()
        one[k] = "source"
        settings[f"single_{k}"] = (one, 0)
    for n in (2, 3):
        tags = original.copy()
        for j in chosen[:n]:
            tags[j] = "source"
        settings[f"on_route_{n}"] = (tags, 0)
    rng = np.random.default_rng(931672)
    out = {}
    for name, (tags, added) in settings.items():
        # Adding a double leaves the single count, head class and gate
        # response fixed. Every collection has enough empty resident sites.
        prices = [1 / x for x in source.chances(board.n, tags, added)]
        assert max(prices) < 4
        local = [prices[k] for k in chosen]
        read = geom.measure(selected, local)
        d, stars, paths = geom.measure(selected, local, [X, V])
        assert read.equal_radial[0] and d == L - 1 and min(s.clearance[[X, V]]) >= 2
        comparison = geom.measure(whole, prices)
        row = dict(
            cost=str(sum(local)),
            chord=str(d),
            tide=str(read.curvature[0]),
            unrestricted_tide=str(comparison.curvature[0]),
            unrestricted_chord=str(comparison.length[0]),
            cross_ranges=[
                [int(a), int(b), str(tx)]
                for (a, b), tx in zip(paths, stars, strict=True)
            ],
        )
        if name in ("original", "all", "ordinary_equivalent") or name.startswith(
            "single_"
        ):
            row["local_measurement"] = measured_sum(
                board.n, tags, added, chosen, rng, 4 * (3 + m)
            )
            lo, hi = row["local_measurement"]["combined_cost"]
            if name in ("all", "ordinary_equivalent"):
                assert hi < 13
            else:
                assert lo > 13
        if name in ("on_route_2", "on_route_3"):
            row["local_measurement"] = measured_sum(
                board.n, tags, added, chosen, rng, 8, 2**32
            )
            lo, hi = row["local_measurement"]["combined_cost"]
            assert hi < 13 if name == "on_route_3" else lo > 13
        if name in ("original", "all", "on_route_3"):
            row["held_routes"] = tide.interval_acquisition(selected, local, [X, V], rng)
        out[name] = row
    zero = Q(out["original"]["tide"])
    for name, row in out.items():
        row["change"] = str(Q(row["tide"]) - zero)
        row["unrestricted_change"] = str(
            Q(row["unrestricted_tide"]) - Q(out["original"]["unrestricted_tide"])
        )
    assert all(Q(out[f"single_{j}"]["change"]) == 0 for j in range(m))
    assert (
        Q(out["all"]["change"]) < 0 and Q(out["ordinary_one_per_patch"]["change"]) == 0
    )
    assert (
        Q(out["ordinary_below"]["change"]) == 0
        and Q(out["ordinary_first"]["change"]) < 0
    )
    assert out["all"]["change"] == out["ordinary_equivalent"]["change"]
    ref, src = (out[k]["held_routes"] for k in ("original", "all"))
    held = [src["lower"] - ref["upper"], src["upper"] - ref["lower"]]
    assert held[1] < 0
    triple = out["on_route_3"]["held_routes"]
    held_three = [triple["lower"] - ref["upper"], triple["upper"] - ref["lower"]]
    assert held_three[1] < 0
    # Test all 256 preparations on the actual release stars, preserving
    # the source locations even when the local counts are exchangeable.
    dose = {}
    for bits in range(1 << m):
        tags = ["source" if bits >> i & 1 else "original" for i in range(m)]
        on = sum(tags[i] == "source" for i in chosen)
        off = tags.count("source") - on
        prices = [1 / x for x in source.chances(board.n, tags)]
        local = [prices[k] for k in chosen]
        read = geom.measure(selected, local)
        correction = read.curvature[0] - zero
        assert (correction < 0) == (sum(local) < 13) and read.length[0] == L - 1
        key = f"{on}/{off}"
        record = dose.setdefault(
            key, dict(subsets=0, change=str(correction), price=str(sum(local)))
        )
        assert record["change"] == str(correction) and record["price"] == str(
            sum(local)
        )
        record["subsets"] += 1
    assert sum(v["subsets"] for v in dose.values()) == 2**m
    volume, chains, offset = source.ordinary_coefficients()
    assert volume * board.n + chains * m + offset == exact_n
    return dict(
        dose=dose,
        ordinary_equivalent_gate=source.pooled_draw(
            board.n, original, chosen[0], exact_n, 477030
        ),
        ordinary_coefficients=list(map(str, (volume, chains, offset))),
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "COLLECTIVE_RELAY_LAB.md", HERE / "collective_tide.json")
        },
        mode=mode,
        owner_cycles_examined=tried,
        source_owner_paths=[p for p, _ in layout],
        cycle=ports,
        start=X,
        probe=V,
        range_target=Y,
        center=C,
        blank_detour=L,
        clearance=s.clearance[[X, V]].tolist(),
        transmitting_sources=chosen,
        first_integer_ordinary_response=threshold,
        equivalent_ordinary_singles=str(exact_n),
        upper_integer_fraction=str(exact_n - math.floor(exact_n)),
        held_change_interval=held,
        held_three_change_interval=held_three,
        preparations=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "collective_relay_lab.json"
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        dst.write_text(raw)
    else:
        assert dst.read_text() == raw
    print(
        result["held_change_interval"],
        result["first_integer_ordinary_response"],
        result["equivalent_ordinary_singles"],
    )
