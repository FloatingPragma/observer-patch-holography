"""An actual five-owner cycle containing the two original conversion ports."""

from functools import cache
import networkx as nx
from model.executable_model import federation
from ....apparatus import fixture


@cache
def layout():
    f = fixture()
    board = federation(0)
    g = nx.Graph(board.seams)
    owners = nx.Graph((u // 12, v // 12) for u, v in board.external)
    local = g.subgraph(range(156, 168))
    gates = {
        (a // 12, b // 12): a for u, v in board.external for a, b in ((u, v), (v, u))
    }
    records = []
    for neighbor in owners[13]:
        for path in nx.all_simple_paths(owners, neighbor, 13, cutoff=4):
            if len(path) != 5 or 0 in path:
                continue
            cycle = [13, *path[:-1]]
            if cycle[1] > cycle[-1]:
                continue
            for first in nx.all_simple_paths(
                local, gates[13, cycle[-1]], gates[13, cycle[1]], cutoff=6
            ):
                if frozenset((160, 167)) not in map(frozenset, zip(first, first[1:])):
                    continue
                marks = [tuple(e) for e in f["marked"] if set(e) & set(first)]
                if any(len(set(e) & set(first)) != 1 for e in marks):
                    continue
                ring = list(first)
                for n, z in enumerate(cycle[1:], 1):
                    a = gates[z, cycle[n - 1]]
                    b = gates[z, cycle[(n + 1) % 5]]
                    ring += min(
                        nx.all_shortest_paths(
                            g.subgraph(range(12 * z, 12 * z + 12)), a, b
                        )
                    )
                records.append(
                    (len(marks), len(ring), tuple(cycle), tuple(ring), tuple(marks))
                )
    m, n, owners, ring, marks = min(records)
    assert m == 2 and n == 14 and len(records) == 35
    assert all(g.has_edge(a, b) for a, b in zip(ring, ring[1:] + ring[:1]))
    leaves = tuple(next(iter(set(e) - set(ring))) for e in marks)
    return ring, marks, leaves, owners, len(records)


def configure(sign, retire=False):
    """Park both signed products, open the two collars and commission one ring."""
    import numpy as np
    from model.experiments.interval_echo_consensus.signed_effects.collar_storage import (
        MEMORIES,
        gates,
    )
    from ..apparatus import stock
    from ....echo_turn.material import chart, transp
    from ....echo_turn.pulse import primitive
    from ....material_interferometer.apparatus import Terminal, inverse

    f = fixture()
    lab = chart()
    original = stock(sign)
    board = federation(0)
    field = nx.Graph(board.seams)
    ring, marks, leaves, _, _ = layout()
    active = set(ring) | set(leaves)

    def edge(a, b):
        return frozenset((a, b))

    stock_pair = set(f["donor"])
    guide = set(f["ports"])
    live = {edge(*z) for z in f["marked"]}
    ordinary = {edge(*z) for z in f["blank"]}
    old_cuts = {
        edge(a, b)
        for a, b in board.seams
        if ({a, b} & (guide | stock_pair))
        and edge(a, b) not in live | ordinary | {frozenset(stock_pair)}
    }
    collars = [(160, 167), (158, 167)]
    for opened in [*collars, (167, 166), lab["service"]]:
        old_cuts.remove(edge(*opened))
    blocked = old_cuts | live | {edge(*z) for z in collars}
    route = field.copy()
    route.remove_edges_from(tuple(e) for e in blocked)
    protected = {edge(*e) for e in marks}
    ring_edges = {edge(a, b) for a, b in zip(ring, ring[1:] + ring[:1])}
    far = lab["storage"] + lab["probes"] + lab["stock"]
    load = np.ones(board.n, int)
    load[160] = 1 + sign
    load[167] = 1 - sign
    load[f["donor"][0]] = 1 - sign
    load[far[5]] = load[far[7]] = 1 + sign
    load[far[8]] = 1 - sign
    assert load.sum() == board.n
    distant = nx.multi_source_dijkstra_path_length(field, active)
    allowed = nx.node_connected_component(route, 160)
    slots = [
        v
        for v in sorted(allowed)
        if distant[v] > 2
        and v not in set(far) | guide | stock_pair
        and load[v] == 1
        and v != lab["reference"]
    ]
    assert len(slots) > 12
    log = []
    parked = {}

    def move(u, v, avoid=()):
        local = route.subgraph(set(route) - set(avoid))
        p = nx.shortest_path(local, u, v)
        initial = load.copy()
        initial[u], initial[v] = initial[v], initial[u]
        for a, b in transp(p):
            load[a], load[b] = load[b], load[a]
        np.testing.assert_array_equal(load, initial)
        log.append([u, v, len(transp(p))])

    for u in (160, 167, far[5], far[7], far[8]):
        v = slots.pop(0)
        assert load[v] == 1
        move(u, v)
        parked[u] = v
    # Opening an old zero head at the exact unit-pair stock restores a blank.
    # Leave remotely parked sharp stock and its original compensator undisturbed.
    for e in old_cuts:
        a, b = sorted(e)
        if load[a] == load[b] == 1:
            route.add_edge(a, b)
    _, pre, echo = gates(1.0)
    from model.experiments.interval_echo_consensus.signed_effects.collar_storage import (
        LAW,
    )

    column = np.kron(MEMORIES[:, 0], np.eye(9)[LAW.SRC.index((1, 1))])
    restored = pre.conj().T @ echo @ column
    assert abs(restored[LAW.SRC.index((2, 0))]) ** 2 > 1 - 1e-12
    for u, v in collars:
        assert load[u] == load[v] == 1 and v == 167 and board.degree[v] == 5
        load[u], load[v] = 2, 0
        route.add_edge(u, v)
        for z in (u, v):
            move(z, slots.pop(0))
    terminal = Terminal("tip")
    if retire:
        # Read both outcomes of each unknown head with a material vacancy,
        # undo its pointer motion, and retain the bit for the known-stock reset.
        for bit, (a, _, g) in enumerate(terminal.arms):
            if (a, g) not in marks:
                continue
            go = [(terminal.hub, a, np.pi), (a, g, np.pi / 8)]
            for mask in range(8):
                final = terminal.run(go, mask)
                pointer = a if (mask >> bit) & 1 else g
                assert (
                    abs(np.sum(abs(final[terminal.ports.index(pointer)]) ** 2) - 1)
                    < 1e-12
                )
                for u, v, dt in inverse(go):
                    terminal.gate(final, u, v, dt)
                np.testing.assert_allclose(final[0, mask], 1, atol=1e-12)
        # The loan is a material vacancy; either sign of the incoming excitation
        # leaves a vacancy among the parked balanced stocks.
        probe = next(v for v in parked.values() if load[v] == 0)
        for a, g in marks:
            for bit in (0, 1):
                # The addressed directional test and inverse distinguish and
                # restore the two stock values. Check each possible reset route.
                move(probe, a)
                assert load[a] == 0 and load[g] == 1
                move(probe, a)
                fixed = (a, g)[bit]
                other = (g, a)[bit]
                move(probe, fixed)
                assert load[fixed] == 0 and load[other] == 1
                p = nx.shortest_path(
                    route.subgraph(set(route) - {fixed}), other, lab["reference"]
                )
                assert p[0] == other and p[-1] == lab["reference"]
                primitive(int(other > fixed))
                move(probe, fixed)
            route.add_edge(a, g)
    # Erase ordinary private references of new degree-five positions with
    # the native (0,0) zero-contrast write. Others are already erased.
    probe = next(v for v in parked.values() if load[v] == 0)
    for v in sorted(active):
        if board.degree[v] == 5 and v not in guide:
            move(probe, v)
            assert load[v] == 0
            move(probe, v)
    target_edges = ring_edges | (set() if retire else protected)
    new_cuts = {
        edge(a, b)
        for a, b in board.seams
        if {a, b} & active and edge(a, b) not in target_edges
    }
    assert all(load[a] == load[b] == 1 for a, b in (tuple(e) for e in new_cuts))
    assert not (new_cuts & live - protected)
    # A single ordinary hatch lets the original two loads return without
    # reopening a gate once the independently timed release has begun.
    for hatch in sorted(tuple(sorted(e)) for e in new_cuts):
        if not route.has_edge(*hatch) or set(hatch) & {160, 167}:
            continue
        trial = route.copy()
        trial.remove_edges_from(tuple(e) for e in new_cuts - {edge(*hatch)})
        if all(nx.has_path(trial, parked[v], v) for v in (160, 167)):
            route = trial
            break
    else:
        raise RuntimeError("no returning stock hatch")
    for v in (160, 167):
        move(parked[v], v)
    assert load[hatch[0]] == load[hatch[1]] == 1
    route.remove_edge(*hatch)
    assert all(
        load[v] == (1 + sign if v == 160 else 1 - sign if v == 167 else 1)
        for v in active
    )
    assert load.sum() == board.n
    return dict(
        sign=sign,
        retired=retire,
        stage_pair=original["target"],
        ring=list(ring),
        active_pendants=[] if retire else list(leaves),
        old_marks_used=marks,
        route_count_before_choice=layout()[-1],
        shutters=len(new_cuts),
        final_hatch=hatch,
        original_pair_parking=[parked[160], parked[167]],
        transpositions=log,
    )
