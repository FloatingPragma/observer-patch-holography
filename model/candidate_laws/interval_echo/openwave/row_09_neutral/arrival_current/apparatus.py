"""Physical capture hatches, collar workstation and balanced receiver stock."""

from functools import cache
from itertools import permutations
import networkx as nx
import numpy as np
from model.executable_model import federation
from model.experiments.interval_echo_consensus.signed_effects.collar_storage import (
    conserved_columns,
    gates,
    MEMORIES,
    LAW,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.apparatus import (
    fixture,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.echo_turn.material import (
    chart,
    transp,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.geometry import (
    layout,
)

SENSORS = {"opposite": (207, 206), "second_owner": (106, 96)}


@cache
def certify():
    h, prep, echo = gates(1.0)
    conserved_columns()
    unit = LAW.SRC.index((1, 1))
    # The pre-existing shutters were closed by a signed zero echo on two units.
    shut = echo[:, unit]
    np.testing.assert_allclose(echo @ shut, np.eye(LAW.DIM)[:, unit], atol=3e-14)
    # Sorted endpoints (2,0), calibrated even magnitude-2 write.
    written = echo @ prep[:, LAW.SRC.index((2, 0))]
    assert (
        abs(np.vdot(np.kron(MEMORIES[:, 0], np.eye(9)[:, unit]), written)) ** 2
        > 1 - 2e-13
    )
    outcomes = {}
    for word in LAW.SRC:
        col = LAW.SRC.index(
            word[::-1]
        )  # Both installed receiver charts list larger address first.
        p = echo[:, col]
        d = word[1] - word[0]
        take = np.repeat([False] + [v == d for _, v in LAW.LABELS], 9)
        np.testing.assert_allclose(
            echo @ (take * p), np.eye(LAW.DIM)[:, col], atol=3e-14
        )
        if abs(d) == 2:
            assert sum(abs(p[unit::9]) ** 2) > 1 - 2e-13
        outcomes[word] = d
    return outcomes


@cache
def station(sensor):
    u, v = SENSORS[sensor]
    b = federation(0)
    board = nx.Graph(b.seams)
    ring, marks, leaves, _, _ = layout()
    active = set(ring) | set(leaves)
    assert (
        u > v
        and u // 12 == v // 12
        and u in ring
        and v in ring
        and board.has_edge(u, v)
    )
    prior = fixture()
    lab = chart()
    # Never cross an old guide shutter, a parked source load, or an active flight seam.
    excluded = (
        active
        | set(prior["ports"])
        | set(prior["donor"])
        | set(lab["storage"] + lab["probes"] + lab["stock"] + [lab["reference"]])
    )
    # Previously parked loads reside here in the frozen release commissioning log.
    # Exclude the whole original owner cells of those parks.
    excluded.update(range(72, 96))
    ext = board.subgraph(set(board) - excluded).copy()
    local = sorted(n for n in ext if n // 12 == u // 12)
    # This blank crossing links the two receiving owners. Only its exterior
    # service hatches need zero-head unwrites while all ring loads are units.
    service = [215, 206, 106, 103] if sensor == "opposite" else []
    if service:
        assert all(board.has_edge(i, j) for i, j in zip(service, service[1:]))
        assert service[1] in active and service[2] in active
        assert all(i in ext for i in (service[0], service[-1]))
        ext.add_edges_from(zip(service, service[1:]))
        assert board.has_edge(208, 206) and 208 in ext
        ext.add_edge(208, 206)
    for x, a, g, d in permutations(local, 4):
        if (
            max(x, g) >= a
            or b.degree[a] != 5
            or not all(
                board.has_edge(i, j)
                for i, j in [(u, x), (v, a), (x, a), (g, a), (x, g), (a, d)]
            )
        ):
            continue
        occupied = {x, a, g, d}
        for flag in sorted(
            tuple(sorted(e)) for e in ext.subgraph(set(local) - occupied).edges
        ):
            # A physically independent read on each output rail against zero.
            for read in sorted(
                (
                    tuple(sorted(e))
                    for e in ext.edges
                    if set(e).isdisjoint(occupied | set(flag))
                    and set(e).isdisjoint(active)
                ),
                key=lambda edge: (
                    min(nx.shortest_path_length(board, x, z) for z in edge),
                    edge,
                ),
            ):
                graph = ext.copy()
                graph.remove_edges_from([(x, a), (g, a)])
                routable = graph.subgraph(set(graph) - {read[1]})
                if not all(
                    nx.has_path(routable, z, read[0]) for z in [x, a, g, d, *flag]
                ):
                    continue
                far = nx.multi_source_dijkstra_path_length(board, active)
                slots = sorted(
                    (
                        i
                        for i in nx.node_connected_component(routable, x)
                        if far[i] > 2 and i not in occupied | set(flag) | set(read)
                    ),
                    key=lambda i: (nx.shortest_path_length(routable, x, i), i),
                )[:12]
                if len(slots) == 12:
                    return (u, v), (x, a, g, d), flag, read, slots, ext
    raise ValueError("No local contact receiver")


def install(sensor, guard, donor):
    (u, v), (x, a, g, d), flag, read, slots, graph = station(sensor)
    graph = graph.copy()
    b = federation(0)
    load = np.ones(b.n, dtype=int)
    load[slots] = [2, 0] * 6
    assert load.sum() == b.n
    log = []

    def move(i, j, avoid=()):
        path = nx.shortest_path(graph.subgraph(set(graph) - set(avoid)), i, j)
        before = load.copy()
        for z, w in transp(path):
            load[z], load[w] = load[w], load[z]
        before[i], before[j] = before[j], before[i]
        np.testing.assert_array_equal(load, before)
        log.append(path)

    for dest, subset in ((x, slots[:2]), (g, slots[2:4])):
        move(subset[0], dest)
        move(subset[1], a)
        assert (load[dest], load[a]) == (2, 0) and dest < a
        load[dest] = load[a] = (
            1  # certified native even collar plus local private pulse
        )
        graph.remove_edge(dest, a)
    move(slots[4 if guard == 2 else 5], g)
    move(slots[6], flag[0])
    move(slots[7], flag[1])
    move(slots[9], read[1])
    if donor != 1:
        move(slots[10 if donor == 2 else 11], d)
    assert (load[x], load[a], load[g], load[d]) == (1, 1, guard, donor)
    assert list(load[list(flag)]) == [2, 0] and tuple(load[list(read)]) == (1, 0)
    assert load.sum() == b.n
    # Both hatches join the active ring to the excluded exterior. Their old
    # zero heads open on two units while the pair is held in the upstream echo.
    for z, w in [(u, x), (v, a)]:
        assert load[z] == load[w] == 1
    # An independent ternary read on each post-gate rail or retained flag can
    # route through other rails: all internal values are restored by transp.
    for z in [x, a, d, *flag]:
        move(z, read[0], (read[1],))
        move(z, read[0], (read[1],))
    return dict(
        receiving=[u, v],
        service_edges=[[215, 206], [206, 106], [106, 103], [206, 208]]
        if sensor == "opposite"
        else [],
        workstation=[x, a, g, d],
        hatches=[[u, x], [v, a]],
        persistent_flag=list(flag),
        charge_reader=list(read),
        balanced_stock=slots,
        guard=guard,
        donor=donor,
        stock_routes=log[:-10],
        read_routes=log[-10::2],
        mounted_nonunit=[[int(i), int(load[i])] for i in np.flatnonzero(load != 1)],
    )


@cache
def standards(sensor):
    # On a reset calibration board the receiver opens service hatches on unit
    # stock and accepts an independently stocked extreme pair in either order.
    (u, v), (x, a, _, _), _, _, slots, graph = station(sensor)
    graph = graph.copy()
    for i, j in [(u, x), (v, a)]:
        graph.add_edge(i, j)
    paths = [nx.shortest_path(graph, slots[0], u), nx.shortest_path(graph, slots[1], v)]
    words = []
    for order in ((0, 2), (2, 0)):
        load = np.ones(federation(0).n, dtype=int)
        load[slots[:2]] = order
        for path in paths:
            for i, j in transp(path):
                load[i], load[j] = load[j], load[i]
        assert tuple(load[[u, v]]) == order and load.sum() == len(load)
        words.append([int(load[u]), int(load[v]), certify()[order]])
    return dict(routes=paths, stocked_words_and_signed_head=words)
