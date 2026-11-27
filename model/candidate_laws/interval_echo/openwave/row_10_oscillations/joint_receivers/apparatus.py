"""Material source read and head-preserving paths on the held federation."""

from functools import cache
from itertools import pairwise
from unittest.mock import patch

import networkx as nx
import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin import (
    apparatus,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.echo_turn import (
    material,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.material_interferometer.apparatus import (
    PHASE_PULSES,
    Terminal,
    inverse,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder import (
    apparatus as stocks,
    measure as conversion,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    detectors,
    geometry,
)


def recipe(sense, omit=False):
    lab = Terminal("tip")
    split = [
        (lab.hub, arm[0], dwell)
        for arm, dwell in zip(lab.arms[:2], (np.pi / 2, np.pi), strict=True)
    ]
    route = split.copy()
    for i, (a, c, g) in enumerate(lab.arms[:2]):
        if i:
            blanks = ((a, c), (c, g))
            phase = [(*blanks[j % 2], t) for j, t in enumerate(PHASE_PULSES)]
            route.extend(phase if sense == 1 else inverse(phase))
        if not (i and omit):
            route.append((a, g, np.pi / 8))
        route.extend(((g, c, np.pi), (c, a, np.pi)))
    return route + inverse(split)


@cache
def pointer(sense, omit=False):
    """Effects for the sharp hub return followed by both vacancy flag tests."""
    lab = Terminal("tip")
    effects = []
    for spectator in (0, 4):
        columns = []
        for m in range(4):
            returned = lab.run(inverse(recipe(sense, omit)), spectator | m)[0]
            columns.append(returned)
        hub = np.column_stack(columns)
        hit = np.outer(hub[spectator].conj(), hub[spectator])
        hub_effect = hub.conj().T @ hub
        effects.append(np.array([hit, hub_effect - hit, np.eye(4) - hub_effect]))
        # The vacancy at the hub tests each active flag, without reading bit 2.
        for m in range(8):
            for i, (a, _, g) in enumerate(lab.arms[:2]):
                go = [(lab.hub, a, np.pi), (a, g, np.pi / 8)]
                flag = lab.run(go, m)
                passed = not ((m >> i) & 1)
                np.testing.assert_allclose(
                    np.sum(abs(flag[lab.ports.index(g)]) ** 2), int(passed), atol=2e-13
                )
                for u, v, dt in inverse(go):
                    lab.gate(flag, u, v, dt)
                np.testing.assert_allclose(flag[0, m], 1, atol=3e-13)
    np.testing.assert_allclose(effects[0], effects[1], atol=3e-13)
    np.testing.assert_allclose(effects[0].sum(axis=0), np.eye(4), atol=2e-13)
    return effects[0]


@cache
def access():
    """Audit actual stock transpositions and open a route for an existing loan."""
    lab = Terminal("tip")
    guide = apparatus.fixture()
    board = federation(0)
    ring, marks, leaves, _, _ = geometry.layout()
    active = set(ring) | set(leaves)
    original_marks = set(map(frozenset, guide["marked"]))
    assert tuple((a, g) for a, _, g in lab.arms[:2]) == marks
    walks = []
    transp = material.transp

    def audited(path):
        steps = transp(path)
        assert all(frozenset(edge) not in original_marks for edge in steps)
        walks.append(steps)
        return steps

    # Both compiler modules reference the blank transposition. Trace their
    # actual steps, including stock before formation of the neutral pair.
    with (
        patch.object(material, "transp", audited),
        patch.object(stocks, "transp", audited),
    ):
        commissioning = geometry.configure(-1)
    raw = commissioning["transpositions"]
    far = (
        material.chart()["storage"]
        + material.chart()["probes"]
        + material.chart()["stock"]
    )
    assert [row[0] for row in raw[:5]] == [160, 167, far[5], far[7], far[8]]
    loan = raw[2][1]
    # After the two full conversions the restored collars give two distant
    # pairs of sharp storage stocks. The returned 160,167 stocks are remote.
    assert [r[0] for r in raw[5:9]] == [160, 167, 158, 167]
    assert [r[1] for r in raw[-2:]] == [160, 167]
    private_resets = raw[9:-2]
    assert private_resets[::2] == private_resets[1::2]
    blocked = {row[1] for row in raw[3:9]} | {guide["donor"][0]}
    graph = nx.Graph(board.seams)
    graph.remove_edges_from(tuple(sorted(e)) for e in original_marks)
    # Avoid both distal owners in their entirety, including their counters.
    graph.remove_nodes_from(blocked | {v for v in graph if v // 12 in (12, 17)})
    # At the first hop the loan is nonunit. Demand an originally ordinary
    # stock edge which no old guide shutter or final ring shutter touches.
    sites = active | set(guide["ports"]) | set(guide["donor"])
    assert loan not in sites
    graph.remove_edges_from((loan, v) for v in list(graph[loan]) if v in sites)
    path = min(nx.all_shortest_paths(graph, loan, lab.hub))
    # Every remaining endpoint on the route is unit when the complete far
    # pair is seen. Open its known unit shutters before moving the vacancy.
    route = transp(path)
    load = np.ones(board.n, int)
    load[loan] = 0
    for a, b in route:
        load[a], load[b] = load[b], load[a]
    assert load[lab.hub] == 0 and load[loan] == 1 and np.count_nonzero(load != 1) == 1
    assert set(sum(lab.arms[:2], []) + [lab.hub]).isdisjoint(blocked)
    assert all(frozenset((u, v)) not in original_marks for u, v in route)
    detectors.verify()  # Native blank swap is a permutation on all nine stocks.
    # The converter acts on hub/target/donor and two collars distinct from
    # the retained directions. Its selected stock branch has a common phase.
    gate = conversion.gate(np.pi, 2)
    col = gate[:, conversion.WORDS.index((0, 1, 2))]
    selected = conversion.WORDS.index((0, 2, 1))
    np.testing.assert_allclose(abs(col[selected]), 1, atol=3e-13)
    np.testing.assert_allclose(np.delete(col, selected), 0, atol=6e-13)
    assert all(
        frozenset(e) not in original_marks for e in ((160, 167), (158, 167), (166, 167))
    )
    assert walks and all(
        all(board.degree[v] >= 5 for v in e) for walk in walks for e in walk
    )
    return dict(
        loan=loan,
        loan_route=path,
        source=lab.hub,
        arms=lab.arms[:2],
        fixed_stocks=sorted(blocked),
        audited_blank_edges=len({frozenset(e) for path in walks for e in path}),
        conversion_amplitude=[float(col[selected].real), float(col[selected].imag)],
    )
