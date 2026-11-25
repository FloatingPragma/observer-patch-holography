"""Equal, physical detours in the three blank arms, with the original terminals."""

from itertools import pairwise
import networkx as nx
import numpy as np
from scipy.sparse import block_diag, coo_matrix, csr_matrix

from .. import apparatus
from ..measure import reduce


def detour(count):
    f = apparatus.fixture()
    board = apparatus.winding.flight.federation(0)
    old = f["ports"]
    used = set(old) | set(f["donor"])
    replace = {}
    for j in range(3, len(f["legs"][0]) - 4):
        candidates = [
            sorted(
                (set(board.neighbours[p[j]]) & set(board.neighbours[p[j + 1]])) - used
            )
            for p in f["legs"]
        ]
        if not all(candidates):
            continue
        inserts = [int(group[0]) for group in candidates]
        if len(set(inserts)) != 3:
            continue
        if len(replace) == count:
            break
        used.update(inserts)
        replace[j] = inserts
    assert len(replace) == count
    port_ids = old + [p for trip in replace.values() for p in trip]
    port_index = {p: i for i, p in enumerate(port_ids)}
    legs = [
        [
            site
            for j, p in enumerate(leg)
            for site in ([p, replace[j][a]] if j in replace else [p])
        ]
        for a, leg in enumerate(f["legs"])
    ]
    blanks = [(a, b) for leg in legs for a, b in pairwise([f["hub"], *leg, f["tip"]])]
    assert nx.is_connected(nx.Graph(blanks))
    edges = {tuple(sorted(edge)) for edge in blanks + f["marked"]}
    assert edges <= set(map(tuple, map(sorted, board.seams)))
    old_blank = {tuple(sorted(e)) for e in f["blank"]}
    new_blank = {tuple(sorted(e)) for e in blanks}
    delta = [(1, e) for e in new_blank - old_blank] + [
        (-1, e) for e in old_blank - new_blank
    ]
    # Full 99-column reduction supplies the blank one-vacancy term on its two rails.
    raw = apparatus.blocks()["blank"] / 8
    rails = [apparatus.LAW.SRC.index(w) for w in ((0, 1), (1, 0))]
    term = raw[np.ix_(rails, rails)] - np.eye(2) * raw[4, 4]
    np.testing.assert_array_equal(term, [[8.5, -0.5], [-0.5, 8.5]])
    base = apparatus.native_graph()
    size = 64 * len(port_ids)
    h = block_diag([base, csr_matrix((size - base.shape[0],) * 2)], format="csr")
    rows, cols, vals = [], [], []
    for sign, (a, b) in delta:
        inds = [64 * port_index[p] + np.arange(64) for p in (a, b)]
        for i, j in np.ndindex(2, 2):
            rows.extend(inds[i])
            cols.extend(inds[j])
            vals.extend(np.repeat(sign * term[i, j], 64))
    h += coo_matrix((vals, (rows, cols)), shape=h.shape).tocsr()
    original_actions = apparatus.permutations_of_chart()
    actions = {}
    for perm, p in original_actions.items():
        extras = [
            64 * port_index[trip[perm[i]]] + p[:64]
            for trip in replace.values()
            for i in range(3)
        ]
        actions[perm] = np.concatenate([p, *extras])
        action = actions[perm]
        assert (h[action, :][:, action] != h).nnz == 0
    # On this material the original hub vacancy can erase the private rail of
    # every added degree-five site before all unused exterior seams are shut.
    # Both stock loads can be parked off each shutter during its (1,1) write.
    blank_graph = nx.Graph(blanks)
    assert not set(f["donor"]) & set(port_ids)
    exterior = [
        e
        for e in board.seams
        if set(e) & (set(port_ids) | set(f["donor"]))
        and tuple(sorted(e)) not in edges | {tuple(sorted(f["donor"]))}
    ]
    for e in exterior:
        if f["hub"] in e:
            assert set(blank_graph[f["hub"]]) - set(e)
        assert not set(f["donor"]) <= set(e)
    lift, small = reduce(h, actions)
    return replace, port_ids, legs, h, actions, lift, small
