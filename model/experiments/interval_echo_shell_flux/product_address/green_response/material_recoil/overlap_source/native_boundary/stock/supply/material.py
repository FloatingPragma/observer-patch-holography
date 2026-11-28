"""Allocate a uniform ordinary resource and capture fixture on the three native boards."""

from itertools import combinations
from functools import lru_cache

import networkx as nx
import numpy as np
from scipy import sparse

from ..measure import observations


@lru_cache(None)
def layouts(local_edges):
    cell = nx.Graph(local_edges)

    def wings(remaining, used):
        if not remaining:
            return ()
        start, *tail = remaining
        for mid in sorted(set(cell[start]) - used):
            for end in sorted(set(cell[mid]) - used - {mid}):
                result = wings(tail, used | {mid, end})
                if result is not None:
                    return ((start, mid, end), *result)
        return None

    out = []
    for core in combinations(range(12), 4):
        interior = cell.subgraph(core)
        if not nx.is_connected(interior):
            continue
        pairs = sorted(tuple(sorted(e)) for e in interior.edges)
        matches = [p for p in combinations(pairs, 2) if len(set(sum(p, ()))) == 4]
        if matches and (parking := wings(core, set(core))) is not None:
            out.append((core, min(matches), parking))
    assert out
    return out


def fixture(level):
    board = observations.federation(level)
    external = np.array(sorted(board.external))
    gateways = [set() for _ in range(board.cells)]
    for v in external.ravel():
        gateways[v // 12].add(v % 12)
    local_edges = [[] for _ in gateways]
    for a, b in board.seams:
        if a // 12 == b // 12:
            local_edges[a // 12].append((a % 12, b % 12))
    donors, wings, ports = [], [], []
    for owner, (edges, gates) in enumerate(zip(local_edges, gateways, strict=True)):
        choices = layouts(tuple(sorted(edges)))
        core, pairs, parking = min(choices, key=lambda c: (-len(set(c[0]) & gates), c))
        offset = 12 * owner
        donors.extend([(a + offset, b + offset) for a, b in pairs])
        wings.extend([(a + offset, b + offset, c + offset) for a, b, c in parking])
        ports.extend(v + offset for v in core)
    pindex = {v: i for i, v in enumerate(ports)}
    active = [(a, b) for a, b in board.seams if a in pindex and b in pindex]
    donor_set = set(donors)
    others = [(pindex[a], pindex[b]) for a, b in active if (a, b) not in donor_set]
    pairs = np.array([(pindex[a], pindex[b]) for a, b in donors])
    core = nx.Graph(active)
    owner_graph = nx.Graph((a // 12, b // 12) for a, b in external)
    distance = nx.single_source_shortest_path_length(owner_graph, ports[0] // 12)
    assert len(core) == 4 * board.cells and len(wings) == len(core)
    assert len(set(np.asarray(wings).ravel())) == board.n
    adjacency = nx.to_scipy_sparse_array(
        core, nodelist=ports, format="csr", dtype=float
    )
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    h = sparse.diags(8.5 * degree, format="csr") - adjacency / 2
    edgeset = set(board.seams)
    assert all(
        tuple(sorted(edge)) in edgeset for wing in wings for edge in zip(wing, wing[1:])
    )
    # Each donor can time an ordinary endpoint phase with an existing,
    # initially empty private-reference wing port after inverse capture.
    wing_at = {a: (b, c) for a, b, c in wings}
    meters = sum(
        any(board.degree[v] == 5 for endpoint in pair for v in wing_at[endpoint])
        for pair in donors
    )
    assert meters == len(donors)
    audit = dict(
        private_reference_analyzers=meters,
        owners=board.cells,
        core_ports=len(core),
        core_links=len(active),
        native_links=len(board.seams),
        donor_links=len(pairs),
        parking_links=2 * len(wings),
        private_dark=int(sum(6 - board.degree)),
        full_native_degree_max=int(max(board.degree)),
        external_retained=sum(a // 12 != b // 12 for a, b in active),
        external_total=len(external),
        components=nx.number_connected_components(core),
        reachable=len(nx.node_connected_component(core, ports[0])),
        root_port=ports[0],
        root_core_degree=int(degree[0]),
    )
    return (
        h,
        pairs,
        np.array(others),
        np.array(donors),
        np.array([distance[a // 12] for a, b in donors]),
        audit,
    )
