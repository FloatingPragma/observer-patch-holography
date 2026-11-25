"""Three equivalent physical arms between local source and receiving triangles."""

from functools import cache
from itertools import combinations, permutations, product
import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix

from ..winding_memory import measure as winding
from ..live_trimer.measure import LAW, blocks
from ..two_winding_loops.measure import stock


def terminal(board, owner):
    ports = set(range(12 * owner, 12 * owner + 12))
    gates = sorted(v for v in ports if board.degree[v] == 6)
    setups = []
    interior = ports - set(gates)
    for h in sorted(interior):
        candidates = []
        for g in gates:
            paths = []
            for a in sorted(
                set(board.neighbours[h]) & set(board.neighbours[g]) & (interior - {h})
            ):
                for c in sorted(
                    set(board.neighbours[a])
                    & set(board.neighbours[g])
                    & (interior - {h, a})
                ):
                    paths.append((int(a), int(c), g))
            candidates.append(paths)
        for trip in product(*candidates):
            if len(set(sum((list(z) for z in trip), []))) == 9:
                setups.append((h, {g: (a, c, g) for a, c, g in trip}))
    return setups


@cache
def fixture():
    board = winding.flight.federation(0)
    g = nx.Graph(board.seams)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    dist = nx.single_source_shortest_path_length(owners, 0)
    far = min(z for z, t in dist.items() if t == max(dist.values()))
    tours = next(
        paths
        for paths in combinations(sorted(nx.all_shortest_paths(owners, 0, far)), 3)
        if len(set(sum((p[1:-1] for p in paths), []))) == 3 * (dist[far] - 1)
    )
    gateway = {
        (a // 12, b // 12): a for x, y in board.external for a, b in [(x, y), (y, x)]
    }
    lefts, rights = terminal(board, 0), terminal(board, far)
    hub, left = lefts[0]
    tip, right = rights[0]
    routes = []
    for tour in tours:
        route = [gateway[tour[0], tour[1]]]
        for u, v, w in zip(tour, tour[1:], tour[2:]):
            path = min(
                nx.all_shortest_paths(
                    g.subgraph(range(12 * v, 12 * v + 12)), gateway[v, u], gateway[v, w]
                )
            )
            route += path
        routes.append(route + [gateway[tour[-1], tour[-2]]])
    length = max(map(len, routes))
    for path in routes:
        while len(path) < length:
            j, d = next(
                (j, min(q))
                for j, (a, b) in enumerate(zip(path, path[1:]))
                if a // 12 == b // 12
                and a // 12 not in (0, far)
                and (
                    q := (set(board.neighbours[a]) & set(board.neighbours[b]))
                    - set(path)
                )
            )
            path.insert(j + 1, int(d))
    legs = []
    blank = []
    marked = []
    for path in routes:
        a, c, g0 = left[path[0]]
        z, y, g1 = right[path[-1]]
        legs.append([a, c] + path + [y, z])
        blank += [
            (hub, a),
            (a, c),
            (c, g0),
            *zip(path, path[1:]),
            (g1, y),
            (y, z),
            (z, tip),
        ]
        marked.extend(
            [(a, g0), (z, g1)]
        )  # direction at local address a/z is initially enabled
    marked = marked[::2] + marked[1::2]
    ports = [hub, tip] + sum(legs, [])
    assert len(ports) == len(set(ports))
    edges = blank + marked
    assert all(g.has_edge(a, b) for a, b in edges)
    active = {tuple(sorted(e)) for e in edges}
    # During writes, terminals and all six marked edges can use blank routes,
    # including final guide edges except deposited seams.
    writing = g.copy()
    writing.remove_edges_from(marked)
    assert nx.is_biconnected(writing)
    ds = nx.multi_source_dijkstra_path_length(g, ports)
    donor, mate = next((p, q) for p, q in board.seams if min(ds[p], ds[q]) > 2)
    cuts = [
        e
        for e in board.seams
        if (set(e) & (set(ports) | {donor, mate}))
        and tuple(sorted(e)) not in active | {tuple(sorted((donor, mate)))}
    ]
    blanks = nx.Graph(blank)
    assert nx.is_connected(blanks)
    # Prepare the six direction heads, erase all private references and restore
    # the sole vacancy to the hub after each write. A blank-path transposition
    # restores the loads at every intermediate port, including a spectator.
    loads = np.ones(board.n, dtype=int)
    loads[hub], loads[donor] = 0, 2
    initial = loads.copy()

    def transpose(route, a, b):
        for u, v in winding.flight.closed_walk(nx.shortest_path(route, a, b)):
            loads[u], loads[v] = loads[v], loads[u]

    for a, b in marked:
        transpose(writing, hub, a)
        assert loads[a] == 0 and loads[b] == 1
        transpose(writing, a, hub)
    for port in ports:
        if board.degree[port] == 5:
            transpose(blanks, hub, port)
            assert loads[port] == 0
            transpose(blanks, port, hub)
    for ends in cuts:
        parked = []
        for src, alt in ((hub, min(blanks[hub])), (donor, mate)):
            if src in ends:
                assert alt not in ends
                loads[src], loads[alt] = loads[alt], loads[src]
                parked.append((src, alt))
        assert all(loads[k] == 1 for k in ends)
        for src, alt in parked:
            loads[src], loads[alt] = loads[alt], loads[src]
    np.testing.assert_array_equal(loads, initial)
    stock(False)
    return dict(
        ports=ports,
        legs=legs,
        marked=marked,
        blank=blank,
        hub=hub,
        tip=tip,
        tours=tours,
        candidate_terminal_counts=[len(lefts), len(rights)],
        donor=[donor, mate],
        private_writes=sum(int(board.degree[p] == 5) for p in ports),
        exterior_writes=len(cuts),
    )


def native_graph():
    f = fixture()
    ports = f["ports"]
    # Heads are addressed in route-relative direction, accounting for sorted
    # physical edge orientation. A bit toggles exactly when the charge crosses
    # that head. A coherent blank detour can leave all other heads unchanged.
    n = len(ports)
    size = 64 * n
    diag = np.zeros(size)
    r = []
    c = []
    v = []
    for bit, inds in [(False, f["blank"]), (True, f["marked"])]:
        for j, (a, b) in enumerate(inds):
            s, t = ports.index(a), ports.index(b)
            raw = blocks()[(1, 1) if bit else "blank"] / 8
            filled = LAW.SRC.index((1, 1))
            host = raw[filled::9, filled::9]
            np.testing.assert_array_equal(host, host[0, 0] * np.eye(len(host)))
            for h in range(64):
                coin = (((h >> j) & 1) ^ int(a > b)) if bit else 0
                for src, dst, forward in [(s, t, True), (t, s, False)]:
                    pair = (0, 1) if forward == (a < b) else (1, 0)
                    loc = 9 * coin + LAW.SRC.index(pair)
                    col = 64 * src + h
                    diag[col] += raw[loc, loc] - host[0, 0]
                    for o in np.flatnonzero(raw[:, loc]):
                        if o == loc:
                            continue
                        assert LAW.SRC[o % 9] == pair[::-1]
                        assert int(o) // 9 == (1 - coin if bit else coin)
                        h1 = h ^ (1 << j) if bit else h
                        r.append(64 * dst + h1)
                        c.append(col)
                        v.append(raw[o, loc])
    h = coo_matrix(
        (v + diag.tolist(), (r + list(range(size)), c + list(range(size)))),
        shape=(size, size),
    ).tocsr()
    assert (h - h.T).nnz == 0
    return h


def permutations_of_chart():
    f = fixture()
    out = {}
    for perm in permutations(range(3)):
        pm = {f["hub"]: f["hub"], f["tip"]: f["tip"]}
        for i, j in enumerate(perm):
            pm.update(zip(f["legs"][i], f["legs"][j]))
        assert {tuple(sorted((pm[a], pm[b]))) for a, b in f["blank"]} == {
            tuple(sorted(z)) for z in f["blank"]
        }
        assert {(pm[a], pm[b]) for a, b in f["marked"]} == set(f["marked"])
        head = [
            sum(
                ((h >> (3 * k + i)) & 1) << (3 * k + perm[i])
                for k in (0, 1)
                for i in range(3)
            )
            for h in range(64)
        ]
        out[perm] = np.array(
            [
                64 * f["ports"].index(pm[v]) + head[h]
                for v in f["ports"]
                for h in range(64)
            ]
        )
    return out
