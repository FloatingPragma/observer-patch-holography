"""Find six-port dipole guides disjoint from each regional input."""

from itertools import combinations
from hashlib import sha256
import json

import networkx as nx
from model.executable_model.carrier import federation


def arms(board, gateway):
    owner = gateway // 12
    for b in sorted(board.neighbours[gateway]):
        if b // 12 == owner and board.degree[b] == 5:
            for a in sorted(board.neighbours[b]):
                if a // 12 == owner and board.degree[a] == 5:
                    yield a, b


def guides():
    board = federation(0)
    net = nx.Graph(board.seams)
    candidates = []
    for u, v in sorted(board.external):
        start = next(arms(board, u))
        end = tuple(reversed(next(arms(board, v))))
        path = (*start, u, v, *end)
        blocked = set(path) | set(range(v // 12 * 12, v // 12 * 12 + 12))
        # The written operand visits a private head contact; the active guide and
        # receiving owner never serve as transit registers for the source write.
        graph = net.subgraph(set(net) - blocked)
        j = min(
            i for i in board.neighbours[start[1]] if i in graph and board.degree[i] == 5
        )
        routes = nx.single_source_shortest_path(graph, j)
        candidates.append((set((u // 12, v // 12)), path, j, routes))
    return candidates


def fixtures(profiles):
    candidates = guides()
    result = []
    for profile in profiles:
        owners = set(profile["owners"])
        ports = [12 * owner + p for owner in sorted(owners) for p in range(12)]
        opts = [
            c for c in candidates if not owners & c[0] and all(i in c[3] for i in ports)
        ]
        chosen = next(
            t for t in combinations(opts, 3) if len(set.union(*(x[0] for x in t))) == 6
        )
        sites = []
        for _, guide, contact, table in chosen:
            routes = [table[p][::-1] for p in ports]
            sites.append(
                dict(
                    guide=list(guide),
                    contact=contact,
                    shortest_write_route_sha256=sha256(
                        json.dumps(routes, separators=(",", ":")).encode()
                    ).hexdigest(),
                    maximum_write_hops=max(len(r) - 1 for r in routes),
                )
            )
        result.append(
            dict(
                owners=sorted(owners),
                center=profile["center"],
                radius=profile["radius"],
                sites=sites,
            )
        )
    return result
