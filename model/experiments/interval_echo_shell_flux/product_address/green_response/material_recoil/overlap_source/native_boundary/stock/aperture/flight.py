"""One-load propagation bounds and a full native source/probe exercise."""

import argparse
from fractions import Fraction
import json
from math import factorial, isqrt
import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm

from model.executable_model import Carrier, CELL_SEAMS
from . import measure as source
from ...free_probe import measure as flight

HERE = source.HERE


def propagation_bound(board, path, start, times, order=96):
    # Integer matrix |2(H_blank-51)|, on the six-neighbor-completed board.
    # All original entries had the same (negative) sign. Powers have one sign.
    near = set(path)
    v = [0] * board.n
    v[start] = 1
    norms = []
    for k in range(order + 1):
        n = sum(v[i] ** 2 for i in near)
        a = isqrt(n)
        norms.append(a + int(a * a < n))
        v = [
            (6 - board.degree[i].item()) * x + sum(v[j] for j in board.neighbours[i])
            for i, x in enumerate(v)
        ]
    first = next((i for i, x in enumerate(norms) if x), None)
    assert first == min(
        nx.shortest_path_length(nx.Graph(board.seams), start, p) for p in path
    )
    bounds = []
    for t in times:
        time = Fraction(t)
        # Use integer upper square roots on all computed coefficients and the
        # norm-3 unitary Taylor integral remainder on the rest.
        integral = sum(
            Fraction(m, factorial(k + 1)) * time * (time / 2) ** k
            for k, m in enumerate(norms)
        )
        rest = time * (3 * time) ** (order + 1) / factorial(order + 2)
        refined = min(Fraction(1), 36 * (integral + rest))
        elementary = min(
            Fraction(1), 36 * time * (3 * time) ** first / factorial(first + 1)
        )
        bounds.append(
            dict(
                time=t,
                uniform_probability_bound=float(min(refined, elementary)),
                integer_bound=str(min(refined, elementary)),
                distance_only_bound=float(elementary),
            )
        )
    return dict(port=start, distance=first, taylor_degree=order, bounds=bounds)


def cell_operators():
    board = Carrier("one cell", 12, CELL_SEAMS)
    graph = nx.Graph(board.seams)
    dist = nx.single_source_shortest_path_length(graph, 0)
    end = min(i for i, k in dist.items() if k == max(dist.values()))
    path = nx.shortest_path(nx.bfs_tree(graph, 0, sort_neighbors=sorted), 0, end)
    edges = [tuple(sorted(pair)) for pair in zip(path, path[1:])]
    n = len(edges)
    eye = sparse.eye(4**n, format="csr")

    def embed(matrix, j):
        return sparse.kron(
            sparse.eye(4**j),
            sparse.kron(matrix, sparse.eye(4 ** (n - 1 - j))),
            format="csr",
        )

    def block(a, b, m):
        return sparse.kron(
            sparse.csr_matrix(([1], ([a], [b])), shape=(board.n, board.n)),
            m,
            format="csr",
        )

    h0 = 48 * sparse.eye(board.n, format="csr") + board.laplacian() / 2
    ref = sparse.kron(h0, eye, format="csr")
    h = sparse.kron(
        sparse.diags(8 * (6 - board.degree), dtype=float), eye, format="csr"
    )
    local, _ = flight.restricted(8)
    tags = [(c, d) for c, d in flight.LAW.LABELS if abs(d) == 1]
    factors = []
    costs = []
    for a, b in board.seams:
        if (a, b) not in edges:
            for p, q, v in ((a, a, 8.5), (b, b, 8.5), (a, b, -0.5), (b, a, -0.5)):
                h += v * block(p, q, eye)
            continue
        j = edges.index((a, b))
        g = embed(local[0, :, 0, :], j)
        elsewhere = sparse.eye(board.n, format="csr")
        elsewhere = elsewhere - sparse.csr_matrix(
            ([1, 1], ([a, b], [a, b])), shape=(board.n, board.n)
        )
        term = sparse.kron(elsewhere, g, format="csr")
        cost = 0 * ref
        for p, word in ((a, 2), (b, 1)):
            for q, other in ((a, 2), (b, 1)):
                elem = embed(local[word, :, other, :], j)
                term += block(p, q, elem)
                cost += block(
                    p,
                    q,
                    (g if p == q else 0 * g) + (8.5 if p == q else -0.5) * eye - elem,
                )
        h += term
        ref += sparse.kron(sparse.eye(board.n), g, format="csr")
        factors.append(j)
        costs.append((j, cost))
    costs = [m for _, m in sorted(costs, key=lambda p: p[0])]
    assert (ref - h - sum(costs)).nnz == 0
    assert (h != h.T).nnz == 0
    assert len(factors) == n
    memory = np.ones(1)
    for u, v in zip(path, path[1:]):
        # Success coin and written input difference from the frozen frontier writer.
        label = (0, 1) if u < v else (1, -1)
        memory = np.kron(memory, np.eye(4)[tags.index(label)])
    return board, path, h, h0, local, tags, memory, costs


def small_flight():
    board, path, h, h0, local, tags, memory, costs = cell_operators()
    graph = nx.Graph(board.seams)
    n = len(path) - 1
    remote = nx.multi_source_dijkstra_path_length(graph, path, weight=None)
    far = min(i for i, k in remote.items() if k == max(remote.values()))
    times = (0, 1, 2, 4, 8)
    rows = []
    for start in (far, path[-1]):
        for age in (0, 1):
            # Vacuum evolution before launch factors across actual stored words.
            held = memory.astype(complex)
            for j in range(n):
                held = np.moveaxis(held.reshape([4] * n), j, 0).reshape(4, -1)
                held = expm(-1j * age * local[0, :, 0, :]) @ held
                held = np.moveaxis(held.reshape([4] * n), 0, j).ravel()
            init = np.kron(np.eye(board.n)[start], held)
            for time in times:
                state = expm_multiply(
                    -1j * time * (h - 72 * sparse.eye(h.shape[0])), init
                )
                standard = expm_multiply(
                    -1j * time * (h0 - 51 * sparse.eye(board.n)), np.eye(board.n)[start]
                )
                wave = state.reshape(board.n, -1)
                ordinary = wave @ wave.conj().T
                difference = ordinary - np.outer(standard, standard.conj())
                tv = np.abs(np.linalg.eigvalsh(difference)).sum() / 2
                loss = np.array([np.vdot(state, c @ state).real for c in costs])
                assert (loss >= -1e-10).all()
                if start == path[-1]:
                    # This one-cell formation ends on a degree-five port like the
                    # larger trails. Aging the empty stock before re-inserting the
                    # terminal load can change its formation energy.
                    if age == 0:
                        assert loss.sum() <= 16.25 + 1e-10
                rows.append(
                    dict(
                        port=start,
                        age=age,
                        time=time,
                        trace_distance=float(tv),
                        maximum_port_shift=float(max(abs(np.diag(difference)))),
                        ordinary_port_contrast=np.diag(difference).real.tolist(),
                        cumulative_excess=np.cumsum(8.5 - loss).tolist(),
                    )
                )
    return dict(path=path, far=far, dimension=h.shape[0], flights=rows)


def run():
    data = json.loads((source.trail.HERE / "trail.json").read_text())
    full = []
    for level, row in enumerate(data["boards"]):
        board = source.trail.measure.observations.federation(level)
        distances = nx.multi_source_dijkstra_path_length(
            nx.Graph(board.seams), row["path"], weight=None
        )
        far = min(i for i, k in distances.items() if k == max(distances.values()))
        full.append(
            dict(
                owners=board.cells,
                far=propagation_bound(board, row["path"], far, (1, 4, 8)),
                on_source=propagation_bound(
                    board, row["path"], row["path"][0], (1, 4, 8)
                ),
            )
        )
    return dict(federations=full, one_cell=small_flight())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    path = HERE / "flight.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        source.check(json.loads(path.read_text()), result)
