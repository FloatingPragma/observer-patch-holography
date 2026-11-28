"""A paid native write-and-read service with passive, spatially persistent deposited energy."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import networkx as nx
import numpy as np

from . import measure

HERE = measure.HERE
LOCK = "b90cb41f5a0bcc401f914265fd859f72d512972266934b065e2de9474d231cf5"


def instrument(forward, W, V, H):
    pair = (1, 0) if forward else (0, 1)
    d = pair[0] - pair[1]
    initial = measure.LAW.SRC.index(pair)

    def written(c):
        return 9 * (1 + measure.LAW.LABELS.index((c, d))) + measure.LAW.SRC.index(
            measure.LAW.native_event(*pair, c)
        )

    success = written(0 if forward else 1)
    retry = written(1 if forward else 0)
    basis = [initial, retry, success]
    small = W[np.ix_(basis, basis)]
    np.testing.assert_array_equal(small, [[0, 32, 32], [32, 16, 16], [32, 16, 16]])
    assert not W[np.ix_(basis, np.setdiff1d(np.arange(99), basis))].any()
    # For the half-exchange control the arriving load still has a blank head.
    assert not V[initial, 9:].any()
    assert V[initial, measure.LAW.SRC.index(pair[::-1])] == 32
    before = [F(int(H[x, x]), 4) for x in basis]
    # H=U G U*, so the coherent pre-read exposure U|x> has the G diagonal.
    g0, g1 = measure.LAW.generator_parts()
    g = g0 + 8 * g1
    work_read = []
    for col, x in enumerate(basis[:2]):
        target = sum(F(int(small[col, j]), 64) * before[j] for j in range(3))
        work_read.append(target - F(float(g[x, x])))
    visits = [F(3, 2), F(1)]
    assert [
        sum(visits[j] * F(int(small[j, k]), 64) for j in range(2)) + int(k == 0)
        for k in range(2)
    ] == visits
    assert sum(visits[j] * F(int(small[j, 2]), 64) for j in range(2)) == 1
    return dict(
        basis=basis,
        kernel_times64=small.tolist(),
        visits=visits,
        mean_echoes=sum(visits),
        one_echo_success=F(1, 2),
        eventual_success=1,
        read_work_per_visit=work_read,
        mean_read_work=sum(x * y for x, y in zip(visits, work_read)),
        geometric_failure_bound_per_attempt=F(3, 4),
    ), success // 9


def energy(board, port, mem, H):
    load = np.zeros(board.n, dtype=int)
    load[port] = 1
    edges = np.array(board.seams)
    a, b = edges.T
    E = np.diag(H)[9 * mem + 3 * load[a] + load[b]]
    degree = 32 * int((6 - board.degree) @ load)
    return F(int(E.sum()) + degree, 4), E


def run():
    assert sha256((HERE / "TRAIL.md").read_bytes()).hexdigest() == LOCK
    W, V, H, h0 = measure.kernels()
    forward, memf = instrument(True, W, V, H)
    reverse, memr = instrument(False, W, V, H)
    assert forward["mean_read_work"] == reverse["mean_read_work"] == F(17, 4)
    sharp = [9 * (i + 1) for i, (_, d) in enumerate(measure.LAW.LABELS) if abs(d) == 1]
    gram = H[np.ix_(sharp, sharp)] - 32 * np.eye(4, dtype=int)
    assert np.array_equal(gram, gram.T) and np.array_equal(gram @ gram, 4 * gram)
    assert not H[np.ix_(sharp, np.setdiff1d(np.arange(99), sharp))].any()
    assert not np.any((H - h0)[np.ix_(sharp, sharp)] - 32 * np.eye(4, dtype=int))
    empty = measure.LAW.SRC.index((0, 0))
    one = [measure.LAW.SRC.index(pair) for pair in [(0, 1), (1, 0)]]
    assert not H[:, empty].any()
    np.testing.assert_array_equal(H[np.ix_(one, one)], [[34, -2], [-2, 34]])
    out = []
    for level in range(3):
        board = measure.observations.federation(level)
        graph = nx.Graph(board.seams)
        tree = nx.bfs_tree(graph, 0, sort_neighbors=sorted)
        dist = nx.single_source_shortest_path_length(tree, 0)
        end = min(p for p in dist if dist[p] == max(dist.values()))
        route = nx.shortest_path(tree, 0, end)
        n = len(route) - 1
        ids = {e: i for i, e in enumerate(board.seams)}
        trail = [ids[tuple(sorted((a, b)))] for a, b in zip(route, route[1:])]
        blank = np.zeros(len(ids), dtype=int)
        filled = blank.copy()
        filled[trail] = [memf if a < b else memr for a, b in zip(route, route[1:])]
        before, _ = energy(board, route[0], blank, H)
        after, terms = energy(board, route[-1], filled, H)
        bare, _ = energy(board, route[-1], blank, H)
        stored = F(int(terms[trail].sum()), 4)
        excess = after - bare
        assert stored == F(17, 2) * (n - 1) + F(17, 4)
        assert excess == F(17, 2) * n - F(51, 4)
        # A sharp branch is unique on success. Any other native term is
        # diagonal between different frontier memory labels, so the local read
        # work is its exact full-board work at every failed or successful try.
        read = forward["mean_read_work"] * n
        out.append(
            dict(
                owners=board.cells,
                ports=board.n,
                path=route,
                live_seams=n,
                supplied_mean_echoes=n * forward["mean_echoes"],
                initial_energy=before,
                simultaneous_terminal_energy=after,
                terminal_path_energy=stored,
                terminal_matched_blank=bare,
                excess_vs_terminal_ordinary=excess,
                net_work=after - before,
                work_from_sharp_reads=read,
                work_from_addressed_echoes=after - before - read,
                passive_path_lower_bound=8 * (n - 2),
                all_blank_one_load_path_upper=18,
                passive_excess_vs_any_one_load_blank=8 * (n - 2) - 18,
                path_floor_over_initial=F(8 * (n - 2), 1) / before,
                replica_memory_swap_gain=0,
                one_unretried_success_probability=F(1, 2**n),
            )
        )
    return dict(
        protocol_sha256=LOCK,
        forward=forward,
        reverse=reverse,
        empty_unit_block_4H=H[np.ix_(sharp, sharp)].tolist(),
        boards=out,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    w = p.parse_args().write

    def encode(x):
        return str(x) if isinstance(x, F) else x.item()

    out = json.dumps(run(), indent=2, sort_keys=True, default=encode) + "\n"
    path = HERE / "trail.json"
    if w:
        path.write_text(out)
    elif path.read_text() != out:
        raise AssertionError("Paid trail replay mismatch")
    print("Verified paid native trail and passive all-time energy floor.")
