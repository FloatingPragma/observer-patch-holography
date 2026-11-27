"""Optional closed-board load reservoir for the sharp injection and removal."""

import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.sparse import eye

from ..measure import source_and_reader, defect_generator
from . import ring


def certificate(level, n):
    board = ring.federation(level)
    ring_ports = ring.geometry(level, n)
    network = nx.Graph((u // 12, v // 12) for u, v in board.external)
    used = {p // 12 for p in ring_ports}
    distances = nx.multi_source_dijkstra_path_length(network, used, weight=None)
    owner = min(k for k, d in distances.items() if d >= 2)
    packet, mark = source_and_reader(board, owner)
    source = min(np.flatnonzero(packet))
    ring_cut = {q for p in ring_ports for q in board.neighbours[p]} - set(ring_ports)
    assert all(packet[p] == 0 for p in set(ring_ports) | ring_cut)
    route = nx.shortest_path(nx.Graph(board.seams), int(source), ring_ports[0])
    swaps = list(zip(route, route[1:]))
    swaps += swaps[-2::-1]
    outcomes = []
    for amount in (1, 2):
        h = (
            board.laplacian() / 2 + 6 * ring.K * eye(board.n, format="csr")
            if amount == 1
            else defect_generator(board, amount, ring.K, 0, background=0)
        )
        energy = np.vdot(packet, h @ packet).real
        np.testing.assert_allclose(h @ packet, energy * packet, atol=1e-12)
        for sign in (-1, 1):
            loads = np.zeros(board.n, int)
            checker = amount * ((np.arange(n) + (sign < 0)) % 2)
            loads[ring_ports] = checker
            loads[source] = amount if sign > 0 else 0
            initial = loads.copy()
            for i, j in swaps:
                loads[i], loads[j] = loads[j], loads[i]
            target = initial.copy()
            p = ring_ports[0]
            target[source], target[p] = target[p], target[source]
            np.testing.assert_equal(target, loads)
            assert np.all(loads[list(ring_cut)] == 0)
            outcomes.append(
                dict(
                    amount=amount,
                    sign=sign,
                    conserved=int(loads.sum()),
                    ring_before=int(initial[ring_ports].sum()),
                    ring_after=int(loads[ring_ports].sum()),
                    reservoir_after=int(loads[source]),
                    reservoir_mode_energy=float(energy),
                )
            )
    return dict(
        level=level,
        n=n,
        reservoir_owner=owner,
        reservoir_port=int(source),
        path=route,
        compact_mode=mark["integer_amplitudes_in_cell"],
        stocks=outcomes,
    )


def report():
    return [certificate(level, n) for level, n in ((0, 14), (1, 14), (2, 14), (0, 16))]


if __name__ == "__main__":
    dst = Path(__file__).with_name("stock.json")
    data = report()
    dst.write_text(json.dumps(data, indent=1) + "\n")
