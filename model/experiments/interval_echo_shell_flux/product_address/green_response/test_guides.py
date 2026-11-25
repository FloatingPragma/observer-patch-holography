"""Check the native endpoints and the underlying midpoint-strip construction."""

from itertools import pairwise, product
import json

import pytest

from model.executable_model import federation
from . import guides


def test_bundled_boards_reconstruct_the_original_deposits():
    stored = json.loads(guides.RECEIPT.read_text())
    observed = guides.run()
    assert observed["ordered_cell_endpoint_pairs"] == 132
    for level, row in observed["levels"].items():
        assert row == stored["levels"][level]
        board = federation(int(level))
        paths = [*row["ports"], row["source_flag"]]
        middle = row["addresses_per_guide"] // 2
        central = [p[middle - 1 : middle + 1] for p in row["ports"]]
        assert row["released_port_triples"] == [list(c) for c in product(*central)]
        used = [v for path in paths for v in path]
        assert len(used) == len(set(used)) == 3 * row["addresses_per_guide"] + 2
        # The bundled full-port table is independent of the walker's split
        # internal-template/external-seam verifier.
        for path in paths:
            assert all(
                v in board.neighbours[u] and u in board.neighbours[v]
                for u, v in pairwise(path)
            )
        assert row["directed_sharp_transfers"] == 2 * (
            3 * (row["addresses_per_guide"] - 1) + 1
        )


def test_missing_crossing_and_flag_collision_fail():
    row = json.loads(guides.RECEIPT.read_text())["levels"]["1"]
    external = set(federation(1).external)
    ports = row["ports"]
    seam = next(
        tuple(sorted((a, b))) for a, b in pairwise(ports[0]) if a // 12 != b // 12
    )
    with pytest.raises(ValueError, match="Missing deposited seam"):
        guides.verify(1, external - {seam}, ports, row["source_flag"])
    with pytest.raises(ValueError, match="fourth coarse face"):
        guides.verify(1, external, ports, ports[0][:2])


def test_all_native_gateway_pairs_and_a_relabelled_matching():
    for a in range(12):
        for b in range(a + 1, 12):
            for i, j in ((a, b), (b, a)):
                path = guides.visit_cell(i, j)
                assert path[0] == i and path[-1] == j
                assert sorted(path) == list(range(12))
                assert all(
                    tuple(sorted(e)) in guides.CELL_SEAMS for e in pairwise(path)
                )

    # Change only the assignment of physical gateways on every owner's ports.
    board = federation(2)

    def relabel(v):
        owner, local = divmod(v, 12)
        return 12 * owner + (local + 5) % 12

    changed = guides.construct(2, [(relabel(u), relabel(v)) for u, v in board.external])
    native = json.loads(guides.RECEIPT.read_text())["levels"]["2"]
    assert changed["owners"] == native["owners"]
    assert changed["addresses_per_guide"] == native["addresses_per_guide"]
    assert changed["ports"] != native["ports"]


def test_continuation_side_row_incidence():
    # A side of the subdivided triangle x+y<=m in integer coordinates.
    # Alternating the up- and down-triangles gives a constructive dual path.
    for m in (1, 2, 4, 8, 16, 32):
        strip = []
        for x in range(m):
            strip.append(frozenset(((x, 0), (x + 1, 0), (x, 1))))
            if x + 1 < m:
                strip.append(frozenset(((x + 1, 0), (x, 1), (x + 1, 1))))
        assert len(set(strip)) == len(strip) == 2 * m - 1
        assert all(0 <= x and 0 <= y and x + y <= m for t in strip for x, y in t)
        assert all(len(a & b) == 2 for a, b in pairwise(strip))
