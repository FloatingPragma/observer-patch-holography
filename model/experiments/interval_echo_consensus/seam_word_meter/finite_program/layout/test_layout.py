"""Independent packing, route, and existing physical service checks."""

import json
from itertools import combinations
from types import SimpleNamespace

import networkx as nx
import pytest

from model.executable_model import carrier
from model.executable_model.refinement import DEPOSIT_SHA256
from ..machine import Bench
from . import measure
from .mount import Mount, PORTS, station


def test_every_station_chart_keeps_transit_and_a_perfect_shelf():
    graph = nx.Graph(carrier.CELL_SEAMS)
    cases = 0
    for fixed in range(12):
        for other in combinations(sorted(PORTS - {fixed}), 2):
            gateways = (fixed, *other)
            head, refs, match = station(fixed, gateways)
            frozen = {fixed, head, *refs}
            remainder = PORTS - frozen
            assert nx.is_connected(graph.subgraph(remainder))
            assert set(other) <= remainder and not (set((head, *refs)) & set(gateways))
            assert all(set(graph[n]) & remainder for n in frozen)
            assert all(graph.has_edge(head, r) for r in refs)
            assert sorted(j for e in match for j in e) == sorted(remainder)
            assert all(graph.has_edge(*e) for e in match)
            cases += 1
    assert cases == 660


def test_packed_tree_agrees_with_original_routing_machine():
    board = carrier.federation(0)
    mount = Mount(board)
    packed = mount.bench(registers=19, axes=2, width=3)
    original = Bench(
        board=board,
        words=packed.words,
        writer=mount.writer,
        heads=mount.heads,
        references=mount.references,
        flags=mount.flags,
    )
    for bench in (original, packed):
        inputs = [3] * 38
        inputs[0], inputs[1], inputs[36], inputs[37] = 2, 0, 6, 4
        bench.start(inputs)
        assert bench.mean_axes([(0, 36), (1, 37)]) == "accepted"
        assert bench.mean(0, 2) == "off_grid"
        bench.public("high")
    assert packed.loads == original.loads and packed.tape == original.tape
    assert packed.graph.capacity == 114
    # This read can use data ports during transit; compare each route with a
    # different graph search that knows nothing about the selected owner tree.
    graph = nx.Graph(board.seams)
    for source in range(board.n):
        for target in (*mount.writer, *mount.heads):
            if source == target:
                continue
            path = mount.path(source, target)
            mask = mount.protected - {source, target}
            assert not set(path) & mask
            assert all(
                graph.has_edge(x, y) and {x, y} != set(mount.writer)
                for x, y in zip(path, path[1:])
            )
            assert len(path) == len(set(path))


def test_limits_and_exact_replay():
    data = json.loads((measure.HERE / "records.json").read_text())
    assert json.loads(json.dumps(measure.bundled())) == data
    external = json.loads((measure.HERE / "level3.json").read_text())
    assert external["level"] == 3 and external["sha256"] == DEPOSIT_SHA256[3]
    assert (
        external["layout"]["cells"] == 1280 and external["layout"]["capacity"] == 7674
    )
    board = SimpleNamespace(
        cells=2, cell_seams=carrier.CELL_SEAMS, external=[(0, 12), (1, 13), (2, 14)]
    )
    small = Mount(board)
    assert small.capacity == 6
    with pytest.raises(ValueError, match="Insufficient"):
        small.bench(registers=7, axes=1, width=1)
    bare = small.bench(registers=0, axes=1, width=1)
    bare.start([])
    assert bare.mean_axes([]) == "accepted" and len(bare.public("low")) == 2


def test_historical_fixed_sites_and_matching_witnesses():
    from .legacy import replay

    replay()
