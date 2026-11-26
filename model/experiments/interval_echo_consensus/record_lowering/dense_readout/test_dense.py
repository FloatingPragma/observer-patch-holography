import gzip
import json
from collections import Counter, defaultdict

import networkx as nx

from model.executable_model.carrier import federation
from .apparatus import SITES, squared_distance
from .measure import HERE, read_paths, report


def test_complete_native_tapes_order_and_metadata():
    stored = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    assert json.loads(json.dumps(report())) == stored
    for m in stored["mounts"]:
        event = m["native_plan"]
        assert len(event) == 1192 and m["boundary_pairs"] == 192
        assert m["lost_remote_single_bit_interventions_per_control"] == 1192 - 64
        for w, answers in zip(m["inputs"], m["native"], strict=True):
            for word, tape in zip(w, answers, strict=True):
                wanted = int(word, 16)
                last = int(tape, 16)
                for j, e in enumerate(event):
                    assert (wanted >> e["source"]) & 1 == (last >> j) & 1
        assert Counter(Counter(e["receiver"] for e in event).values()) == {
            11: 8,
            15: 24,
            21: 24,
            30: 8,
        }


def test_shortest_path_comparison_and_parallel_occupancy():
    stored = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    for m in stored["mounts"]:
        graph = nx.Graph(federation(m["level"]).seams)
        protected = {p for s in m["layout"] for p in s}
        for calendar in [m["native_plan"], m["stale_plan"]]:
            busy = defaultdict(list)
            for e in calendar:
                path = e["path"]
                assert path[-1] == m["layout"][e["receiver"]][4]
                assert protected.isdisjoint(path[1:-1])
                assert all(graph.has_edge(a, b) for a, b in zip(path, path[1:]))
                for v in [*path, e["ref"]]:
                    busy[v].append((e["start"], e["end"]))
            for slots in busy.values():
                slots.sort()
                assert all(b <= c for (_, b), (c, _) in zip(slots, slots[1:]))
        times = m["clock"]["ordered_pair_slots"]
        for r, s in [(30, 9), (52, 20), (3, 0), (44, 44)]:
            a, h = m["layout"][s][0], m["layout"][r][4]
            keep = (set(graph) - protected) | {a, h}
            assert (
                times[r][s] == nx.shortest_path_length(graph.subgraph(keep), a, h) + 2
            )
        for r in range(len(SITES)):
            inside = max(
                t for s, t in enumerate(times[r]) if squared_distance(r, s) <= 4
            )
            outside = min(
                t for s, t in enumerate(times[r]) if squared_distance(r, s) > 4
            )
            assert inside >= outside
        assert (
            m["clock"]["any_placement"]["best"]["lower_bound"]
            == {80: 130, 320: 138}[m["owners"]]
        )


def test_parked_head_still_reads_its_own_physical_rail():
    appointments = [dict(path=[0, 1, 2], ref=3)]
    loads = [0, 0, 1, 0]
    assert read_paths(appointments, loads, parked=True, zero=0) == [1]
    assert read_paths(appointments, loads, zero=0) == [0]
    assert loads == [0, 0, 1, 0]
