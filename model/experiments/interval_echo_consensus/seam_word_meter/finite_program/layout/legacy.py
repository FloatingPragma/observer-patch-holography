"""Matching witnesses for the original 240-port addresses, with or without shelf moves."""

import json
from itertools import chain
from pathlib import Path

import networkx as nx
from ..machine import Bench

HERE = Path(__file__).resolve().parent


def witness():
    bench = Bench()
    old = set(chain.from_iterable(chain.from_iterable(bench.words)))
    readout = {*bench.heads, *bench.writer, *chain.from_iterable(bench.references)}
    service = readout | set(chain.from_iterable(bench.flags))
    graph = bench.graph.subgraph(set(bench.graph) - readout)
    assert nx.is_connected(graph) and all(
        set(bench.graph[p]) & set(graph) for p in readout
    )
    output = {}
    for name, fixed in (
        ("fixed_services", service),
        ("fixed_services_and_84", service | old),
    ):
        available = bench.graph.subgraph(set(bench.graph) - fixed)
        pairs = sorted(
            sorted(e) for e in nx.max_weight_matching(available, maxcardinality=True)
        )
        sizes = sorted(map(len, nx.connected_components(available)))
        # Any data dimer stays inside one remaining connected component.
        upper = sum(s // 2 for s in sizes)
        assert len(pairs) == upper
        output[name] = dict(
            component_sizes=sizes, maximum_new_dimers=upper, matching=pairs
        )
    return output


def replay():
    data = json.loads((HERE / "original.json").read_text())
    # A different maximum matching also attains the same bound. Replay the
    # stored witness directly against the original, physically labelled graph.
    book = Bench()
    old = set(chain.from_iterable(chain.from_iterable(book.words)))
    service = {
        *book.heads,
        *book.writer,
        *chain.from_iterable(book.references),
        *chain.from_iterable(book.flags),
    }
    for name, fixed in (
        ("fixed_services", service),
        ("fixed_services_and_84", service | old),
    ):
        row = data[name]
        edges = row["matching"]
        use = list(chain.from_iterable(edges))
        assert len(set(use)) == len(use) and set(use).isdisjoint(fixed)
        assert all(book.graph.has_edge(*e) for e in edges)
        sizes = sorted(
            map(
                len,
                nx.connected_components(book.graph.subgraph(set(book.graph) - fixed)),
            )
        )
        assert sizes == row["component_sizes"] and sum(s // 2 for s in sizes) == row[
            "maximum_new_dimers"
        ] == len(edges)


if __name__ == "__main__":
    replay()
    print("original shelf bounds and both matching witnesses passed")
