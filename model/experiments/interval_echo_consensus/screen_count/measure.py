"""Acquire sharp digit responses, dual-read joint histograms and targeted preparations."""

import argparse
from functools import reduce
import gzip
import json
from operator import xor
from pathlib import Path
import networkx as nx
import numpy as np

from .apparatus import (
    READ_CODES,
    calibrate,
    decode,
    layout,
    rail_swaps,
    raw_word,
    read_words,
)

HERE = Path(__file__).parent
SIZE = 2**18
SEED = 1242603


def edges(config):
    sites = config["owners"]
    return [(sites.index(a // 12), sites.index(b // 12)) for a, b in config["inside"]]


def orientation(config, target, deposit):
    """Choose tree exchanges leaf-first; reject an incompatible component charge."""
    es = edges(config)
    n = len(config["owners"])
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    graph.add_edges_from(es)
    start = reduce(xor, (1 << u for u, _ in es), 0)
    source = config["owners"].index(config["boundary"][0] // 12)
    diff = target ^ start ^ (deposit << source)
    moves = 0
    for comp in nx.connected_components(graph):
        tree = list(nx.bfs_edges(graph, min(comp)))
        for parent, child in reversed(tree):
            if diff & (1 << child):
                j = next(j for j, e in enumerate(es) if {parent, child} == set(e))
                moves |= 1 << j
                diff ^= (1 << parent) | (1 << child)
    if diff:
        return None
    return moves


def bulk_decode(config, stations):
    """Push each sharp two-station rail outcome through the physical read menu."""
    n = len(config["owners"])
    e = config["inside"]
    k = 2 ** len(e)
    ids = np.arange(k, dtype=int)
    port_col = {
        12 * owner + j: i * 12 + j
        for i, owner in enumerate(config["owners"])
        for j in range(12)
    }
    base = np.zeros((k, 12 * n + 1), dtype=np.uint8)
    for j, (a, b) in enumerate(e):
        move = (ids >> j) & 1
        base[:, port_col[a]], base[:, port_col[b]] = 1 - move, move
    reports = []
    expected_ports = [*port_col, config["boundary"][1]]
    for dep in (0, 1):
        base[:, port_col[config["boundary"][0]]] = dep
        base[:, -1] = 1 - dep
        outputs = []
        for station in stations:
            assert [r["port"] for r in station["paths"]] == expected_ports
            flags = np.take(station["sharp"], base)
            ternary = np.take([READ_CODES.index(i) for i in range(3)], flags)
            parity = ternary[:, :-1].reshape(k, n, 12).sum(axis=2) % 2
            outputs.append(parity @ (1 << np.arange(n)) + (ternary[:, -1] << n))
        np.testing.assert_array_equal(*outputs)
        reports.append(outputs[0])
    return np.array(reports)


def acquire():
    rng = np.random.Generator(np.random.PCG64(SEED))
    transfer = rail_swaps()
    archive = []
    for level in range(3):
        config = layout(level)
        stations = calibrate(config)
        assert all(len(set(s["shuttered"])) < 3 for s in stations)
        n = len(config["owners"])
        m = len(config["inside"])
        rows = bulk_decode(config, stations)
        panels = {}
        for mode in (*transfer, "cut", "diagonal"):
            mask = 0
            if mode == "cut":
                mask = sum(1 << i for i in config["skip"])
            elif mode == "diagonal":
                mask = 2**m - 1
            valid = (np.arange(2**m) & mask) == 0
            dep = transfer.get(mode, 0) / 4
            probs = np.array([1 - dep, dep])[:, None] * valid / valid.sum()
            raw = rng.multinomial(SIZE, probs.ravel()).reshape(2, -1)
            joint = np.zeros(2 ** (n + 1), dtype=int)
            np.add.at(joint, rows.ravel(), raw.ravel())
            panels[mode] = dict(
                rail_counts=raw.tolist(), joint=joint.tolist(), refused=0
            )
        clocks = {}
        for name, num in transfer.items():
            clocks[name] = rng.multinomial(2**16, [1 - num / 4, num / 4]).tolist()
        clocks["diagonal"] = rng.multinomial(2**16, [1, 0]).tolist()
        targets = []
        if level == 0:
            for word in range(2**n):
                for deposit in (0, 1):
                    move = orientation(config, word, deposit)
                    raw = (
                        None
                        if move is None
                        else read_words(raw_word(config, move, deposit), stations)
                    )
                    if raw is not None:
                        assert decode(raw, n) == (word, 1 - deposit)
                    targets.append(
                        dict(target=word, deposit=deposit, swap=move, raw=raw)
                    )
        cycle = None
        graph = nx.Graph(edges(config))
        cycles = nx.cycle_basis(graph)
        if cycles:
            path = cycles[0]
            es = edges(config)
            mask = sum(
                1 << next(i for i, e in enumerate(es) if {u, v} == set(e))
                for u, v in zip(path, path[1:] + path[:1], strict=True)
            )
            cycle = [read_words(raw_word(config, mv, 0), stations) for mv in (0, mask)]
            assert cycle[0] != cycle[1] and decode(cycle[0], n) == decode(cycle[1], n)
        archive.append(
            dict(
                config=config,
                stations=stations,
                clocks=clocks,
                targets=targets,
                cycle=cycle,
                panels=panels,
            )
        )
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    f = HERE / "records.json.gz"
    if args.write:
        f.write_bytes(gzip.compress(raw, mtime=0))
    elif gzip.decompress(f.read_bytes()) != raw:
        raise ValueError("Sharp screen bank differs from its receipt")


if __name__ == "__main__":
    main()
