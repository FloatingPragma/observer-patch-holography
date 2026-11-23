"""Compose independently prepared native repair triples in one connected parent patch."""

import argparse
from itertools import combinations, product
import json
from math import log
from pathlib import Path

import networkx as nx

from model.executable_model import federation
from .measure import branch


def coarse_patch(level):
    car = federation(level)
    graph = nx.Graph((a//12, b//12) for a, b in car.external)
    eligible = []
    for home in sorted(graph):
        private = [c for c in sorted(graph[home]) if all(car.degree[12*c+p] == 5 for p in (6, 7))]
        eligible.extend((home, left, right) for left, right in combinations(private, 2))
    used = set()
    groups = []
    for triple in eligible:
        if used.isdisjoint(triple):
            groups.append(triple)
            used.update(triple)
    connected = max(nx.connected_components(graph.subgraph(used)), key=len)
    groups = [g for g in groups if set(g) <= connected]
    if len(groups)*3 != len(connected):
        raise ValueError("A repair preparation crosses the selected connected region")
    cut = [(a, b) for a, b in car.external if (a//12 in connected) != (b//12 in connected)]
    if not cut:
        raise ValueError("The coarse region requires an exterior")
    for home, left, right in groups:
        local = {home, left, right}
        local_cut = sorted(p for a, b in car.external if (a//12 in local) != (b//12 in local)
                           for p in (a, b))
        by_shared, words = {}, {}
        for x, a, b in product((0, 1), repeat=3):
            for mode in ("independent", "common", "disjoint"):
                # Independent private pairs choose either matched or mismatched outputs.
                bits = (a, b, x, a, a, b, b)
                result = branch(car, (left, right), local_cut, bits, (0, 1), mode, home=home)
                row, disagreement, word, boundary, loads = result
                expected_errors = {"independent": (a, b), "common": (a & b, a & b),
                                   "disjoint": (a & b, a & (1-b))}[mode]
                if row != (x, x ^ expected_errors[0], x ^ expected_errors[1]):
                    raise ValueError("A translated preparation changed its reported word")
                if loads != (6, 6, 6) or disagreement != (row[1] ^ row[2]):
                    raise ValueError("A translated volume proxy or material pointer changed")
                if by_shared.setdefault(x, boundary) != boundary:
                    raise ValueError("A private error reaches this group's instantaneous cut")
                if words.setdefault(row, word) != word:
                    raise ValueError("A retained input digit survives in the group's material state")
        if len(set(words.values())) != 8:
            raise ValueError("Region does not preserve all three-bit records")
    return dict(board_cells=car.cells, region_cells=len(connected), groups=groups,
                cut_contacts=len(cut), regional_load=6*len(connected),
                entropy_slope_at_independence=0., entropy_slope_at_one_eighth=-len(groups)*log(5),
                mean_routed_mismatch_slope=-2., total_routed_mismatch_slope=-2*len(groups))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps([coarse_patch(i) for i in range(3)], indent=2)+"\n")
