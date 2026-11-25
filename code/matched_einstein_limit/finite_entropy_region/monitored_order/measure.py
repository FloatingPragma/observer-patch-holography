"""Recover distances from native monitoring and test the registered regions as centered balls."""

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from model.experiments.interval_echo_einstein_station.equilibrium_refinement import (
    measure as monitored,
)
from model.experiments.interval_echo_einstein_station.proper_ball_ruler import (
    measure as ruler,
)

RECORD = Path(__file__).with_name("results.json")


def word_order(board, words, seams):
    h = 2 * monitored.generator(board, words, seams)
    np.testing.assert_array_equal(h, h.astype(np.int64))
    doubled = sparse.csr_matrix(h.astype(np.int64))
    graph = doubled.copy()
    graph.setdiag(0)
    graph.eliminate_zeros()
    np.testing.assert_array_equal(graph.data, -np.ones_like(graph.data))
    successors = -graph
    # Each native port link has all M-2 occupied spectator calibrations.
    recovered = Counter()
    for i, j in zip(*sparse.triu(successors).nonzero(), strict=True):
        a, b = set(words[i]), set(words[j])
        assert len(a & b) == 1
        recovered[tuple(sorted(a ^ b))] += 1
    ports = sorted(set(v for w in words for v in w))
    assert recovered == Counter({edge: len(ports) - 2 for edge in seams})
    d = shortest_path(successors, directed=False, unweighted=True)
    assert np.isfinite(d).all()
    distances = d.astype(np.int64)
    np.testing.assert_array_equal(distances, d)
    src, dest = np.unravel_index(np.argmax(distances), distances.shape)
    depth = int(distances[src, dest])
    # At distance d the exact number of length-d paths is a positive integer.
    counts = [0] * len(words)
    counts[src] = 1
    for step in range(1, depth + 1):
        for j in np.flatnonzero(distances[src] == step):
            adjacent = successors.indices[
                successors.indptr[j] : successors.indptr[j + 1]
            ]
            counts[j] = sum(
                counts[k] for k in adjacent if distances[src, k] == step - 1
            )
    ways = counts[dest]
    assert ways > 0
    # Independent exact modular check of the native matrix powers at the diameter pair.
    prime = 65521
    vector = np.zeros(len(words), dtype=np.int64)
    vector[src] = 1
    for _ in range(depth):
        assert vector[dest] == 0
        vector = (doubled @ vector) % prime
    assert vector[dest] == (-1) ** depth * ways % prime
    return dict(
        ports=len(ports),
        native_undirected_links=len(seams),
        complete_words=len(words),
        nonzero_conditional_transitions_at_algebraic_dwell=len(words) ** 2,
        word_hop_diameter=depth,
        from_word=list(words[src]),
        to_word=list(words[dest]),
        shortest_word_walks=ways,
        leading_probability_coefficient=str(
            Fraction(ways**2, 4**depth * factorial(depth) ** 2)
        ),
        quadratic_contact_spectators=len(ports) - 2,
    )


def integer_distances(local, blocked):
    active, d = ruler.distances(local, blocked)
    assert np.isfinite(d).all()
    z = d.astype(np.int64)
    np.testing.assert_array_equal(z, d)
    return active, z


def gap_test(active, distances, ports, cut):
    labels = ports[active]
    gates = np.flatnonzero(np.isin(labels, list(cut)))
    inside_radius = distances.max(axis=1)
    outside_radius = distances[:, gates].min(axis=1) + 1
    gaps = outside_radius - inside_radius
    best = int(gaps.max())
    i = int(np.argmax(gaps))
    j = gates[np.argmin(distances[i, gates])]
    example = dict(
        center=int(labels[i]),
        farthest_inside=int(labels[np.argmax(distances[i])]),
        last_inside=int(inside_radius[i]),
        closest_outside=int(cut[int(labels[j])]),
        first_outside=int(outside_radius[i]),
    )
    balls = []
    for i in np.flatnonzero(gaps > 0):
        radial = distances[i]
        balls.append(
            dict(
                center=int(labels[i]),
                radius=int(inside_radius[i]),
                first_outside=int(outside_radius[i]),
                shell_counts=np.bincount(radial).tolist(),
                exposed_rim_radii=sorted(map(int, radial[gates])),
                transmitting_contacts=len(gates),
                volume_ports=len(labels),
            )
        )
    return best, example, balls


def verify_ambient(board, ports, residents, active, reduced):
    # Test the eliminated exterior distances against direct full-board integer hop counts.
    rows = np.setdiff1d(np.arange(board.n), residents)
    inside = np.searchsorted(rows, ports[active])
    graph = board.adjacency().tocsr()[rows][:, rows]
    direct = shortest_path(
        graph, method="D", indices=inside, directed=False, unweighted=True
    )
    np.testing.assert_array_equal(direct[:, inside], reduced)


def local_relabelings(ports, cut, distances, positives, words, h):
    # These transports must preserve the *blank* cut before inspecting either occupied stock.
    graph = nx.Graph()
    for i, p in enumerate(ports):
        graph.add_node(int(p), cut=int(p) in cut)
        for j in range(i):
            if distances[i, j] == 1:
                graph.add_edge(int(p), int(ports[j]))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph, node_match=lambda a, b: a["cut"] == b["cut"]
    )
    permutations = list(matcher.isomorphisms_iter())
    word_id = {word: i for i, word in enumerate(words)}
    result = []
    for origin in positives:
        mapping = []
        for other in positives:
            checked = None
            for perm in permutations:
                if {perm[k] for k in origin["resident_word"]} != set(
                    other["resident_word"]
                ) or perm[origin["center"]] != other["center"]:
                    continue
                relabel = np.searchsorted(ports, [perm[p] for p in ports])
                if not np.array_equal(distances, distances[np.ix_(relabel, relabel)]):
                    continue
                ix = np.searchsorted(
                    other["_active"], [perm[k] for k in origin["_active"]]
                )
                if np.array_equal(
                    origin["_distance"], other["_distance"][np.ix_(ix, ix)]
                ):
                    order = [
                        word_id[tuple(sorted((perm[a], perm[b])))] for a, b in words
                    ]
                    np.testing.assert_array_equal(h, h[np.ix_(order, order)])
                    checked = [int(perm[k]) for k in ports]
                    break
            if checked is not None:
                mapping.append(dict(word=other["resident_word"], permutation=checked))
        assert len(mapping) == 3
        result.append(mapping)
    for row, maps in zip(positives, result, strict=True):
        row["local_rooted_cut_isometries"] = maps
        del row["_active"], row["_distance"]


def region(board, record):
    ports, seams, cut, words = monitored.apparatus(board, record["cells"])
    assert sorted(cut) == record["cut_ports"] and len(words) == record["enrolled_words"]
    ports = np.asarray(ports)
    local = ruler.ambient_reduction(board, ports, cut)
    all_ports, blank_d = integer_distances(local, [])
    blank_gap, witnesses, _ = gap_test(all_ports, blank_d, ports, cut)
    unrestricted = dict(best_gap=blank_gap, best_gap_witness=witnesses)
    positives = []
    counts = Counter()
    for j, word in enumerate(words):
        blocked = np.searchsorted(ports, word)
        active, d = integer_distances(local, blocked)
        gap, _, balls = gap_test(active, d, ports, cut)
        counts[gap] += 1
        if j in (0, len(words) // 2, len(words) - 1) or balls:
            verify_ambient(board, ports, word, active, d)
        for ball in balls:
            center = ball["center"]
            gate_indices = np.searchsorted(ports, list(cut))
            k = len(set(word) & set(cut))
            positives.append(
                dict(
                    resident_word=list(word),
                    cut_occupancy=k,
                    equilibrium_interior_multiplicity=comb(len(cut), k),
                    unoccupied_center_to_gates=sorted(
                        map(int, blank_d[np.searchsorted(ports, center), gate_indices])
                    ),
                    _active=ports[active],
                    _distance=d,
                    **ball,
                )
            )
    if positives:
        local_relabelings(
            ports,
            cut,
            blank_d,
            positives,
            words,
            monitored.generator(board, words, seams),
        )
    total_biased = sum(Fraction(1, 2 ** len(set(w) & set(cut))) for w in words)
    valid_words = {tuple(r["resident_word"]) for r in positives}
    hit_biased = sum(Fraction(1, 2 ** len(set(w) & set(cut))) for w in valid_words)
    return dict(
        board_cells=board.cells,
        owners=record["cells"],
        cut_ports=list(map(int, sorted(cut))),
        monitoring=word_order(board, words, seams),
        unoccupied=unrestricted,
        possible_protected_balls=positives,
        protected_word_gap_counts={str(k): v for k, v in sorted(counts.items())},
        any_center_equilibrium_probability=str(Fraction(len(valid_words), len(words))),
        any_center_outward_rejection_equilibrium_probability=str(
            hit_biased / total_biased
        ),
    )


def run():
    data = json.loads(Path(monitored.__file__).with_name("results.json").read_text())
    rows = []
    for level in (0, 1, 2):
        board = monitored.federation(level)
        for row in data["observations"]:
            if row["board_cells"] == board.cells:
                rows.append(region(board, row))
    return dict(monitored_dwell="1/2", rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    if args.write:
        RECORD.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        assert output == json.loads(RECORD.read_text())
        print(
            "Verified the nine native regions, 2,916 resident words, and centered cut isometries."
        )
