"""Addressed preparations converting three winding bits into physical probe vacancies."""

from functools import cache
from itertools import chain, pairwise
import networkx as nx
import numpy as np

from model.executable_model import federation
from ..apparatus import fixture
from .pulse import primitive


def edge(a, b):
    return tuple(sorted((a, b)))


def transp(path):
    links = list(pairwise(path))
    return links + links[-2::-1] if links else []


@cache
def chart():
    raw = fixture()
    board = federation(0)
    graph = nx.Graph(board.seams)
    tip = raw["tip"]
    cell = set(range(12 * (tip // 12), 12 * (tip // 12 + 1)))
    written = raw["marked"][3:]
    original_opening = {edge(*p) for p in raw["blank"] + written}
    quiet = {
        edge(a, b)
        for a, b in graph.subgraph(cell).edges
        if edge(a, b) not in original_opening
        and (a in raw["ports"] or b in raw["ports"])
    }
    assert len(quiet) == 18
    physical = set(raw["ports"])
    exterior = set(graph) - physical - set(raw["donor"])
    g = nx.Graph(raw["blank"])
    # Open one real equal-load shutter well outside the receiving owner. The
    # staging area contains no emitting or receiving head and connects to all
    # other ordinary guide ports by blank bonds.
    region = max(nx.connected_components(graph.subgraph(exterior)), key=len)
    distances = nx.single_source_shortest_path_length(g, tip)
    path_edge = min(
        ((a, b) for a in physical - cell for b in graph[a] if b in region),
        key=lambda e: (distances[e[0]], e),
    )
    g.add_edge(*path_edge)
    g.add_edges_from(graph.subgraph(exterior).edges)
    g.add_edges_from(graph.subgraph(cell).edges)
    distance_to_arms = nx.multi_source_dijkstra_path_length(graph, physical)
    stage = sorted(v for v in region if v != path_edge[1] and distance_to_arms[v] > 1)
    holds = stage[:6]
    probes = stage[6:9]
    stock = stage[9:12]
    reference = next(p for p in stage[12:] if board.degree[p] == 5)
    assert all(nx.has_path(g, tip, p) for p in [*holds, *probes, *stock, reference])
    return dict(
        g=g,
        all_edges={edge(*p) for p in graph.subgraph(cell).edges},
        tip=tip,
        marks=written,
        quiet=quiet,
        cell=sorted(cell),
        base=12 * (tip // 12),
        storage=holds,
        probes=probes,
        stock=stock,
        reference=reference,
        service=path_edge,
        donor=raw["donor"][0],
        private={p for p in cell & physical if board.degree[p] == 5},
        gateways={p for p in cell if board.degree[p] == 6},
    )


@cache
def mount(rotation):
    data = chart()
    base = data["base"]

    def turn(v):
        return base + rotation[v - base]

    marks = [(turn(a), turn(b)) for a, b in data["marks"]]
    hub = turn(data["tip"])
    quiet = {edge(turn(a), turn(b)) for a, b in data["quiet"]}
    return hub, marks, quiet


@cache
def path(active, a, b, excluded=()):
    graph = chart()["g"].copy()
    graph.remove_edges_from(active)
    graph.remove_nodes_from(excluded)
    return tuple(transp(nx.shortest_path(graph, a, b)))


def swap(word, links):
    for a, b in links:
        word[a], word[b] = word[b], word[a]


def clear_or_restore(word, hub, quiet, marks):
    # The vacancy at the sharp hub can temporarily vacate it when a private
    # read or equal-load shutter write needs a different endpoint condition.
    blank = [
        b if a == hub else a
        for a, b in chart()["all_edges"] - quiet - {edge(*e) for e in marks}
        if hub in (a, b)
    ]
    assert len(blank) >= 3
    for a, b in sorted(quiet):
        relay = (
            [(hub, next(p for p in blank if p not in (a, b)))] if hub in (a, b) else []
        )
        swap(word, relay)
        assert word[a] == word[b] == 1
        swap(word, relay[::-1])


@cache
def encode(rotation):
    """Return the complete material word and phase on all eight input basis states."""
    lab = chart()
    hub, marks, quiet = mount(rotation)
    edges = tuple(sorted(edge(a, b) for a, b in marks))
    stages = lab["storage"]
    stocks = lab["probes"]
    data = []
    phase_paths = []
    for bits in range(8):
        word = np.ones(240, dtype=int)
        word[hub] = 0
        word[stocks] = 0
        word[lab["stock"]] = 2
        word[lab["donor"]] = 2
        initial = word.copy()
        assert all(word[p] == 1 for p in lab["service"])
        clear_or_restore(word, hub, quiet, marks)
        base = lab["base"]
        rotated_private = {base + rotation[p - base] for p in lab["private"]}
        if {base + rotation[p - base] for p in lab["gateways"]} == lab["gateways"]:
            assert rotated_private.isdisjoint(lab["gateways"])
            for p in sorted(rotated_private):
                route = path(edges, hub, p)
                swap(word, route)
                assert word[p] == 0
                swap(word, route[::-1])
        # Each zero head clears at its stock pair (checked in primitive).
        # Only the three direction heads remain written inside the owner.
        amp = 1.0 + 0j
        for j, (a, g) in enumerate(marks):
            swap(word, path(edges, stocks[j], a))
            assert (word[a], word[g]) == (0, 1)
            bit = (bits >> j) & 1
            phases, _ = primitive(int(a > g))
            amp *= phases[bit]
            if not bit:
                word[a], word[g] = word[g], word[a]
        # Three material dual rails carry the whole unknown code. Park them
        # while the original charge acts as a reset probe on the fixed heads.
        for j, p in enumerate(chain.from_iterable(marks)):
            swap(word, path(edges, p, stages[j]))
        for a, g in marks:
            swap(word, path(edges, hub, g))
            assert word[g] == 0 and word[a] == 1 and word[lab["reference"]] == 1
            # In both equal-weight splitter branches, the phase-rail
            # transposition restores every intermediate particle and head.
            rail = path(edges, a, lab["reference"], (g,))
            w = word.copy()
            w[a], w[g] = w[g], w[a]
            swap(w, rail)
            assert w[lab["reference"]] == 0 and w[g] == 1
            swap(w, rail[::-1])
            assert w[a] == 0 and w[g] == 1
            amp *= primitive(int(a > g))[1]
            swap(word, path(edges, hub, g)[::-1])
            if bits == 0:
                phase_paths.append(len(rail))
        # The three heads are blank. Put the data back onto the six real
        # receiving ports; the upcoming spatial turn operates on these loads.
        for j, p in enumerate(chain.from_iterable(marks)):
            swap(word, path((), p, stages[j])[::-1])
        assert word[hub] == 0
        assert all(
            word[a] == (1 - ((bits >> j) & 1)) and word[b] == ((bits >> j) & 1)
            for j, (a, b) in enumerate(marks)
        )
        assert np.array_equal(word[lab["storage"]], initial[lab["storage"]])
        assert all(word[p] == 1 for p in stocks)
        assert (word == 0).sum() == 4 and (word == 2).sum() == 4
        data.append((tuple(word[lab["cell"]].tolist()), amp))
    return data, phase_paths


@cache
def rotate_program(rotation):
    lab = chart()
    base = lab["base"]
    # Factor the actual icosahedral vertex permutation into transpositions,
    # each compiled along blank nearest-neighbor links inside the cell.
    graph = nx.Graph(lab["all_edges"])
    perm = [base + p for p in rotation]
    moved = list(lab["cell"])
    dest = [None] * 12
    for k, p in enumerate(perm):
        dest[p - base] = moved[k]
    steps = []
    for target, value in enumerate(dest):
        j = moved.index(value)
        if j == target:
            continue
        bonds = transp(nx.shortest_path(graph, base + j, base + target))
        steps.extend(bonds)
        moved[j], moved[target] = moved[target], moved[j]
    assert moved == dest
    return steps


def compose(x, y):
    return tuple(x[y[i]] for i in range(12))


def dressed(rotation, turn, rail_phase=False):
    # E_{R D}† P_R E_D. Use independent complete basis replays for the two
    # mounts; neither a representation sign nor an internal rotation is inserted.
    lab = chart()
    start, _ = encode(rotation)
    end = compose(turn, rotation)
    out, _ = encode(end)
    output = {w: (j, z) for j, (w, z) in enumerate(out)}
    t = np.zeros((8, 8), complex)
    base = lab["base"]
    for i, (word, z) in enumerate(start):
        w = list(word)
        if rail_phase:
            # Each of the three g rails visits a real unfrozen private zero
            # reference for pi/8 and returns with its original load.
            for _, g in mount(rotation)[1]:
                rail = path((), g, lab["reference"])
                assert rail and g != lab["reference"]
                material = np.ones(240, dtype=int)
                material[lab["cell"]] = w
                material[lab["stock"]] = 2
                material[lab["donor"]] = 2
                before = material.copy()
                swap(material, rail)
                assert material[lab["reference"]] == before[g]
                z *= np.exp(-1j * np.pi * material[lab["reference"]] ** 2)
                swap(material, rail[::-1])
                assert np.array_equal(material, before)
        for a, b in rotate_program(turn):
            w[a - base], w[b - base] = w[b - base], w[a - base]
        j, other = output[tuple(w)]
        t[j, i] = z * other.conjugate()
    return end, t
