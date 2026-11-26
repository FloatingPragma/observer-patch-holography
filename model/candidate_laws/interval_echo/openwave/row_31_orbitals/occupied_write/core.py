"""Native signed-memory blocks and their one-load spectra."""

import networkx as nx
import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import EMITTER_SEAM, federation
from ..written_transitions.measure import generator

LAW = candidate.native_law


def integer_blocks():
    """Eight times the written Hamiltonian, from the actual integer echo pattern."""
    W = np.array(LAW._write_pattern(), dtype=np.int64)
    R, D2 = [np.array(m, dtype=np.int64) for m in LAW._label_ops()]
    complement = 2 * np.eye(90, dtype=np.int64) - W @ W.T
    identity = np.eye(6, dtype=np.int64)
    # coin 0 / second port, coin 1 / first port
    seed = np.array([0, 0, 1, 0, 1, 0], dtype=np.int64)
    words = [LAW.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1))]
    for k in (0, 8):
        G2 = np.kron(
            np.eye(10, dtype=np.int64) - R + 2 * k * D2, np.eye(9, dtype=np.int64)
        )
        H8 = complement @ G2 @ complement
        if np.any(complement @ G2 @ W):
            raise ArithmeticError("Blank and written blocks couple")
        for coin in (0, 1):
            for word in words:
                index = 9 * LAW.LABELS.index((coin, 0)) + word
                if np.any(H8[:, index]):
                    raise ArithmeticError("Nonzero zero-direction action")
        for d in (1, 2):
            for sign in (1, -1):
                embedding = np.zeros((90, 6), dtype=np.int64)
                for coin in (0, 1):
                    for j, word in enumerate(words):
                        for direction, factor in ((d, 1), (-d, sign)):
                            embedding[
                                9 * LAW.LABELS.index((coin, direction)) + word,
                                3 * coin + j,
                            ] = factor
                mu = k * d * d + (sign < 0)
                expected = 8 * mu * identity
                if d == 1:
                    expected = expected - 4 * mu * np.outer(seed, seed)
                if not np.array_equal(H8 @ embedding, embedding @ expected):
                    raise ArithmeticError("Occupied-core reduction failed")
    return dict(
        memory_dimensions=[1, 2, 4, 4],
        sector_order=["blank", "zero", "unit", "double"],
        occupied_seed_coin_and_position=[[0, "second"], [1, "first"]],
        global_subthreshold_multiplicities=[0, 4, 2, 0],
    )


def chart(level):
    board = federation(level)
    a, b = EMITTER_SEAM[level]
    graph = nx.Graph(board.seams)
    graph.remove_edge(a, b)

    def route(p, forbidden, predicate):
        good = graph.subgraph(set(graph) - set(forbidden))
        paths = nx.single_source_shortest_path(good, p)
        return min(
            (v for i, v in paths.items() if predicate(i)), key=lambda v: (len(v), v)
        )

    phase = route(a, [b], lambda i: board.degree[i] == 5)
    # Transpose to two distinct empty sites and restore all intermediate loads.
    park_a = route(a, [b], lambda i: i not in (a, b))
    park_b = route(b, [a, park_a[-1]], lambda i: i not in (a, b, park_a[-1]))
    if max(board.degree) > 6 or (a, b) not in board.seams:
        raise ValueError("Unlicensed board")
    return (
        board,
        (a, b),
        dict(
            core=[a, b],
            endpoint_degrees=[int(board.degree[v]) for v in (a, b)],
            phase_route=phase,
            reference_parks=[park_a, park_b],
        ),
    )


def standing_modes(board, seam, k, kind, odd):
    """Yield every eigenfrequency and two launches in a common spectral frame."""
    n = board.n
    a, b = seam
    A = (
        board.laplacian().toarray() / 2
        if kind == "blank"
        else generator(board, seam, k)
    )
    e, u = eigh(A, driver="evd")
    orientation = -1 if odd else 1
    if kind != "occupied":
        src = np.array([(u[a] + orientation * s * u[b]) / np.sqrt(2) for s in (1, -1)])
        return e, src, lambda vec: abs(u @ vec) ** 2

    # Bright and dark coin lines at each eigenvector of the cut graph.
    lengths = np.hypot(u[b], u[a])
    c = np.divide(u[b], lengths, out=np.ones(n), where=lengths != 0)
    s = np.divide(u[a], lengths, out=np.zeros(n), where=lengths != 0)
    z = lengths / np.sqrt(2)
    mu = k + int(odd)
    if mu:
        bright, v = eigh(np.diag(e) - mu * np.outer(z, z), driver="evd")
    else:
        bright, v = e, np.eye(n)
    src = []
    for phase in (1, -1):
        first, second = u[b] / np.sqrt(2), phase * u[a] / np.sqrt(2)
        src.append(np.r_[v.T @ (c * first + s * second), -s * first + c * second])

    def positions(coeff):
        bc = v @ coeff[:n]
        first = u @ (c * bc - s * coeff[n:])
        second = u @ (s * bc + c * coeff[n:])
        return abs(first) ** 2 + abs(second) ** 2

    return np.r_[bright, e], np.array(src), positions


def local_index(label, word):
    m = 0 if label is None else 1 + LAW.LABELS.index(label)
    return 9 * m + LAW.SRC.index(word)


def occupied_source(odd):
    """Build the occupied standing write with a blank swap, private phase and echo."""
    from scipy.linalg import expm

    a = local_index(None, (1, 0))
    b = local_index(None, (0, 1))
    h0, h1 = LAW.local_term_parts()
    before = expm(-1j * np.pi / 2 * (h0 + 8 * h1))[:, a]
    # A private degree-reference contact at cost 8 advances the opposite relative phase.
    pulse_time = (1 if odd else 3) * np.pi / 16
    before[a] *= np.exp(-8j * pulse_time)
    sign = -1 if odd else 1
    expected = np.zeros(99, dtype=complex)
    for c, word in ((0, (0, 1)), (1, (1, 0))):
        expected[local_index((c, 1), word)] = 0.5
        expected[local_index((c, -1), word)] = 0.5 * sign
    source = LAW.echo() @ before
    align = np.vdot(expected, source)
    if abs(abs(align) - 1) > 1e-12 or abs(before[b] / before[a] - sign) > 1e-12:
        raise ArithmeticError("The occupied source did not produce the standing write")
    return source / align


def separated_standards(wait):
    """Evolve actual separated copies to obtain the reference-arm overlaps."""
    from scipy.linalg import expm

    standards = {}
    h0, h1 = LAW.local_term_parts()
    for k in (0, 8):
        unitary = expm(-1j * wait * (h0 + k * h1))
        for odd in (False, True):
            prepared = occupied_source(odd)
            # Two orthogonal parked positions tag the coin. The original endpoints
            # have word (0,0) in both branches. All other contacts are shuttered.
            overlap = 0j
            for coin, word in ((0, (0, 1)), (1, (1, 0))):
                vacuum = np.zeros(99, dtype=complex)
                for direction in (1, -1):
                    label = (coin, direction)
                    vacuum[local_index(label, (0, 0))] = prepared[
                        local_index(label, word)
                    ]
                overlap += np.vdot(vacuum, unitary @ vacuum)
            standards[f"core_{k}_{'odd' if odd else 'even'}"] = overlap.conjugate()

    board = federation(0)
    uniform = np.ones(board.n) / np.sqrt(board.n)
    # Sum the physical one-load edge terms and the blank degree references.
    for k in (8, 0):
        action = k * (6 - board.degree) * uniform
        words = [LAW.SRC.index(w) for w in ((1, 0), (0, 1))]
        blank = (h0 + k * h1)[np.ix_(words, words)]
        for a, b in board.seams:
            action[[a, b]] += blank @ uniform[[a, b]]
        energy = float(np.vdot(uniform, action))
        if max(abs(action - energy * uniform)) > 1e-12:
            raise ArithmeticError("Separated uniform load has no sharp reference")
        standards[f"uniform_blank_{k}"] = np.exp(1j * wait * energy)
    return standards


def terminal_reader(board, seam):
    """Reconstruct every local single-load or vacant word from fresh signed writes."""
    # The joint sharp-word effects of separate echo heads commute. Sequential
    # adjacent echoes would move the load, so each head restores the material by
    # inverse echo after its signed label is copied to a public register.
    W = LAW._write_pattern()
    for col, (first, second) in enumerate(LAW.SRC):
        directions = {LAW.LABELS[i // 9][1] for i in np.flatnonzero(W[:, col])}
        if directions != {first - second} or int(W[:, col] @ W[:, col]) != 2:
            raise ArithmeticError("The signed echo head mixes differences")
    trees = {}
    trees_checked = 0
    for owner in range(board.cells):
        ports = set(range(12 * owner, 12 * (owner + 1)))
        graph = nx.Graph()
        graph.add_nodes_from(sorted(ports))
        graph.add_edges_from(
            (a, b)
            for a in sorted(ports)
            for b in board.neighbours[a]
            if b in ports and (min(a, b), max(a, b)) != seam
        )
        tree = list(nx.bfs_edges(graph, 12 * owner))
        if len(tree) != 11:
            raise ArithmeticError("Terminal tree is disconnected")
        local = [(a % 12, b % 12) for a, b in tree]
        key = str(local)
        if key not in trees:
            for word in range(-1, 12):
                digits = [int(p == word) for p in range(12)]
                relative = {0: 0}
                for a, b in local:
                    relative[b] = relative[a] - (digits[a] - digits[b])
                low = min(relative.values())
                if [relative[i] - low for i in range(12)] != digits:
                    raise ArithmeticError("Terminal tree misreads occupancy")
            trees[key] = local
        trees_checked += 1
    return dict(
        owners=trees_checked,
        local_tree_geometries=list(trees.values()),
        observable_words_per_owner=13,
        heads_per_owner=11,
    )
