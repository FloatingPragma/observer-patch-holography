"""Full vacancy-sector hierarchy bounds and the commutant of charge-position readout."""

from collections import Counter
from itertools import product

import networkx as nx
import numpy as np

from ..apparatus import fixture, native_graph
from ...live_trimer.measure import blocks, LAW


def certify():
    lab = fixture()
    h = native_graph()
    h2 = h.copy()
    h2.data = np.rint(2 * h.data).astype(np.int64)
    np.testing.assert_array_equal(h2.data, 2 * h.data)
    diagonal = h2.diagonal()
    radius = np.asarray(abs(h2).sum(axis=1)).ravel() - abs(diagonal)
    lower, upper = int(np.min(diagonal - radius)), int(np.max(diagonal + radius))
    assert (lower, upper) == (16, 54)
    # Scale sqrt(energies) by their minimum. Each coordinate is in [1,2).
    # 8 sum(x_i^2)-3(sum x_i)^2 is convex in each coordinate separately.
    vertices = {str(k): 8 * (3 - k + 4 * k) - 3 * (3 + k) ** 2 for k in range(4)}
    assert max(vertices.values()) == 0
    ix = {v: i for i, v in enumerate(lab["ports"])}

    def block(a, b):
        i, j = ix[b], ix[a]
        return h2[64 * i : 64 * (i + 1), 64 * j : 64 * (j + 1)].toarray()

    graph = nx.Graph(lab["blank"])
    assert set(graph) == set(ix) and nx.is_connected(graph)
    for a, b in graph.edges:
        np.testing.assert_array_equal(block(a, b), -np.eye(64))
    # For any F commuting with all port projectors, F has one 64x64 block
    # per vacancy address. Blank hops identify all blocks with a common X.
    # The six directed marked hops and their adjoints constrain that X.
    generators = []
    for bit, (a, b) in enumerate(lab["marked"]):
        up = np.zeros((64, 64), int)
        for s in range(64):
            if not (s & (1 << bit)):
                up[s ^ (1 << bit), s] = 1
        np.testing.assert_array_equal(block(a, b), -8 * up)
        np.testing.assert_array_equal(block(b, a), -8 * up.T)
        generators.append(up)
    # Products of number projectors isolate every computational head word.
    # Raising/lowering thereafter supplies |u><v| for every pair.
    id64 = np.eye(64, dtype=int)
    ones = [np.diag(g @ g.T) for g in generators]
    for v in range(64):
        p = np.prod([d if v & (1 << j) else 1 - d for j, d in enumerate(ones)], axis=0)
        np.testing.assert_array_equal(p, id64[v])
    mutable = Counter()
    for u, v in product(range(64), repeat=2):
        tip = v
        for j, g in enumerate(generators):
            if (u ^ v) & (1 << j):
                next_tip = tip ^ (1 << j)
                assert (g if u & (1 << j) else g.T)[next_tip, tip] == 1
                tip = next_tip
        assert tip == u
        mutable[tip] += 1
    assert list(mutable.values()) == [64] * 64
    # Globally there is one vacancy on the entire connected bench, no excess.
    # Whenever both endpoints have unit load, every active physical seam
    # on this preparation has its sharp matched energy, including either coin.
    f = LAW.SRC.index((1, 1))
    filled = {}
    for name, sector, level in (
        ("blank", "blank", 0),
        ("silent", (0, 1), 0),
        ("unit_even", (1, 1), 8),
    ):
        raw = blocks()[sector]
        b = raw[f::9, f::9]
        np.testing.assert_array_equal(b, 8 * level * np.eye(len(b)))
        mask = np.tile(np.arange(9) == f, len(b))
        assert not raw[np.ix_(~mask, mask)].any()
        filled[name] = level
    return dict(
        dimension=h.shape[0],
        twice_spectral_interval=[lower, upper],
        filled_seam_energy=filled,
        upper_mass_ratio="27/8",
        koide_upper="3/8",
        corner_numerators=vertices,
        connected_blank_ports=len(graph),
        blank_hops=graph.number_of_edges(),
        marked_head_bits=len(generators),
        head_matrix_units=len(mutable) * 64,
        common_commutant_dimension=1,
    )
