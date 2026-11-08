"""Source-specific and independent small-graph completeness checks."""
from itertools import combinations

import numpy as np
import pytest

from . import pauli


def test_actual_source_complete_census():
    result = pauli.census()
    assert result['selected_field_charge_words_modulo_phase'] == [[0, 0], [255, 0]]
    assert result['zero_field_charge_words_modulo_phase'] == [[0, 0], [0, 255], [255, 0], [255, 255]]
    assert result['nonidentity_strictly_local_Pauli_charge_count'] == 0


@pytest.mark.parametrize('sites,edges,fields', [
    (1, [], []), (1, [], [0]), (2, [(0, 1)], []),
    (2, [(0, 1)], [1]), (3, [(0, 1)], [0]),
    (3, [(0, 1), (1, 2)], []), (4, [(0, 1), (2, 3)], [0]),
])
def test_independent_matrix_completeness(sites, edges, fields):
    # Build an ordinary Heisenberg action in an independent computational
    # frame. Distinct rational edge strengths prevent accidental edge swaps.
    h = np.zeros((1 << sites, 1 << sites), complex)
    for edge_number, (a, b) in enumerate(edges, 1):
        support = (1 << a) | (1 << b)
        for word in ((support, 0), (0, support), (support, support)):
            h += edge_number * pauli.word_matrix(sites, *word)
    for a in fields:
        h += (a + 1) * pauli.word_matrix(sites, 1 << a, 0)
    actual = []
    for x in range(1 << sites):
        for z in range(1 << sites):
            p = pauli.word_matrix(sites, x, z)
            if np.array_equal(p @ h, h @ p):
                actual.append((x, z))
    assert tuple(actual) == pauli.charges(sites, edges, fields)


def test_every_three_vertex_graph_and_field_support():
    possible = list(combinations(range(3), 2))
    for graph in range(8):
        edges = [e for i, e in enumerate(possible) if graph >> i & 1]
        for field in range(8):
            fields = [i for i in range(3) if field >> i & 1]
            assert pauli.charges(3, edges, fields) == pauli.exhaustive(3, edges, fields)


def test_local_creation_changes_only_the_global_charge():
    global_x = pauli.word_matrix(8, 255, 0)
    local_z = pauli.word_matrix(8, 0, 1)
    assert np.array_equal(global_x @ local_z, -local_z @ global_x)
    assert (0, 1) not in pauli.charges(8, pauli.source.cube()[0], range(8))


@pytest.mark.parametrize('sites,edges,fields', [
    (0, [], []), (2, [(0, 0)], []), (2, [(0, 2)], []),
    (2, [(0, 1), (1, 0)], []), (2, [], [2]), (2, [], [0, 0]),
])
def test_invalid_domains_rejected(sites, edges, fields):
    with pytest.raises(ValueError):
        pauli.charges(sites, edges, fields)
