from fractions import Fraction as F

import pytest
import sympy as sp

from . import contact as c, leakage as l


def test_every_boundary_triplet_has_exact_cubic_leakage():
    result = l.certificate()
    assert result['distinct_boundary_triplets'] == 220
    assert result['positive_principal_minors_per_triplet'] == 15
    assert result['outside_action_norm_squared_lower'] == F(1, 100)
    assert sum(row['cells'] for row in result['evaluated_federations']) == 420


def test_outside_gram_against_actual_connected_source_action():
    graph = c.source.source.federation_family(0, (1,)*12)
    n = l.cubic_basis(); pairs = c.data()['pairs']
    output = []
    for column in range(n.cols):
        initial = {pair: int(n[i, column]) for i, pair in enumerate(pairs) if n[i, column]}
        # Off-cell output of G is minus the original transfer counts.
        counts = c.source.counts_step(graph, initial, 2)
        output.append({pair: -value for pair, value in counts.items() if pair[1] >= 12 and value})
    outside = tuple(sorted({pair for col in output for pair in col}))
    actual = sp.Matrix([[col.get(pair, 0) for col in output] for pair in outside])
    ports = tuple(sorted(a for a, b in graph.seams if a < 12 <= b))
    grams, norm = l.port_grams()
    expected = sp.Matrix([[sum(grams[p][i][j] for p in ports) for j in range(15)] for i in range(15)])
    assert actual.T*actual == expected
    assert all(x > 0 for x in l.leading_minors(l.boundary_matrix(ports)))


def test_sylvester_control_rejects_an_overstated_leakage_bound():
    assert min(l.leading_minors(l.boundary_matrix((0, 1, 2), F(2)))) <= 0
    assert l.leading_minors(((2, 1), (1, 2))) == (2, 3)
    with pytest.raises(ValueError):
        l.leading_minors(((2, 0), (1, 2)))


@pytest.mark.parametrize('ports', [(False, 1, 2), (0, 0, 2), (2, 1, 0), (0, 1, 12), [0, 1, 2]])
def test_boundary_ports_are_typed_and_canonical(ports):
    with pytest.raises(ValueError):
        l.boundary_matrix(ports)


def test_fixed_count_endpoint_lower_bound_and_zero_case():
    t = F(1, 5760)
    assert l.endpoint_lower(t) == F(1, 115200)
    assert l.endpoint_lower(t, 1000000) > 0
    assert l.endpoint_lower(t, 1) == 0
    assert l.endpoint_lower(F()) == 0
    with pytest.raises(ValueError):
        l.endpoint_lower(t, False)
