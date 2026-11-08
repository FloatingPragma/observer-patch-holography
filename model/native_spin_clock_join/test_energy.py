"""Exact current orientation, torque and locality controls."""
import pytest
import sympy as sp

from . import energy


def test_complete_current_checks():
    result = energy.exact_checks()
    assert len(result['finite_graphs']) == 4
    assert result['oriented_bond_chirality_exact']
    assert not result['continuum_stress_or_momentum_tensor_constructed']


def test_staggered_field_changes_density_even_when_current_agrees():
    edges = ((0, 1), (1, 2))
    zero = energy.source(3, edges)
    selected = energy.source(3, edges, field=sp.Rational(3, 8), staggered=(1, -1, 1))
    # Opposite field signs cancel in the pair current. This does not
    # identify the two Hamiltonians or their local energy densities.
    assert zero['currents'] == selected['currents']
    assert zero['density'][0] != selected['density'][0]
    assert energy.commutator(zero['hamiltonian'], selected['hamiltonian']) != sp.zeros(8)
    assert energy.commutator(selected['hamiltonian'] - zero['hamiltonian'],
                             zero['density'][0]) != sp.zeros(8)


def test_distant_currents_vanish():
    data = energy.source(4, ((0, 1), (1, 2), (2, 3)))
    assert data['currents'][0, 3] == sp.zeros(16)
    assert data['currents'][0, 2] != sp.zeros(16)


@pytest.mark.parametrize('size,edges,options', [
    (True, (), {}), (7, (), {}), (2, ((0, 1), (0, 1)), {}),
    (2, ((1, 0),), {}), (2, ((0, 2),), {}),
    (2, ((0, 1),), {'coupling': 1.0}),
    (2, ((0, 1),), {'field': -1}),
    (2, ((0, 1),), {'field': 1}),
    (2, ((0, 1),), {'field': 1, 'staggered': (1, 1)}),
])
def test_invalid_source_is_rejected(size, edges, options):
    with pytest.raises(ValueError):
        energy.source(size, edges, **options)
