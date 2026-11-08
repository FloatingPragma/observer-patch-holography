from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
import sympy as sp

from . import contact as c


def graph(edges, ports=4):
    return c.source.source.Family((1,)*ports, edges, (1,)*len(edges))


def test_source_compression_independently_on_full_tensor_square():
    d = c.data(); pairs = d['pairs']; one = np.array(d['one'], dtype=np.int64)
    insertion = np.zeros((144, 66), dtype=np.int64)
    for k, (i, j) in enumerate(pairs):
        insertion[12*i+j, k] = insertion[12*j+i, k] = 1
    assert np.array_equal(insertion.T@insertion, 2*np.eye(66, dtype=np.int64))
    independent = np.kron(one, np.eye(12, dtype=np.int64))+np.kron(np.eye(12, dtype=np.int64), one)
    compressed_twice = insertion.T@independent@insertion
    expected = 2*np.array(d['two']+2*d['adjacent'], dtype=np.int64)
    assert np.array_equal(compressed_twice, expected)
    assert not np.array_equal(compressed_twice, 2*np.array(d['two'], dtype=np.int64))
    assert d['incidence']*d['probe'] == sp.zeros(12, 1)


def test_complete_native_generator_against_all_actual_one_attempt_columns():
    from model.native_program import interpreter
    value = interpreter.program(biased=False)
    native_pairs = [tuple(i for i, x in enumerate(row) if x) for row in c.source.source.local.basis()]
    where = {p: i for i, p in enumerate(c.data()['pairs'])}
    marker = {(i, i, i, i): F(1) for i in range(66)}
    actual = sp.zeros(66)
    output = interpreter._outcomes(value, value['commands'][0], marker)
    assert len(output) == 721
    for block in output:
        for (i, a, j, b), weight in block.density.items():
            assert i == j and a == b
            actual[where[native_pairs[i]], where[native_pairs[a]]] += sp.Rational(weight.numerator, weight.denominator)
    assert actual == sp.eye(66)-c.data()['two']/60


def test_full_spectral_and_minimal_response_identity():
    result = c.spectral_data()
    assert result['primitive_characteristic'].degree() == 54
    assert result['denominator'].degree() == 12
    assert result['numerator'].degree() == 11
    assert result['moments'][:4] == (F(1), F(12), F(148), F(1860))
    q = c.data()['primitive']; h = c.data()['probe']; g = c.data()['two']
    for z in (F(-1), F(0), F(6), F(15)):
        # Invert only the primitive block, retaining an identity on its complement.
        matrix = sp.Rational(z.numerator, z.denominator)*q-q*g*q+(sp.eye(66)-q)
        answer = (h.T*matrix.inv()*h)[0]/4
        assert c.response(z) == c.fraction(answer)


def test_positive_observed_pole_against_independent_dense_spectral_measure():
    cert = c.pole_certificate(); lower, upper = cert['pole_interval']
    assert F(5) < lower < upper < F(51, 10)
    assert upper-lower < F(1, 10**20)
    assert cert['subthreshold_separation_lower'] == F(2, 5)
    values, vectors = np.linalg.eigh(np.array(c.data()['two'], dtype=float))
    h = np.array(c.data()['probe'], dtype=float).ravel()/2
    selected = (values > 5) & (values < 5.1)
    assert np.count_nonzero(selected) == 5
    weight = np.linalg.norm(vectors[:, selected].T@h)**2
    assert abs(values[selected].mean()-float((lower+upper)/2)) < 1e-12
    lo, hi = cert['normalized_residue_interval']
    assert abs(weight-float((lo+hi)/2)) < 1e-13
    # This test diagnoses a complete matrix calculation; the rational proof
    # and interval, not floating diagonalization, certify the pole.


@pytest.mark.parametrize('z', [8, 10, 12])
def test_observed_poles_are_rejected(z):
    with pytest.raises(ValueError, match='pole'):
        c.response(z)


@pytest.mark.parametrize('z', [True, 0.0, '0', sp.Rational(1, 2)])
def test_response_does_not_alias_rational_input_types(z):
    with pytest.raises(ValueError):
        c.response(z)


@pytest.mark.parametrize('edges', [((0, 1), (1, 2), (2, 3)), tuple(combinations(range(4), 2))])
@pytest.mark.parametrize('z', [F(-1), F(), F(9, 2), F(7)])
def test_contact_resolvent_including_signed_denominator(edges, z):
    result = c.contact_identity(graph(edges), z)
    assert result['primitive_dimension'] == 2
    assert result['identity_verified']


@pytest.mark.parametrize('z', [F(6), F(8)])
def test_contact_inverse_rejects_degenerate_and_reference_poles(z):
    with pytest.raises(ValueError, match='pole'):
        c.contact_identity(graph(tuple(combinations(range(4), 2))), z)


def test_full_matrix_caps_are_not_raised():
    with pytest.raises(ValueError, match='twelve'):
        c.generators(c.source.source.federation_family(0, (1,)*12))
    with pytest.raises(ValueError, match='eight'):
        c.contact_identity(c.family(), F(-1))


def test_persistence_budget_requires_proved_inputs_without_supplying_them():
    assert c.persistence_budget(F(1, 2), F(1, 100), F(1, 1000))['outside_spectral_weight_upper'] == F(1, 2500)
    with pytest.raises(ValueError):
        c.persistence_budget(F(), F(1), F())
