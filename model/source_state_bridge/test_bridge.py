"""Exact native state bounds and explicit controls on their hypotheses."""
from fractions import Fraction as F
import json

import pytest
import sympy as s

from . import bridge as b


def test_actual_native_reference_and_modular_spectrum():
    source = b.native_spectrum()
    assert source['normalizer'] == 7505
    assert sum(source['probabilities']) == 1
    assert len(source['levels']) == 66
    assert sum(source['gaps'].values()) == 4356
    assert source['gaps'][0] == 2426
    # Check the annihilating polynomial on every native matrix unit.
    assert all(k*(k*k-1)*(k*k-4) == 0
               for a in source['levels'] for c in source['levels'] for k in (a-c,))


def test_isometric_pure_and_subspace_bounds_are_exact():
    assert b.support_state_distance_bound(1) == F(1376, 1501)
    assert b.support_state_distance_bound(20) == F(616, 1501)
    assert b.support_state_distance_bound(21) == F(576, 1501)
    assert b.support_state_distance_bound(66) == b.support_state_distance_bound(80) == 0
    source = sorted(b.native_spectrum()['probabilities'], reverse=True)
    # The top eigenvector attains the pure-target trace bound.
    assert (abs(source[0]-1)+sum(source[1:]))/2 == b.support_state_distance_bound(1)


@pytest.mark.parametrize('dimension,expected', [(438, F(62, 73)), (1596, F(255, 266)), (66, F()), (1, F())])
def test_complete_carrier_bounds(dimension, expected):
    assert b.carrier_channel_distance_bound(dimension) == expected


@pytest.mark.parametrize('dimension,carrier', [(2, 1), (3, 2), (4, 2), (3, 3)])
def test_independent_flagged_encode_decode_is_complete(dimension, carrier):
    rows = b.exact_flagged_code(dimension, carrier)
    assert sum((e.T*e for e, _ in rows), s.zeros(dimension)) == s.eye(dimension)
    effective = []
    for encode, decoders in rows:
        assert sum((d.T*d for d in decoders), s.zeros(carrier)) == s.eye(carrier)
        effective.extend(d*encode for d in decoders)
    assert sum((k.T*k for k in effective), s.zeros(dimension)) == s.eye(dimension)
    # Direct Choi/Bell overlap, independently of the bound function.
    bell = s.zeros(dimension**2, 1)
    for j in range(dimension): bell[j*dimension+j] = 1/s.sqrt(dimension)
    choi = s.zeros(dimension**2)
    for k in effective:
        v = s.kronecker_product(k, s.eye(dimension))*bell
        choi += v*v.T
    fidelity = (bell.T*choi*bell)[0]
    assert s.trace(choi) == 1
    assert fidelity <= s.Rational(carrier, dimension)
    if dimension % carrier == 0:
        assert fidelity == s.Rational(carrier, dimension)
    # The flag carries classical distinctions even with a one-state carrier.
    for j in range(dimension):
        state = s.zeros(dimension); state[j, j] = 1
        assert sum((k*state*k.T for k in effective), s.zeros(dimension)) == state


def test_quantum_residual_memory_can_leave_field_pure():
    rho = s.diag(s.Rational(3, 4), s.Rational(1, 4))
    # V|i> = |0>_field |i>_memory; the mixedness remains accessible in memory.
    v = s.Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])
    assert v.T*v == s.eye(2)
    joint = v*rho*v.T
    field = s.Matrix(2, 2, lambda i, j: sum(joint[2*i+r, 2*j+r] for r in range(2)))
    assert field == s.diag(1, 0)
    assert joint.rank() == 2
    assert joint != s.diag(1, 0, 0, 0)
    # A reset likewise prepares a pure field after its environment is traced.
    kraus = (s.Matrix([[1, 0], [0, 0]]), s.Matrix([[0, 1], [0, 0]]))
    assert sum((k.T*k for k in kraus), s.zeros(2)) == s.eye(2)
    assert sum((k*rho*k.T for k in kraus), s.zeros(2)) == field


@pytest.mark.parametrize('bad', [0, -1, True, 2.0, '2'])
def test_exact_domain_rejects_invalid_dimensions(bad):
    with pytest.raises(ValueError): b.support_state_distance_bound(bad)
    with pytest.raises(ValueError): b.carrier_channel_distance_bound(bad)
    with pytest.raises(ValueError): b.carrier_channel_distance_bound(3, bad)


def test_certificate_roundtrip_and_boundaries():
    result = b.certificate()
    assert json.loads(json.dumps(result)) == result
    assert result['all_input_channel_half_diamond_lower'] == {'438': '62/73', '1596': '255/266'}
    assert not result['scope']['source_selected_interacting_Higgs_action_or_state_derived']
    with pytest.raises(ValueError): b.exact_flagged_code(17, 2)
