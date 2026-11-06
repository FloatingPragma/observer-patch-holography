"""Independent algebra and fail-closed controls for phase selection."""
import json

import pytest
import sympy as s

from . import selection as p


def test_complete_native_law_and_twenty_rotations():
    result = p.source_certificate()
    assert result['event_factors'] == {'p': '1/12', 'q': '1/6', 'joint': '1/72'}
    assert result['all_raw_labels_retained'] == 128
    assert result['regular_cyclic_actions_checked'] == 20
    assert not result['native_kernel_commutes_with_any_tested_C3']
    for row in result['cyclic_actions']:
        assert s.sympify(row['commutator_norm_squared']).is_positive


def test_exact_classical_instruments_do_not_select_coherent_spectrum():
    result = p.phase_lift_certificate()
    assert result['raw_complete_record_instruments_identical']
    assert result['distinct_spectra']
    assert result['positive_square_root_operators']


def test_phase_invariance_on_an_entangled_input():
    original, changed = (p.coherent_lift(z) for z in (1, (99+20*s.I)/101))
    # A rank-one maximally entangled state checks retained system/reference blocks.
    v = s.zeros(9, 1)
    for i in range(3):
        v[4*i] = 1/s.sqrt(3)
    rho = v*v.H
    for K, L in zip(original['raw'], changed['raw'], strict=True):
        extended_K, extended_L = s.kronecker_product(K, s.eye(3)), s.kronecker_product(L, s.eye(3))
        assert p.clean(extended_K*rho*extended_K.H-extended_L*rho*extended_L.H) == s.zeros(9)


def test_erasure_outcomes_are_not_the_original_classical_outcomes():
    first, second = (p.coherent_lift(z) for z in (1, (99+20*s.I)/101))
    assert p.superoperator(first['erased'][0]) != p.superoperator(second['erased'][0])
    # Discarding failures would make the claimed unconditional channel false.
    gram = p.clean(sum((K.H*K for K in first['erased'][:3]), s.zeros(3)))
    assert gram == s.eye(3)/3
    assert gram != s.eye(3)


def test_general_symmetric_kernel_cyclic_twirl():
    a, b, c, d, e, f = s.symbols('a b c d e f', real=True)
    X = s.Matrix([[a, d, e], [d, b, f], [e, f, c]])
    twirl, singlet, doublet = p.real_cyclic_twirl(X, p.cyclic_rotation())
    expected = s.Matrix([[a+b+c, d+e+f, d+e+f],
                         [d+e+f, a+b+c, d+e+f],
                         [d+e+f, d+e+f, a+b+c]])/3
    assert p.clean(twirl-expected) == s.zeros(3)
    assert s.simplify(singlet-doublet-d-e-f) == 0


def test_complex_hermitian_countercontrol():
    C = p.coherent_lift((99+20*s.I)/101)['C']
    assert C == C.H
    assert C != C.T
    R = p.cyclic_rotation()
    assert p.clean(C*R-R*C) == s.zeros(3)
    assert s.factor(C.charpoly().discriminant()) != 0
    with pytest.raises(ValueError, match='real symmetric'):
        p.real_cyclic_twirl(C, R)


@pytest.mark.parametrize('phase', [0, 2, 1.0, s.Float('.6')+s.Float('.8')*s.I,
                                 s.Rational(3, 4)+s.I/4, s.Symbol('z')])
def test_nonunit_or_uncertified_phase_rejected(phase):
    with pytest.raises(ValueError, match='unit complex'):
        p.coherent_lift(phase)


def test_nonregular_rotation_rejected():
    with pytest.raises(ValueError, match='regular order-three'):
        p.real_cyclic_twirl(s.eye(3), s.eye(3))


def test_receipt_replays_without_experimental_inputs():
    stored = json.loads((p.HERE/'results.json').read_text())
    assert stored == p.receipt()
    assert stored['empirical_inputs'] == []
    assert stored['measured_spectrum_claimed'] is False
