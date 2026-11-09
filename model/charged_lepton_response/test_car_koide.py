"""Independent occupied-sector, conjugate-frame and archive controls."""
from fractions import Fraction as F
import itertools
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import car_koide as k


def _jw(modes, index):
    factors = [s.Matrix([[0, 0], [1, 0]]) if j == index else
               s.diag(1, -1) if j < index else s.eye(2) for j in reversed(range(modes))]
    return s.kronecker_product(*factors)


def test_pairing_uses_all_occupied_states_and_two_creator_terms():
    w = s.Matrix([[s.Rational(1, 3), s.I/4], [-s.I/4, s.Rational(2, 3)]])
    cs = [_jw(4, j) for j in range(4)]
    forward = sum((w[i, j]*cs[i]*cs[j+2] for i, j in itertools.product(range(2), repeat=2)), s.zeros(16))
    expected = forward+forward.H
    assert k._pairing(w) == expected
    parity = s.diag(*[(-1)**m.bit_count() for m in range(16)])
    assert expected*parity == parity*expected
    assert any(expected[15, j] != 0 for j in range(16))
    number = sum((a*a.H for a in cs), s.zeros(16))
    assert number*expected != expected*number


def test_full_pairing_norm_has_spectator_multiplicity():
    g = k._pairing(s.diag(1, 1))
    assert g.eigenvals() == {-2: 1, -1: 4, 0: 6, 1: 4, 2: 1}
    assert max(abs(x) for x in g.eigenvals()) == 2


def test_pairing_requires_conjugate_leg_frames_not_two_equal_frames():
    delta = s.Matrix([[0, 1], [-1, 0]])
    correct = s.diag(s.I, -s.I)
    wrong = s.I*s.eye(2)
    assert correct.H*delta*correct.conjugate() == delta
    assert wrong.H*delta*wrong.conjugate() == -delta


@pytest.fixture(scope='module')
def finite():
    return k.finite_control()


def test_complex_full_CAR_mean_is_not_a_scalar(finite):
    y = finite['family_coefficient']
    assert y == s.Matrix([[s.Rational(197, 250), s.Rational(24, 125)*(-1+s.I)],
                          [s.Rational(24, 125)*(-1-s.I), s.Rational(1, 2)]])
    assert y == y.H and y[0, 1].as_real_imag()[1] != 0
    assert y.det() > 0


def test_sixteen_original_ancillary_outputs_are_complete(finite):
    ks = finite['kraus']
    assert len(ks) == 16
    assert sum((a.H*a for a in ks), s.zeros(16)).applyfunc(s.expand) == s.eye(16)
    probabilities = [s.expand(sum(s.conjugate(x)*x for x in a)/16) for a in ks]
    assert sum(probabilities) == 1 and all(p >= 0 for p in probabilities)
    assert probabilities[0] < 1
    assert any(probabilities[j] > 0 for j in range(1, 16))


def test_all_vacuum_block_and_complete_reference_archive_bound(finite):
    h = 2*math.atan(1/64)
    H = np.asarray(finite['physical_pairing'], complex)
    ks = [np.asarray(a, complex) for a in finite['kraus']]
    ideal = expm(-2j*h*H)
    a = ks[0]@ks[0]
    # Maximally entangled physical/reference state, sixteen dimensions each.
    v = a.reshape(256)/4
    w = ideal.reshape(256)/4
    norm = np.vdot(v, v).real
    overlap = abs(np.vdot(w, v))**2
    zero_error = math.sqrt(max(0., (1+norm)**2-4*overlap))
    complete_error = zero_error+1-norm
    total = sum(np.linalg.norm(b@c, 'fro')**2/16 for b in ks for c in ks)
    assert abs(total-1) < 1e-13
    assert 0 < 1-norm <= 2*4*h*h
    assert 0 < complete_error <= 3*4*(2*h)**2/2
    assert np.linalg.norm(ks[0]-expm(-1j*h*H), 2) <= 4*h*h


def test_fine_records_and_coherent_quantum_archive_are_different(finite):
    # Nonzero quantum amplitudes survive. Merely omitting their records
    # gives a strict contraction instead of a complete channel.
    ks = finite['kraus']
    vacuum_weight = s.expand(sum(s.conjugate(x)*x for x in ks[0])/16)
    assert vacuum_weight < 1
    assert any(a != s.zeros(16) for a in ks[1:])
    # Trace of the tracial output, computed without constructing irrelevant
    # off-diagonal entries of sixteen symbolic products.
    total_weight = s.expand(sum(s.conjugate(x)*x for a in ks for x in a)/16)
    assert total_weight == 1


@pytest.mark.parametrize('phase', [s.S.One, (s.Integer(99)+20*s.I)/101])
def test_actual_analyzer_compression_is_exact_existing_source(phase):
    data = k.source_lift(phase)
    row = data['analyzer_flag_zero_row']; y = data['physical_family_coefficient']
    assert row.shape == (3, 24)
    assert k._clean(row*row.H) == s.eye(3)
    assert k._clean(row[:, :3].H*row[:, :3]) == y
    assert y == k.source.exact_source(phase)['action_coefficient']
    assert s.trace(y) == s.Rational(5, 18)


def test_actual_source_supports_preserve_full_doublet_and_singlet():
    row = k.higgs_covariance()
    assert row['left_internal_rank'] == 2 and row['right_internal_rank'] == 1
    assert row['source_current_support_checks'] == 24
    assert row['both_Higgs_components_checked'] == 2
    assert row['fixed_neutral_Higgs_commutes_all_twelve_currents'] is False


@pytest.mark.parametrize('phase', [1.0, 0, 2])
def test_source_phase_admission_is_not_bypassed_by_cached_equal_values(phase):
    k.source_lift(1)
    with pytest.raises(ValueError):
        k.source_lift(phase)


def test_exact_720_mode_bound_and_implementation_error():
    row = k.archive_bound(65536, 1)
    assert F(row['pairing_operator_norm']) == 6
    assert F(row['ideal_any_nonvacuum_probability_upper']) == F(36, 65536)
    assert F(row['complete_unhalved_output_error_upper']) == F(108, 65536)
    eps = F(1, 10**12)
    changed = k.archive_bound(65536, 1, round_dilation_error=eps)
    assert F(changed['complete_unhalved_output_error_upper']) == F(108, 65536)+2*65536*eps
    assert F(changed['implemented_any_nonvacuum_probability_upper']) == F(36, 65536)+65536*eps
    assert not row['vacuum_record_postselection']


@pytest.mark.parametrize('args', [(0, 1), (True, 1), (10, -1), (10, .1)])
def test_invalid_clock_inputs_rejected(args):
    with pytest.raises(ValueError):
        k.archive_bound(*args)


def test_bound_is_reference_safe_not_a_ground_or_native_claim():
    row = k.archive_bound(16, F(1, 3), higgs_norm=F(2, 3))
    assert row['arbitrary_physical_occupation_and_reference_inputs']
    assert row['all_nonvacuum_quantum_archives_retained']
    assert F(row['pairing_operator_norm']) == 4


def test_certificate_keeps_exact_original_weights_and_scope():
    c = k.certificate()
    json.dumps(c, allow_nan=False)
    assert c['physical_target_modes']+c['fresh_ancillary_modes_per_round'] == c['enlarged_modes'] == 720
    assert sum(map(F, c['finite_control']['original_tracial_probabilities_exact'])) == 1
    for flag in ('nonvacuum_quantum_archives_discarded',
                 'coherent_ancillary_archive_identified_with_classical_flags',
                 'native_full_720_mode_tape_executed',
                 'source_gate_overhead_equals_target_pulse_clock',
                 'source_global_physical_clock_selected',
                 'phase_or_physical_pole_masses_derived', 'same_as_four_current_gated_actions'):
        assert c[flag] is False
    assert c['physical_particle_assignments'] == []
    c['source_controls'][0]['Y_exact'][0][0] = 'corrupt'
    assert k.certificate()['source_controls'][0]['Y_exact'][0][0] != 'corrupt'
