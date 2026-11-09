from fractions import Fraction as F
import pytest
import sympy as s
from model.particle_atlas_horizontal import measurement as m


def test_same_preparation_two_complete_analyzers():
    c = m.certificate()
    assert c['position_probabilities'] == ['1/2', '1/2', '0', '0']
    assert c['momentum_probabilities'] == ['1/2', '1/4', '0', '1/4']
    assert c['after_position_x0']['position_probabilities'] == ['1', '0', '0', '0']
    assert c['after_position_x0']['momentum_probabilities'] == ['1/4']*4
    assert c['after_momentum_k0']['momentum_probabilities'] == ['1', '0', '0', '0']
    assert c['after_momentum_k0']['position_probabilities'] == ['1/4']*4


def test_every_history_and_zero_weight_branch_is_retained():
    c = m.certificate()
    for name in ('position_then_momentum_records', 'momentum_then_position_records', 'public_copies'):
        rows = c[name]
        assert len(rows) == 16
        assert sum(F(x['probability']) for x in rows) == 1
        assert any(F(x['probability']) == 0 for x in rows)
    assert len(c['all_conditional_branches']) == 8
    for row in c['all_conditional_branches']:
        assert row['normalized_conditional_state_exists'] == (F(row['probability']) > 0)
        if row['zero_weight_record_retained']:
            assert row['conditional_other_distribution'] is None


def test_phase_changes_momentum_without_changing_position():
    f, _, px, pk = m.construction()
    minus = s.Matrix([1, -1, 0, 0])/s.sqrt(2)
    rho = minus*minus.H
    assert m.outcome_probabilities(rho, px) == (s.Rational(1,2), s.Rational(1,2), 0, 0)
    assert m.outcome_probabilities(rho, pk) == (0, s.Rational(1,4), s.Rational(1,2), s.Rational(1,4))
    assert m.clean(f.H*f) == s.eye(4)


def test_unread_record_changes_future_interference_and_copying_does_not_choose_it():
    c = m.certificate()
    assert c['unread_position_then_momentum'] == ['1/4']*4
    assert c['unread_position_then_momentum'] != c['momentum_probabilities']
    assert all(F(x['probability']) == 0 for x in c['public_copies']
               if x['observer_A'] != x['observer_B'])
    assert c['record_sharing_changes_Born_weights'] is False


def test_mixed_state_and_spectator_reference_are_not_discarded():
    _, _, px, pk = m.construction()
    # Maximally entangled reference: each retained conditional branch updates
    # the reference, while the complete channel preserves its marginal.
    omega = s.Matrix([int(i//4 == i%4) for i in range(16)])/2
    joint = omega*omega.H
    out = s.zeros(16)
    for k in pk:
        kr = s.kronecker_product(k, s.eye(4))
        out += kr*joint*kr.H
    reference = s.Matrix(4, 4, lambda r, t: sum(out[4*j+r, 4*j+t] for j in range(4)))
    assert m.clean(reference) == s.eye(4)/4
    assert m.outcome_probabilities(s.eye(4)/4, px) == (s.Rational(1,4),)*4


def test_finite_bound_does_not_masquerade_as_canonical_ccr():
    c = m.certificate()
    assert F(c['finite_uncertainty']['variance_product']) >= F(c['finite_uncertainty']['Robertson_rhs_squared'])
    assert c['finite_uncertainty']['trace_commutator'] == '0'
    zero = c['finite_position_eigenstate_uncertainty']
    assert zero['variance_position_bins'] == '0'
    assert zero['Robertson_rhs_squared'] == '0'
    assert zero['canonical_CCR_claimed'] is False
    assert c['physical_lengths_momenta_or_hbar_calibrated'] is False
    assert c['native_physical_detector_derived'] is False


def test_double_path_marker_changes_interference_without_a_human_reader():
    c = m.double_path_certificate()
    assert c['unmarked']['screen_probabilities'] == ['1/2', '1/4', '0', '1/4']
    assert c['orthogonal_path_tags']['screen_probabilities'] == ['1/4']*4
    assert c['unmarked_phase_pi_screen_probabilities'] == ['0', '1/4', '1/2', '1/4']
    records = c['orthogonal_path_tags']['path_screen_records']
    assert len(records) == 8
    assert all(F(row['probability']) == F(1, 8) for row in records)
    assert sum(F(row['probability']) for row in records) == 1
    assert c['orthogonal_path_tags']['unread_path_record_suffices']
    assert not c['orthogonal_path_tags']['human_observation_required']
    assert not c['postselection_used']
    assert not c['source_selected_apparatus_or_particle_identity']


def test_double_path_isometry_and_partial_trace_retain_complete_inputs():
    for overlap in (0, 1, s.Rational(1, 2), (1+s.I)/2):
        c = m.double_path_construction(0, overlap)
        assert m.clean(c['marker_isometry'].H*c['marker_isometry']) == s.eye(4)
        assert len(c['joint_marker_screen_records']) == 8
        rho = c['path_density_after_unread_marker']
        assert s.trace(rho) == 1
        assert rho[0, 1] == s.conjugate(overlap)/2
        assert m.clean(rho[1, 0]-s.sympify(overlap)/2) == 0
        f, _, _, _ = m.construction()
        screen = m.clean(f*rho*f.H)
        assert tuple(screen[k, k] for k in range(4)) == c['screen_probabilities']
    assert m.double_path_construction(0, 0)['path_density_after_unread_marker'] == s.diag(s.Rational(1, 2), s.Rational(1, 2), 0, 0)


def test_double_path_complex_overlap_fixes_the_fringe_phase_convention():
    c = m.double_path_construction(0, s.I/2)
    assert c['screen_probabilities'] == (s.Rational(1, 4), s.Rational(3, 8), s.Rational(1, 4), s.Rational(1, 8))
    assert c['phase_scan_visibility'] == s.Rational(1, 2)
    zero_phase = m.double_path_construction(0, s.Rational(1, 2))
    opposite_phase = m.double_path_construction(s.pi, s.Rational(1, 2))
    hi, lo = zero_phase['screen_probabilities'][0], opposite_phase['screen_probabilities'][0]
    assert (hi-lo)/(hi+lo) == s.Rational(1, 2)
    # A common phase changes both marker-overlap and path conventions only
    # through gamma exp(i phi); swapping its sign would move these fringes.
    assert m.double_path_construction(s.pi/2, s.Rational(1, 2))['screen_probabilities'] == c['screen_probabilities']


def test_double_path_rejects_nonphysical_or_inexact_marker_inputs():
    for overlap in (2, 1+s.I, 0.5):
        with pytest.raises(ValueError):
            m.double_path_construction(0, overlap)
    for phase in (s.I, 0.1):
        with pytest.raises(ValueError):
            m.double_path_construction(phase, 1)
