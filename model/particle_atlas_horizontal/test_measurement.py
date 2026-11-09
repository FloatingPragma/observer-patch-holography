from fractions import Fraction as F
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
