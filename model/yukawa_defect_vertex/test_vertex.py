"""Source, complete-record, scalar-charge and normalization controls."""
import sympy as s

from . import vertex as v


def test_literal_source_and_every_CAR_column():
    data = v.source_data()
    actual, currents, _, _, _ = v.source.quantum_lepton_block()
    assert data['H'] == actual
    assert data['currents'] == tuple(currents)
    for i, ci in enumerate(data['c']):
        for j, cj in enumerate(data['c']):
            assert v.zero(ci*cj.H+cj.H*ci-(v.eye(320) if i == j else s.zeros(320)))
            assert v.zero(ci*cj+cj*ci)


def test_full_two_point_records_fix_both_leg_normalizations():
    result = v.two_point_data()
    assert result['fermion_CAR_sum'] == s.eye(6)
    assert result['scalar_commutator_sum'] == s.zeros(2)
    assert result['scalar_first_moment'] == s.Rational(14, 5)*s.eye(2)
    assert all(probability == s.Rational(rank, 320)
               for _, probability, rank in v.spectral_weights(v.eye(320)/320))
    for field in ('fermion', 'scalar'):
        assert len(result[field]) == 14
        assert sum(row['rank'] for row in result[field]) == 320
        for row in result[field]:
            assert row['addition'] == row['addition'].H
            assert row['removal'] == row['removal'].H


def test_equation_of_motion_three_points_reconstruct_actual_scalar_operator():
    result = v.vertex_data()
    assert result['scalar_covariance'] == s.Rational(7, 10)*s.eye(2)
    assert result['normalized_unit_coefficient'] == s.sqrt(s.Rational(14, 5))
    assert result['all_fermion_component_pairs_checked'] == 36
    assert len(result['nonzero_components']) == 8
    negative = {(row['left'], row['right']): tuple(row['source_coefficients'])
                for row in result['nonzero_components'] if row['left'] < row['right']}
    assert negative == {(0, 4): (-1, 0), (0, 5): (0, -1),
                        (1, 3): (-1, 0), (2, 3): (0, -1)}


def test_complete_charged_probes_retain_no_operations_and_preparation_failures():
    result = v.complete_probe_records()
    assert len(result['all_preparation_records']) == 14
    assert len(result['ground_conditional_records']) == 70
    all_records = result['unconditional_preparation_and_probe_records']
    assert len(all_records) == 96
    for field in ('fermion', 'scalar'):
        records = [r for r in all_records if r['probe'] == field]
        assert s.simplify(sum(r['probability'] for r in records)) == 1
        assert s.simplify(sum(r['probability'] for r in records
                              if r['stage'] == 'ground_then_probe')) == s.Rational(1, 160)
    scalar_noop = sum(r['probability'] for r in result['ground_conditional_records']
                     if r['branch'] == 'no_operation')
    assert scalar_noop == s.Rational(3, 10)


def test_scalar_pulse_is_an_actual_gauge_compensated_intervention():
    result = v.scalar_intervention()
    assert result['total_Gauss_operator_identities'] == 24
    assert result['derivative_of_pair_expectation_at_zero_pulse'] == -s.Rational(3, 10)
    assert result['pair_expectation'] == -3*s.sqrt(2)/40
    assert result['conditional_pointer_probabilities'] == {
        'minus': s.Rational(1, 2)+3*s.sqrt(2)/80,
        'plus': s.Rational(1, 2)-3*s.sqrt(2)/80}
    assert result['scalar_cost_or_interaction_coefficient_differentiated'] is False
    assert result['quantum_reference_coherences_retained'] is True
    assert len(result['complete_preparation_and_pointer_records']) == 15
    assert s.simplify(sum(r['probability'] for r in result['complete_preparation_and_pointer_records'])) == 1
