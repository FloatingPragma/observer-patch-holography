"""All-charge positivity, full-source perturbation and hostile controls."""
import pytest
import sympy as s
from . import local_bound,lattice,perturbation,ground_exchange,ground_records,resolved_exchange


def test_complete_unrestricted_weight_bound_and_every_other_Casimir():
    result=local_bound.certificate();block=result['maximal_quark_module']
    assert len(block['unrestricted_weight_blocks'])==1686
    assert sum(row['complete_scalar_blank_dimension'] for row in block['unrestricted_weight_blocks'])==8960
    assert any(row['twice_hypercharge_weight']!=0 for row in block['unrestricted_weight_blocks'])
    assert result['every_other_quark_module']['largest_nonmaximal_Casimir']=='95/2'
    assert result['full_unconstrained_local_dimension']==5*2**30


def test_wrong_metric_or_lower_bound_fails_exact_positivity():
    with pytest.raises(ValueError):local_bound.positive_pivots(s.Matrix([[1,2],[2,1]]))
    with pytest.raises(ValueError):local_bound.positive_pivots(s.Matrix([[1,1],[0,1]]))
    label=max(local_bound.blank_blocks(),key=lambda label:len(local_bound.blank_blocks()[label]))
    _,g,k,h=local_bound.exact_block(label)
    with pytest.raises(ValueError):local_bound.positive_pivots(-g-k)


def test_unit_penalty_atomic_gap_and_scalable_complete_link_incidence():
    a=lattice.atomic_certificate();b=lattice.atomic_certificate(side=4)
    assert a['full_physical_atomic_gap_lower_bound']==b['full_physical_atomic_gap_lower_bound']
    assert a['exact_ground_rank']==256 and b['exact_ground_rank']==2**64
    assert s.Rational(a['full_physical_atomic_gap_lower_bound'])>s.Rational(4,100)
    with pytest.raises(ValueError):lattice.gap_data(kappa=s.Rational(1,2))


def test_full_hopping_norm_zero_first_order_and_admitted_graph_errors():
    result=lattice.hopping_certificate()
    assert result['full_hopping_norm_upper_bound']=='360' and result['exact_atomic_P0_T_P0']=='0'
    p=perturbation.parameters();assert p['radius']==s.Rational(1,2) and p['contraction']==s.Rational(1,4)
    with pytest.raises(ValueError):perturbation.parameters(coupling=1)
    with pytest.raises(ValueError):perturbation.parameters(coupling=0.0)
    a=perturbation.iteration_bounds(4);b=perturbation.iteration_bounds(8)
    assert s.Rational(b['full_H_residual_operator_upper_bound'])==s.Rational(a['full_H_residual_operator_upper_bound'])/256
    assert perturbation.iteration_bounds(0,coupling=0)['full_H_residual_operator_upper_bound']=='0'


def test_actual_ground_CAR_second_moment_and_complete_unit_quench_records():
    exchange=ground_exchange.second_moment_certificate()
    assert len(exchange['axis_checks'])==3
    assert s.Rational(exchange['actual_local_ground_spin_coefficients']['beta_ec_interval'][1])<0
    assert s.Rational(exchange['z_edge_antialigned_minus_aligned_second_moment_interval'][0])>s.Rational(1,24)
    records=ground_records.certificate()
    assert len(records['all_eight_reference_and_detector_outcomes'])==4
    for row in records['all_eight_reference_and_detector_outcomes']:
        for name in ('success','failure'):
            lo,hi=map(s.Rational,row['original_reference_and_link_'+name+'_probability_interval'])
            assert 0<lo<hi<s.Rational(1,4)
    lo,hi=map(s.Rational,records['complete_probability_enclosure']);assert lo<=1<=hi
    assert s.Rational(records['antialigned_minus_aligned_original_joint_success_lower_bound'])>0


def test_actual_full_denominator_exchange_and_nonconstant_controlled_band():
    result=resolved_exchange.certificate()
    lo,hi=map(s.Rational,result['certified_J_edge_interval'])
    assert 0<lo<hi
    assert result['denominator_is_evaluated_at_actual_atomic_ground_energy']
    assert result['tensor_and_source_checks']['actual_CAR_charge_conjugation_mode_identities']==30
    assert s.Rational(result['exact_band_energy_expectation_contrast_lower_bound'])>0
    assert s.Rational(result['admitted_positive_lambda_for_nonconstant_exact_band'])<perturbation.parameters()['lambda_max']
