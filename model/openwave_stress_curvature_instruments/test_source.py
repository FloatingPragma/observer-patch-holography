from fractions import Fraction as F
import sympy as s
import pytest
from .source import operators,instrument_records,backaction_controls,acquisition_budget,tensor_certificate

def test_full_carrier_all_ten_settings_and_every_outcome():
    H,ops,rho,_,_=operators()
    assert H.shape==(60,60) and len(ops)==10 and rho.rank()==60
    rows=instrument_records()
    assert len(rows)==10 and all(len(r['all_three_probabilities_exact'])==3 for r in rows)
    assert all(sum(map(F,r['all_three_probabilities_exact']))==1 for r in rows)

def test_incompatibility_reference_and_missing_outcome_are_detected():
    c=backaction_controls()
    assert c['X_mean_after_prior_Z_instrument_exact']=='3/4'
    assert c['dropped_zero_event_total_effect_exact']=='5I/8'
    assert c['reference_joint_trace_exact']=='1'

def test_sampling_resources_are_explicit_and_precision_quadratic():
    a=acquisition_budget();b=acquisition_budget(F(1,2000))
    assert b['total_prepared_copies']>=4*a['total_prepared_copies']-30
    assert a['total_prepared_copies']>10**6
    with pytest.raises(ValueError):acquisition_budget(0)

def test_complete_tensor_decoder_retains_nonzero_residual_and_error():
    c=tensor_certificate()
    assert c['same_state_mean_tensor_decoder_identity'] and c['complete_scalar_support_count']>1
    assert s.sympify(c['actual_mean_maximum_null_residual_exact'])>0
    assert F(c['acquisition_error']['every_null_charge_error_upper_exact'])>0
    assert not c['Einstein_equation_enforced_by_acquisition']
