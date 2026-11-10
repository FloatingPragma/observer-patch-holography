from copy import deepcopy
from fractions import Fraction as F
import math
import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm
from . import apparatus,model
from .__main__ import compare


def test_complete80_state_projectors_and_inverse_decoder():
    W,f,projectors=apparatus.current_algebra()
    assert sum(projectors,s.zeros(80))==s.eye(80)
    for i,P in enumerate(projectors):
        assert P*P==P and f*P==apparatus.LABELS[i]*P
        for j,Q in enumerate(projectors):
            if i!=j:assert P*Q==s.zeros(80)
    x=s.symbols('x')
    inverse=s.Poly(s.interpolate(list(zip(apparatus.LABELS,apparatus.VALUES)),x),x)
    restored=sum((coefficient*f**power[0]for power,coefficient in inverse.terms()),s.zeros(80))
    assert restored==W


def test_rational_pointer_and_all400_coupling_columns():
    P,error,norm,coefficients=apparatus.pointer()
    assert P.H==P and norm<3 and 0<error<F(1,10**12)
    assert all((x/s.I).is_Rational for x in P)
    for q in coefficients.values():assert (q*(1<<apparatus.BITS)).denominator==1
    row=apparatus.local_execution()
    assert row['complete_local_coupling_dimension']==400
    assert row['raw_numerical_execution']['all400_columns_executed']
    assert row['raw_numerical_execution']['complete_five_Kraus_residual']<1e-12
    assert 0<F(row['rational_Taylor_operator_error_upper'])<F(1,10**35)
    assert not row['raw_numerical_execution']['numerical_residual_is_exact_error_bound']


def test_five_pointer_slots_are_bijective_and_reversal_changes_negative_label():
    _,_,projectors=apparatus.current_algebra()
    shift=s.zeros(5)
    for j in range(5):shift[(j+1)%5,j]=1
    blank=s.eye(5)[:,0]
    slots=tuple(j%5 for j in apparatus.LABELS)
    assert set(slots)==set(range(5))
    for j,slot in zip(apparatus.LABELS,slots):assert shift**(j%5)*blank==s.eye(5)[:,slot]
    assert shift**2*blank!=shift**3*blank
    assert s.trace(projectors[0])==4  # Four actual W=-1 states reject a reversed decoder.


def test_full_source_bound_includes_actual_quark_scalar_terms_and_six_paths():
    row=model.local_source_bound()
    assert row['nonzero_actual_unit_Yukawa_W_commutator_terms']=={'up':144,'down':144,'charged_lepton':40}
    assert row['touched_directed_terms']==6 and F(row['touched_coefficient_absolute_sum'])==F(3,2)
    assert F(row['touched_full_Fock_hopping_norm_upper'])==90
    assert F(row['full_source_commutator_H_fW_norm_upper'])==224+F(360,2**34)
    assert row['scalar_penalty_and_link_energy_commutators_zero']


def test_noncommuting_finite_clock_control_exercises_the_Duhamel_bound():
    # Independent2-state mathematical control; no replacement of the actual
    # eight-cell source by this diagnostic example.
    H=np.diag([1.,-1.]);G=.7*np.array([[0.,1.],[1.,0.]])
    for tau in (.1,.03):
        actual=expm(-1j*(tau*H+G));sequential=expm(-1j*tau*H)@expm(-1j*G)
        difference=np.linalg.norm(actual-sequential,2)
        bound=tau*np.linalg.norm(H@G-G@H,2)/2
        assert .001<difference<bound


def test_complete_history_bound_and_both_negative_pointer_means():
    value=model.calculate();error=F(value['complete_clock_and_error']['all_retained_history_total_variation_error_upper'])
    assert 0<error<F(27,100000)
    assert F(value['complete_clock_and_error']['conditional_current_mean_absolute_error_upper'])==4*error
    records=value['complete_prepared_records']
    assert all(F(pair[1])<F('-.0043')for pair in records['current_means'].values())
    for prefix in('actual_elementary','actual_composite'):
        intervals=[tuple(map(F,row[prefix+'_bin_conditional_probability_interval']))for row in records['five_complete_current_records']]
        assert sum(lo for lo,_ in intervals)<=1<=sum(hi for _,hi in intervals)
    minus=next(row for row in records['five_complete_current_records']if row['W_value']=='-1/2')
    for kind in('elementary','composite'):
        p=F(minus['original_'+kind+'_attempt_joint_probability_lower'])
        assert p>0 and F(minus['fair_setting_original_'+kind+'_joint_probability_lower'])==p/2
    assert records['complete_setting_pump_energy_current_alphabet_size']==80
    assert records['disabled_pointer_zero_is_not_a_W_measurement']


def test_source_clock_reference_work_and_parent_services_stay_distinct():
    value=model.calculate();clock=value['complete_clock_and_error'];resources=value['control_resources']
    assert F(clock['total_elapsed_source_time'])==3+F(1,2**20)!=4
    assert not value['source']['full_Hext_cancellation_used_for_final_readout']
    assert not resources['full_matrix_preparation_and_exact_spectral_bin_premises_removed']
    assert not resources['a_separate_work_measurement_instrument_added']
    assert resources['work_is_unconditioned_mean_bare_source_energy_accounting']
    assert F(resources['final_pulse_gross_switching_work_upper'])==2*F(resources['local_final_drive_norm_upper'])
    assert F(resources['final_pulse_net_bare_source_work_absolute_upper'])<=F(clock['full_source_drive_commutator_norm_upper'])
    assert not value['scope']['Standard_Model_species_or_GeV_mass_identified']


def test_unlisted_float_nan_scope_and_error_mutations_rejected():
    value=model.calculate()
    for mutation in('nan','large','bound','status','clock'):
        changed=deepcopy(value)
        if mutation in('nan','large'):
            changed['local_apparatus']['raw_numerical_execution']['direct_vs_exact_cyclic_operator_residual']=float('nan')if mutation=='nan'else .1
        elif mutation=='bound':changed['complete_clock_and_error']['all_retained_history_total_variation_error_upper']='0'
        elif mutation=='clock':changed['complete_clock_and_error']['total_elapsed_source_time']='4'
        else:changed['scope']['autonomous_native_repair_compiler']=True
        with pytest.raises(ValueError):compare(value,changed)


def test_nonpositive_or_inexact_duration_rejected(monkeypatch):
    try:
        for duration in(F(),F(-1),1e-6):
            monkeypatch.setattr(model,'TAU',duration);model.calculate.cache_clear()
            with pytest.raises(ValueError,match='exact rational'):model.calculate()
    finally:
        monkeypatch.setattr(model,'TAU',F(1,2**20));model.calculate.cache_clear()
