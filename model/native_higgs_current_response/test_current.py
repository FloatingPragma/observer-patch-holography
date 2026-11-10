from fractions import Fraction as F
import sympy as s
from . import local,model
from model.native_weak_current_interface.model import dressed_weak_polynomial
from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial


def test_actual_conjugation_changes_current():
    W,*_=dressed_weak_polynomial(); Wbar=local.polynomial()
    assert Wbar!=W and conjugate_polynomial(Wbar)==W


def test_complete_atomic_projectors_and_two_outcomes():
    result=local.calculate()
    assert [x['complete_target_dimension'] for x in result['complete_target_controls']]==[7,17,7]
    values=result['atomic_addition_branch_W'];x=s.Symbol('x')
    poly=lambda name:s.Poly.from_list(list(map(s.Rational,values[name]['coefficients_in_charged_energy'])),x)
    assert poly('W2') == -poly('W')/2
    assert poly('probability_-1/2') == -2*poly('W')
    assert poly('probability_0') == 1+2*poly('W')
    for value in ('-1','3/2','2'):assert poly('probability_'+value).is_zero
    assert F(values['W']['outward_interval'][1])<0


def test_both_actual_current_means_are_negative():
    result=model.calculate()['prepared_current']
    for key in ('elementary_W_mean_interval','Higgs_composite_W_mean_interval'):
        low,high=map(F,result[key]); assert F('-.071')<low<high<F('-.0053')
    assert F(result['same_history_TV_upper'])<F('0.000004')


def test_complete_five_record_laws_keep_original_weights():
    result=model.calculate()['prepared_current']; rows=result['five_complete_current_records']
    assert len(rows)==5
    for field in ('actual_elementary_bin_conditional_probability_interval','actual_composite_bin_conditional_probability_interval'):
        assert sum(F(x[field][0]) for x in rows)<=1<=sum(F(x[field][1]) for x in rows)
    negative=next(x for x in rows if x['W_value']=='-1/2')
    for pump in ('elementary','composite'):
        p=F(negative[f'original_{pump}_attempt_joint_probability_lower'])
        assert 0<p<F(result[f'{pump}_original_attempt_bin_probability_lower'])<1
        assert F(negative[f'fair_setting_original_{pump}_joint_probability_lower'])==p/2


def test_complete_source_clock_and_complement_are_retained():
    result=model.calculate()
    assert result['source']['full_source_reference_complements_retained']
    assert result['apparatus']['complete_setting_pump_energy_current_alphabet_size']==80
    assert result['apparatus']['failures_and_outside_energy_bin_retained']
    assert result['apparatus']['current_measurement_can_change_admitted_energy']
    assert result['apparatus']['additional_zero_energy_current_pointer_capacity_nats']=='ln5'
    assert result['apparatus']['complete_register_capacity_nats_formula'].endswith('+ln5')
    assert result['apparatus']['complete_elapsed_source_time']=='4 source units'
    assert F(result['apparatus']['gross_four_constant_pulse_switching_work_upper'])==8*F(result['apparatus']['constant_control_norm_upper'])
    assert result['energy_admission']['common_elapsed_time_before_current_service']=='3 source units including preparation'
    assert not result['preparation']['exact_interacting_ground_projection_assumed']


def test_no_charge_clock_or_particle_promotion():
    result=model.calculate();conj=result['conjugation'];scope=result['scope']
    assert conj['source_time_reversed'] and not conj['old_reference_4N_particle_hole_symmetry_used']
    assert not conj['same_untransformed_current_on_opposite_charge_branches_equal']
    assert not scope['equal_unconditional_A_and_B_probabilities']
    assert not scope['weak_boson_chiral_Lorentz_interaction_or_SM_identification']
    assert not scope['GeV_rest_mass_or_measured_lepton_ratio']
