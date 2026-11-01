from fractions import Fraction as F
import math
import numpy as np

from . import experiment as x,circuits as c


def test_full_native_clock_and_all_outcome_accounting(native_packet):
    value,arrays=native_packet;result=value['executed_native_detector']
    assert result['invariant_ground_preparation_executed']
    assert result['complete_binary_phase_clock_source_executed']
    assert arrays['phase_clock_complete_success_state'].shape==(262144,30)
    assert len(result['all_selector_pointer_outcomes'])==1296
    assert len(result['all_joint_selector_charge_clock_outcomes'])==3240
    assert len(result['all_charge_clock_marginal_outcomes'])==216
    assert abs(result['complete_probability']-1)<2e-8
    assert result['first_failure_and_invalid_selector_stop_mass']>=0
    assert not result['success_charge_or_bin_postselection_used']


def test_actual_native_results_obey_the_independent_full_H_and_interval_bounds(native_packet):
    value,_=native_packet;result=value['executed_native_detector']
    bound=float(F(value['complete_native_error_certificate']['complete_state_and_total_variation_error']))
    assert result['complete_output_error_diagnostic']<bound
    assert result['total_variation_error_diagnostic']<bound
    for tag,target in ((0,.75),(1,.25)):
        mass=math.fsum(row['probability'] for row in result['all_charge_clock_marginal_outcomes'] if row['charge_record']==tag)
        assert abs(mass-target)<bound
    assert not result['all78_lines_separately_resolved_by_native_bins']


def test_discarding_the_ground_energy_phase_is_detectably_wrong(native_packet):
    value,_=native_packet
    ideal=x.ideal_clock().reshape(c.M,65536,30)
    fourier=np.exp(2j*np.pi*np.arange(c.M)[:,None]*np.arange(c.M)[None,:]/c.M)/math.sqrt(c.M)
    wrong=fourier@np.diag(np.exp(1j*np.arange(c.M)*math.sqrt(3)/4))@fourier.conj().T
    wrong_state=np.tensordot(wrong,ideal,axes=(1,0)).reshape(c.M*65536,30)
    wrong_prob=x.pointer_probabilities(wrong_state)
    right=value['executed_native_detector']['independent_original_H_marginal_outcomes']
    variation=math.fsum(abs(a['probability']-b['probability']) for a,b in zip(wrong_prob,right,strict=True))/2
    bound=float(F(value['complete_native_error_certificate']['complete_state_and_total_variation_error']))
    assert variation>2*bound
