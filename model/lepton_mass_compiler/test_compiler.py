"""Native factorization, every CAR occupation and complete mass records."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import analyzer,compiler,family_map


def test_actual_complete_analyzer_and_exact_source_partition():
    result=analyzer.certificate()
    assert result['native_mask_counts']==[500,250,250]
    assert result['all_inverse_columns_checked']==result['all_projector_columns_checked']==10800
    assert result['maximum_inverse_column_error']<2e-12
    assert result['maximum_projector_column_error']<2e-12
    assert result['wrong_equal_leg_frame_full_projector_difference']>.1
    theta=s.symbols('theta',real=True)
    cs=(s.cos(theta),-s.cos(theta)/2-s.sqrt(3)*s.sin(theta)/2,-s.cos(theta)/2+s.sqrt(3)*s.sin(theta)/2)
    lambdas=[s.sqrt(s.Rational(5,108))+2*s.sqrt(s.Rational(5,216))*c for c in cs]
    assert s.trigsimp(sum(v*v for v in lambdas))==s.Rational(5,18)
    intervals=[tuple(map(F,v)) for v in result['effect_intervals_k0_k1_k2']]
    assert F(1,20000)<intervals[1][0]<intervals[1][1]<intervals[2][0]<intervals[2][1]<intervals[0][0]
    with pytest.raises(ValueError):analyzer.apply(np.zeros((3,3)))


def test_full_occupied_CAR_control_and_original_reference_histories():
    for e in analyzer.certificate()['effect_values_k0_k1_k2']:
        result=compiler.four_mode_control(e)
        assert result['full_Fock_dimension']==16
        assert result['complete_history_count']==256
        assert result['compression_error']<2e-14
        assert result['complete_Kraus_Gram_error']<2e-13
        assert abs(sum(result['complete_original_history_probabilities'].values())-1)<2e-12
        assert len(result['coarse_vacuum_nonvacuum_history_probabilities'])==16
        assert abs(sum(result['coarse_vacuum_nonvacuum_history_probabilities'].values())-1)<2e-12
        assert result['retained_coarse_quantum_ket_probability_error']<2e-14
        assert result['maximum_retained_fine_label_offdiagonal_density_Frobenius_norm']>0
        assert result['any_nonvacuum_probability']>0
        assert result['complete_archive_unhalved_error']<=result['bound']
        assert result['common_H0_commutes_analyzer_but_not_pairing']


def test_complete_mass_instrument_reads_m_not_pair_gap():
    result=compiler.exact_rest_instrument();records=result['complete_records']
    assert len(records)==36 and result['zero_probability_records']==30
    assert sum(F(r['probability_exact']) for r in records)==1
    assert result['full_16_state_CAR_intertwiner_energy_and_electric_charge_checked']
    assert all(r['native_electric_charge']==-r['relative_charge'] for r in records)
    assert {r['energy_in_units_of_same_m'] for r in records if r['success'] and F(r['probability_exact'])>0}=={1}
    assert any(not r['success'] and r['energy_in_units_of_same_m']==2 and F(r['probability_exact'])>0 for r in records)


def test_complete_archive_cost_scales_with_volume_and_keeps_failures():
    result=compiler.archive_bound(108000000000000)
    assert F(result['complete_unhalved_trace_error_exact'])==F(1,10**12)
    assert result['enlarged_CAR_modes']==324000
    assert result['fresh_ancillary_modes_each_round']==323910
    large=compiler.archive_bound(108*64**2*10**12,sites=64)
    assert F(large['complete_unhalved_trace_error_exact'])==F(1,10**12)
    assert not large['vacuum_postselection']
    assert large['all_quantum_archives_retained']
    with pytest.raises(ValueError):compiler.archive_bound(0)
    with pytest.raises(ValueError):compiler.archive_bound(10,time=.1)


def test_finite_ramsey_precision_resolves_small_electron_denominator():
    result=compiler.ramsey_readout()
    assert F(result['complete_channel_error_to_energy_error_exact'])==F(4,10**12)
    for record in result['records']:
        intervals=[tuple(map(F,p)) for p in record['ideal_Xplus_Xminus_Yplus_Yminus_intervals']]
        assert sum(p[0] for p in intervals)<=1<=sum(p[1] for p in intervals)
        assert all(0<=lo<=hi<=1 for lo,hi in intervals)
        assert len(record['compiled_original_probability_intervals'])==4
    for ratio in result['ratio_bounds']:
        low,high=map(F,ratio['ratio_interval']);assert 0<low<high
        assert float(high/low-1)<2e-7
        assert list(map(F,ratio['electron_anchored_GeV_interval']))==[F(result['electron_anchor_GeV_exact'])*low,F(result['electron_anchor_GeV_exact'])*high]
    with pytest.raises(ValueError):compiler.ramsey_readout(steps=10000)


def test_same_native_dirac_symbol_keeps_every_valley_and_zero_spectator():
    result=compiler.full_native_dirac()
    assert len(result['all_eight_valleys'])==8
    assert all(r['complete_Nambu_dimension']==180 and r['spectator_internal_zero_coefficients']==39 for r in result['all_eight_valleys'])
    assert all(len(r['complete_180_energy_rows'])==180 for r in result['all_eight_valleys'])
    assert result['complete_45_internal_rest_coefficients'].count(0)==39
    assert result['extended_kinetic_vacuum_reducing'] and result['extended_kinetic_commutes_full_analyzer']
    assert not result['ancillary_vacuum_is_extended_ground']
    assert np.allclose(result['same_effect_rest_energies_k0_k1_k2'],analyzer.certificate()['effect_values_k0_k1_k2'],atol=0,rtol=0)


def test_one_anchor_does_not_replace_the_compiled_source_coefficients():
    result=compiler.electron_calibration()
    values=[float(row['electron_anchored_GeV']) for row in result['rows']]
    assert np.allclose(values,[.00051099895069,.1056594144959,1.7769849708134],atol=5e-14,rtol=0)
    assert result['same_compiled_free_source_rest_coefficients']
    assert not result['observed_interacting_lepton_identification_proved']


def test_receipt_comparison_rejects_missing_outcomes_and_false_scope():
    from .packet import equal
    expected={'complete_records':[{'probability_exact':'0'},{'probability_exact':'1'}],
              'physical_identification':False,'error':1e-12}
    assert equal(expected,expected)
    assert not equal({**expected,'complete_records':expected['complete_records'][1:]},expected)
    assert not equal({**expected,'physical_identification':0},expected)
    assert not equal({**expected,'error':float('nan')},expected)
    tiny={'complete_original_history_probabilities':{'0,1':1e-18}}
    assert not equal({'complete_original_history_probabilities':{'0,1':0.}},tiny)


def test_literal_map_to_existing_charged_family_source_and_negative_controls():
    result=family_map.certificate()
    assert result['actual_native_effect_numeric_residual']<2e-13
    assert result['prepared_to_native_pair_factor_exact']=='5/4'
    assert result['prepared_double_trace_exact']=='5/18'
    assert result['distinct_native_pair_trace_exact']=='2/9'
    assert min(result['negative_control_operator_distances'].values())>1e-3
    assert result['admission']['literal_effect_to_existing_configured_action_frame_map']
    assert not result['admission']['configured_eight_flag_and_native_3600_address_complete_source_laws_identical']
    assert not result['admission']['observed_interacting_lepton_dictionary_proved']
