"""Exact original-weight and phase-access controls."""
from fractions import Fraction as F
import json

import pytest
import sympy as s

from . import phase_readout as lane
from model.standard_model_source_map import koide_source as source


Z = (s.Integer(99)+20*s.I)/101


@pytest.fixture(scope='module')
def result():
    return lane.certificate()


def test_fine_instrument_all_inputs_and_zeros(result):
    raw = result['raw_phase_equivalence']
    assert raw['fine_records_checked'] == 128
    assert raw['zero_weight_records_retained'] == 104
    assert raw['all_input_matrix_units_per_record'] == 9
    assert raw['complete_CP_maps_identical_for_every_unit_phase']
    assert not raw['coherent_flag_access_included_in_identity']


def test_literal_two_event_reference_entangled_countercontrol():
    # Two adaptive branches, with a noncommuting fixed intervention.
    # Equality on matrix units extends to arbitrary reference blocks.
    p = source.P; u = s.diag(1, s.I, -1)
    for i in range(3):
        for j in range(3):
            e = s.zeros(3); e[i, j] = 1
            for sign in (-1, 1):
                a = Z*p; b = Z**sign*p**(1 if i == 0 else 2)
                phased = b*u*a
                plain = p**(1 if i == 0 else 2)*u*p
                assert lane._clean(phased*e*phased.H-plain*e*plain.H) == s.zeros(3)


def test_all_phase_laurent_identity(result):
    assert result['all_phase_identity']['symbolic_all_phase_identity']


def test_tracial_analyzer_blind_but_pure_readout_sensitive(result):
    assert result['tracial_analyzer_records_identical_for_distinct_phases']
    invariant = s.ones(3, 1)/s.sqrt(3)
    probabilities = []
    for z in (1, Z):
        k = source.analyzer_block(0, 0, z)[:3, :3]
        probabilities.append(s.simplify((invariant.H*k.H*k*invariant)[0]))
    assert s.simplify(probabilities[0]-probabilities[1]) != 0


@pytest.mark.parametrize('z',[s.Integer(1)])
def test_complete_score_against_original_six_dimensional_kraus(z):
    row = lane.complete_score(z)
    ks = source.weak_kraus(z, s.Integer(1))
    assert lane._clean(sum((k.H*k for k in ks),s.zeros(6))) == s.eye(6)
    actual = [s.simplify(s.trace(k.H*k)/6) for k in ks]
    assert all(s.simplify(p-s.sympify(r['probability'])) == 0
               for p,r in zip(actual,row['records'],strict=True))
    assert len(actual) == 8 and s.simplify(sum(actual)) == 1
    assert all(p.is_nonnegative is True for p in actual)


def test_phase_counterfamily_is_not_just_coordinate_relabeling(result):
    assert s.sympify(result['distinct_action_second_moment_difference']) != 0
    controls = result['complete_weak_score_controls']
    assert controls[0]['nonzero_flag_probability'] != controls[1]['nonzero_flag_probability']
    assert controls[1]['nonzero_flag_probability'] == controls[2]['nonzero_flag_probability']


def test_exact_inverse_covers_known_program_phase_and_action_ratios(result):
    inv = result['certified_action_ratio_readout']
    lo,hi = map(F, inv['cos_three_theta_interval'])
    expected = F(s.re(Z**3))
    assert lo <= expected <= hi
    assert inv['action_ratios_certified']
    # Independent 70-digit direct diagonalization of the exact circulant.
    a,b = s.sqrt(source.WEIGHTS[0]),s.sqrt(source.WEIGHTS[1])
    omega = -s.Rational(1,2)+s.I*s.sqrt(3)/2
    values = sorted(float(s.re((a+b*(Z*omega**k+s.conjugate(Z)*omega**(-k))).evalf(70)))
                    for k in range(3))
    for (lo,hi), value in zip(inv['action_ratios_to_smallest'],
                              [(values[1]/values[0])**2,(values[2]/values[0])**2],strict=True):
        assert float(F(lo))-1e-9 <= value <= float(F(hi))+1e-9


def test_interval_error_is_charged_and_coverage_is_conditional(result):
    p = s.sympify(result['complete_weak_score_controls'][1]['nonzero_flag_probability'])
    interval = lane._rational_enclosure(p)
    tight = lane.cosine_interval(interval)
    loose = lane.cosine_interval(interval,event_error=F(1,10**6))
    a,b = map(F,tight['cos_three_theta_interval']); c,d = map(F,loose['cos_three_theta_interval'])
    assert c < a <= b < d
    assert not loose['score_coverage_or_native_error_proved_here']


def test_source_inconsistent_record_is_refused():
    out = lane.cosine_interval((F(0), F(0)))
    assert not out['certified']
    assert not out['positive_root_chamber_certified']
    assert not lane.action_ratio_intervals((F(0),F(0)))['action_ratios_certified']


def test_zero_phase_default_keeps_degeneracy(result):
    assert result['default_twofold_degeneracy']
    assert result['default_action_ratios_to_smallest'] == ['1','34+24*sqrt(2)']
    c = source.exact_source(1)['C_exclusive']
    a,b = s.sqrt(source.WEIGHTS[0]),s.sqrt(source.WEIGHTS[1])
    assert lane._clean(c-(a-b)*s.eye(3)) == b*s.ones(3)
    assert (c-(a-b)*s.eye(3)).rank() == 1
    assert (c-(a+2*b)*s.eye(3)).rank() == 2
    p = s.sympify(result['complete_weak_score_controls'][0]['nonzero_flag_probability'])
    interval = lane._rational_enclosure(p)
    assert not lane.action_ratio_intervals(interval)['action_ratios_certified']


def test_full_score_noise_never_asserts_positive_spectrum():
    out = lane.cosine_interval((F(0), F(1)))
    assert out['certified'] and not out['positive_root_chamber_certified']
    assert not lane.action_ratio_intervals((F(0),F(1)))['action_ratios_certified']


@pytest.mark.parametrize('bad',[0,-1,True,0.5,'1'])
def test_invalid_kick_is_rejected(bad):
    with pytest.raises(ValueError): lane.kick_factor(bad)


@pytest.mark.parametrize('bounds',[(F(-1),F(0)),(F(1),F(0)),(F(0),F(2)),(0.1,0.2)])
def test_invalid_probability_intervals_are_rejected(bounds):
    with pytest.raises(ValueError): lane.cosine_interval(bounds)


def test_json_output_and_scope(result):
    json.dumps(result)
    assert result['source_parameter_readable_without_inference']
    assert not result['coherent_phase_selected_by_raw_repair_records']
    assert not result['coherent_phase_inaccessible_to_all_programs']
    assert not result['physical_common_clock_selected']
    assert not result['physical_pole_mass_or_nature_assignment']
    assert not result['empirical_inputs_used']


def test_fresh_public_certificate_is_detached(result):
    result['complete_weak_score_controls'][0]['records'][0]['probability'] = 'forged'
    assert lane.complete_score()['records'][0]['probability'] != 'forged'


def test_known_program_phase_has_direct_ratio_enclosure():
    row = lane.phase_ratio_enclosure(Z)
    assert row['action_ratios_certified']
    assert row['exact_cos_three_theta'] == str(s.re(Z**3))
    assert row['phase_is_readable_specification_input']
    assert not row['phase_source_selection_or_calibration_provenance_supplied']
    assert lane.phase_ratio_enclosure(1)['default_degenerate_spectrum']


def test_complete_tracial_law_matches_literal_source_analyzer():
    assert lane.tracial_analyzer_law() == source.exact_analyzer_records(1)


def test_nondefault_clock_is_reported_exactly():
    result = lane.complete_score(1, F(1,64))
    assert result['half_angle_tangent'] == '1/64'
    assert result['kick_factor'] == '4/4097'


def test_complex_source_full_padding_numeric_crosscontrol():
    import numpy as np
    prep,fourier,blocks = source.analyzer_factors(Z)
    block = np.zeros((32,32),complex)
    for flag,b in enumerate(blocks):
        block[4*flag:4*flag+4,4*flag:4*flag+4] = np.array(b,complex)
    unitary = np.kron(np.array(fourier,complex),np.eye(4))@block@np.kron(np.array(prep,complex),np.eye(4))
    q = np.zeros((32,32),complex);q[:4,:4]=np.eye(4)
    projection = unitary.conj().T@q@unitary
    x = np.array([[0,1],[1,0]],complex)
    full = np.eye(64)-np.kron(np.eye(2),projection)-1j*np.kron(x,projection)
    assert np.linalg.norm(full.conj().T@full-np.eye(64)) < 1e-12
    columns = [0,1,2,32,33,34]
    probabilities=[]
    for flag in range(8):
        rows=[32*n+4*flag+j for n in range(2) for j in range(3)]
        k=full[np.ix_(rows,columns)]
        probabilities.append(np.linalg.norm(k,'fro')**2/6)
    exact=lane.complete_score(Z)['records']
    assert max(abs(p-float(s.sympify(e['probability']))) for p,e in zip(probabilities,exact,strict=True)) < 1e-12
    assert abs(sum(probabilities)-1) < 1e-12
