from fractions import Fraction as F
import numpy as np
from .source import lepton_source,local_norm_certificate,exact_baryon_window,full_source_terms,parameters,G


def test_complete_gaussian_reference_and_ward_identities():
    d=lepton_source()
    assert d['H'].shape==(96,96)
    assert min(d['E'])>0
    assert np.max(abs(d['current_covariance']-.75*np.eye(24)))<1e-15
    assert np.max(abs(d['current_means']))<1e-15


def test_full_local_contact_bound():
    d=local_norm_certificate()
    assert d['exact_operator_norm']=='6'
    assert d['three_quark_all_family_norm']=='54'
    p=parameters()
    assert p['mbar']-p['delta']/2-18*G>1


def test_actual_trial_complete_probability_bounds():
    d=exact_baryon_window()
    assert F(d['contact_second_moment_exact'])<37
    for r in d['complete_energy_window_instrument_bounds']:
        inside=F(r['original_inside_probability_lower'])
        outside=F(r['original_outside_probability_upper'])
        failure=F(r['preparation_failure_probability'])
        assert inside+outside+failure==1
        assert F(r['conditional_inside_probability_lower'])>F(96,100)


def test_entire_carrier_polynomial_rebuilds_native_baryon_action():
    d=full_source_terms()
    assert d['native_CAR_modes']==192
    assert d['colourful_modes']+d['Dirac_lepton_modes']==192
    assert d['complete_first_family_baryon_intertwiner_error']<1e-9
    assert len([x for x in d['terms']if x['kind']=='positive_linear_fermion_square'])==48
