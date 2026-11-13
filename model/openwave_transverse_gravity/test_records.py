import numpy as np
from .field import field,detector
from .records import execute,comparisons
from .units import comparison


def test_full_transverse_gauss_and_absorber_completeness():
    d=field();mask,K,c=detector();F=d['frame']
    assert F.shape==(882,540)
    assert np.linalg.norm(d['gradient'].T@F)<1e-11
    assert np.linalg.norm(K@K+.5*F[mask].T@F[mask]-np.eye(540))<1e-11
    assert c['no_click_effect_lower_exact']=='1/2'


def test_actual_photon_clock_arrival_and_work_share_one_source():
    a,b=execute(0),execute(1);c=comparisons(a,b)
    assert c['conditional_arrival_delay']>c['arrival_delay_error_upper']
    assert c['transverse_position_difference']+c['position_difference_error_upper']<0
    assert c['separate_Fourier_direction_deflection']+c['direction_error_upper']<0
    assert all(0<r<1 for r in c['same_joint_native_clock_phase_ratios'])
    for x in(a,b):
        assert abs(x['arrival_success_probability']+x['retained_no_arrival_and_photon_vacuum']['total_probability']-1)<1e-9
        assert abs(x['complete_history_work']['residual'])<1e-9
        account=x['complete_preparation_and_work']
        assert account['input_energy_exact']=='0'
        assert account['mean_preparation_work']==x['complete_history_work']['initial']
        assert abs(account['mean_total_work']-account['final_mean_energy'])<1e-9


def test_independent_units_retain_both_physical_failures():
    r=comparison()
    assert r['source_stiffness_or_q_tuned_to_G']is False
    assert r['physical_G_strength_match']is False
    assert r['physical_common_light_speed_match']is False
    assert r['conditional_to_observed_strength_ratio']>1e39
