from fractions import Fraction as F
from .records import clock_ratios, weak_clock_bound, execute, field_controls

def test_all_native_clocks_redshift_with_retained_universality_residual():
    r = clock_ratios(); w = weak_clock_bound()
    assert all(F(x['same_source_dressed_clock_ratio_interval_exact'][1]) < 1 for x in r)
    assert F(r[-1]['same_source_dressed_clock_ratio_interval_exact'][0]) > F(r[0]['same_source_dressed_clock_ratio_interval_exact'][1])
    assert F(w['at_probe_energy_one_thousandth_error_upper_exact']) < F(1, 100000)

def test_joint_histories_include_failure_vacuum_and_work():
    r = execute(1)
    assert abs(sum(x['total_probability'] for x in r['first_arrival_joint_histories'])+r['retained_no_arrival_including_field_vacuum']['total_probability']-1) < 1e-9
    assert r['retained_no_arrival_including_field_vacuum']['total_probability'] >= .5
    assert abs(r['complete_history_energy_account']['residual']) < 1e-9
    assert r['full_history_TV_error_upper'] < 1e-6

def test_metric_cancellation_is_an_actual_complete_wave_control():
    r = field_controls(-1)
    assert abs(r['conditional_arrival_delay']) < 1e-11
    assert abs(r['transverse_position_difference']) < 1e-11
