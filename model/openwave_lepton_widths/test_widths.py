import numpy as np
from . import widths as w

def test_massive_phase_space_recovers_exact_massless_integral_and_reduces_width():
    for r in (.0048,.0595):
        zero=w.massive_integral(r,0,0)
        assert abs(12*zero-w.phase_factor(r))<2e-12
        assert 0<w.massive_integral(r,.01,.02)<zero
    assert w.mass_error_certificate()['relative_width_loss_upper_decimal']<1.7e-8

def test_complete_family_sum_units_and_held_out_ratio():
    e=w.channel(2,0);mu=w.channel(2,1)
    assert len(e['all9_massive_neutral_channels'])==9
    assert abs(sum(x['native_current_family_weight'] for x in e['all9_massive_neutral_channels'])-1)<2e-15
    ratio=mu['massive_neutral_width_GeV_diagnostic']/e['massive_neutral_width_GeV_diagnostic']
    assert .97256<ratio<.97257
    assert 6.582119568e-25<w.units()['hbar_GeV_s']<6.582119570e-25

def test_complete_finite_weak_source_records_and_uniform_kernel_bounds():
    w.kernel_certificate()
    for parent,daughter in ((1,0),(2,0),(2,1)):
        r=w.finite_records(parent,daughter,bins=32)
        assert abs(sum(r['all_survival_and_emission_probabilities'])-1)<2e-13
        assert min(r['all_survival_and_emission_probabilities'])>=0
        assert abs(r['energy_after_free_readout']-1-r['measurement_work_source_units'])<1e-15

def test_threshold_weighted_spectral_bound_and_actual_weak_lifetime_transfer():
    # Non-negligible test masses probe the inequality independently of tiny data.
    for r,a,b in ((.05,.001,.002),(.4,.03,.02),(.9,.04,.03),(.99,.006,.008)):
        I0=w.phase_factor(r)/12;Im=w.massive_integral(r,a,b)
        assert -1e-13<=I0-Im<=21*max(a,b)**2*(1-r*r)**3+1e-13
    w.spectral_transfer_certificate()
    for row in w.channel(1,0)['calibrated_massive_finite_lifetime_windows']:
        assert row['massive_exact_survival_additional_probability_error_upper']<5e-7
        assert row['massive_exact_survival_total_probability_error_upper']<5e-7
