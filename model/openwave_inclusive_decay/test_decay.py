from fractions import Fraction as F
import numpy as np
from . import decay as d

def test_native_current_trace_and_identical_counts():
    assert d.source_data()['coefficient']=='1'
    rows=d.spin_certificate();assert len(rows)==9
    assert all(abs(x['omitted_parent_spin_average']-2*x['trace'])<2e-14 for x in rows)
    phase=d.exact_phase_space();assert phase['massless_neutral_helicities_each']==1
    assert phase['identical_particle_factor'].startswith('1:')

def test_full_phase_space_and_stable_threshold():
    assert d.phase_factor(0)==1 and d.phase_factor(1)==0
    assert 0<d.phase_factor(1-1e-8)<1e-36
    assert all(d.density(e)==0 for e in(0,2,3))
    assert d.density(1)>0
    assert d.exact_phase_space()['rate']=='G^2*M^5*f(m/M)/(1536*pi^3)'

def test_rational_decay_window_and_finite_error():
    value=d.survival_bounds();b=value['rational_certificate']
    assert F(b['B0_squared_margin'])>0 and F(b['B1_cubed_margin'])>0
    assert value['rows'][0]['survival_interval_diagnostic'][0]>.3678
    assert value['rows'][0]['survival_interval_diagnostic'][1]<.3680
    assert value['rows'][-1]['inclusive_emission_interval_diagnostic'][0]>.9814
    assert F(value['rows'][-1]['finite_bin_survival_error_upper'])<F(1,10**9)

def test_complete_finite_source_energy_and_control():
    value=d.finite_records(64,1000.)
    assert abs(value['parent_survival']+sum(value['all_emission_bin_probabilities'])-1)<2e-13
    assert value['inclusive_emission']>0
    assert abs(value['full_energy_before_readout']-1)<2e-13
    assert d.finite_records(32,1000.,0.)['inclusive_emission']==0
    full=d.finite_records(32,.1);half=d.finite_records(32,.1,float(d.G)/2)
    assert np.isclose(full['inclusive_emission']/half['inclusive_emission'],4,atol=1e-6)
