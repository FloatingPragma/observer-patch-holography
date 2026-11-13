import numpy as np
from . import matter as m

def test_exact_native_forward_contraction_and_massive_all_channel_norm():
    d=m.native_forward();assert d['exact_polynomial']=='x**9*(x - 1)*(x + 1)**2'
    s=m.source(1e7);bar=m.source(1e7,True)
    assert np.linalg.norm(s['M'].reshape(6,2,6,2).transpose(2,1,0,3).reshape(12,12)+bar['M'])<1e-14
    assert np.linalg.norm(s['mean']+bar['mean'].conj())<1e-14
    assert np.linalg.norm(s['M'],2)<=1+1e-14
    assert np.linalg.norm(s['mean'],2)<=.25+1e-14

def test_complete_electron_medium_channel_and_refinement():
    a=m.run(.6,295.,2.6,N=128);b=m.run(.6,295.,2.6,N=1024)
    for r in (a,b):
        assert abs(sum(r['complete_finite_medium_records'].values())-1)<1e-10
        assert r['finite_medium_vs_mean_TV_diagnostic']<r['uniform_complete_channel_TV_bound']
        assert r['complete_finite_medium_records']['sterile_quantum_output']>=0
    assert b['finite_medium_vs_mean_TV_diagnostic']<a['finite_medium_vs_mean_TV_diagnostic']

def test_vacuum_limit_and_physical_potential_dictionary():
    r=m.run(.6,295.,rho_g_cm3=0,N=16)
    assert max(abs(r['complete_finite_medium_records'][k]-r['complete_vacuum_records'][k]) for k in r['complete_vacuum_records'])<1e-13
    u=m.unit_dictionary();assert 1e-13<u['charged_current_potential_eV']<1.1e-13
    assert u['native_common_neutral_current_shift_eV']==-.5*u['charged_current_potential_eV']
