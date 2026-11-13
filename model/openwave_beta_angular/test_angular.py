import numpy as np
from . import angular as a

def test_actual_native_current_polynomial_and_noncollinear_spin_records():
    assert 'ax*bx' in a.spin_polynomial()
    d=a.source()
    for x,y in [(np.array([.5,.2,-.3]),np.array([-.2,.7,.1])),(np.array([0.,0.,.9]),np.array([.3,.4,-.5]))]:
        assert abs(a.explicit_spin_weight(x,y)-a.native_spin_weight(x[None],y[None],d['ga'])[0,0])<2e-11

def test_complete_3d_record_negative_parity_and_full_source_covariance():
    for L in (4,8):
        r=a.records(L);assert r['parity_odd_first_order_probability']<r['exact_full_Fock_parity_upper']<0
        assert 30<r['full_current_Wick_column_squared_norm']<40
        assert abs(sum(r['selected_probability_by_neutral_mass'])-r['first_order_selected_probability'])<1e-25
        assert len(r['complete_original_record_intervals'])==11
        assert r['retained_hadron_state_count']==120

def test_group_velocity_and_momentum_are_distinct_instruments():
    theta,p=a.grid(8);v=2*p*np.cos(theta)/np.sqrt(1+(p*p).sum(axis=1))[:,None]
    assert np.any(p*v<-.1)
    d,E=a.densities(p,.2,1);assert np.max(abs(np.trace(d,axis1=1,axis2=2)-1))<1e-14
    assert np.min(np.linalg.eigvalsh(d))>=-1e-14
