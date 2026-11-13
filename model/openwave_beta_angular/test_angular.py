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

def test_actual_positive_spin_gradient_quartet_transfer():
    from .transfer import certificate
    r=certificate();assert abs(r['actual_native_vector']-1)<2e-10
    assert abs(r['actual_native_axial']-1.2753)<1e-8
    assert r['actual_strong_parameters']['spin_gradient']=='5/2'
    assert all(x['exact_full_Fock_parity_upper']<0 for x in r['finite_cube_records'])


def test_exact_energy_census_retains_every_spatial_momentum():
    import sympy as sy
    from collections import Counter
    from itertools import product
    for L,levels,counts in [(4,[1,5,9,13],[8,24,24,8]),(8,[1,3,5,7,9,11,13],[8,48,120,160,120,48,8])]:
        one=[sy.simplify(4*sy.sin(2*sy.pi*n/L)**2)for n in range(L)]
        actual=Counter(1+sum(v)for v in product(one,repeat=3))
        c=a.exact_free_electron_energy_census(L)
        assert c['squared_energy_source_units']==levels
        assert c['spatial_momentum_multiplicities']==counts
        assert dict(zip(levels,counts))==actual and sum(counts)==L**3
        theta,p=a.grid(L)
        numerical=1+(p*p).sum(axis=1)
        exact=np.array([1+sum(v)for v in product(one,repeat=3)],dtype=float)
        assert np.max(abs(numerical-exact))<1e-14
        assert c['electron_energy_GeV']==[float(np.sqrt(x)*.00051099895069)for x in levels]


def test_complete_energy_records_use_exact_classes_without_roundoff_duplicates():
    for L,n in [(4,4),(8,7)]:
        r=a.records(L)
        assert len(r['electron_energy_GeV'])==n
        assert sum(r['exact_free_electron_energy_census']['spatial_momentum_multiplicities'])==r['spatial_sites']
        assert len(r['complete_original_record_intervals'])==11
        assert len(r['first_order_joint_z_sign_probabilities'])==3
        assert len(r['selected_probability_by_neutral_mass'])==3
