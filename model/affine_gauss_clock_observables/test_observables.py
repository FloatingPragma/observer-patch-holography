from fractions import Fraction as F
import numpy as np
import pytest
from .certificate import roots,projector_bounds,certificate


def test_integer_root_enclosures():
    for x in [F(0),F(1),F(2),F(1,10**80),F(12345,789)]:
        lo,hi=roots(x);assert lo*lo<=x<=hi*hi
    with pytest.raises(ValueError):roots(-1)


def test_sharp_bound_saturates_two_dimensional_quantum_states():
    # These pure states have the declared fidelity and extremize a fixed
    # projector. The test exercises interference, not a classical mixture.
    for p in [F(1,10000),F(1,5),F(1,2),F(4,5)]:
        for f in [F(999999,1000000),F(9,10)]:
            g=np.array([np.sqrt(float(p)),np.sqrt(float(1-p))])
            orth=np.array([g[1],-g[0]])
            lo,hi=projector_bounds((p,p),f)
            for sign in [-1,1]:
                v=np.sqrt(float(f))*g+sign*np.sqrt(float(1-f))*orth
                q=abs(v[0])**2
                assert float(lo)-3e-15<=q<=float(hi)+3e-15
            if p>1-f:assert abs(float(lo)-abs((np.sqrt(float(f))*g-np.sqrt(float(1-f))*orth)[0])**2)<3e-15
            if p<f:assert abs(float(hi)-abs((np.sqrt(float(f))*g+np.sqrt(float(1-f))*orth)[0])**2)<3e-15


def test_fidelity_boundary_and_invalid_input():
    assert projector_bounds((F(0),F(0)),F(1))==(0,0)
    assert projector_bounds((F(0),F(1)),F(1))==(0,1)
    assert projector_bounds((F(1,2),F(1,2)),F(0))==(0,1)
    for interval,f in [((1,0),1),((-1,0),1),((0,1),F(11,10))]:
        with pytest.raises(ValueError):projector_bounds(interval,f)


def test_actual_source_complete_joint_records():
    d=certificate();rows=d['joint_post_clock_outcomes']
    assert {(r['source_vertex'],r['transverse_nonvacuum']) for r in rows}=={(x,b) for x in range(13) for b in [False,True]}
    assert sum(F(r['exact_normalized_trial_probability']) for r in rows)==1
    assert d['positive_joint_record_lower_bounds']==26
    assert F(d['conditional_post_clock_nonvacuum_probability_interval'][0])>F(7,100000)
    assert F(d['unconditional_accepted_nonvacuum_probability_interval'][0])>F(17,1000000)
    assert not d['exact_ground_projection_or_stationary_post_clock_state_claimed']
    assert d['all_energy_labels_retained']==2**20
    assert d['scope']['mass_GeV'] is None


def test_noncommuting_readout_preserves_the_autonomous_clock():
    from scipy.linalg import expm
    from scipy.integrate import solve_ivp
    x=np.array([[0,1],[1,0]],dtype=complex);eye=np.eye(2)
    h=np.kron(x,eye)
    p0=np.diag([1,0]);p1=np.diag([0,1])
    uread=np.kron(p0,eye)+np.kron(p1,x)
    log=np.pi*(np.eye(4)-uread)/2
    assert np.linalg.norm(h@log-log@h)>1
    def derivative(t,u):
        e=expm(-1j*t*h);v=e@log@e.conj().T
        return (-1j*(h+v)@u.reshape(4,4)).ravel()
    solved=solve_ivp(derivative,(0,1),np.eye(4,dtype=complex).ravel(),rtol=1e-11,atol=1e-13)
    assert solved.success
    assert np.linalg.norm(solved.y[:,-1].reshape(4,4)-expm(-1j*h)@uread)<2e-10
