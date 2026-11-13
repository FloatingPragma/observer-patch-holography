from fractions import Fraction as F
import numpy as np
import sympy as sy
from model.openwave_pole_weak_current import source as s
from model.openwave_pole_weak_current.native import contraction


def test_native_full_current_completeness_retains_extra_channels():
    d=contraction();M=d['axial_full_fixed_weight_square'];I=sy.eye(M.rows)
    P1=(4*I-M)/3;P4=(M-I)/3
    assert P1*P1==P1 and P4*P4==P4 and P1*P4==sy.zeros(M.rows)
    K=sy.sqrt(8)*P1/3+sy.sqrt(5)*P4/3
    assert sy.simplify(K*K+M/9)==I
    B=d['native_neutron_embedding'];C=d['native_proton_embedding']
    A=d['axial_numerator'];G=sy.diag(6,1,1,1,1)
    D=sy.diag(sy.sqrt(6),1,1,1,1);a=D.inv()*A*D.inv()
    leakage=sy.simplify(B.T*M*B-a.T*a)
    assert leakage!=sy.zeros(5)
    assert all(x>=-1e-12 for x in np.linalg.eigvalsh(np.array(leakage,float)))


def test_full_poles_have_positive_native_current_response():
    r=s.run();assert len(r['momentum_records'])==6
    for row in r['momentum_records']:
        for current in row['currents'].values():
            lo,hi=map(F,current['actual_source_pole_amplitude_magnitude_interval'])
            assert 0<lo<=hi<=3
            p,q=map(F,current['current_success_then_proton_pole_probability_interval'])
            assert p==lo**2/9 and q==hi**2/9
    rest=r['momentum_records'][0]['currents']
    assert rest['vector']['actual_source_pole_amplitude_magnitude_interval']==['1','1']
    assert r['autonomous_weak_decay_or_neutron_lifetime_proved']is False
    assert r['new_empirical_current_or_mass_input']is False


def test_current_instrument_independent_charge_and_energy_check():
    # Entire two-charge state space, including failed and already-raised inputs.
    J=np.array([[0.,0.],[1.,0.]])
    K=J/3;Fop=np.diag([np.sqrt(8/9),1.])
    H=np.diag([.9395654,.9382721]);rho=np.diag([1.,0.])
    post=K@rho@K.T+Fop@rho@Fop.T
    assert np.max(abs(K.T@K+Fop.T@Fop-np.eye(2)))<1e-15
    assert abs(np.trace(H@(post-rho))-(H[1,1]-H[0,0])/9)<1e-15
    assert abs(np.trace(post)-1)<1e-15
    # Tensor the success with a periodic reference shift; charge accounting
    # is exact on an input away from this finite test reference's wrap edge.
    U=np.roll(np.eye(5),-1,axis=0);full=np.kron(K,U)
    charge=np.kron(np.diag([0,1]),np.eye(5))+np.kron(np.eye(2),np.diag(range(5)))
    psi=np.zeros(10);psi[2]=1
    assert np.linalg.norm((charge@full-full@charge)@psi)<1e-15


def test_absolute_amplitude_enclosure_uses_complex_phase_and_error():
    square,interval=s.absolute_interval(3+4*sy.I,F(25),F(1,10),F(3))
    assert square==1 and interval==(F(9,10),F(11,10))
    assert s.absolute_interval(0,F(1),F(1,10),F(3))[1][0]==0
