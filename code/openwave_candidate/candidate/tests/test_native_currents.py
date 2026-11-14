import numpy as np
import pytest
from scipy.linalg import expm
from oph_openwave.native_currents import run,regional_current,complete_readout

def test_actual_native_outcomes_and_continuity():
    r=run();assert all(r['checks'].values())
    assert not r['full_criterion_met']
    assert len(r['native_complete_records'])>4

def test_direction_and_time_derivative_reject_sign_flip():
    h=np.array([[0,1j],[-1j,0]],dtype=complex);rho=np.ones((2,2))/2
    p,j=regional_current(h,[0]);dt=1e-6;u=expm(-1j*h*dt)
    finite=np.trace(p@(u@rho@u.conj().T-rho)).real/dt
    expected=np.trace(j@rho).real
    assert abs(finite-expected)<2e-6
    assert abs(finite+expected)>1
    _,all_current=regional_current(h,[0,1]);assert np.max(abs(all_current))==0

def test_reference_uniform_complete_PVM_and_validation():
    h=np.array([[0,1],[1,0]],dtype=complex);_,j=regional_current(h,[0])
    rho=np.array([[.5,.2j],[-.2j,.5]])
    records,post=complete_readout(j,rho)
    assert abs(np.trace(post)-1)<1e-12
    assert abs(sum(x['value']*x['probability'] for x in records)-np.trace(j@rho))<1e-12
    with pytest.raises(ValueError):regional_current(h,[True])
