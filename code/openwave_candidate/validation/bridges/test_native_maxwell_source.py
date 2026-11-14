import numpy as np
import pytest
from scipy.linalg import expm
from native_maxwell_source import Construction,UniformFiberReduction,source,execute,run


def test_fourier_continuity_independent_finite_difference_and_wrong_charge_control():
    c=Construction(3,((0,1),(1,2),(2,0)),.23);m=UniformFiberReduction(c,6)
    psi=np.array([1,2j,0,1-1j,0,0],complex);psi/=np.linalg.norm(psi)
    xyz=np.array([[0,0,0],[1,0,0],[.5,.8,0]]);k=np.array([.7,.4,.3]);dt=1e-5
    s=source(m,psi,xyz,k)
    plus=source(m,expm(-1j*dt*m.hamiltonian)@psi,xyz,k)
    minus=source(m,expm(1j*dt*m.hamiltonian)@psi,xyz,k)
    assert abs((plus['charge_fourier']-minus['charge_fourier'])/(2*dt)+1j*k@s['current_fourier'])<1e-9
    # Dropping transit fractions destroys this SAME continuity identity.
    rho=np.outer(psi,psi.conj());derivative=-1j*(m.hamiltonian@rho-rho@m.hamiltonian)
    wrong=np.diag(derivative).real[:3]@np.exp(-1j*xyz@k)
    assert abs(wrong+1j*k@s['current_fourier'])>.01


def test_complete_bridge_accounting_and_hold_refinement():
    r=run();assert all(r['checks'].values());assert not r['full_physical_criterion_met']
    assert len(r['fine']['records'])==40
    assert r['fine']['source_hold_error']<r['coarse_source_hold_error']


def test_invalid_embedding_state_and_steps():
    m=UniformFiberReduction(Construction(2,((0,1),),.23),6)
    with pytest.raises(ValueError):source(m,np.ones(3),np.zeros((2,3)),np.ones(3))
    with pytest.raises(ValueError):source(m,np.array([1,0,0]),np.zeros((2,3)),np.array([np.nan,0,0]))
    with pytest.raises(ValueError):execute(0)
