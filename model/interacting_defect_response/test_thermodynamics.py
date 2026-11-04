import numpy as np
import pytest
from scipy.linalg import expm

from .finite import numerical_block
from .thermodynamics import gibbs,response,control


def test_noncommuting_BKM_response_and_energy_work_first_law():
    result = control()
    assert float(result['BKM_variance']) > 0
    assert float(result['free_energy_second_derivative']) < 0


def test_state_first_law_for_an_independent_trace_zero_perturbation():
    h,_,_ = numerical_block(1,.5)
    value = gibbs(h,.7)
    tangent = np.diag(np.linspace(-1,1,16))
    u,p = value['basis'],value['probability']
    log_rho = (u*np.log(p))@u.conj().T
    assert abs(-np.trace(tangent@log_rho)-.7*np.trace(tangent@h)) < 1e-12


def test_noncommuting_response_differs_from_classical_variance():
    h,x,_ = numerical_block(1,.5)
    value = response(h,x,.7)
    mean = value['mean_perturbation']
    classical = np.trace(value['rho']@x@x).real-mean*mean
    assert classical-value['BKM_variance'] > .05


def test_spectral_reset_is_trace_preserving_and_prepares_same_Gibbs_state():
    h,_,_ = numerical_block(1,.5)
    value = gibbs(h,.7)
    rng = np.random.default_rng(20261104)
    psi = rng.normal(size=16)+1j*rng.normal(size=16)
    psi /= np.linalg.norm(psi)
    initial = np.outer(psi,psi.conj())
    out = np.zeros((16,16),complex)
    gram = np.zeros_like(out)
    u,p = value['basis'],value['probability']
    for i in range(16):
        for j in range(16):
            k = np.sqrt(p[i])*np.outer(u[:,i],u[:,j].conj())
            gram += k.conj().T@k
            out += k@initial@k.conj().T
    assert np.linalg.norm(gram-np.eye(16)) < 1e-12
    assert np.linalg.norm(out-value['rho']) < 1e-12


@pytest.mark.parametrize('beta',[0,-1,float('nan')])
def test_invalid_Gibbs_domain(beta):
    with pytest.raises(ValueError):
        gibbs(np.eye(2),beta)
