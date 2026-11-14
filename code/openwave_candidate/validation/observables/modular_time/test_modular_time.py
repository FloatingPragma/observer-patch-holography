import json
import numpy as np
import pytest
from run_modular_time import modular_flow,run


def test_modular_automorphism_group_star_and_fixed_state():
    rho=np.array([[.7,.1],[.1,.3]],complex)
    a=np.array([[1,2+1j],[3-2j,4]],complex)
    assert np.allclose(modular_flow(rho,modular_flow(rho,a,.21),.37),modular_flow(rho,a,.58))
    assert np.allclose(modular_flow(rho,a.conj().T,.37),modular_flow(rho,a,.37).conj().T)
    assert np.allclose(modular_flow(rho,rho,.37),rho)
    assert np.isclose(np.trace(rho@modular_flow(rho,a,.37)),np.trace(rho@a))


def test_singular_state_requires_explicit_algebra_change_no_hidden_floor():
    with pytest.raises(ValueError):modular_flow(np.diag([1.,0.]),np.array([[0,1],[1,0]]),.37)
    # Explicit support algebra is M1: its modular flow is trivial.
    assert np.array_equal(modular_flow(np.array([[1.]]),np.array([[2.]]),.37),np.array([[2.]]))
    with pytest.raises(ValueError):modular_flow(np.eye(2)/2,np.eye(2),float('nan'))


def test_actual_native_observer_corner_and_receipt():
    result=run()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    assert result['observer']['depolarizing_floor_added'] is False
    assert result['observer']['thermal_or_Gibbs_state_supplied'] is False
    assert result['full_physical_criterion_met'] is False
    assert np.isclose(result['observer']['corner_probability']+result['observer']['complement_probability'],1)
    json.dumps(result,allow_nan=False)
