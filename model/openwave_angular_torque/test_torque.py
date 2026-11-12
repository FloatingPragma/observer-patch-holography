import numpy as np
import pytest
from .source import oscillator,operators,evaluate
from .experiment import algebra
def test_actual_differential_rotation():
    for shell in (2,3):algebra(oscillator(shell),operators(oscillator(shell)))
def test_signed_current_in_both_energy_sectors():
    data=oscillator();ops=operators(data);rng=np.random.default_rng(31)
    for _ in range(3):
        psi=rng.normal(size=len(ops[0]))+1j*rng.normal(size=len(ops[0]));psi/=np.linalg.norm(psi)
        b=rng.normal(size=3)/10;h,g,_=evaluate(psi,b,ops)
        for j in range(3):
            v=np.eye(3)[j]*1e-5
            actual=np.vdot(psi,(evaluate(psi,b+v,ops)[0]-evaluate(psi,b-v,ops)[0])@psi).real/2e-5
            assert abs(actual-g[j])<4e-9
def test_incomplete_shell_rejected():
    with pytest.raises(ValueError):oscillator(1)
