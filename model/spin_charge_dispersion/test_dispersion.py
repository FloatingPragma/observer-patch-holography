import copy
from fractions import Fraction as F
import json
import subprocess
import sys

import numpy as np
import pytest
import sympy as s
from scipy.integrate import quad

from . import bounds, finite, packet


@pytest.fixture(scope='module')
def result():
    return packet.build()


def test_angular_projection_orientation_and_radial_coefficient():
    for x, y in ((.1, .2), (1., -2.), (3., 4.)):
        def integrand(theta):
            a = x*np.cos(theta)+y*np.sin(theta)
            return np.exp(-1j*theta)*a/(1+a*a)/(2*np.pi)
        integral = quad(lambda t: integrand(t).real, 0, 2*np.pi, epsabs=1e-12)[0]
        integral += 1j*quad(lambda t: integrand(t).imag, 0, 2*np.pi, epsabs=1e-12)[0]
        radius = x*x+y*y
        expected = (x-1j*y)*(1-1/np.sqrt(1+radius))/radius
        assert abs(integral-expected) < 1e-12


def test_radial_derivative_and_jensen_convexity_exact():
    u = s.symbols('u', positive=True)
    h = (1-(1+u)**(-s.Rational(1, 2)))/u
    assert s.simplify(h+u*s.diff(h,u)-1/(2*(1+u)**s.Rational(3, 2))) == 0
    assert s.simplify(s.diff((1+u)**(-s.Rational(3, 2)),u,2)-s.Rational(15,4)*(1+u)**(-s.Rational(7,2))) == 0
    assert s.limit(h,u,0) == s.Rational(1,2)


def test_same_source_exact_charge_projection_independent_fourier():
    op, values, numerator, minus, divisor = finite.operators()
    w = op['walsh']/16
    charges = np.array(op['charges'])
    f = np.diag(np.array(values, float))
    projected = np.zeros((256,256), complex)
    for j in range(16):
        theta = 2*np.pi*j/16
        rotation = (w*np.exp(1j*theta*charges))@w.T
        projected += np.exp(-1j*theta)*(rotation@f@rotation.conj().T)/16
    np.testing.assert_allclose(projected,numerator/divisor,atol=4e-15)
    np.testing.assert_allclose(projected.conj().T,minus/divisor,atol=4e-15)


def test_finite_complete_instrument_retains_arbitrary_reference():
    _, _, plus, minus, divisor = finite.operators()
    a,b=plus/divisor,minus/divisor
    effect=np.eye(256)-a.conj().T@a-b.conj().T@b
    values,vectors=np.linalg.eigh(effect)
    assert values.min() >= .5-1e-14
    complement=(vectors*np.sqrt(values))@vectors.conj().T
    # Finite noncommuting configurations need no normal-operator assumption.
    assert np.linalg.norm(a.conj().T@a-a@a.conj().T)>1e-3
    rng=np.random.default_rng(13217)
    state=rng.normal(size=(256,3))+1j*rng.normal(size=(256,3))
    state/=np.linalg.norm(state)
    output=np.concatenate([a@state,b@state,complement@state])
    np.testing.assert_allclose(output.conj().T@output,state.conj().T@state,atol=3e-15)


def test_ground_probabilities_use_actual_complete_source(result):
    op, _, plus, minus, divisor=finite.operators()
    _, vectors=np.linalg.eigh(op['H']/op['scale'])
    ground=vectors[:,0]
    rows=result['exact_finite_source_creator']['creation_records']
    for matrix,row in zip((plus,minus),rows):
        probability=np.linalg.norm(matrix@ground/divisor)**2
        lower,upper=map(F,row['original_probability_interval'])
        assert float(lower)<=probability<=float(upper)
    assert rows[0]['original_probability_interval']==rows[1]['original_probability_interval']


def test_band_scaling_and_original_probability():
    for m0 in (F(1,1000),F(1,4),F(1,2)):
        for radius in (F(1,100),F(1),F(7)):
            row=bounds.band(m0,radius)
            d2=F(row['squared_commutator_expectation_lower'])
            a,b=map(F,row['dimensionless_energy_over_radius_window'])
            assert a/12+1/(4*b)==d2/4
            assert F(row['each_charge_original_joint_band_probability_lower'])==m0*m0/17576
            assert list(map(F,row['energy_window_in_same_scaled_clock']))==[radius*a,radius*b]
            assert a<b
    assert F(bounds.event_budget()['observed_original_probability_lower'])>0


@pytest.mark.parametrize('m0,radius', [(0,1),(F(3,4),1),(F(1,4),0),(F(1,4),-1)])
def test_invalid_source_band_inputs_rejected(m0,radius):
    with pytest.raises(ValueError):
        bounds.band(m0,radius)


def test_scope_and_complete_outcomes_cannot_be_promoted(result):
    altered=copy.deepcopy(result)
    altered['scope']['sharp_nonzero_residue_dispersion_or_Lorentz_covariance_proved']=True
    with pytest.raises(ValueError):
        packet.validate(altered)
    altered=copy.deepcopy(result)
    altered['exact_finite_source_creator']['creation_records'].pop()
    with pytest.raises(ValueError):
        packet.validate(altered)


def test_actual_runtime_source_imports_bound():
    script='''
import json,sys
from pathlib import Path
from model.spin_charge_dispersion import packet
packet.build()
root=packet.ROOT.resolve()
local=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None: continue
    path=Path(path).resolve()
    if path.suffix=='.py' and path.is_relative_to(root):
        local.append(path.relative_to(root).as_posix())
print(json.dumps(sorted(set(local)-set(packet.bindings()))))
'''
    completed=subprocess.run([sys.executable,'-W','error','-c',script],cwd=packet.ROOT,
        text=True,capture_output=True,check=True)
    assert json.loads(completed.stdout)==[]
