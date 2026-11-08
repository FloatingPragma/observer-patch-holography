from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import algebra, readout


def test_actual_normalized_color_and_full_gauss_source():
    row=algebra.source_certificate()
    assert row['full_source_Gauss_checks']==240
    assert row['full_source_action_checks']==30
    assert row['complete_link_dimension']==60
    assert row['committed_source_defaults']=={'hopping_J':'1','electric_epsilon':'1',
                                               'epsilon_over_abs_J':'1'}
    assert row['alpha_s_MSbar_prediction'] is None


def test_exact_electric_witness_annihilates_every_quadratic_casimir():
    row=algebra.source_certificate()
    matrix=s.Matrix([[s.Rational(v[k]) for k in ('color_C2','weak_C2','Y_squared')]
                     for v in row['blocks']])
    vector=s.Matrix(list(map(s.Rational,row['electric_incompatibility_left_witness'])))
    assert vector.T*matrix==s.zeros(1,3)
    assert (vector.T*s.ones(5,1))[0]!=0
    # Positive control: an actual specified Casimir action is in the range.
    chosen=matrix*s.Matrix([2,3,5])
    assert vector.dot(chosen)==0


@pytest.mark.parametrize('d',[1,2,3,6])
def test_exact_native_default_readout_and_reference_backaction(d):
    h=np.array([[0,np.sqrt(d)],[np.sqrt(d),1]],complex)
    t=0.137; unitary=expm(-1j*t*h)
    projectors=(np.diag([1,0]),np.diag([0,1]))
    branches=[p@unitary for p in projectors]
    assert np.linalg.norm(sum(k.conj().T@k for k in branches)-np.eye(2))<2e-15
    # Untouched reference has coherent support in both outcomes.
    psi=np.array([[1,2j],[3+2j,-1j]],complex);psi/=np.linalg.norm(psi)
    complete=np.stack([k@psi for k in branches])
    assert abs(np.linalg.norm(complete)-1)<2e-15
    assert np.linalg.norm(complete[0]-projectors[0]@unitary@psi)<2e-15
    omega2=1+4*d
    p=4*d/omega2*np.sin(np.sqrt(omega2)*t/2)**2
    assert abs(abs(unitary[1,0])**2-p)<2e-15


@pytest.mark.parametrize('d',[1,2,3,6])
def test_rational_two_time_inference_contracts_and_recovers_source(d):
    widths=[]
    for time in (F(1,256),F(1,512)):
        p1=readout.probability(d,F(1),F(1),time)
        p2=readout.probability(d,F(1),F(1),2*time)
        answer=readout.infer_coefficients(d,time,p1,p2,F(1),F(1))
        for name in ('hopping_squared','electric_squared','electric_over_hopping_squared'):
            assert answer[name].contains(F(1))
        widths.append(answer['electric_squared'].high-answer['electric_squared'].low)
    assert widths[1]<widths[0]/3


@pytest.mark.parametrize('j,energy',[(F(1,2),F(1,3)),(F(2),F(1)),(F(1),F(0))])
def test_inference_is_not_fitted_to_default(j,energy):
    d=3;time=F(1,1024)
    p1=readout.probability(d,j,energy,time)
    p2=readout.probability(d,j,energy,2*time)
    answer=readout.infer_coefficients(d,time,p1,p2,F(2),F(2))
    assert answer['hopping_squared'].contains(j*j)
    assert answer['electric_squared'].contains(energy*energy)


def test_probability_error_is_retained_and_low_precision_rejected():
    d=3;t=F(1,256);eta=F(1,10**14)
    p1=readout.probability(d,F(1),F(1),t)
    p2=readout.probability(d,F(1),F(1),2*t)
    fine=readout.infer_coefficients(d,t,p1,p2,F(1),F(1))
    wide=readout.infer_coefficients(d,t,p1.expanded(eta),p2.expanded(eta),F(1),F(1))
    assert wide['electric_squared'].low<fine['electric_squared'].low
    assert wide['electric_squared'].high>fine['electric_squared'].high
    with pytest.raises(ValueError):
        readout.infer_coefficients(d,t,readout.Interval(F(),F(1)),
                                  readout.Interval(F(),F(1)),F(1),F(1))


def test_common_clock_rescaling_preserves_complete_probability():
    d=3;j=F(1);e=F(2);t=F(1,1024);scale=F(7,3)
    assert readout.probability(d,j,e,t)==readout.probability(d,scale*j,scale*e,t/scale)


def test_invalid_or_aliasing_inputs_are_rejected():
    with pytest.raises(ValueError):readout.probability(3,F(1),F(1),F(1))
    with pytest.raises(ValueError):readout.probability(True,F(1),F(1),F(1,1024))
    with pytest.raises(ValueError):readout.probability(3,1.0,F(1),F(1,1024))
    with pytest.raises(ValueError):readout.Interval(F(2),F(1))
