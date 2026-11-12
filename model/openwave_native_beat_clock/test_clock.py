"""Independent source-space unitary and reference controls for the beat."""
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import scipy.linalg as la
import sympy as s
import pytest

from .clock import (partial_isometry_algebra,interchannel_algebra,trig_box,
                    parameters,interchannel_parameters,observation_bounds,complete_records)


def test_native_operator_algebras_and_reference_cost_cancellation():
    assert not partial_isometry_algebra()['ongoing_phase_rotations_or_periodic_drive']
    result=interchannel_algebra()
    assert result['exact_reference_cost_cancellation']
    assert result['reference_family_dephasing_kills_both_quadratures_for_all_time']


@pytest.mark.parametrize('seed',range(8))
def test_full_H_free_beat_on_arbitrary_retained_reference(seed):
    rng=np.random.default_rng(seed)
    omega=4.7;nu=.001;radius=.1
    h0=np.diag([0.,0.,3.,4.,omega,omega,omega+3,omega+4])
    perturb=np.zeros((8,8),complex)
    for start in (0,4):
        raw=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        block=raw+raw.conj().T;block/=np.linalg.norm(block,2)
        perturb[start:start+4,start:start+4]=block
    h=h0+nu*perturb
    A=np.zeros((8,8),complex);A[4,0]=A[5,1]=1
    P=A.conj().T@A;Q=A@A.conj().T;K=P+Q
    X=A+A.conj().T;Y=1j*(A-A.conj().T)
    preparation=np.eye(8)-K+(K+A-A.conj().T)/np.sqrt(2)
    # A nontrivial untouched two-state reference is retained throughout.
    initial=np.zeros(16,complex)
    initial[0]=np.sqrt(.37);initial[3]=np.exp(.7j)*np.sqrt(.63)
    prepared=np.kron(preparation,np.eye(2))@initial
    values,vectors=np.linalg.eigh(h)
    selected=(abs(values)<radius)|(abs(values-omega)<radius)
    herald=vectors[:,selected]@vectors[:,selected].conj().T
    post=np.kron(herald,np.eye(2))@prepared
    probability=np.vdot(post,post).real
    assert probability>=1-(nu/radius)**2-1e-13
    post/=np.sqrt(probability)
    # The full-H cross-reference frequency support is tested directly.
    low=values[abs(values)<radius];high=values[abs(values-omega)<radius]
    assert np.min(high)-np.max(low)>=omega-2*radius
    assert np.max(high)-np.min(low)<=omega+2*radius
    effects=[(np.eye(8)+X)/4,(np.eye(8)-X)/4,(np.eye(8)+Y)/4,(np.eye(8)-Y)/4]
    r0=np.diag([1]*4+[0]*4);r1=np.eye(8)-r0
    density=np.outer(post,post.conj())
    dephased=sum((np.kron(r,np.eye(2))@density@np.kron(r,np.eye(2)) for r in (r0,r1)))
    for t in (0,.2,1,10,100):
        unitary=la.expm(-1j*t*h)
        evolved=np.kron(unitary,np.eye(2))@post
        measured_x=np.vdot(evolved,np.kron(X,np.eye(2))@evolved).real
        measured_y=np.vdot(evolved,np.kron(Y,np.eye(2))@evolved).real
        error=2*nu/radius+2*nu*t
        assert abs(measured_x-np.cos(omega*t))<=error+1e-12
        assert abs(measured_y+np.sin(omega*t))<=error+1e-12
        weights=[np.vdot(evolved,np.kron(effect,np.eye(2))@evolved).real for effect in effects]
        assert min(weights)>-1e-13 and abs(sum(weights)-1)<1e-12
        joined=np.kron(unitary,np.eye(2))
        dephased_out=joined@dephased@joined.conj().T
        assert abs(np.trace(dephased_out@np.kron(X,np.eye(2))))<1e-12
        assert abs(np.trace(dephased_out@np.kron(Y,np.eye(2))))<1e-12
        unprepared=np.kron(unitary,np.eye(2))@initial
        assert abs(np.vdot(unprepared,np.kron(X,np.eye(2))@unprepared))<1e-12
        assert abs(np.vdot(unprepared,np.kron(Y,np.eye(2))@unprepared))<1e-12
    # An exact stationary input cannot self-start a fixed observable.
    ground=vectors[:,0]
    for t in (0,1,10):
        evolved=la.expm(-1j*t*h)@ground
        assert abs(np.vdot(evolved,X@evolved)-np.vdot(ground,X@ground))<1e-12


def test_integer_phase_reduction_covers_long_times_and_negative_angles():
    mp.mp.dps=90
    for center in (F(0),F(-23,10),F(12345,1000),F(10000001,10)):
        half=F(1,10**22)
        sine,cosine=trig_box((center-half,center+half))
        for point in (center-half,center,center+half):
            angle=mp.mpf(point.numerator)/point.denominator
            for value,interval in ((mp.sin(angle),sine),(mp.cos(angle),cosine)):
                lo=mp.mpf(interval[0].numerator)/interval[0].denominator
                hi=mp.mpf(interval[1].numerator)/interval[1].denominator
                assert lo<=value<=hi


def test_source_bins_and_nonzero_native_interchannel_frequency():
    base=parameters();other=interchannel_parameters()
    assert base['omega'][0]>4
    assert s.Rational(59,100)<other['omega'][0]<s.Rational(3,5)
    assert base['herald_success_lower']>s.Rational(999999999,10**9)
    row=observation_bounds(10**6)
    assert s.Rational(row['expectation_error_upper'])<s.Rational(43,1000)
    assert complete_records()['coarse_record_count']==7
    with pytest.raises(ValueError):observation_bounds(-1)
    with pytest.raises(ValueError):observation_bounds(1,channel='physical_neutrino')
