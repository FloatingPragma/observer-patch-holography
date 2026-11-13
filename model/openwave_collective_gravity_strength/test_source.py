from fractions import Fraction as F
import numpy as np
import pytest
from . import source as s


def test_complete_small_hypercubes_give_the_binomial_response():
    X=np.array([[0.,1.],[1.,0.]]); Z=np.diag([1.,-1.]); I=np.eye(2)
    for n in (1,2,3,4):
        def embed(op,i):
            out=np.ones((1,1))
            for j in range(n):out=np.kron(out,op if i==j else I)
            return out
        H=float(s.OMEGA)/2*(n*np.eye(2**n)-sum(embed(X,i)for i in range(n)))
        Q=float(s.Q)/n*sum(embed(Z,i)for i in range(n))
        eigenvectors=np.linalg.eigh(H); e,U=eigenvectors; u=U[:,0]
        A=np.diag(np.exp(np.diag(Q))); B=np.diag(np.exp(np.diag(Q)/2))
        actual=-2*sum((u@A@U[:,r])*(U[:,r]@B@u)/e[r]for r in range(1,2**n))/(u@B@u)
        interval,weights=s.susceptibility(n,F(1,2))
        assert abs(actual-float((interval.lo+interval.hi)/2))<1e-12
        assert len(weights)==n


def test_integer_calibration_has_explicit_ancestry_and_keeps_electron():
    u,hc,target,_,_=s.inherited_units(); n,_=s.calibrated_count(u,target)
    lo,hi=s.coefficient_bounds(n)
    assert max(abs(lo/(target*u*u)-1),abs(hi/(target*u*u)-1))<F(1,10**40)
    assert u*(hc/u)/hc==1
    other,_=s.calibrated_count(u,2*target)
    assert other<n  # Hardware explicitly depends on the empirical target.


def test_free_gap_does_not_close_and_strength_dilutes():
    intervals=[s.coefficient_bounds(n)for n in (1,4,16,64)]
    assert all(right[1]<left[0]for left,right in zip(intervals,intervals[1:]))
    assert s.OMEGA==1


@pytest.mark.parametrize('n',[0,-1,True,F(1,2)])
def test_non_native_register_counts_rejected(n):
    with pytest.raises(ValueError):s.coefficient_bounds(n)
