from fractions import Fraction as F
import numpy as np
from . import source as s
from model.openwave_source_ward.source import clock_slope

def test_original_register_is_recovered_and_full_multiplicities_retained():
    for f in (F(0),F(1,4),F(1)):
        a,_=s.susceptibility(1,f);lo,hi=clock_slope(f)
        assert a.lo<=hi and a.hi>=lo
    rows=s.vacuum_records(16,F(1))['coordinate_outcomes']
    assert sum(F(x['probability_exact'])for x in rows)==1

def test_independent_full_hypercube_spectral_response():
    X=np.array([[0.,1.],[1.,0.]]);Z=np.diag([1.,-1.]);I=np.eye(2)
    for n in (2,3,4):
        def embed(op,i):
            out=np.ones((1,1))
            for j in range(n):out=np.kron(out,op if i==j else I)
            return out
        omega=s.frequency(n);w=float((omega.lo+omega.hi)/2)
        H=w/(2*n)*(n*np.eye(2**n)-sum(embed(X,i)for i in range(n)))
        q=float(s.Q)/n*sum(embed(Z,i)for i in range(n))
        e,U=np.linalg.eigh(H);u=U[:,0];f=.5
        A=np.diag(np.exp(np.diag(q)));B=np.diag(np.exp(f*np.diag(q)))
        total=-2*sum((u@A@U[:,r])*(U[:,r]@B@u)/e[r]for r in range(1,2**n))/(u@B@u)
        exact,_=s.susceptibility(n,F(1,2))
        assert abs(total-float((exact.lo+exact.hi)/2))<1e-12

def test_noise_reduces_but_retained_clock_gap_also_closes():
    variances=[];gaps=[]
    for n in (1,4,16,64):
        v=s.vacuum_records(n,F(1));variances.append(F(v['relative_lapse_variance_interval'][1]))
        gaps.append(s.frequency(n).hi/n)
    assert all(b<a for a,b in zip(variances,variances[1:]))
    assert all(b<a for a,b in zip(gaps,gaps[1:]))

def test_complete_small_register_ground_lies_in_classical_enclosure():
    X=np.array([[0.,1.],[1.,0.]]);Z=np.diag([1.,-1.]);I=np.eye(2)
    for n in (2,4):
        def embed(op,i):
            out=np.ones((1,1))
            for j in range(n):out=np.kron(out,op if i==j else I)
            return out
        omega=s.frequency(n);w=float((omega.lo+omega.hi)/2)
        H=w/(2*n)*(n*np.eye(2**n)-sum(embed(X,i)for i in range(n)))
        q=float(s.Q)/n*sum(embed(Z,i)for i in range(n))
        H+=np.diag(np.exp(np.diag(q))+F(1,2)*np.exp(np.diag(q)/2)).astype(float)
        energy=np.linalg.eigvalsh(H)[0]
        row=s.classical_limit_record(n,s.Interval.point(1),s.Interval.point(F(1,2)),F(1,2))
        lo,hi=map(F,row['complete_quantum_ground_energy_enclosure'])
        assert float(lo)<=energy<=float(hi)
