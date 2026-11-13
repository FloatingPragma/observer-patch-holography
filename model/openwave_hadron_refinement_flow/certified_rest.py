"""Exact rational calibration-band enclosures for the three-state rest block."""
from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from .analytic import exact_rest,constants
from .parameters import degree_coefficients


def inertia(M):
    A=[[F(v)for v in row]for row in M.tolist()];count=0
    while A:
        n=len(A);j=next((i for i in range(n)if A[i][i]),None)
        if j is None:
            pair=next(((i,k)for i in range(n)for k in range(i)if A[i][k]),None)
            if pair is None:break
            i,k=pair;old=[row[:]for row in A]
            for r in range(n):A[i][r]=old[i][r]+old[k][r]
            for r in range(n):A[r][i]=old[r][i]+old[r][k]
            A[i][i]=old[i][i]+old[k][k]+2*old[i][k]
            j=i
        A[0],A[j]=A[j],A[0]
        for row in A:row[0],row[j]=row[j],row[0]
        d=A[0][0];count+=d<0
        A=[[A[i][k]-A[i][0]*A[0][k]/d for k in range(1,n)]for i in range(1,n)]
    return count


def sqrt_upper(q,den=10**24):
    q=F(q);assert q>=0
    n=isqrt(q.numerator*den*den//q.denominator)
    if n*n*q.denominator<q.numerator*den*den:n+=1
    return F(n,den)


def eigen_intervals(A,G):
    """Deterministic rational bisection; floating seeds do not enter receipts."""
    assert inertia(A)==0 and inertia(A-16384*G)==A.rows
    result=[]
    for j in range(A.rows):
        lo,hi=F(0),F(16384)
        while hi-lo>F(1,10**8):
            middle=(lo+hi)/2
            if inertia(A-sy.Rational(middle)*G)<=j:lo=middle
            else:hi=middle
        result.append((lo,hi))
    return result


def rest(dU=F(0),dm=F(0)):
    H,O,J=exact_rest(dU,dm);D=sy.diag(sy.sqrt(3),1,1)
    A=sy.simplify(D*H*D);G=sy.diag(3,1,1);AJ=sy.simplify(D*J*D)
    assert all(x.is_Rational for x in A) and all(x.is_Rational for x in AJ)
    endpoints=eigen_intervals(A,G);oddlo=eigen_intervals(O,sy.eye(2))[0][0]
    witnesses=json.loads((Path(__file__).parent/'witnesses.json').read_text())
    assert (dU,dm)==(F(0),F(0)) or (dU,dm)==degree_coefficients()
    key='fixed' if (dU,dm)==(F(0),F(0)) else 'flow'
    vector=sy.Matrix([sy.Rational(x)for x in witnesses[key]])
    n=F((vector.T*G*vector)[0]);mean=F((vector.T*A*vector)[0])/n
    residual=A*vector-sy.Rational(mean)*G*vector
    variance=F((residual.T*G.inv()*residual)[0])/n
    gap=min(endpoints[1][0],oddlo)-endpoints[0][1]
    sine=sqrt_upper(variance)/(endpoints[1][0]-mean)
    assert gap>0 and sine<F(1,10**6)
    axial=F((vector.T*AJ*vector)[0])/n
    axerror=6*sine
    centered=AJ*vector-sy.Rational(axial)*G*vector
    axvar=F((centered.T*G.inv()*centered)[0])/n
    return {'energy_interval':endpoints[0],'next_energy_lower':min(endpoints[1][0],oddlo),
            'gap_lower':gap,'axial_interval':(axial-axerror,axial+axerror),
            'axial_variance_sqrt_upper':sqrt_upper(axvar)+6*sine,
            'exact_trial_variance':variance,'exact_trial_sine_bound':sine,
            'rational_trial':[str(v)for v in vector],
            'Gram':[3,1,1]}


def full(dU=F(0),dm=F(0)):
    from .analytic import full_bounds
    base=full_bounds(dU,dm);r=rest(dU,dm);C=base['Q_floor'];eps=base['uniform_energy_error'];_,_,_,p,_,unit=constants(dU,dm)
    sine=eps/(r['gap_lower']-eps)
    tail=5*p['t']**2/(C-r['energy_interval'][1])**2
    gaerr=2*r['axial_variance_sqrt_upper']*sine+6*sine*sine+6*tail
    return {'rest_certificate':r,'Q_floor':C,'full_Q_probability_upper':tail,
            'P_eigenvector_sine_upper':sine,'full_source_GeV_interval':[(r['energy_interval'][0]-eps)*unit,r['energy_interval'][1]*unit],
            'full_source_axial_interval':[r['axial_interval'][0]-gaerr,r['axial_interval'][1]+gaerr]}


if __name__=='__main__':
    import json
    print(json.dumps({'fixed':full(),'flow':full(*degree_coefficients())},default=str),flush=True)
