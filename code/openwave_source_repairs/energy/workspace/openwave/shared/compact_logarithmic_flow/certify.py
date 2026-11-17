"""Exact, gap-free full-source spectral and logarithm certificates."""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
GRAM_GATE=F(1,10**9)
RESIDUAL_GATE=F(1,10**9)

def load():
    proposal=json.loads((HERE/'proposal.json').read_text())
    source=np.load(HERE/'source_columns.npz')
    return proposal,source

def certify(law):
    proposal,source=load();Q=proposal['denominator'];K=1260
    rows=proposal['laws'][law];V=np.asarray(rows['columns'],dtype=object).T
    values=np.asarray(rows['eigenvalues'],dtype=object);dimension=len(values)
    C=np.asarray(source['counts' if law=='native' else 'changed_counts'],dtype=object)
    assert V.shape==C.shape==(dimension,dimension)
    assert np.array_equal(C,C.T) and all(x==K for x in C.sum(axis=0))
    degree=max(K-int(C[i,i]) for i in range(dimension))
    floor=F(K-2*degree,K)
    assert all(floor<=F(int(x),Q)<=1 for x in values)
    gram=V.T@V-Q*Q*np.eye(dimension,dtype=object)
    gram_max=max(sum(abs(int(x)) for x in row) for row in gram)
    residual=Q*(C@V)-K*V*values[None,:]
    residual_max=max(max(sum(abs(int(x)) for x in row) for row in residual),
                     max(sum(abs(int(x)) for x in column) for column in residual.T))
    g=F(gram_max,Q*Q);r=F(residual_max,K*Q*Q)
    assert g<=GRAM_GATE and r<=RESIDUAL_GATE
    energy=K*(r+2*g)/floor
    return {'dimension':dimension,'count_denominator':K,'proposal_denominator':Q,
            'source_positive_floor':str(floor),'gram_operator_upper':str(g),
            'residual_operator_upper':str(r),'energy_operator_upper':str(energy),
            'phase_projector_product_upper':str(8*dimension*g),
            'gram_gate':str(GRAM_GATE),'residual_gate':str(RESIDUAL_GATE),
            'gap_required':False,'comparison_polar_basis_executed':False}

def energy_interval(value,denominator=1<<44,count_denominator=1260,degree=20):
    value=F(int(value),denominator);d=1-value
    assert 0<=d<1
    lo=sum((d**k/k for k in range(1,degree+1)),F())
    hi=lo+d**(degree+1)/((degree+1)*(1-d))
    return count_denominator*lo,count_denominator*hi

def phase(value,t):
    lo,hi=energy_interval(value);a,b=-F(t)*hi,-F(t)*lo
    angle=float((a+b)/2)
    error=max(abs(F(angle)-a),abs(F(angle)-b))
    assert error<=F(1,10**12)
    return angle,error
