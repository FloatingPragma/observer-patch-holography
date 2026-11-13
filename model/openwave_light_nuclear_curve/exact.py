"""Rational generalized-matrix inertia certifies every compressed Ritz level."""
from fractions import Fraction as F
import numpy as np
import sympy as sy
from .source import projected_source,UNIT
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals


def rational_matrix(a):
    # Every spin matrix entry is a dyadic rational, checked before admission.
    z=np.rint(a*4).astype(np.int64)
    assert np.array_equal(z/4,a)
    return sy.Matrix(z.tolist())/4


def certified_levels(cells,B):
    d=projected_source(cells,B,True);rows=[]
    Iz=F(-1,2) if cells%2 else F(0)
    for twice_s in range(cells%2,cells+1,2):
        for twice_i in range(cells%2,cells+1,2):
            s,i=F(twice_s,2),F(twice_i,2)
            ids=np.flatnonzero((np.diag(d['Sz'])==float(s))&(np.diag(d['Iz'])==float(Iz)))
            S=rational_matrix(d['S2'][ids][:,ids]);I=rational_matrix(d['I2'][ids][:,ids]);eye=sy.eye(len(ids))
            selector=(S-sy.Rational(s*(s+1))*eye)**2+(I-sy.Rational(i*(i+1))*eye)**2
            cols=selector.nullspace()
            if not cols:continue
            P=sy.Matrix.hstack(*cols);G=P.T*P
            H=sy.zeros(len(ids))
            for coefficient,a in d['rational_terms']:
                H+=sy.Rational(coefficient)*rational_matrix(a[ids][:,ids])
            A=P.T*H*P
            intervals=eigen_intervals(A,G)
            rows.append({'spin':str(s),'isospin':str(i),'Iz':str(Iz),'rank':len(cols),
                         'source_energy_intervals':intervals,
                         'GeV_intervals':[(lo*UNIT,hi*UNIT) for lo,hi in intervals]})
    return rows

