"""Exact source-bond and bounded configuration-instrument identities."""
from functools import lru_cache
import copy
import sympy as s

from model.spin_canonical_energy.algebra import exact_checks as inherited_checks


@lru_cache(maxsize=1)
def _checks():
    inherited=inherited_checks()
    pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    left=[s.kronecker_product(a/2,s.eye(2))for a in pauli]
    right=[s.kronecker_product(s.eye(2),a/2)for a in pauli]
    h=sum((a*b for a,b in zip(left,right)),s.zeros(4))-s.eye(4)/4
    gx,gy,fx,fy=s.symbols('gx gy fx fy',real=True)
    comm=lambda a,b:a*b-b*a
    controls=[]
    for axis in range(3):
        b,c=(axis+1)%3,(axis+2)%3
        current=left[b]*right[c]-left[c]*right[b]
        af=fx*left[axis]+fy*right[axis]
        agf=gx*fx*left[axis]+gy*fy*right[axis]
        residual=s.I*(gy-gx)*(fx+fy)*current/2
        assert s.simplify(comm((gx+gy)*h/2,af)-comm(h,agf)-residual)==s.zeros(4)
        assert max(abs(v)for v in current.eigenvals())==s.Rational(1,2)
        controls.append({'axis':axis,'exact_weighted_current_commutator':True,'current_norm':'1/2'})
    # A two-edge finite source keeps the cutoff plateau on all touched bonds.
    a,b,d,e=s.symbols('a b d e',real=True)
    assert s.expand((a+d)/2*(fx-fy)-(a*fx-d*fy)-(d-a)*(fx+fy)/2)==0
    # Common configuration calculus: |B+-B+'|=|B--B-'|=d,
    # |B0-B0'|<=sqrt(2)d, hence full isometry squared error<=4d².
    x,y=s.symbols('x y',real=True)
    assert s.simplify(s.diff(x/(1+x*x),x)-(1-x*x)/(1+x*x)**2)==0
    assert s.expand((1+x*x)**2-(1-x*x))==x**4+3*x*x
    assert s.expand((1+x*x)**2+(1-x*x))==x**4+x*x+2
    return {'inherited_source_twist_algebra':inherited,'weighted_current_identities':controls,
            'scalar_cutoff_edge_identity':True,'bounded_creator_function_global_Lipschitz_one':True,
            'complete_creator_isometry_error_coefficient':4}


def exact_checks():return copy.deepcopy(_checks())
