"""Exact local source identities and finite selected-ground concentration controls."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
import sympy as sp

from model.native_ordered_ground import source
from model.spin_twist_closure import field_gap
from . import bounds


def comm(a,b):
    return a*b-b*a


@lru_cache(maxsize=1)
def local_certificate():
    one = sp.eye(2)
    pauli = (sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),
             sp.diag(1,-1))
    left = tuple(sp.kronecker_product(x,one)/2 for x in pauli)
    right = tuple(sp.kronecker_product(one,x)/2 for x in pauli)
    a,b,j,h = sp.symbols('a b J B',real=True)
    exchange = j*(sum((left[k]*right[k] for k in range(3)),sp.zeros(4))-sp.eye(4)/4)
    selection = -h*(left[0]-right[0])
    checks=[]
    for axis in range(3):
        observable = a*left[axis]+b*right[axis]
        transverse = sum((left[k]*right[k] for k in range(3) if k!=axis),sp.zeros(4))
        expected = -j*(a-b)**2*transverse
        if axis:
            expected += h*(a*a*left[0]-b*b*right[0])
        residual = comm(observable,comm(exchange+selection,observable))-expected
        if sp.simplify(residual)!=sp.zeros(4):
            raise ArithmeticError('positive-field source double commutator differs')
        checks.append({'axis':axis+1,'all_symbolic_matrix_entries_zero':True,
                       'selection_commutes_exactly':axis==0})
    # A compact bond observable may overlap two onsite selection terms.
    # Every nonzero nested onsite term is bounded by B/2; each edge has 12 choices.
    return {'all_three_axes':checks,'local_spin_norm':'1/2','transverse_bond_norm':'1/2',
            'exchange_overlap_coefficient':198,'selection_overlap_coefficient':6,
            'selection_kick_coefficient':'1/4',
            'finite_gap_source_substitution':field_gap.source_certificate()}


@lru_cache(maxsize=1)
def finite_ground_controls():
    edges = ((0,1),(1,2),(2,3),(0,3))
    x = np.array([[0,1],[1,0]],float)/2
    ident = np.eye(2)
    def at(site):
        value = np.array([[1.]])
        for k in range(4):
            value = np.kron(value,x if k==site else ident)
        return value
    # Native Marshall coordinates turn the staggered physical S1 into uniform X/2.
    axes = [at(k) for k in range(4)]
    rows=[]
    for n in (2,3,4):
        family = bounds.family(n)
        integer,scale = source.matrix(4,edges,F(1),family['field'])
        h = np.asarray(integer,dtype=float)/scale
        values,vectors = np.linalg.eigh(h)
        omega = vectors[:,0]; generator = h-values[0]*np.eye(16)
        observable = family['p']**3*sum((i+1)*axes[i] for i in range(4))
        observable = np.asarray(observable,dtype=float)
        centered = observable@omega-float(omega@observable@omega)*omega
        variance = float(centered@centered)
        energy = float(centered@generator@centered)
        if values[1]-values[0]+1e-12 < float(family['gap']):
            raise ArithmeticError('finite source violates imported gap lower bound')
        if variance > energy/float(family['gap'])+1e-12:
            raise ArithmeticError('finite source variance-gap inequality differs')
        rows.append({'scale':n,'sites':4,'full_dimension':16,
                     'full_gap_B_over_8_control_passed':True,
                     'centered_variance_gap_control_passed':True,
                     'floating_diagnostic_tolerance':'1/1000000000000',
                     'numerical_diagonalization_is_a_theorem_proof':False})
    return rows
