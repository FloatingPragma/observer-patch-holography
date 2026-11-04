"""Construct the complete C-duty memory coefficient from pixel-source data.

The inverse design uses the conditional MCPR law, never measured lepton masses.
Its target-informed architecture remains an explicit design input.
"""
from functools import lru_cache
import importlib.util
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as s

from model.native_flavor_orbit import algebra
from model.common_source import native_family_record_memory as memory
from model.common_source import native_gauge_family_records as gates
from model.common_source import native_repair_mass as repair
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
FORWARD = 'code/particles/hierarchy/certificates/R_P_source_audit_pixel_certificate.json'
MEASURED = 'code/particles/hierarchy/certificates/R_EW_tick_projection_certificate.json'
MCPR = 'code/particles/leptons/derive_charged_mcpr_completion_conditional.py'
BRANCHES = ('forward', 'measured_endpoint')


def source(branch='forward'):
    if branch == 'forward':
        row = strict_json.load(ROOT/FORWARD)
        return row['P_cand'], row['alpha_U_P_cand']
    if branch == 'measured_endpoint':
        row = strict_json.load(ROOT/MEASURED)['source_values']
        return row['P_star'], row['alpha_U']
    raise ValueError('unknown pixel-source branch')


@lru_cache(maxsize=1)
def _mcpr():
    spec = importlib.util.spec_from_file_location('_memory_bridge_mcpr', ROOT/MCPR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target(branch='forward'):
    module = _mcpr()
    with mp.workdps(100):
        pixel, alpha = map(mp.mpf, source(branch))
        geom = module.icosahedral_incidence_solve()
        fixed = module.response_fixed_point(pixel, alpha, module.register_dimensions(geom))
        return module.charged_triple(pixel, alpha, fixed, geom)


def _bisect(f, lo, hi):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    left = f(lo)
    if left*f(hi) >= 0:
        raise ArithmeticError('root bracket does not change sign')
    for _ in range(300):
        middle = (lo+hi)/2
        if left*f(middle) <= 0:
            hi = middle
        else:
            lo, left = middle, f(middle)
    return (lo+hi)/2


def design(branch='forward', inverse_branch=0):
    """Return finite numerical matrices; existence is certified separately.

    ``mass_matrix_over_v`` has the conditional MCPR masses/v as eigenvalues.
    ``coefficient_matrix`` is the unconditional C-duty matrix of trace p.
    The feedback amplitude ``gain_over_v`` multiplies every enabled pulse;
    disabled pulses are identity and their probabilities are retained.
    """
    if inverse_branch not in (0, 1):
        raise ValueError('the inverse has exactly two admitted branches')
    with mp.workdps(100):
        row = target(branch)
        masses = row['masses_over_v']
        total = mp.fsum(masses)
        lam = [m/total for m in masses]
        e2 = sum(lam[i]*lam[j] for i in range(3) for j in range(i))
        e3 = mp.fprod(lam)
        determinant = mp.mpf(2)/5-2*mp.sqrt(5)/25
        polynomial = lambda z: z**3-z*z+5*e2*z/4-e3/determinant
        w = _bisect(polynomial, '.0009', '.0011')
        u = _bisect(polynomial, '.07', '.071')
        v = _bisect(polynomial, '.92', '.94')
        A = v*(v+w)
        B = v*(1-v+w)+(v+w)*w
        discriminant = B*B-4*A*w
        p = (B+(-1 if inverse_branch == 0 else 1)*mp.sqrt(discriminant))/(2*A)
        q = 1-v/(1-p*(v+w))
        if not 0 < p < mp.mpf(1)/12 or not 0 < q < mp.mpf(1)/6:
            raise ArithmeticError('inverse design leaves the native event domain')
        theta_left = mp.asin(mp.sqrt(12*p))/mp.sqrt(mp.mpf(9)/8)
        theta_right = mp.asin(mp.sqrt(6*q))/mp.sqrt(3)
        projectors = memory._projectors(False)
        normalized = sum((float(x)*np.asarray(projectors[i], float)
                          for x, i in zip((u,v,w), (2,3,10))), np.zeros((3,3)))
        gain = total/p
        return {'source_branch': branch, 'inverse_branch': inverse_branch,
                'P': source(branch)[0], 'alpha_U': source(branch)[1],
                'weights': tuple(float(x) for x in (u,v,w)),
                'p': float(p), 'q': float(q),
                'theta_left': float(theta_left), 'theta_right': float(theta_right),
                'normalized_matrix': normalized, 'coefficient_matrix': float(p)*normalized,
                'mass_matrix_over_v': float(total)*normalized,
                'gain_over_v': float(gain), 'masses_over_v': tuple(float(x) for x in masses),
                'masses_over_E_star': tuple(float(x) for x in row['masses_over_e_star']),
                'v_over_E_star': float(row['v_over_e_star'])}


@lru_cache(maxsize=1)
def exact_identities():
    """Check the generic invariant and inverse equations over Q(sqrt(5))."""
    P = memory._projectors(False)
    u,v,w,p = s.symbols('u v w p', real=True)
    Y = u*P[2]+v*P[3]+w*P[10]
    determinant = s.Rational(2,5)-2*s.sqrt(5)/25
    if any(s.simplify(s.trace(P[i]*P[j])-s.Rational(1,5)) != 0
           for i,j in ((2,3),(2,10),(3,10))):
        raise ArithmeticError('projector overlap changed')
    expected_e2 = s.Rational(4,5)*(u*v+u*w+v*w)
    if s.simplify((s.trace(Y)**2-s.trace(Y*Y))/2-expected_e2) != 0:
        raise ArithmeticError('second invariant identity failed')
    if s.simplify(Y.det()-determinant*u*v*w) != 0:
        raise ArithmeticError('determinant identity failed')
    q = 1-v/(1-p*(v+w))
    a,b = p*(1-q),p*q
    D = 1-a*a+a*b
    remainder = -v*(v+w)*p*p+(v*(1-v+w)+(v+w)*w)*p-w
    # Equality of all three coefficient weights follows after u=1-v-w
    # and reduction by the inverse quadratic. No numerical eigenvalue fit.
    numerators = [s.cancel((1-a)/D-(1-p)-p*(1-v-w)).as_numer_denom()[0],
                  s.cancel(a*(1-a)/D-p*v).as_numer_denom()[0],
                  s.cancel(a*b/D-p*w).as_numer_denom()[0]]
    if any(s.rem(n, remainder, p).factor() != 0 for n in numerators):
        raise ArithmeticError('inverse memory coefficient identity failed')
    return {'trace': 'u+v+w', 'second_elementary_invariant': '4*(u*v+u*w+v*w)/5',
            'determinant': str(determinant)+' times u*v*w',
            'inverse_polynomial': '-v*(v+w)*p^2+(v*(1-v+w)+(v+w)*w)*p-w',
            'all_three_inverse_weight_identities_exact': True}


def direct_memory(design_row):
    """Independently sum all 512 retained binary/native records on 12 inputs."""
    p,q = design_row['p'],design_row['q']
    binary = ((1-12*p)*(1-6*q), (1-12*p)*6*q, 12*p*(1-6*q), 72*p*q)
    _, (gC,gE) = memory.geometry()
    transition = np.zeros((12,12))
    enabled = np.zeros((12,12))
    for old in range(12):
        for cc,cr,ec,er,_,raw in repair.native_records():
            for left in (False,True):
                for right in (False,True):
                    C = cc == 0 and cr == 1 and left
                    E = ec == 0 and er in (4,5) and right
                    new = memory._destination(old,C,E,2,gC,gE)
                    probability = float(raw)*binary[2*int(left)+int(right)]
                    transition[new,old] += probability
                    if C:
                        enabled[new,old] += probability
    equation = transition-np.eye(12)
    equation[-1,:] = 1
    rhs = np.zeros(12); rhs[-1] = 1
    stationary = np.linalg.solve(equation,rhs)
    weights = enabled@stationary
    projectors = memory._projectors(False)
    coefficient = sum((x*np.asarray(P,float) for x,P in zip(weights,projectors)),np.zeros((3,3)))
    return {'transition': transition, 'stationary': stationary,
            'enabled_transition':enabled,
            'coefficient_matrix': coefficient, 'retained_records_per_input': 512,
            'retained_old_inputs': 12}


def instrument_errors(row):
    left = gates.instrument(s.Float(row['theta_left'],17))
    right = gates.instrument(s.Float(row['theta_right'],17))
    return [left['completeness_error'],right['completeness_error']]


def averaging_bound(steps,time=1,*,p,gain,stationary_product=False):
    """Unhalved trace-norm bound including an untouched finite reference."""
    if type(steps) is not int or steps <= 0 or type(stationary_product) is not bool:
        raise ValueError('positive integral steps and explicit preparation required')
    if any(not math.isfinite(x) for x in (time,p,gain)) or not time >= 0 or not 0 <= p < 1 or gain < 0:
        raise ValueError('finite nonnegative time/gain and 0<=p<1 required')
    transient = 0 if stationary_product else 4*gain*time/(steps*(1-p))
    return transient+(8/(1-p)+4)*gain*gain*time*time/steps


def evolve(row,blocks,steps,time=1):
    """Complete C-duty channel with the actual correlated memory retained."""
    averaging_bound(steps,time,p=row['p'],gain=row['gain_over_v'])
    blocks = np.asarray(blocks,complex)
    if blocks.shape != (12,6,6) or not np.isfinite(blocks).all():
        raise ValueError('twelve finite six-dimensional CQ blocks required')
    law = direct_memory(row)
    enabled = law['enabled_transition']; disabled = law['transition']-enabled
    h = row['gain_over_v']*time/steps
    X = np.array([[0,1],[1,0]])
    U = [np.eye(6)+(np.cos(h)-1)*np.kron(np.eye(2),np.asarray(P,float))
         -1j*np.sin(h)*np.kron(X,np.asarray(P,float)) for P in memory._projectors(False)]
    for _ in range(steps):
        on = np.einsum('ij,jab->iab',enabled,blocks)
        blocks = np.einsum('ij,jab->iab',disabled,blocks)+np.asarray([
            unitary@block@unitary.conj().T for unitary,block in zip(U,on)])
    return blocks


source_design = design
