"""Wilson source clock on the literal twenty-component flux/Clifford code.

`symbol(a,p)` uses physical momentum p. `symbol_lattice(a,k)` uses
dimensionless lattice momentum k=a*p. This Hamiltonian is explicitly
declared; it is not a logarithm of the ordered routing tick.
"""
from itertools import product
from functools import lru_cache
import sympy as s

from ..spatial.dirac_route import clifford
from ..mobile_flux.model import observer
from ..mobile_flux.routing import odd_words

PROGRAM_ID='coherent-graded-flux-wilson-v1'


def exact_real(value,positive=False):
    value=s.sympify(value)
    if value.has(s.Float) or value.is_number is not True or value.is_real is not True or value.is_finite is not True:
        raise ValueError('finite exact real number required')
    if positive and value.is_positive is not True:raise ValueError('positive value required')
    return value


@lru_cache(maxsize=1)
def source_mass():
    counts={observer().base.law.syndrome(g).bit_count() for g in odd_words()}
    if counts!={2}:raise ArithmeticError('literal paired-flux face count changed')
    return s.Integer(2)


def symbol_lattice(spacing,momentum):
    a=exact_real(spacing,True)
    if len(momentum)!=3:raise ValueError('three momentum coordinates required')
    k=tuple(exact_real(x) for x in momentum);alpha,beta=clifford()
    mass=source_mass()+sum((1-s.cos(x) for x in k),s.Integer(0))/a
    return (sum((s.sin(x)*matrix/a for x,matrix in zip(k,alpha,strict=True)),s.zeros(4))+mass*beta).applyfunc(s.simplify)


def symbol(spacing,momentum):
    """Four by four symbol at physical momentum p; identical in all five faces."""
    a=exact_real(spacing,True)
    return symbol_lattice(a,tuple(a*exact_real(p) for p in momentum))


def energy(spacing,momentum):
    a=exact_real(spacing,True);p=tuple(exact_real(x) for x in momentum)
    if len(p)!=3:raise ValueError('three physical momentum coordinates required')
    mass=source_mass()+sum((1-s.cos(a*x) for x in p),s.Integer(0))/a
    return s.simplify(s.sqrt(sum((s.sin(a*x)**2/a**2 for x in p),s.Integer(0))+mass**2))


def projectors(spacing,momentum):
    """Return (positive, negative) exact spectral projectors."""
    h=symbol(spacing,momentum);e=energy(spacing,momentum)
    return tuple(((s.eye(4)+sign*h/e)/2).applyfunc(s.simplify) for sign in (1,-1))


def ground_covariance(spacing,momentum):
    return projectors(spacing,momentum)[1]


def continuum_bound(spacing,window):
    a=exact_real(spacing,True);p=exact_real(window)
    if p.is_nonnegative is not True:raise ValueError('nonnegative momentum window required')
    return a*a*p**3/2+3*a*p*p/2


def corners(spacing):
    a=exact_real(spacing,True)
    return [{'corner':list(bits),'positive_energy':str(source_mass()+2*sum(bits)/a)}
            for bits in product((0,1),repeat=3)]
