"""Exact two-body exchange, contact resolvents and finite spectral witnesses."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import sympy as sp

from model.native_correlation_sectors import sectors as source

X = sp.Symbol('x')
CUBIC = sp.Poly(X**3-30*X**2+280*X-776, X)
FACTORS = ((sp.Poly(X-12, X), 5), (sp.Poly(X-10, X), 11),
           (sp.Poly(X-8, X), 5), (sp.Poly(X**2-20*X+88, X), 1),
           (sp.Poly(X**2-20*X+90, X), 4),
           (sp.Poly(X**2-20*X+98, X), 4), (CUBIC, 5))
PROBE = ((0, 1), (2, 4))


def rational(value):
    if type(value) not in (int, F):
        raise ValueError('an exact rational, excluding bool and float, is required')
    return F(value)


def integer(value, lower=0, upper=None):
    if type(value) is not int or value < lower or (upper is not None and value > upper):
        raise ValueError('an exact integer in the declared range is required')
    return value


def fraction(value):
    return F(int(sp.numer(value)), int(sp.denom(value)))


def family():
    return source.source.local_family(source.source.local.Parameters(port_weights=(1,)*12))


def polynomial(matrix, coefficients):
    result = sp.zeros(matrix.rows)
    identity = sp.eye(matrix.rows)
    for coefficient in coefficients:
        result = result*matrix+coefficient*identity
    return result


def generators(graph):
    source.fair_graph(graph)
    v = len(graph.weights)
    if not 4 <= v <= 12:
        raise ValueError('full exact matrix construction requires four through twelve ports')
    pairs = tuple(combinations(range(v), 2))
    where = {pair: i for i, pair in enumerate(pairs)}
    one = sp.zeros(v)
    for a, b in graph.seams:
        one[a, a] += 1; one[b, b] += 1
        one[a, b] -= 1; one[b, a] -= 1
    two = sp.zeros(len(pairs))
    for j, pair in enumerate(pairs):
        two[j, j] = 2*len(graph.seams)
        for image, count in source.counts_step(graph, {pair: 1}, 2).items():
            two[where[image], j] -= count
    incidence = sp.Matrix([[int(i in pair) for pair in pairs] for i in range(v)])
    q = sp.eye(len(pairs))-incidence.T*(incidence*incidence.T).inv()*incidence
    adjacent = sp.diag(*(int(pair in graph.seams) for pair in pairs))
    return {'pairs': pairs, 'one': sp.ImmutableMatrix(one),
            'two': sp.ImmutableMatrix(two), 'incidence': sp.ImmutableMatrix(incidence),
            'primitive': sp.ImmutableMatrix(q), 'adjacent': sp.ImmutableMatrix(adjacent)}


@lru_cache(maxsize=1)
def data():
    result = generators(family())
    values = source.contrast(*PROBE, 12)
    h = sp.ImmutableMatrix([values.get(pair, 0) for pair in result['pairs']])
    result['probe'] = h
    return result


def spectral_data():
    d = data(); g = d['two']; q = d['primitive']; h = d['probe']
    one_poly = sp.Poly(X*(X-6)**5*(X**2-10*X+20)**3, X)
    primitive_poly = sp.Poly(sp.prod(p.as_expr()**m for p, m in FACTORS), X)
    minimal = sp.Poly(sp.prod(p.as_expr() for p, _ in FACTORS), X)
    if sp.Poly(d['one'].charpoly(X).as_expr(), X) != one_poly:
        raise ArithmeticError('native one-body factorization changed')
    if sp.Poly(g.charpoly(X).as_expr(), X) != one_poly*primitive_poly:
        raise ArithmeticError('complete native two-body factorization changed')
    if q*q != q or q*g != g*q or sp.trace(q) != 54 or q*h != h or (h.T*h)[0] != 4:
        raise ArithmeticError('primitive sector or actual seam probe changed')
    current = h; moments = []; columns = []
    for _ in range(12):
        columns.append(current)
        moments.append((h.T*current)[0]/4)
        current = g*current
    krylov = sp.Matrix.hstack(*columns)
    if krylov.rank() != 12 or polynomial(g, minimal.all_coeffs())*h != sp.zeros(66, 1):
        raise ArithmeticError('the twelve observable poles are not a minimal native response')
    coefficients = minimal.all_coeffs()
    numerator = sp.Poly.from_list([sum(coefficients[k]*moments[j-k] for k in range(j+1))
                                  for j in range(12)], X)
    if sp.gcd(numerator, minimal).degree() != 0:
        raise ArithmeticError('a claimed native pole has zero residue')
    return {'denominator': minimal, 'numerator': numerator,
            'moments': tuple(fraction(x) for x in moments),
            'primitive_characteristic': primitive_poly,
            'one_characteristic': one_poly}


@lru_cache(maxsize=1)
def response_polynomials():
    d = spectral_data()
    return d['numerator'], d['denominator']


def response(z):
    """<h,(z-G2)^(-1)h>/4 on the cyclic primitive sector.

    This projected response can be finite at an eigenvalue outside that sector.
    """
    z = rational(z); numerator, denominator = response_polynomials()
    den = denominator.eval(sp.Rational(z.numerator, z.denominator))
    if den == 0:
        raise ValueError('the specified value is an observed spectral pole')
    return fraction(numerator.eval(sp.Rational(z.numerator, z.denominator))/den)


def _interval_multiply(a, b):
    values = tuple(x*y for x in a for y in b)
    return min(values), max(values)


def polynomial_interval(poly, interval):
    lower, upper = interval
    if type(lower) is not F or type(upper) is not F or lower > upper:
        raise ValueError('an ordered exact rational interval is required')
    result = (F(), F())
    for coefficient in poly.all_coeffs():
        result = _interval_multiply(result, interval)
        coefficient = fraction(coefficient)
        result = (result[0]+coefficient, result[1]+coefficient)
    return result


def pole_interval(bits=80):
    integer(bits, 12, 160)
    lower, upper = F(5), F(51, 10)
    if CUBIC.count_roots(5, sp.Rational(51, 10)) != 1:
        raise ArithmeticError('the lowest cubic root is not uniquely isolated')
    # C'(x)>0 on [5,51/10], and C has opposite endpoint signs.
    for _ in range(bits):
        middle = (lower+upper)/2
        if CUBIC.eval(sp.Rational(middle.numerator, middle.denominator)) < 0:
            lower = middle
        else:
            upper = middle
    return lower, upper


def pole_certificate():
    lower, upper = pole_interval()
    numerator, denominator = response_polynomials()
    n = polynomial_interval(numerator, (lower, upper))
    d = polynomial_interval(denominator.diff(), (lower, upper))
    if d[0] <= 0 <= d[1]:
        raise ArithmeticError('the simple pole derivative interval contains zero')
    residue = _interval_multiply(n, (1/d[1], 1/d[0]))
    if not F(1, 100) < residue[0] < residue[1] < F(1, 50):
        raise ArithmeticError('the visible lowest-pole weight is not certified')
    # sqrt(5)<9/4, hence 2 gamma_1>11/2, while E<51/10.
    if upper >= F(51, 10) or F(9, 4)**2 <= 5:
        raise ArithmeticError('subthreshold separation failed')
    return {'pole_interval': (lower, upper), 'normalized_residue_interval': residue,
            'independent_centered_threshold': '2*(5-sqrt(5))',
            'subthreshold_separation_lower': F(2, 5),
            'pole_multiplicity_in_primitive_sector': 5,
            'native_attempt_transfer_interval': (1-upper/60, 1-lower/60),
            'raw_stationary_score_variance': F(2, 33),
            'raw_pole_weight_interval': tuple(F(2, 33)*x for x in residue)}


def contact_identity(graph, z):
    """Exact reduced Birman--Schwinger/Woodbury check on bounded controls.

    The theorem is dimension independent; this inverse evaluator is not.
    """
    z = rational(z)
    if len(graph.weights) > 8:
        raise ValueError('the direct contact inverse control is bounded to eight ports')
    d = generators(graph); q = d['primitive']; g = d['two']; a = d['adjacent']
    n = g.rows; eye = sp.eye(n)
    e = sp.Matrix.hstack(*(eye[:, i] for i, pair in enumerate(d['pairs']) if pair in graph.seams))
    u = q*e
    h = q*g*q
    h0 = h+2*u*u.T
    z = sp.Rational(z.numerator, z.denominator)
    left_matrix = h-z*q+(eye-q)
    free_matrix = h0-z*q+(eye-q)
    if free_matrix.det() == 0:
        raise ValueError('the compressed independent reference resolvent has a pole')
    r0 = free_matrix.inv()*q
    contact = sp.eye(e.cols)-2*u.T*r0*u
    if left_matrix.det() == 0 or contact.det() == 0:
        raise ValueError('the correlated resolvent has a spectral pole')
    actual = left_matrix.inv()*q
    expected = r0+2*r0*u*contact.inv()*u.T*r0
    if actual != expected or h0-h != 2*q*a*q:
        raise ArithmeticError('the source contact resolvent identity failed')
    return {'primitive_dimension': int(sp.trace(q)), 'contact_dimension': e.cols,
            'contact_denominator': fraction(contact.det()), 'identity_verified': True}


def persistence_budget(gap, residual, observation_error):
    """Sufficient finite spectral-measure concentration, given proved inputs.

    A residual alone does not prove a volume-uniform gap, localization or units.
    """
    gap = rational(gap); residual = rational(residual); observation_error = rational(observation_error)
    if gap <= 0 or residual < 0 or observation_error < 0:
        raise ValueError('positive spectral window and nonnegative proved errors required')
    return {'outside_spectral_weight_upper': min(F(1), (residual/gap)**2),
            'bounded_observation_error_upper': min(F(2), 2*residual/gap)+observation_error}
