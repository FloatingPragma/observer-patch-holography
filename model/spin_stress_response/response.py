"""Exact finite identities on the unchanged spin carrier and energy origin."""
from fractions import Fraction as F

import sympy as s

from model.native_spin_clock_join import algebra, energy


def weighted_source(size, edges, weights, rates, **kwargs):
    """Actual weighted local action, including onsite and scalar terms."""
    if (len(weights) != size or len(rates) != size
            or any(isinstance(v, bool) or not isinstance(v, (int, s.Rational))
                   for v in tuple(weights) + tuple(rates))):
        raise ValueError('one exact rational weight and rate per site required')
    data = energy.source(size, edges, **kwargs)
    zero = s.zeros(2**size)
    e = data['density']
    local = [a*x for a, x in zip(weights, e, strict=True)]
    h = sum(local, zero)
    work = [a*x for a, x in zip(rates, e, strict=True)]
    current = {(x, y): weights[x]*weights[y]*j
               for (x, y), j in data['currents'].items()}
    for x in range(size):
        outgoing = sum((current[x, y] if x < y else -current[y, x]
                        for y in range(size) if y != x), zero)
        if s.I*energy.commutator(h, local[x]) + outgoing != zero:
            raise ArithmeticError('weighted energy continuity failed')
    phi = tuple(s.Rational((x+2)**2, x+1) for x in range(size))
    observable = sum((p*x for p, x in zip(phi, local, strict=True)), zero)
    explicit = sum((p*x for p, x in zip(phi, work, strict=True)), zero)
    flux = sum(((phi[y]-phi[x])*j for (x, y), j in current.items()), zero)
    if s.I*energy.commutator(h, observable) + explicit != flux + explicit:
        raise ArithmeticError('weighted work/flux identity failed')
    # The derivative in an arbitrary independent source direction is checked.
    f = tuple(s.Rational(1-x, x+2) for x in range(size))
    perturbation = sum((a*x for a, x in zip(f, e, strict=True)), zero)
    lam = s.Symbol('lambda', real=True)
    hl = h + lam*perturbation
    for (x, y), base in data['currents'].items():
        jl = s.I*energy.commutator((weights[x]+lam*f[x])*e[x],
                                   (weights[y]+lam*f[y])*e[y])
        if s.diff(jl, lam).subs(lam, 0) != (f[x]*weights[y]+weights[x]*f[y])*base:
            raise ArithmeticError('current source derivative failed')
    return {'hamiltonian': h, 'density': local, 'work': work,
            'currents': current, 'perturbation': perturbation,
            'deformed_hamiltonian': hl}


def first_law_bound(floor, halftrace):
    """Relative entropy <= epsilon^2/(mu-epsilon) on a faithful common carrier."""
    if (type(floor) is not F or type(halftrace) is not F
            or not 0 < floor <= 1 or not 0 <= halftrace < floor):
        raise ValueError('require exact 0 <= epsilon < mu <= 1')
    return halftrace**2/(floor-halftrace)


def modular_stability(floor, source_error):
    """Fixed complete decoder: actual spectral floor and logarithm error."""
    first_law_bound(floor, source_error)
    return {'actual_regional_floor': floor-source_error,
            'modular_generator_error': source_error/(floor-source_error)}


def response_bounds(duration, perturbation_norm, observable_norm, parameter,
                    source_error=F(0)):
    """Duhamel distance and a uniform second-order response remainder."""
    values = duration, perturbation_norm, observable_norm, parameter, source_error
    if (any(type(x) is not F for x in values) or duration < 0
            or perturbation_norm < 0 or observable_norm < 0
            or not 0 <= source_error <= 1):
        raise ValueError('exact nonnegative norms, duration and source error required')
    scale = abs(parameter)*duration*perturbation_norm
    # Duhamel twice at each real lambda, using unitarity rather than a series.
    return {'ideal_halftrace': min(F(1), scale),
            'complete_halftrace_with_source_error': min(F(1), scale+source_error),
            'linear_expectation_remainder': 2*observable_norm*scale**2,
            'source_expectation_error': 2*observable_norm*source_error}


def partial_first(rho):
    if rho.shape != (4, 4):
        raise ValueError('two-spin density required')
    return s.Matrix(2, 2, lambda a, b: sum(rho[2*a+c, 2*b+c] for c in range(2)))


def exact_response():
    owner = algebra.controls()
    h = energy.source(2, ((0, 1),))['hamiltonian']
    ground = -h
    omega = s.Matrix([-1, 3, -3, 1])/s.sqrt(20)
    hb = energy.source(2, ((0, 1),), field=s.Rational(3, 8),
                       staggered=(1, -1))['hamiltonian']
    if hb*omega != -s.Rational(9, 8)*omega:
        raise ArithmeticError('selected source ground changed')
    spectrum = {-s.Rational(9, 8): 1, s.Integer(0): 2, s.Rational(1, 8): 1}
    if hb.eigenvals() != spectrum:
        raise ArithmeticError('complete selection spectrum or unique ground changed')
    t = s.Symbol('t', real=True)
    u = s.eye(4)+(s.cos(t)+s.I*s.sin(t)-1)*ground
    rho = s.simplify(u*omega*omega.H*u.H)
    regional = s.simplify(partial_first(rho))
    x = energy.PAULI[0]
    target = (s.eye(2)+s.Rational(3, 5)*s.cos(t)*x)/2
    if s.simplify(regional-target) != s.zeros(2):
        raise ArithmeticError('actual regional marginal changed')
    if s.simplify(s.trace(h*rho)) != -s.Rational(9, 10):
        raise ArithmeticError('common-clock energy changed')
    y = s.diff(regional, t)*t
    # H_lambda=(1+lambda)H_J, with the same input and observation time.
    tangent = s.simplify(y.subs(t, s.pi/4))
    r = 3*s.sqrt(2)/10
    k = -(s.log((1+r)/2)+s.log((1-r)/2))/2*s.eye(2)
    k -= s.log((1+r)/(1-r))/2*x
    slope = s.simplify(s.trace(tangent*k))
    desired = 3*s.sqrt(2)*s.pi*s.log((10+3*s.sqrt(2))/(10-3*s.sqrt(2)))/80
    if s.simplify(slope-desired) != 0:
        raise ArithmeticError('entropy/modular response slope changed')
    # A physical regional energy alone cannot reproduce the modular tangent:
    # e_0=H_J/2 has constant expectation throughout this experiment.
    if s.simplify(s.diff(s.trace(h*rho)/2, t)) != 0 or tangent == s.zeros(2):
        raise ArithmeticError('stationarity countercontrol lost')
    return {'selection_field': owner['field'], 'HJ_energy_origin': '-1/4 per bond',
            'prepared_vector_big_endian': ['-1/sqrt(20)', '3/sqrt(20)',
                                          '-3/sqrt(20)', '1/sqrt(20)'],
            'selected_HB_ground_energy': '-9/8',
            'complete_HB_spectrum_with_multiplicity': {'-9/8': 1, '0': 2, '1/8': 1},
            'evolution_clock': 'H_J, hbar=1, common time t',
            'regional_density': '(I+(3/5)*cos(t)*X)/2',
            'uniform_regional_floor': '1/5',
            'complete_X_probabilities': ['(1+(3/5)*cos(t))/2', '(1-(3/5)*cos(t))/2'],
            'source_perturbation': 'H_lambda=(1+lambda)*H_J; initial vector and t fixed',
            'reference_time': 'pi/4', 'regional_tangent': '-3*sqrt(2)*pi*X/80',
            'entropy_and_modular_slope': str(desired),
            'HJ_mean_energy_all_times': '-9/10',
            'physical_energy_response': '0',
            'entropy_slope_strictly_positive': True,
            'source_selected_state_is_HJ_stationary': False,
            'modular_generator_equals_physical_local_energy': False}


def exact_currents():
    rows = []
    for size, edges, field, signs in (
            (3, ((0, 1), (1, 2)), 0, (1, -1, 1)),
            (3, ((0, 1), (0, 2), (1, 2)), 0, (1, 1, 1)),
            (4, ((0, 1), (0, 2), (1, 3), (2, 3)), s.Rational(3, 8), (1, -1, -1, 1))):
        weights = tuple(s.Rational(x+2, x+1) for x in range(size))
        rates = tuple(s.Rational(1-x, x+1) for x in range(size))
        data = weighted_source(size, edges, weights, rates, field=field, staggered=signs)
        nonzero = sum(j != s.zeros(2**size) for j in data['currents'].values())
        if not nonzero:
            raise ArithmeticError('noncommuting current control vanished')
        rows.append({'sites': size, 'edges': [list(e) for e in edges],
                     'field': str(field), 'weights': [str(a) for a in weights],
                     'weight_rates': [str(a) for a in rates],
                     'nonzero_current_pairs': nonzero,
                     'driven_continuity_and_current_response_exact': True})
    return rows
