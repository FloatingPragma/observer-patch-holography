"""Exact polynomial moments for first variations on three-dimensional balls."""
from math import factorial

import sympy as s

X = s.symbols('x y z', real=True)


def exact(value, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, float, complex)):
        raise ValueError('finite exact real input required')
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols or value.is_real is not True or value.is_finite is not True:
        raise ValueError('finite exact real input required')
    if positive and value.is_positive is not True:
        raise ValueError('strictly positive input required')
    if nonnegative and value.is_nonnegative is not True:
        raise ValueError('nonnegative input required')
    return value


def polynomial(value):
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols-set(X):
        raise ValueError('an exact spatial polynomial is required')
    try:
        poly = s.Poly(value, *X)
    except s.PolynomialError as error:
        raise ValueError('an exact spatial polynomial is required') from error
    if any(c.is_real is not True or c.is_finite is not True for c in poly.coeffs()):
        raise ValueError('finite real polynomial coefficients required')
    return poly


def angular_average(powers):
    """Integral of n^powers over S2 divided by4pi, with odd moments retained."""
    if len(powers) != 3 or any(type(p) is not int or p < 0 for p in powers):
        raise ValueError('three nonnegative integer powers required')
    if any(p % 2 for p in powers):
        return s.S.Zero
    halves = [p//2 for p in powers]
    numerator = 1
    for k in halves:
        numerator *= factorial(2*k)//(2**k*factorial(k))
    n = sum(halves)
    denominator = factorial(2*n+1)//(2**n*factorial(n))
    return s.Rational(numerator, denominator)


def integrate(value, radius, *, boundary=False):
    """Full polynomial surface/volume integral, without truncating its degree."""
    radius = exact(radius, positive=True)
    if type(boundary) is not bool:
        raise ValueError('Boolean boundary selector required')
    total = s.S.Zero
    for powers, coefficient in polynomial(value).terms():
        degree = sum(powers)
        scale = radius**(degree+2) if boundary else radius**(degree+3)/(degree+3)
        total += 4*s.pi*coefficient*angular_average(powers)*scale
    return s.expand(total)


def laplacian_at_origin(value):
    value = polynomial(value).as_expr()
    return s.simplify(sum(s.diff(value, x, 2) for x in X).subs(dict.fromkeys(X, 0)))


def radial_metric(seed, factor=1):
    """A symmetric radial-gauge metric tangent for any exact symmetric seed."""
    seed = s.ImmutableMatrix(seed)
    if seed.shape != (3, 3) or seed != seed.T:
        raise ValueError('an exact symmetric3x3 seed required')
    for value in seed:
        exact(value)
    factor = polynomial(factor).as_expr()
    x = s.Matrix(X)
    r2 = (x.T*x)[0]
    result = factor*(r2*seed+(x.T*seed*x)[0]*s.eye(3)-seed*x*x.T-x*x.T*seed)
    return s.ImmutableMatrix(result.applyfunc(s.expand))


def metric_tangent(value):
    h = s.ImmutableMatrix(value)
    if h.shape != (3, 3) or h != h.T:
        raise ValueError('symmetric3x3 spatial metric tangent required')
    for entry in h:
        polynomial(entry)
    origin = dict.fromkeys(X, 0)
    if h.subs(origin) != s.zeros(3) or any(h.diff(x).subs(origin) != s.zeros(3) for x in X):
        raise ValueError('normal-coordinate metric tangent needs zero value and first jet')
    if any(s.expand(v) != 0 for v in h*s.Matrix(X)):
        raise ValueError('complete radial gauge h_ij x^j=0 required')
    return h


def variation(h, f, radius, Fstar):
    """First variation only; f and h are derivatives at one flat reference."""
    h = metric_tangent(h)
    f = polynomial(f).as_expr()
    radius, Fstar = exact(radius, positive=True), exact(Fstar, positive=True)
    trace = s.trace(h)
    delta_volume = integrate(trace, radius)/2
    delta_area = integrate(trace, radius, boundary=True)/2
    scalar_surface = integrate(f, radius, boundary=True)
    scalar_volume = integrate(f, radius)
    delta_wald = 4*s.pi*(Fstar*delta_area+scalar_surface)
    delta_generalized_volume = delta_volume+3*scalar_volume/(2*Fstar)
    fixed_w = s.expand(delta_wald-8*s.pi*Fstar*delta_generalized_volume/radius)
    fixed_v = s.expand(delta_wald-8*s.pi*Fstar*delta_volume/radius)
    origin = dict.fromkeys(X, 0)
    ricci_scalar = sum(s.diff(h[i, j], X[i], X[j]) for i in range(3) for j in range(3))
    ricci_scalar -= sum(s.diff(trace, x, 2) for x in X)
    delta_Guu = s.simplify(ricci_scalar.subs(origin)/2)
    if s.simplify(laplacian_at_origin(trace)+4*delta_Guu/3) != 0:
        raise ArithmeticError('radial-gauge curvature identity failed')
    delta_F_laplacian = laplacian_at_origin(f)
    C = 4*s.pi*radius**4/15
    metric_euler_uu = 2*Fstar*delta_Guu-2*delta_F_laplacian
    leading = -2*s.pi*C*metric_euler_uu
    return {
        'delta_volume': s.simplify(delta_volume),
        'delta_area': s.simplify(delta_area),
        'delta_wald': s.simplify(delta_wald),
        'delta_generalized_volume': s.simplify(delta_generalized_volume),
        'fixed_generalized_volume_variation': s.simplify(fixed_w),
        'fixed_proper_volume_variation': s.simplify(fixed_v),
        'delta_Guu': delta_Guu,
        'delta_F_laplacian': delta_F_laplacian,
        'metric_euler_uu': metric_euler_uu,
        'leading_variation': s.simplify(leading),
        'remainder': s.simplify(fixed_w-leading),
        'C_ell': C,
    }


def fourth_derivative_bound(value, radius):
    """A coefficient bound for every fourth directional derivative on the ball.

    Expand(direction.d)^4 using|direction_i|<=1 and|x_i|<=radius.
    This is a conservative exact upper bound, independent of the moment code.
    """
    radius = exact(radius, positive=True)
    value = polynomial(value).as_expr()
    bound = s.S.Zero
    for a in range(5):
        for b in range(5-a):
            c = 4-a-b
            derivative = s.diff(value, X[0], a, X[1], b, X[2], c)
            coefficient = s.Rational(factorial(4), factorial(a)*factorial(b)*factorial(c))
            bound += coefficient*sum(abs(v)*radius**sum(p) for p, v in polynomial(derivative).terms())
    return s.simplify(bound)


def remainder_bound(radius, Fstar, Mh, Mf):
    radius, Fstar = exact(radius, positive=True), exact(Fstar, positive=True)
    Mh, Mf = exact(Mh, nonnegative=True), exact(Mf, nonnegative=True)
    return s.pi**2*radius**6*(s.Rational(3, 7)*Fstar*Mh+s.Rational(20, 21)*Mf)


def certificate():
    x, y, z = X
    seed = s.Matrix([[1, s.Rational(1, 3), 0], [s.Rational(1, 3), -2, 1], [0, 1, 3]])
    h2 = radial_metric(seed)
    f3 = 2+x+2*y*z+x*x-s.Rational(3, 5)*z*z+x*y*z
    exact_row = variation(h2, f3, s.Rational(1, 3), s.Rational(5, 2))
    if exact_row['remainder'] != 0:
        raise ArithmeticError('quadratic/cubic first-variation identity failed')
    rows = []
    for radius in (s.Rational(1, 2), s.Rational(1, 4), s.Rational(1, 8)):
        h = radial_metric(seed, 1+x*x+2*y*z)
        f = f3+x**4+2*x*x*y*y+z**6/7
        row = variation(h, f, radius, s.Rational(5, 2))
        Mh = fourth_derivative_bound(s.trace(h), radius)
        Mf = fourth_derivative_bound(f, radius)
        bound = remainder_bound(radius, s.Rational(5, 2), Mh, Mf)
        if s.simplify(bound-abs(row['remainder'])).is_nonnegative is not True:
            raise ArithmeticError('independent fourth-derivative remainder bound failed')
        rows.append({'radius': str(radius), 'Mh': str(Mh), 'Mf': str(Mf),
                     **{k: str(v) for k, v in row.items()}, 'remainder_upper': str(bound)})
    alpha, F = s.symbols('alpha F', positive=True)
    cancellation = s.solve(4-8*F*alpha/3, alpha)
    if cancellation != [3/(2*F)]:
        raise ArithmeticError('constant-scalar cancellation did not fix generalized volume')
    return {'first_variation_only': True,
            'reference': 'flat metric, constant positive Fstar, radial normal coordinates',
            'generalized_volume_coefficient': '3/(2*Fstar)',
            'coefficient_forced_by_constant_scalar_variations': True,
            'exact_low_degree_control': {k: str(v) for k, v in exact_row.items()},
            'higher_degree_controls': rows,
            'ordinary_fixed_volume_identified_with_generalized_volume': False,
            'finite_variation_or_nonlinear_equation_claimed': False}
