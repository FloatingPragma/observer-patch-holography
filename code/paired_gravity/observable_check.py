"""Producer-independent action, curvature, detector and data-chart checks."""
from fractions import Fraction as F
from functools import lru_cache
import mpmath as mp
import sympy as sp
from .format import need, keys, equal, rational


@lru_cache(maxsize=1)
def algebra():
    g = sp.symbols('g', real=True)
    k = sp.symbols('k', positive=True)
    u, v, j = sp.symbols('u v j')
    E = (u*u+k*(v-g*u)**2)/2-j*u
    S = sp.hessian(E, (u, v))
    need(sp.simplify(S.det()-k) == 0, 'positive constitutive determinant')
    need(sp.simplify(S[0, 0]-S[0, 1]**2/S[1, 1]-1) == 0, 'fixed Schur complement')
    need(sp.simplify(S.inv()[0, 0]-1) == 0, 'complete Gaussian marginal')
    need(sp.simplify(S.inv()[1, 0]-g) == 0, 'cross response')
    need(sp.expand(E+j*j/2-((u-j)**2+k*(v-g*u)**2)/2) == 0,
         'sourced energy is positive above its common minimum')
    # Reconstruct linearized curvature from h rather than checking copied coefficients.
    t, x, y, z = sp.symbols('t x y z', real=True)
    coordinates = (t, x, y, z)
    U, V = sp.Function('U')(x, y, z), sp.Function('V')(x, y, z)
    eta = (-1, 1, 1, 1)
    h = sp.diag(2*U, 2*V, 2*V, 2*V)
    trace = sum(eta[a]*h[a, a] for a in range(4))
    ricci = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            ricci[a, b] = sp.simplify(sum(eta[c]*(sp.diff(h[c, b], coordinates[c], coordinates[a])+
                sp.diff(h[c, a], coordinates[c], coordinates[b])-sp.diff(h[a, b], coordinates[c], 2))
                for c in range(4))/2-sp.diff(trace, coordinates[a], coordinates[b])/2)
    scalar = sum(eta[a]*ricci[a, a] for a in range(4))
    einstein = ricci-sp.diag(*eta)*scalar/2
    lap = lambda f: sum(sp.diff(f, q, 2) for q in (x, y, z))
    need(sp.simplify(scalar-2*lap(U)+4*lap(V)) == 0, 'linear scalar curvature')
    need(sp.simplify(einstein[0, 0]+2*lap(V)) == 0, 'temporal curvature')
    for a in range(1, 4):
        for b in range(1, 4):
            expected = sp.diff(U-V, coordinates[a], coordinates[b])-(lap(U-V) if a == b else 0)
            need(sp.simplify(einstein[a, b]-expected) == 0, 'spatial tensor selector')
    # Ward balance alone does not select g: the divergence is an identity.
    for b in range(1, 4):
        need(sp.simplify(sum(sp.diff(einstein[a, b], coordinates[a]) for a in range(1, 4))) == 0,
             'linear Ward identity for every constitutive member')
    r, mu, R = sp.symbols('r mu R', positive=True)
    exterior = mu/sp.sqrt(x*x+y*y+z*z)
    exterior_scalar = scalar.subs({U: exterior, V: g*exterior}).doit()
    need(sp.simplify(exterior_scalar) == 0, 'scalar vacuum does not select gamma')
    for a in range(1, 4):
        for b in range(1, 4):
            axis = einstein[a, b].subs({U: exterior, V: g*exterior}).doit().subs({x: r, y: 0, z: 0})
            expected = (1-g)*mu/r**3*((2 if a == 1 else -1) if a == b else 0)
            need(sp.simplify(axis-expected) == 0, 'nonzero exterior tensor selects gamma one')
    rho = 105*mu/(8*R**3)*(1-r*r/R**2)**2
    interior = mu/(16*R)*(35-35*(r/R)**2+21*(r/R)**4-5*(r/R)**6)
    need(sp.simplify(sp.integrate(r*r*rho, (r, 0, R))-mu) == 0, 'independently normalized source mass')
    need(sp.simplify(-sp.diff(r*r*sp.diff(interior, r), r)/r**2-rho) == 0, 'compact source produces its potential')
    for n in range(4):
        need(sp.simplify(sp.diff(interior, r, n).subs(r, R)-sp.diff(mu/r, r, n).subs(r, R)) == 0,
             'source boundary matching through third derivative')
    # Exact optical integrals complement the separate numerical-coordinate replay.
    need(sp.simplify(sp.diff(z/sp.sqrt(1+z*z), z)-(1+z*z)**sp.Rational(-3, 2)) == 0,
         'exact straight-ray optical primitive')
    need(sp.limit(z/sp.sqrt(1+z*z), z, sp.oo)-sp.limit(z/sp.sqrt(1+z*z), z, -sp.oo) == 2,
         'exact complete linear deflection integral')
    need(sp.simplify(sp.diff(t/sp.sqrt(2-t*t), t)-2/(2-t*t)**sp.Rational(3, 2)) == 0,
         'exact derivative of nonlinear ray at zero strength')
    return True


def exp_encloses(interval, exponent):
    need(type(interval) is list and len(interval) == 2, 'clock interval')
    a, b = map(rational, interval)
    # Independent alternating-series bounds for exp(-x), exact rational arithmetic.
    x = abs(exponent)
    need(0 <= x <= F(1, 100), 'small clock exponent domain')
    term = value = F(1)
    for n in range(1, 29):
        term *= -x/n
        value += term
    lower = value+term*(-x)/29
    upper = value
    if exponent > 0:
        lower, upper = 1/upper, 1/lower
    need(0 < a <= lower <= upper <= b and b-a < F(1, 10**45), 'outward clock interval')


def verify(data, row, sources):
    algebra()
    keys(row, 'branches finite_reads nonlinear_rays redshift outcome physical_establishment raw_measurements_reanalysed')
    equal([row['outcome'], row['physical_establishment'], row['raw_measurements_reanalysed']],
          ['SOURCE_CLASS_DOES_NOT_SELECT_PAIRED_GRAVITY', False, False], 'scientific outcome boundary')
    need(type(row['branches']) is list and len(row['branches']) == 4, 'all constitutive controls retained')
    for entry, g in zip(row['branches'], (-1, 0, 1, 2)):
        keys(entry, 'gamma stiffness schur marginal_covariance_coefficient clock_exponent clock_ratio optical_exponent asymptotic_linear_deflection exterior_tensor_coefficient comparisons')
        equal(entry['gamma'], g, 'constitutive member')
        need(type(entry['stiffness']) is list and len(entry['stiffness']) == 2 and
             all(type(a) is list and len(a) == 2 for a in entry['stiffness']), 'bounded action shape')
        S = sp.Matrix([[rational(x) for x in a] for a in entry['stiffness']])
        need(S.shape == (2, 2) and S == S.T and S.det() == 1 and S[1, 1] == 1, 'positive normalized two-field action')
        need(-S[1, 0]/S[1, 1] == g, 'source-derived cross response')
        equal([entry['schur'], entry['marginal_covariance_coefficient']], ['1', '1'], 'same full clock marginal')
        equal(entry['clock_exponent'], '-1/20000', 'same source clock calibration')
        exp_encloses(entry['clock_ratio'], F(-1, 20000))
        equal(entry['optical_exponent'], str(F(g+1, 10000)), 'both metric potentials in optics')
        equal(entry['exterior_tensor_coefficient'], str(1-g), 'vacuum tensor selector')
        # Integrate the transverse force independently; producer uses the primitive.
        with mp.workdps(60):
            bend = (g+1)*mp.mpf('0.0001')*mp.quad(lambda z: (1+z*z)**mp.mpf('-1.5'), [-mp.inf, 0, mp.inf])
            value = rational(entry['asymptotic_linear_deflection'])
            need(abs(bend-mp.mpf(value.numerator)/value.denominator) < mp.mpf('1e-55'), 'integrated optical force')
        need(type(entry['comparisons']) is list and len(entry['comparisons']) == 4, 'all overlapping fits retained')
        for actual, measured in zip(entry['comparisons'], data['deflection']['rows']):
            keys(actual, 'name measured sd prediction residual')
            d, sd = F(measured['gamma_minus_one']), F(measured['sd_gamma'])
            prediction = F(g+1, 2)
            equal(actual, dict(name=measured['name'], measured=str(1+d/2), sd=str(sd/2),
                              prediction=str(prediction), residual=str((g-1-d)/2)), 'measurement residual chart')
    red = data['redshift']
    equal(row['redshift'], dict(prediction='1', residual=str(1-F(red['ratio'])),
          profiled_statistical_residual=str(max(F(0), abs(1-F(red['ratio']))-F(red['systematic_bound']))/F(red['statistical_sd'])),
          reference_fraction=red['reference_fraction']), 'bounded redshift systematic retained')
    need(type(row['finite_reads']) is list and len(row['finite_reads']) == len(sources) == 16,
         'same-history finite read census')
    for actual, state in zip(row['finite_reads'], sources):
        keys(actual, 'projection emit receive clock_exponent clock_ratio log_stationary_clock_error optical')
        equal([actual['projection'], actual['emit'], actual['receive']],
              [state['projection'], [21, 33], [42, 33]], 'actual source writers consumed by detectors')
        need(all(type(state[key]) is list and len(state[key]) == 64 for key in ('u', 'w')),
             'bounded source field dimensions')
        u, w = [list(map(rational, state[key])) for key in ('u', 'w')]
        need(len(u) == len(w) == 64, 'complete source field dimensions')
        exponent = u[42]-u[21]
        equal(actual['clock_exponent'], str(exponent), 'same-history clock difference')
        exp_encloses(actual['clock_ratio'], exponent)
        eps_u, eps_w = [rational(a[1]) for a in state['error_intervals']]
        equal(actual['log_stationary_clock_error'], str(2*eps_u), 'source error propagated into clock')
        need(type(actual['optical']) is list and len(actual['optical']) == 4, 'all source optical reads')
        for optical, g in zip(actual['optical'], (-1, 0, 1, 2)):
            keys(optical, 'gamma exponent index log_stationary_error')
            # Reconstruct the two separate tetrad factors, not a copied index formula.
            lapse_log = -u[21]
            length_log = g*u[21]+w[21]
            equal([optical['gamma'], optical['exponent'], optical['log_stationary_error']],
                  [g, str(length_log-lapse_log), str(abs(g+1)*eps_u+eps_w)], 'common source metric and error')
            exp_encloses(optical['index'], length_log-lapse_log)
        # The actual finite read resolves the conformal and Einstein coefficients;
        # a nominal branch difference hidden inside solver errors is insufficient.
        conformal, einstein = actual['optical'][0], actual['optical'][2]
        need(rational(conformal['exponent'])+rational(conformal['log_stationary_error']) <
             rational(einstein['exponent'])-rational(einstein['log_stationary_error']),
             'stationary optical intervals must be disjoint')
    need(4*F(500, 499)**2+6*F(500, 499)**3 < 11, 'uniform nonlinear second-derivative bound')
    need(type(row['nonlinear_rays']) is list and len(row['nonlinear_rays']) == 4, 'nonlinear ray census')
    for actual, p in zip(row['nonlinear_rays'], (F(0), F(1, 10**6), F(1, 10**4), F(1, 1000))):
        keys(actual, 'p numerical_angle linear_angle rigorous_remainder')
        equal([actual['p'], actual['linear_angle'], actual['rigorous_remainder']],
              [str(p), str(2*p), str(11*p*p)], 'nonlinear ray domain and analytic bound')
        numerical = rational(actual['numerical_angle'])
        # Independent unbounded radial coordinate r/r_min=1+x^2.
        with mp.workdps(65):
            z = mp.mpf(p.numerator)/p.denominator
            def radial(x):
                if not x:
                    return 2/mp.sqrt(2-2*z)
                y = x*x
                d = mp.expm1(2*mp.log1p(y)-2*z*y/(1+y))
                return 2*x/((1+y)*mp.sqrt(d))
            expected = 2*mp.quad(radial, [0, 1, mp.inf])-mp.pi
            need(abs(expected-mp.mpf(numerical.numerator)/numerical.denominator) < mp.mpf('1e-44'),
                 'independent nonlinear null-ray integral')
        need(abs(numerical-2*p) <= 11*p*p, 'nonlinear read lies inside proved weak-field remainder')
