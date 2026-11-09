"""Exact source marginal obstruction to a longitudinal mean-spin force."""
from fractions import Fraction as F

import sympy as s

from model.native_ordered_ground import source
from model.spin_charge_dispersion import dimer


def _rational(value):
    if type(value) not in (int, F):
        raise TypeError('exact rational required')
    return F(value.numerator, value.denominator) if type(value) is F else F(value)


def marginal(a, b, c, z, delta=F(0)):
    """Check the complete axial two-site density and its exact inequality."""
    a, b, c, z, delta = map(_rational, (a, b, c, z, delta))
    determinant = b*c-z*z
    if min(a, b, c, delta) < 0 or 2*a+b+c != 1 or determinant < 0:
        raise ValueError('positive normalized axial marginal required')
    m = (b-c)/2
    kappa = F(1, 4)-a-z/2
    d = -a-z
    if d < delta:
        raise ValueError('the declared energy improvement is unsupported')
    slack = kappa-m*m
    terms = (determinant, a*(2*d+F(1, 2)),
             (d-delta)*(d+delta+F(1, 2)))
    floor = delta*delta+delta/2
    if any(x < 0 for x in terms) or slack-floor != sum(terms):
        raise ArithmeticError('exact source marginal decomposition differs')
    return {'magnetization': m, 'bond_coefficient': kappa,
            'spin_dot_energy_coefficient': 3*(d+F(1, 4)),
            'energy_improvement': d, 'force_coefficient_excess': slack,
            'excess_lower': floor, 'nonnegative_remainder_terms': terms,
            'ground_selection_inferred_from_marginal_alone': False}


def joint_floor(n):
    if type(n) is not int or n < 2:
        raise ValueError('joint-family integer n>=2 required')
    field = F(1, n**5)
    delta = F(1, 120)-field/150
    if delta < F(13, 1600):
        raise ArithmeticError('joint source energy floor differs')
    return {'n': n, 'p': F(1, n*n), 'B_over_J': field,
            'dimer_trial_order_lower': F(12, 25),
            'exchange_energy_e_lower': F(31, 40)-field/50,
            'marginal_energy_improvement_lower': delta,
            'kappa_minus_m_squared_lower': delta*delta+delta/2}


def _spin_control():
    x = s.Matrix([[0, 1], [1, 0]])/2
    y = s.Matrix([[0, -s.I], [s.I, 0]])/2
    z = s.diag(1, -1)/2
    one = s.eye(2)
    pair = lambda a, b: s.kronecker_product(a, b)
    integer, scale = source.matrix(2, ((0, 1),), F(1), F(2, 3))
    marshall = s.diag(1, -1, 1, -1)
    hadamard = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    frame = pair(hadamard, hadamard)*marshall
    hb = (frame*s.Matrix(integer.tolist())*frame.H/scale).applyfunc(s.simplify)
    expected = sum((pair(a, a) for a in (x, y, z)), s.zeros(4))-s.eye(4)/4
    expected -= s.Rational(2, 3)*(pair(z, one)-pair(one, z))
    if hb != expected:
        raise ArithmeticError('the actual source frame changed')
    rho = s.Matrix([[0, 0, 0, 0], [0, s.Rational(9, 10), -s.Rational(3, 10), 0],
                    [0, -s.Rational(3, 10), s.Rational(1, 10), 0], [0, 0, 0, 0]])
    k = hb+s.Rational(4, 3)*s.eye(4)
    if k*rho != s.zeros(4) or hb.eigenvals() != {s.Rational(-4, 3): 1, s.Rational(1, 3): 1, s.Integer(0): 2}:
        raise ArithmeticError('actual full dimer ground or spectrum differs')
    expect = lambda a: s.simplify(s.trace(rho*a))
    m = expect(pair(z, one))
    cp, ct = expect(pair(z, z)), expect(pair(x, x))
    if (m, cp, ct) != (s.Rational(2, 5), -s.Rational(1, 4), -s.Rational(3, 20)):
        raise ArithmeticError('physical source marginal conventions differ')
    kap = -cp-ct
    p, beta = s.Rational(1, 4), s.Rational(2, 3)
    profile = (s.Integer(1), s.Integer(2))
    q = p*p*(profile[0]*pair(x, one)-profile[1]*pair(one, x))
    momentum = p*(profile[0]*pair(y, one)+profile[1]*pair(one, y))
    lap = ((profile[0]-profile[1])/p**2, (profile[1]-profile[0])/p**2)
    q_lap = p*p*(lap[0]*pair(x, one)-lap[1]*pair(one, x))
    scaled = k/p
    rq = s.I*scaled*q-2*m*momentum  # degree one, rather than cubic degree six
    rp = s.I*scaled*momentum+m*q_lap
    sf, df = p**3*sum(f*f for f in profile), (profile[0]-profile[1])**2
    iq, ip = s.im(expect(q*rq)), s.im(expect(momentum*rp))
    eq = (kap-m*m)*sf-kap*p**3*df/2+beta*m*sf/2
    ep = (kap-m*m)*p*df/2+beta*m*sf/(2*p*p)
    if s.simplify(iq-eq) or s.simplify(ip-ep):
        raise ArithmeticError('actual ground force-score identity differs')
    return {'full_source_dimension': 4, 'J': '1', 'B': '2/3', 'ground_energy': '-4/3',
            'ground_rank': 1, 'magnetization': str(m), 'kappa': str(kap),
            'force_coefficient_excess': str(kap-m*m),
            'Q_residual_ground_score': str(iq), 'P_residual_ground_score': str(ip),
            'source_degree': 1, 'source_frame_checked': True,
            'finite_dimer_substituted_for_cubic_ground': False}


def certificate():
    trial = dimer.pair_certificate()
    counts = dimer.torus_certificate(4)
    if trial['staggered_single_site_spin_mean_squared'] != '6/25' or F(12, 25)**2 >= F(6, 25):
        raise ArithmeticError('actual variational trial order floor differs')
    witness = marginal(F(0), F(9, 10), F(1, 10), F(-3, 10), F(3, 10))
    return {'source': 'spin-one-half cubic antiferromagnet, J>0',
            'state_and_clock': 'positive-field ground; p=n^-2, B=J n^-5, K_B/(Jp)',
            'source_trial': trial, 'source_torus_count_control': counts,
            'joint_family_floors': [joint_floor(n) for n in (2, 4, 16)],
            'liminf_kappa_minus_m_squared_lower': F(61, 14400),
            'exact_sharp_local_marginal_control': witness,
            'actual_dimer_force_control': _spin_control(),
            'excluded_approximation': 'replace every longitudinal spin by its staggered mean in unfiltered microscopic forces',
            'corrected_kinetic_coefficients_excluded': False,
            'infrared_sharp_component_excluded': False,
            'source_to_history_unbounded_force_norm_preservation_assumed': False,
            'physical_mass_or_species_assigned': False,
            'additional_native_detector_execution': False}
