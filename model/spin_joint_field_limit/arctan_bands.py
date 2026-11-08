"""Bounded charged preparations on the joint positive-field source clock."""
from fractions import Fraction as F
from math import isqrt

import sympy as s

from model.native_ordered_ground import source


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise TypeError('exact rational input required')
    out = F(value)
    if out < 0 or (positive and not out):
        raise ValueError('invalid rational sign')
    return out


def sqrt_interval(value, bits=64):
    x = rational(value)
    if type(bits) is not int or bits < 1:
        raise ValueError('positive integer precision required')
    d = 1 << bits
    k = isqrt(x.numerator*d*d//x.denominator)
    lo = F(k, d)
    hi = lo if lo*lo == x else F(k+1, d)
    assert lo*lo <= x <= hi*hi
    return lo, hi


def atan_interval(value, terms=32):
    x = rational(value)
    if x > 1 or type(terms) is not int or terms < 1:
        raise ValueError('alternating arctangent requires x<=1 and positive terms')
    total = sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)), F())
    next_size = x**(2*terms+1)/F(2*terms+1)
    return (total, total+next_size) if terms % 2 == 0 else (total-next_size, total)


def pi_interval():
    a, b = atan_interval(F(1, 5)), atan_interval(F(1, 239))
    result = (16*a[0]-4*b[1], 16*a[1]-4*b[0])
    assert 3 < result[0] < result[1] < F(22, 7)
    return result


def outward(interval, denominator=10**18):
    if type(denominator) is not int or denominator < 1:
        raise ValueError('positive integer denominator required')
    lo, hi = interval
    if lo > hi:
        raise ValueError('reversed interval')
    a = lo.numerator*denominator//lo.denominator
    b = -((-hi.numerator*denominator)//hi.denominator)
    return [F(a, denominator), F(b, denominator)]


def continuum_band(order_floor=F(1, 10)):
    m = rational(order_floor, positive=True)
    if m > F(1, 2):
        raise ValueError('spin-one-half order floor exceeds its source norm')
    d = m/44
    a, b = 3*d*d/2, 2/(d*d)
    assert a/F(12)+F(1, 4)/b == d*d/4
    return {'order_floor': m, 'configuration_second_moment_upper': F(6),
            'weak_commutator_lower': d,
            'scaled_energy_band': [a, b], 'spatial_momentum_band': [F(1, 2), F(1)],
            'each_original_charge_probability_lower': d*d/2,
            'complete_preparation_records': ['charge_minus', 'charge_plus', 'complement'],
            'full_created_scaled_first_moment_upper': F(1, 3),
            'joint_positive_field_clock_used': True,
            'old_zero_field_clock_substituted': False,
            'original_odd_CAR_identification': False,
            'sharp_shell_or_physical_mass_inferred': False}


def finite_commutator_lower(mean, variance_upper, commutator_norm_upper,
                           configuration_variance_upper):
    mu, var, c, vx = map(rational, (mean, variance_upper,
                                  commutator_norm_upper, configuration_variance_upper))
    sigma = sqrt_interval(var)[1]
    numerator = mu/(1+vx)-sigma-c/2
    return {'weak_commutator_lower': max(F(), numerator)*F(7, 44),
            'longitudinal_norm_error_upper': sigma,
            'numerator_lower': numerator,
            'source_moments_required': True}


def longitudinal_errors(n, bond_square_sum, cubic_sum, momentum_scale=F(1)):
    """Exact coefficient bounds from source profile sums in (9) and (12)."""
    if type(n) is not int or n < 2:
        raise ValueError('joint-family integer n>=2 required')
    edges, cubes = map(rational, (bond_square_sum, cubic_sum))
    r = rational(momentum_scale, positive=True)
    p = F(1, n*n)
    return {'p': p, 'field_over_J': F(1, n**5),
            'longitudinal_variance_upper': 2*p**6*n**5*edges,
            'longitudinal_configuration_commutator_upper':
                sqrt_interval(r)[1]*p**5*cubes/2}


def bounded_projection_weight(commutator_lower, canonical_p_norm_squared_upper,
                              outside_canonical_p_norm_squared_upper):
    """Transfer actual P leakage to an original bounded charge probability."""
    d = rational(commutator_lower)
    v = rational(canonical_p_norm_squared_upper, positive=True)
    outside = rational(outside_canonical_p_norm_squared_upper)
    if 2*d*d > v:
        raise ValueError('commutator floor contradicts the charged creator norm')
    eta = sqrt_interval(min(outside, v))[1]
    return {'each_original_charge_probability_lower': max(F(), d-eta)**2/(2*v),
            'outside_canonical_p_norm_upper': eta,
            'actual_joint_projection_and_charge_symmetries_required': True,
            'growing_Pauli_sum_normalization_used': False,
            'source_leakage_vanishing_proved': False}


def ward_tube_leakage(cutoff, u, b, delta_q, delta_p, tube_radius, outside_window):
    """Source Ward bound for the canonical P vector, with its outer tail."""
    m = rational(cutoff, positive=True)
    dq, dp, tail = map(rational, (delta_q, delta_p, outside_window))
    eps = rational(tube_radius, positive=True)
    if type(u) not in (int, F) or type(b) not in (int, F) or u*b <= 0:
        raise ValueError('source coefficients must have positive frequency squared')
    u, b = F(u), F(b)
    error = abs(b)*sqrt_interval(dq)[1]+m*sqrt_interval(dp)[1]
    return {'outside_canonical_p_norm_squared_upper': tail+m*error*error/(u*b*eps*eps),
            'source_frequency_squared': u*b,
            'original_filtered_cross_moment_required': True,
            'source_slack_saturation_assumed_or_proved': False}


def _tensor(left, right):
    # Native site zero is the least-significant source bit.
    return s.kronecker_product(right, left)


def _zero(matrix):
    return matrix.applyfunc(s.simplify) == s.zeros(*matrix.shape)


def finite_source_control():
    """All-input nine-record instrument on the actual joint n=2 dimer."""
    identity = s.eye(2)
    sx = s.Matrix([[0, 1], [1, 0]])/2
    sy = s.Matrix([[0, -s.I], [s.I, 0]])/2
    sz = s.diag(1, -1)/2
    hadamard = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    rotate = _tensor(hadamard, hadamard)
    integer, scale = source.matrix(2, ((0, 1),), F(1), F(1, 32))
    marshall = s.diag(1, 1, -1, -1)
    hb = rotate*marshall*s.Matrix(integer.tolist())*marshall*rotate/scale
    expected = sum((_tensor(a, a) for a in (sx, sy, sz)), s.zeros(4))
    expected -= s.eye(4)/4+(_tensor(sz, identity)-_tensor(identity, sz))/32
    if not _zero(hb-expected):
        raise ArithmeticError('actual source Hamiltonian or joint frame differs')
    root = s.sqrt(257)
    e0, e1 = -(16+root)/32, (-16+root)/32
    if hb.eigenvals() != {e0: 1, s.Integer(0): 2, e1: 1}:
        raise ArithmeticError('actual source complete spectrum differs')
    ground = (hb*(hb-e1*s.eye(4))/(e0*(e0-e1))).applyfunc(s.simplify)
    upper = (hb*(hb-e0*s.eye(4))/(e1*(e1-e0))).applyfunc(s.simplify)
    middle = (s.eye(4)-ground-upper).applyfunc(s.simplify)
    charge = _tensor(sz, identity)+_tensor(identity, sz)
    plus = _tensor(sx+s.I*sy, identity)-_tensor(identity, sx+s.I*sy)
    assert _zero(hb*charge-charge*hb)
    assert _zero(charge*plus-plus*charge-plus)
    p = s.Rational(1, 4)
    x = p*p*(_tensor(sx, identity)-_tensor(identity, sx))
    y = -p*p*(_tensor(sy, identity)-_tensor(identity, sy))
    canonical_p = p*(_tensor(sy, identity)+_tensor(identity, sy))
    longitudinal = p**3*(_tensor(sz, identity)-_tensor(identity, sz))
    assert _zero(canonical_p*x-x*canonical_p+s.I*longitudinal)
    assert _zero(canonical_p*y-y*canonical_p)
    theta = s.symbols('theta', real=True)
    rotated_x = x*s.cos(theta)+y*s.sin(theta)
    assert _zero((rotated_x**3-p**4*rotated_x).applyfunc(s.trigsimp))
    # f(X_theta) is exactly atan(p^2)/(pi*p^2) X_theta on all four rows.
    assert _zero((x-s.I*y)/p**2-plus)
    singlet = s.Matrix([0, 1, -1, 0])/s.sqrt(2)
    triplet = s.Matrix([0, 1, 1, 0])/s.sqrt(2)
    ps, pt = singlet*singlet.H, triplet*triplet.H
    parallel = s.eye(4)-ps-pt
    k = s.symbols('k', positive=True)
    cs, cp = s.symbols('cs cp', positive=True)
    def complement_reduce(matrix):
        return matrix.applyfunc(lambda value: s.simplify(s.expand(value).subs(
            {cs**2: 1-4*k*k, cp**2: 1-2*k*k})))
    bplus, bminus = k*plus, k*plus.H
    complement = cs*ps+pt+cp*parallel
    assert _zero(complement_reduce(bplus.H*bplus+bminus.H*bminus+
                                   complement.H*complement-s.eye(4)))
    singlet_weight = (1+16/root)/2
    expected_probabilities = {
        ('charge_plus', 'ground'): 0, ('charge_plus', 'middle'): 2*k*k*singlet_weight,
        ('charge_plus', 'upper'): 0, ('charge_minus', 'ground'): 0,
        ('charge_minus', 'middle'): 2*k*k*singlet_weight, ('charge_minus', 'upper'): 0,
        ('complement', 'ground'): (1-singlet_weight*(1-cs))**2,
        ('complement', 'middle'): 0,
        ('complement', 'upper'): singlet_weight*(1-singlet_weight)*(1-cs)**2}
    effect_sum = s.zeros(4)
    for record, kraus in (('charge_plus', bplus), ('charge_minus', bminus),
                          ('complement', complement)):
        for energy, projection in (('ground', ground), ('middle', middle), ('upper', upper)):
            effect = kraus.H*projection*kraus
            effect_sum += effect
            probability = s.trace(ground*effect)
            if s.simplify(probability-expected_probabilities[record, energy]) != 0:
                raise ArithmeticError('complete source probability differs')
    assert _zero(complement_reduce(effect_sum-s.eye(4)))
    assert s.simplify(s.expand(sum(expected_probabilities.values())-1).subs(
        cs**2, 1-4*k*k)) == 0
    weak_commutator = s.simplify(s.I*s.trace(ground*(canonical_p*(bplus+bminus)/2
                                                   -(bplus+bminus)*canonical_p/2)))
    assert s.simplify(weak_commutator-p*k/root) == 0
    p_norm_squared = s.simplify(s.trace(ground*canonical_p**2))
    assert s.simplify(p_norm_squared-p*p*(1-singlet_weight)) == 0
    assert _zero((s.eye(4)-middle)*canonical_p*ground)
    assert s.simplify(weak_commutator**2/(2*p_norm_squared)-
                      expected_probabilities['charge_plus', 'middle']) == 0
    return {'dimension': 4, 'p': F(1, 4), 'field_over_J': F(1, 32),
            'complete_source_spectrum': [str(e0), '0', '0', str(e1)],
            'charged_scaled_energy': str((16+root)/8),
            'charged_preparation_probability_formula':
                'atan(1/16)^2/(4*pi^2)*(1+16/sqrt(257))',
            'finite_weak_commutator_formula': 'atan(1/16)/(8*pi*sqrt(257))',
            'complete_record_count': 9, 'allocated_zero_probability_records': 5,
            'all_input_effect_sum_identity': True,
            'bounded_projection_probability_bound_saturated_on_middle_level': True,
            'canonical_P_outside_middle_norm_squared': F(),
            'whole_source_and_charge_rotated_together': True,
            'finite_control_has_continuum_momentum_assignment': False,
            'additional_native_pulse_tape_executed': False}


def dimer_intervals():
    atan = atan_interval(F(1, 16))
    pi = pi_interval()
    k = (atan[0]/(2*pi[1]), atan[1]/(2*pi[0]))
    k2 = (k[0]*k[0], k[1]*k[1])
    root = sqrt_interval(F(257))
    st = ((1+16/root[1])/2, (1+16/root[0])/2)
    cs = (sqrt_interval(1-4*k2[1])[0], sqrt_interval(1-4*k2[0])[1])
    charge = (2*k2[0]*st[0], 2*k2[1]*st[1])
    fail_ground = ((1-st[1]*(1-cs[0]))**2, (1-st[0]*(1-cs[1]))**2)
    fail_upper = (st[1]*(1-st[1])*(1-cs[1])**2,
                  st[0]*(1-st[0])*(1-cs[0])**2)
    rows = []
    for record in ('charge_minus', 'charge_plus', 'complement'):
        for label in ('ground', 'middle', 'upper'):
            interval = (F(), F())
            if record != 'complement' and label == 'middle':
                interval = charge
            if record == 'complement' and label == 'ground':
                interval = fail_ground
            if record == 'complement' and label == 'upper':
                interval = fail_upper
            rows.append({'preparation': record, 'energy_label': label,
                         'original_probability_interval': outward(interval)})
    assert sum(row['original_probability_interval'][0] for row in rows) <= 1
    assert sum(row['original_probability_interval'][1] for row in rows) >= 1
    return {'records': rows, 'charged_scaled_energy_interval':
            outward(((16+root[0])/8, (16+root[1])/8), 10**12),
            'scalar_bounds_are_exact_rational': True,
            'preparation_success_probability_divided_out': False}


def certificate():
    """Fresh rational data; no cached mutable scalar or source state is exposed."""
    return {'source_family': 'joint p=n^-2, B=J*n^-5; same H_B ground and clock',
            'continuum_band': continuum_band(),
            'actual_finite_source': finite_source_control(),
            'complete_finite_probabilities': dimer_intervals(),
            'finite_source_commutator_error': '(sigma+C/2)/(2*pi)',
            'bounded_Ward_projection_transfer': {
                'each_original_probability_lower': 'max(0,d-eta)^2/(2*||P Omega||^2)',
                'leakage_from_same_clock_Ward_tube_required': True,
                'growing_Pauli_sum_normalization_used': False,
                'shrinking_continuum_width_proved': False},
            'continuum_limit_error_rate_supplied': False,
            'full_history_linear_force_closure_proved': False,
            'physical_species': [], 'physical_mass_prediction': None}
