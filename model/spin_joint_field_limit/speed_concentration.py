"""Uncut source speed concentration from one weighted susceptibility."""
from fractions import Fraction as F

import sympy as s

from model.native_ordered_ground import source


def _q(value, name, positive=False):
    if type(value) not in (int, F):
        raise TypeError(name+' requires an exact rational')
    value = F(value.numerator, value.denominator) if type(value) is F else F(value)
    if value < 0 or (positive and value == 0):
        raise ValueError(name+' has invalid sign')
    return value


def _interval(pair, name, positive=False):
    if type(pair) not in (tuple, list) or len(pair) != 2:
        raise TypeError(name+' requires two endpoints')
    lo, hi = (_q(x, name, positive) for x in pair)
    if lo > hi:
        raise ValueError(name+' interval is reversed')
    return lo, hi


def static_first_moment(kappa, magnetization, field_over_J,
                        s_weight, reciprocal_s_weight):
    """Static interval for <A_f Omega,K/(J s(P_space)) A_f Omega>."""
    kl, ku = _interval(kappa, 'bond coefficient')
    ml, mu = _interval(magnetization, 'staggered order')
    beta = _q(field_over_J, 'B/J', True)
    sl, su = _interval(s_weight, 'integral s |fhat|^2', True)
    il, iu = _interval(reciprocal_s_weight, 'integral |fhat|^2/s', True)
    if ku > F(1, 2) or mu > F(1, 2):
        raise ValueError('spin-one-half source bounds exceeded')
    return {'ratio_first_moment_interval':
            [kl*sl/2+beta*ml*il/2, ku*su/2+beta*mu*iu/2],
            'spatial_symbol': 's(k)=sqrt(2 sum_j(1-cos(k_j)))',
            'requires_same_source_ground_correlations': True,
            'projected_Stone_first_moment_equality_inferred': False}


def concentration(norm, first, inverse, speed, radius, normalization=1):
    """Admit interval data and bound original probability outside a speed tube."""
    wl, wu = _interval(norm, 'source norm squared', True)
    hl, hu = _interval(first, 'ratio first moment')
    il, iu = _interval(inverse, 'ratio inverse moment')
    c, eps, lam = (_q(x, name, True) for x, name in
                   ((speed, 'speed'), (radius, 'radius'),
                    (normalization, 'preparation normalization')))
    if wu > lam*lam:
        raise ValueError('preparation normalization does not bound the source norm')
    if hu == 0 or iu == 0:
        raise ValueError('positive-energy response requires positive moment upper bounds')
    if hu*iu < wl*wl:
        raise ValueError('intervals contradict positive spectral Cauchy-Schwarz')
    lower = max(F(0), hl-2*c*wu+c*c*il)
    upper = hu-2*c*wl+c*c*iu
    if upper < lower:
        raise ValueError('intervals contradict the positive weighted residual')
    outside = min(wu, (c+eps)*upper/(eps*eps))
    return {'weighted_residual_interval': [lower, upper],
            'outside_speed_tube_norm_squared_upper': outside,
            'original_preparation_probability_interval': [wl/(lam*lam), wu/(lam*lam)],
            'original_outside_probability_upper': outside/(lam*lam),
            'original_inside_probability_lower': max(F(0), wl-outside)/(lam*lam),
            'ratio_first_moment_above_twice_speed_upper': min(hu, 4*upper),
            'speed': c, 'radius': eps,
            'upper_energy_cutoff_required': False,
            'success_probability_divided_out': False,
            'actual_source_residual_vanishing_proved': False}


def joint_score_error(regularizer, coupling, s_min, energy_radius,
                      far_label_probability, output_l1_error,
                      input_trace_error=0, sampling_error=0):
    """Complete score error; spatial labels are exact finite translation labels."""
    zeta, j, sm = (_q(x, name, True) for x, name in
                  ((regularizer, 'regularizer'), (coupling, 'J'), (s_min, 's_min')))
    width, far, out, state, sample = (_q(x, name) for x, name in
        ((energy_radius, 'energy radius'), (far_label_probability, 'far-label probability'),
         (output_l1_error, 'output L1 error'), (input_trace_error, 'input trace error'),
         (sampling_error, 'sampling error')))
    if far > 1:
        raise ValueError('far-label probability exceeds one')
    return {'absolute_complete_score_error_upper':
            min(F(1), width/(zeta*j*sm)+far+out+state+sample),
            'score_range': [F(0), F(1)], 'all_failure_scores': F(0),
            'spatial_labels_assumed_exact': True,
            'spatial_labels_require_source_recentered_chart': True,
            'source_energy_origin_error_must_be_in_energy_radius': True}


def inverse_from_score(score, score_error, normalization, regularizer,
                       source_gap, coupling, s_range):
    """Upper susceptibility enclosure from a complete bounded score and true gap."""
    ql, qu = _interval(score, 'measured score')
    err = _q(score_error, 'score error')
    lam, zeta, gap, j = (_q(x, name, True) for x, name in
        ((normalization, 'normalization'), (regularizer, 'regularizer'),
         (source_gap, 'source excitation gap'), (coupling, 'J')))
    sm, sx = _interval(s_range, 'positive spatial symbol range', True)
    if qu > 1 or err > 1:
        raise ValueError('bounded score or error exceeds one')
    ql, qu = max(F(0), ql-err), min(F(1), qu+err)
    xl, xu = lam*lam*ql/zeta, lam*lam*qu/zeta
    gx = gap/(j*sx)
    return {'regularized_ratio_inverse_interval': [xl, xu],
            'ratio_inverse_interval': [xl, (1+zeta/gx)*xu],
            'ratio_gap_lower': gx,
            'regularization_relative_gap': zeta/gx,
            'known_original_score_factor': zeta/(lam*lam),
            'positive_spatial_symbol_range': [sm, sx],
            'successful_branch_has_zero_ground_component_required': True,
            'full_joint_label_and_preparation_certificates_required': True,
            'success_probability_divided_out': False}


def joint_regularization(n, momentum_outer, relative_error):
    """Conservative scale using B/J=n^-5 and s(k)<=|k|<=p*momentum_outer."""
    if type(n) is not int or n < 2:
        raise ValueError('joint-family integer n>=2 required')
    r = _q(momentum_outer, 'outer rescaled momentum', True)
    error = _q(relative_error, 'relative gap error', True)
    if error > 1:
        raise ValueError('relative gap error must not exceed one')
    gx = F(1, 8*n**3)/r
    return {'p': F(1, n*n), 'B_over_J': F(1, n**5),
            'spatial_symbol_upper': r/(n*n), 'ratio_gap_lower': gx,
            'regularizer': error*gx,
            'source_energy_regularizer_at_outer_symbol_over_J': error/(8*n**5),
            'uniform_scaling_limit_convergence_inferred': False}


def finite_source_control():
    """Exact full source dimer; the fixed ratio divisor is one, not momentum."""
    integer, scale = source.matrix(2, ((0, 1),), F(1), F(2, 3))
    x = s.Matrix([[0, 1], [1, 0]])/2
    y = s.Matrix([[0, -s.I], [s.I, 0]])/2
    z = s.diag(1, -1)/2
    one = s.eye(2)
    pair = s.kronecker_product
    h2 = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    frame = pair(h2, h2)*s.diag(1, -1, 1, -1)
    h = (frame*s.Matrix(integer.tolist())*frame.H/scale).applyfunc(s.simplify)
    expected = sum((pair(a, a) for a in (x, y, z)), s.zeros(4))-s.eye(4)/4
    expected -= s.Rational(2, 3)*(pair(z, one)-pair(one, z))
    if h != expected:
        raise ArithmeticError('actual source frame differs')
    rho = s.Matrix([[0, 0, 0, 0], [0, s.Rational(9, 10), -s.Rational(3, 10), 0],
                    [0, -s.Rational(3, 10), s.Rational(1, 10), 0], [0, 0, 0, 0]])
    k = h+s.Rational(4, 3)*s.eye(4)
    if k*rho != s.zeros(4) or k.eigenvals() != {s.Integer(0): 1,
            s.Rational(4, 3): 2, s.Rational(5, 3): 1}:
        raise ArithmeticError('actual full source ground or spectrum differs')
    a = pair(x, one)-pair(one, x)
    failure = s.eye(4)-a*a
    if failure*failure != failure or a.H*a+failure.H*failure != s.eye(4):
        raise ArithmeticError('all-input preparation is incomplete')
    energies = (s.Integer(0), s.Rational(4, 3), s.Rational(5, 3))
    projections = []
    for energy in energies:
        project = s.eye(4)
        for other in energies:
            if energy != other:
                project = project*(k-other*s.eye(4))/(energy-other)
        projections.append(project.applyfunc(s.simplify))
    expect = lambda op: s.simplify(s.trace(rho*op))
    inverse = sum((p/e for p, e in zip(projections[1:], energies[1:])), s.zeros(4))
    w, first, chi = (expect(a*op*a) for op in (s.eye(4), k, inverse))
    kap = -expect(pair(z, z)+pair(x, x))
    order = expect(pair(z, one))
    static = kap*4/2+s.Rational(2, 3)*order*2/2
    if (w, first, chi) != (s.Rational(4, 5), s.Rational(16, 15), s.Rational(3, 5)) or first != static:
        raise ArithmeticError('exact source moment identity differs')
    rows, total = [], s.zeros(4)
    score = s.Integer(0)
    for label, kraus in (('success', a), ('failure', failure)):
        for energy, project in zip(energies, projections):
            effect = kraus.H*project*kraus
            total += effect
            probability = expect(effect)
            rows.append({'preparation': label, 'energy': str(energy),
                         'original_probability': str(probability)})
            if label == 'success':
                score += probability/(1+6*energy)
    if total != s.eye(4) or sum(s.Rational(row['original_probability']) for row in rows) != 1:
        raise ArithmeticError('complete source records do not sum to identity')
    if score != s.Rational(4, 45) or expect(a*projections[0]*a) != 0:
        raise ArithmeticError('inverse score or ground exclusion differs')
    exact_gap = inverse_from_score((F(4, 45), F(4, 45)), 0, 1, F(1, 6),
                                  F(4, 3), 1, (1, 1))
    return {'J': F(1), 'B': F(2, 3), 'source_dimension': 4,
            'source_gap': F(4, 3), 'ratio_divisor': F(1),
            'ratio_divisor_is_spatial_momentum': False,
            'norm_squared': F(4, 5), 'first_moment': F(16, 15),
            'inverse_moment': F(3, 5), 'kappa': F(2, 5), 'order': F(2, 5),
            'regularizer': F(1, 6), 'complete_score': F(4, 45),
            'certified_inverse_upper': exact_gap['ratio_inverse_interval'][1],
            'complete_record_count': len(rows), 'records': rows,
            'all_input_instrument_completeness_exact': True,
            'additional_native_tape_executed': False,
            'finite_dimer_substituted_for_cubic_ground': False}


def certificate():
    source_control = finite_source_control()
    return {'finite_source_control': source_control,
            'exact_source_speed_tube': concentration((F(4, 5),)*2, (F(16, 15),)*2,
                (F(3, 5),)*2, F(4, 3), F(1, 3)),
            'joint_gap_regularization_example': joint_regularization(4, 1, F(1, 16)),
            'dynamic_unknown': 'weighted inverse susceptibility of the same source packet',
            'actual_ordered_source_concentration_proved': False,
            'physical_speed_or_particle_assigned': False}
