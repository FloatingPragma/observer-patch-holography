"""Original complete-record phase score for the configured balanced source."""
from fractions import Fraction as F
from math import isqrt

import sympy as s

from model.standard_model_source_map import koide_source as source


def _clean(a):
    return s.ImmutableMatrix(a.applyfunc(s.simplify))


def _exact_rational(value):
    if type(value) not in (int, F):
        raise ValueError('an exact integer or Fraction is required')
    return F(value)


def _sqrt2(digits=50):
    scale = 10**digits
    floor = isqrt(2*scale*scale)
    lo, hi = F(floor, scale), F(floor+1, scale)
    assert lo*lo < 2 < hi*hi
    return lo, hi


def kick_factor(tangent=F(1)):
    """Four times sin(h/2)^2, with h=2 atan(tangent)."""
    u = _exact_rational(tangent)
    if u <= 0:
        raise ValueError('a positive exact half-angle tangent is required')
    return 4*u*u/(1+u*u)


def raw_instrument_identity():
    """All 128 fine CP maps, including zero weights, for a formal unit phase."""
    rows = source.native_partition()
    z = s.Symbol('z', nonzero=True)
    exponents = (0, 1, -1, 0, 0)
    rotations = (s.eye(3), source.P, source.P**2, s.eye(3), s.eye(3))
    total = s.zeros(3)
    checks = 0
    for row in rows:
        weight = s.sympify(row['probability'])
        group = row['group']; r = rotations[group]; e = exponents[group]
        # Column vectorization, with conjugate(z)=z^-1 on the unit circle.
        phased = weight*s.kronecker_product(z**(-e)*r, z**e*r)
        unphased = weight*s.kronecker_product(r, r)
        if _clean(phased-unphased) != s.zeros(9):
            raise ArithmeticError('a measured fine CP map detects its scalar phase')
        total += weight*r.T*r
        checks += 1
    if total != s.eye(3):
        raise ArithmeticError('the fine instrument lost its complete probability')
    return {'fine_records_checked': checks,
            'zero_weight_records_retained': sum(s.sympify(r['probability']) == 0 for r in rows),
            'all_input_matrix_units_per_record': 9,
            'complete_CP_maps_identical_for_every_unit_phase': True,
            'arbitrary_reference_and_classically_adaptive_histories': True,
            'coherent_flag_access_included_in_identity': False}


def moment_identity():
    """A Laurent-polynomial identity for every phase, not fitted controls."""
    z = s.Symbol('z', nonzero=True)
    a, b = s.sqrt(source.WEIGHTS[0]), s.sqrt(source.WEIGHTS[1])
    c = a*s.eye(3)+b*(z*source.P+source.P**2/z)
    y = c*c/3
    cosine = (z**3+z**(-3))/2
    expected = s.Rational(25, 7776)*(17+4*s.sqrt(2)*cosine)
    if s.simplify(s.trace(y*y)-expected) != 0:
        raise ArithmeticError('the all-phase second action moment differs')
    normalized = c/a
    x = s.Symbol('x')
    poly = x**3-3*x*x+s.Rational(3, 2)*x+s.Rational(1, 2)-cosine/s.sqrt(2)
    if s.simplify(normalized.charpoly(x).as_expr()-poly) != 0:
        raise ArithmeticError('the normalized positive-root cubic differs')
    return {'trace_Y': '5/18',
            'trace_Y_squared': '25*(17+4*sqrt(2)*cos(3*theta))/7776',
            'nonzero_flag_probability': 'kappa*(1735-100*sqrt(2)*cos(3*theta))/23328',
            'normalized_root_cubic': 'x^3-3*x^2+(3/2)*x+1/2-cos(3*theta)/sqrt(2)',
            'action_eigenvalue_from_normalized_root': '(5/108)*x^2',
            'symbolic_all_phase_identity': True}


def complete_score(z=1, tangent=F(1)):
    """Eight original outcomes on I_2/2 tensor I_3/3, one finite kick."""
    u = _exact_rational(tangent)
    kappa = s.Rational(kick_factor(u).numerator, kick_factor(u).denominator)
    data = source.exact_source(z)
    y = data['action_coefficient']
    # A_0f=(u_f I+v_f zP+w_f conjugate(z)P^2)/sqrt(3).
    # Orthogonality of the three prep rows gives sum_f A_0f A_0f*=I.
    # Thus sum_f B_f*B_f=Y for B_f=A_0f* A_00 on every input.
    prep = source.analyzer_factors(1)[0]
    if _clean(prep[:3, :]*prep[:3, :].T) != s.eye(3):
        raise ArithmeticError('the complete source analyzer row Gram differs')
    a, b = s.sqrt(source.WEIGHTS[0]), s.sqrt(source.WEIGHTS[1])
    cosine = s.simplify(s.re(s.sympify(z)**3))
    probabilities = []
    for flag in range(8):
        if flag == 0:
            probability = 1-kappa*s.trace(y-y*y)/3
        else:
            row_u, v, w = prep[0, flag], prep[1, flag], prep[2, flag]
            # The three circulant coefficients of 3B_f are
            # ua+b(v+w), (ub+wa)z+vb/z^2, (ub+va)/z+wb*z^2.
            diagonal = ((row_u*a+b*(v+w))**2+(row_u*b+w*a)**2+(v*b)**2
                        +(row_u*b+v*a)**2+(w*b)**2)
            interference = 2*b*(v*(row_u*b+w*a)+w*(row_u*b+v*a))
            probability = kappa*(diagonal+interference*cosine)/9
        probabilities.append(s.simplify(probability))
    if s.simplify(sum(probabilities)-1) != 0:
        raise ArithmeticError('the complete unconditioned score lost probability')
    actual = s.simplify(sum(probabilities[1:]))
    cosine = s.simplify(s.re(s.sympify(z)**3))
    expected = kappa*(1735-100*s.sqrt(2)*cosine)/23328
    if s.simplify(actual-expected) != 0:
        raise ArithmeticError('the complete score disagrees with the spectral formula')
    decoded = s.simplify((1735-23328*actual/kappa)/(100*s.sqrt(2)))
    if decoded != cosine:
        raise ArithmeticError('original-probability inversion failed')
    return {'phase': str(z), 'half_angle_tangent': str(u), 'kick_factor': str(kappa),
            'input_state': 'I_2/2 tensor I_3/3',
            'records': [{'flag': i, 'probability': str(p)} for i, p in enumerate(probabilities)],
            'nonzero_flag_probability': str(actual),
            'decoded_cos_three_theta': str(decoded),
            'positive_root_chamber_certified': data['positive_root_chamber_certified'],
            'postselection_or_duty_renormalization': False,
            'ideal_finite_matrix_service_executed': True,
            'native_pulse_word_reexecuted_here': False}


def cosine_interval(probability_interval, tangent=F(1), event_error=F(0)):
    """Certified conditional inverse; caller supplies coverage and model error."""
    if type(probability_interval) is not tuple or len(probability_interval) != 2:
        raise ValueError('an original two-endpoint probability interval is required')
    lo, hi = map(_exact_rational, probability_interval)
    error = _exact_rational(event_error)
    if not 0 <= lo <= hi <= 1 or error < 0:
        raise ValueError('ordered probability bounds in [0,1] and nonnegative error required')
    k = kick_factor(tangent)
    probability = max(F(0), lo-error), min(F(1), hi+error)
    numerator = 1735-23328*probability[1]/k, 1735-23328*probability[0]/k
    root = _sqrt2()
    candidates = [n/(100*r) for n in numerator for r in root]
    raw = min(candidates), max(candidates)
    clipped = max(F(-1), raw[0]), min(F(1), raw[1])
    result = {'certified': clipped[0] <= clipped[1], 'source_model_required': True,
              'input_probability_interval': tuple(map(str, (lo, hi))),
              'event_probability_error': str(error),
              'score_coverage_or_native_error_proved_here': False,
              'inverse_Lipschitz_upper': str(F(23328)/(100*root[0]*k)),
              'cos_three_theta_interval': tuple(map(str, clipped)),
              'positive_root_chamber_certified': clipped[0] >= root[1]/2,
              'strictly_positive_and_nondegenerate_certified':
                  clipped[0] > root[1]/2 and clipped[1] < 1}
    if not result['certified']:
        result['positive_root_chamber_certified'] = False
        result['strictly_positive_and_nondegenerate_certified'] = False
    return result


def action_ratio_intervals(probability_interval, tangent=F(1), event_error=F(0), digits=24):
    """Isolate the positive action eigenvalues; no spacetime or pole assertion."""
    if type(digits) is not int or not 8 <= digits <= 60:
        raise ValueError('root precision must be an integer from 8 to 60')
    result = cosine_interval(probability_interval, tangent, event_error)
    result['action_ratios_certified'] = False
    if not result['certified'] or not result['strictly_positive_and_nondegenerate_certified']:
        return result
    lo, hi = map(F, result['cos_three_theta_interval']); root = _sqrt2()
    d = lo/root[1]-F(1, 2), hi/root[0]-F(1, 2)
    if not 0 < d[0] <= d[1] < root[0]/2-F(1, 2):
        return result
    x = s.Symbol('x'); isolated = []
    for endpoint in d:
        p = s.Poly(x**3-3*x*x+s.Rational(3, 2)*x-s.Rational(endpoint.numerator, endpoint.denominator), x)
        rows = p.intervals(eps=s.Rational(1, 10**digits))
        if len(rows) != 3 or any(m != 1 for _, m in rows):
            raise ArithmeticError('three simple real source roots were not isolated')
        isolated.append([(F(a), F(b)) for (a, b), _ in rows])
    low, high = isolated
    roots = [(low[0][0], high[0][1]), (high[1][0], low[1][1]), (low[2][0], high[2][1])]
    if any(not 0 < a <= b for a, b in roots):
        raise ArithmeticError('positive ordered root enclosures failed')
    masses = [(F(5, 108)*a*a, F(5, 108)*b*b) for a, b in roots]
    ratios = [(a/masses[0][1], b/masses[0][0]) for a, b in masses[1:]]
    result.update(action_ratios_certified=True,
                  normalized_positive_root_intervals=[list(map(str, r)) for r in roots],
                  action_eigenvalue_intervals=[list(map(str, r)) for r in masses],
                  action_ratios_to_smallest=[list(map(str, r)) for r in ratios],
                  physical_pole_ratios_identified=False)
    return result


def _rational_enclosure(expression, digits=35):
    """Exact outward enclosure for the affine sqrt(2) score used here."""
    expression = s.expand(expression)
    b = s.simplify(expression.coeff(s.sqrt(2))); a = s.simplify(expression-b*s.sqrt(2))
    if a.is_Rational is not True or b.is_Rational is not True:
        raise ValueError('an affine rational sqrt(2) expression is required')
    a, b = F(a), F(b); roots = _sqrt2(digits)
    return min(a+b*r for r in roots), max(a+b*r for r in roots)


def phase_ratio_enclosure(z=1, digits=24):
    """Ratios for a readable exact program phase; no empirical fit occurs here."""
    z = source._phase(z)
    cosine = s.simplify(s.re(z**3))
    if cosine.is_Rational is not True:
        raise ArithmeticError('the exact rational phase lost its rational triple cosine')
    score = s.Rational(2, 23328)*(1735-100*s.sqrt(2)*cosine)
    out = action_ratio_intervals(_rational_enclosure(score), digits=digits)
    out.update(program_phase=str(z), exact_cos_three_theta=str(cosine),
               phase_is_readable_specification_input=True,
               phase_source_selection_or_calibration_provenance_supplied=False)
    if cosine == 1:
        out.update(default_degenerate_spectrum=True,
                   exact_action_ratios_to_smallest=['1', '34+24*sqrt(2)'])
    return out


def tracial_analyzer_law():
    """Phase-independent complete law from the actual Fourier row and weights."""
    _, fourier, _ = source.analyzer_factors(1)
    rotations = (s.eye(3), source.P, source.P**2)
    if any(s.trace(r.T*t) != (3 if i == j else 0)
           for i, r in enumerate(rotations) for j, t in enumerate(rotations)):
        raise ArithmeticError('the actual C3 rotations are not trace orthogonal')
    weights = list(source.WEIGHTS)+[0]*3
    probabilities = [s.simplify(sum(s.conjugate(fourier[i,j])*fourier[i,j]*weights[j]
                                   for j in range(8))) for i in range(8)]
    expected = [s.Rational(5,54)]*3+[s.Rational(1,36),s.Rational(25,36)]+[0]*3
    if probabilities != expected:
        raise ArithmeticError('the phase-independent analyzer law differs')
    return [{'flag':i,'probability':str(p)} for i,p in enumerate(probabilities)]


def certificate():
    z = (s.Integer(99)+20*s.I)/101
    controls = [complete_score(phase) for phase in (s.Integer(1), z, s.conjugate(z))]
    scores = [s.sympify(row['nonzero_flag_probability']) for row in controls]
    if s.simplify(scores[0]-scores[1]) == 0 or s.simplify(scores[1]-scores[2]) != 0:
        raise ArithmeticError('the phase countercontrols lost their distinction')
    records = tracial_analyzer_law()
    y = [source.exact_source(phase)['action_coefficient'] for phase in (1, z)]
    moment_difference = s.simplify(s.trace(y[0]*y[0])-s.trace(y[1]*y[1]))
    if moment_difference == 0:
        raise ArithmeticError('counterfamily actions have equal second moments')
    interval = _rational_enclosure(scores[1])
    inverse = action_ratio_intervals(interval)
    if not inverse['action_ratios_certified']:
        raise ArithmeticError('the complete-record inverse failed its exact control')
    return {'schema': 'oph.charged_lepton_response.phase_readout.v1',
            'raw_phase_equivalence': raw_instrument_identity(),
            'all_phase_identity': moment_identity(),
            'tracial_analyzer_records_identical_for_distinct_phases': True,
            'complete_tracial_analyzer_records': records,
            'distinct_action_second_moment_difference': str(moment_difference),
            'complete_weak_score_controls': controls,
            'certified_action_ratio_readout': inverse,
            'configured_default_phase': '0',
            'default_action_ratios_to_smallest': ['1', '34+24*sqrt(2)'],
            'default_twofold_degeneracy': True,
            'source_parameter_readable_without_inference': True,
            'coherent_phase_selected_by_raw_repair_records': False,
            'coherent_phase_inaccessible_to_all_programs': False,
            'physical_common_clock_selected': False,
            'physical_pole_mass_or_nature_assignment': False,
            'empirical_inputs_used': False}
