"""Seventh-order unitary clock error and complete rational ring records."""
from functools import lru_cache
import json
from pathlib import Path
import sympy as s

from ..matter_induced_ring.model import physical_split, NORM_BOUND


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def clock_error(time, norm, subdivisions=1):
    time, norm = s.sympify(time), s.sympify(norm)
    if time.is_Rational is not True or norm.is_Rational is not True or norm < 0:
        raise ValueError('rational time and nonnegative rational norm required')
    if type(subdivisions) is not int or subdivisions < 1:
        raise ValueError('positive integer subdivisions required')
    return abs(time * norm)**7 / (100800 * subdivisions**6)


def polynomial_certificate():
    x, z, ez = s.symbols('x z ez', real=True)
    a, b = 1-x*x/10, x/2-x**3/120
    denominator = s.expand(a*a+b*b)
    require(denominator == 1+x*x/20+x**4/600+x**6/14400,
            'unitary denominator identity failed')
    # I_k = integral_0^1 exp(u*z)*u^k du, by integration by parts.
    integrals = [(ez-1)/z]
    for k in range(1, 7):
        integrals.append(ez/z-k*integrals[-1]/z)
    beta_integral = integrals[3]-3*integrals[4]+3*integrals[5]-integrals[6]
    p = lambda w: 1+w/2+w*w/10+w**3/120
    require(s.cancel(ez*p(-z)-p(z)+z**7*beta_integral/720) == 0,
            'exact exponential error identity failed')
    require(s.integrate(x**3*(1-x)**3, (x, 0, 1)) == s.Rational(1,140),
            'beta integral normalization failed')
    series = s.series(s.exp(z)*p(-z)-p(z), z, 0, 8).removeO()
    require(series == -z**7/100800, 'clock order changed')
    return {
        'clock': 'R33(tH)=(a+i*b)^(-1)*(a-i*b)',
        'a': 'I-(tH)^2/10', 'b': 'tH/2-(tH)^3/120',
        'positive_denominator': 'I+(tH)^2/20+(tH)^4/600+(tH)^6/14400',
        'exact_unitary_on_every_finite_self_adjoint_H': True,
        'error_identity': 'exp(z)*p(-z)-p(z)=-z^7/720*integral_0^1 exp(u*z)*u^3*(1-u)^3 du',
        'integral_absolute_bound_on_imaginary_axis': '1/140',
        'one_step_operator_error_bound': 'abs(t*B)^7/100800',
        'n_substep_operator_error_bound': 'abs(t*B)^7/(100800*n^6)',
        'unchanged_H_energy_and_commuting_constraints_preserved': True,
        'eigenphase_log_reinterpreted_as_physical_H': False,
    }


def peter_weyl_transfer():
    path = Path(__file__).parent.parent/'peter_weyl_links/results.json'
    source = json.loads(path.read_text())['complete_neutral_transfer_experiment']
    time = s.Rational(source['controlled_time'])
    bound = s.Rational(source['complete_source_operator_norm_bound'])
    leading = s.Rational(source['leading_amplitude_absolute_value'])
    remainder = s.Rational(source['outward_amplitude_remainder'])
    require(leading == s.Rational(151,3840)*time**5, 'source moment changed')
    require(remainder == leading/2, 'source Taylor enclosure changed')
    error = clock_error(time, bound)
    require(0 < error < leading/4, 'finite clock destroys the certified arrival')
    lower, upper = leading/4, 7*leading/4
    return {
        'source': 'actual all-family Peter-Weyl119 cubic Hamiltonian',
        'source_time': str(time), 'full_H_norm_bound': str(bound),
        'source_target_amplitude_leading_magnitude': str(leading),
        'full_exponential_amplitude_remainder': str(remainder),
        'finite_R33_operator_error': str(error),
        'clock_error_over_leading_amplitude': str(error/leading),
        'clock_error_less_than_quarter_leading': True,
        'exclusive_arrival_probability_interval': [str(lower**2), str(upper**2)],
        'exclusive_other_records_probability_interval': [str(1-upper**2), str(1-lower**2)],
        'local_arrival_probability_interval': [str(lower**2), '1'],
        'local_other_records_probability_interval': ['0', str(1-lower**2)],
        'same_spin_forbidden_amplitude': '0',
        'both_spin_preparations_and_reference_retained': True,
        'complete_channel': 'same physical preparation and two-outcome instrument as the consumed source',
        'small_event_control_uses_amplitude_error': True,
        'full_huge_source_matrix_explicitly_inverted': False,
    }


def family_transfer():
    path = Path(__file__).parent.parent/'family_transport/results.json'
    receipt = json.loads(path.read_text())
    pw = json.loads((path.parent.parent/'peter_weyl_links/results.json').read_text())
    require(receipt['source_and_static_reference']['source'] == pw['source'],
            'family and neutral clocks have different source descriptors')
    source = receipt['same_clock_full_H_probability_records']
    time = s.Rational(source['common_source_time'])
    require(source['common_outward_second_image_norm_bound'] == 377,
            'family second moment certificate changed')
    require(time == s.Rational(1,1508), 'family common clock changed')
    require({row['six_Y'] for row in source['five_prepared_reference_resolved_powers']}
            == {-4,1,6,-3,2}, 'native representation omitted')
    error = clock_error(time,1896,2)
    require(0 < error < time/96, 'finite family clock exceeds contrast margin')
    prepared = (19*time/96)**2
    background = (13*time/96)**2
    contrast = prepared-background
    require(contrast == time*time/48, 'local record contrast changed')
    return {'source':'same full Peter-Weyl H with static dual-R15 reference',
            'source_time':str(time), 'full_H_norm_bound':'1896', 'R33_substeps':2,
            'full_clock_error':str(error), 'clock_error_less_than_t_over_96':True,
            'prepared_local_arrival_probability_interval':[str(prepared),'1'],
            'prepared_local_other_probability_interval':['0',str(1-prepared)],
            'preparation_off_local_arrival_probability_interval':['0',str(background)],
            'preparation_off_local_other_probability_interval':[str(1-background),'1'],
            'strict_prepared_minus_unprepared_record_lower_bound':str(contrast),
            'same_clock_and_bounds_for_all_five_native_blocks':True,
            'all_reference_states_and_complements_retained':True,
            'unprepared_state_substituted_for_stationary_vacuum':False}


@lru_cache(maxsize=1)
def rational_ring_records():
    block = physical_split()
    h, g, p = block['h'], block['g'], block['p']
    size = h.rows
    time = s.Rational(1,64)
    require(h.T*g == g*h, 'physical Hamiltonian is not metric self-adjoint')
    initial = s.zeros(size,1)
    initial[p[0]] = 1
    require((initial.T*g*initial)[0] == 1, 'source preparation is not normalized')
    x = time*h
    x2 = x*x
    a = s.eye(size)-x2/10
    b = x/2-x2*x/120
    denominator = s.ImmutableMatrix(a*a+b*b)
    rhs = (a*a-b*b)*initial
    rhs_im = -2*a*b*initial
    # Exact rational linear algebra on the entire physical block, with no
    # compression to the five loop coordinates and no energy substitution.
    inverse = denominator.inv(method='DM')
    real = inverse*rhs
    imag = inverse*rhs_im
    require(denominator*real == rhs and denominator*imag == rhs_im,
            'finite clock solve residual is nonzero')
    require((real.T*g*real+imag.T*g*imag)[0] == 1,
            'finite source evolution lost norm')
    rows = []
    total = s.Integer(0)
    effects = []
    for name, j in (('singlet_loop',2), ('doublet_loop',4)):
        target = s.zeros(size,1)
        target[p[j]] = block['normalization'][j,j]
        effect = target*target.T*g
        effects.append(effect)
        require(effect*effect == effect and effect.T*g == g*effect,
                'loop effect is not an orthogonal physical projector')
        re, im = (target.T*g*real)[0], (target.T*g*imag)[0]
        probability = re*re+im*im
        require(0 < probability < 1, 'finite loop probability is not positive')
        rows.append({'outcome':name, 'effect_rank':1,
                     'amplitude_real':str(re), 'amplitude_imaginary':str(im),
                     'probability':str(probability)})
        total += probability
    require(effects[0]*effects[1] == s.zeros(size), 'loop effects overlap')
    complement = s.eye(size)-sum(effects, s.zeros(size))
    require(complement*complement == complement and s.trace(complement) == 70,
            'complete other-states outcome changed')
    require(0 < 1-total < 1, 'other states were discarded')
    rows.append({'outcome':'all_other_states', 'effect_rank':70,
                 'probability':str(1-total)})
    return {'full_physical_block_dimension':size, 'retained_loop_effect_ranks':[1,1,70],
            'source_time':str(time), 'finite_clock':'one exact R33 step of the full H',
            'original_H_norm_bound':str(NORM_BOUND),
            'exponential_comparison_error':str(clock_error(time,NORM_BOUND)),
            'exact_initial_norm_and_linear_system_residuals':True,
            'exact_full_final_norm':'1', 'outcomes':rows,
            'positive_loop_records':True, 'all_other_states_and_reference_retained':True,
            'source': '72-state complete compensated periodic lepton/Higgs block',
            'Peter_Weyl_spectrum_transfer':False,
            'unit_conversion_or_empirical_inputs':[],
            'energy_dependent_Schur_kernel_used_as_clock':False}
