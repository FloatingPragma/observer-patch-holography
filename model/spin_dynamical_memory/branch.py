"""Original source spectral windows, finite response bases and readout budgets."""
from fractions import Fraction as F
from math import comb, isqrt
import sympy as s
from . import closure


def rational(x, name='value'):
    if type(x) not in (int, F):
        raise TypeError(name+' must be an exact rational')
    return F(x)


def gram_floor(matrix, entry_error=F(0)):
    """Certified floor with a full Hermitian positivity test before inversion."""
    g = s.Matrix(matrix)
    eta = rational(entry_error, 'entry error')
    if eta < 0 or not g.rows or g.rows != g.cols or g != g.H:
        raise ValueError('nonempty Hermitian Gram matrix and nonnegative error required')
    if g.has(s.Float):
        raise TypeError('exact matrix entries required')
    minors = [s.simplify(g[:j, :j].det()) for j in range(1, g.rows+1)]
    if any(value.is_positive is not True for value in minors):
        raise ValueError('strict Gram positivity is not certified')
    floor = F(g.det()/s.trace(g)**(g.rows-1))-g.rows*eta
    if floor <= 0:
        raise ValueError('entry errors consume the certified Gram floor')
    return floor


def krylov(clock, seeds, depth, resolution=F(1)):
    """Exact finite source compression, retaining dependent-column information."""
    if type(depth) is not int or depth < 1:
        raise ValueError('positive integer response depth required')
    k, seed = s.Matrix(clock), s.Matrix(seeds)
    if seed.rows != k.rows or not seed.cols or seed.has(s.Float):
        raise ValueError('nonempty exact seed columns on the source space required')
    if any(s.simplify((seed[:,j].H*seed[:,j])[0]).is_nonnegative is not True
           or s.simplify(1-(seed[:,j].H*seed[:,j])[0]).is_nonnegative is not True
           for j in range(seed.cols)):
        raise ValueError('seed columns must have original probability at most one')
    bounded = closure.finite_defect(k, s.eye(k.rows), resolution)['resolvent']
    columns = s.Matrix.hstack(*(bounded**j*seed for j in range(depth)))
    basis = columns.columnspace()
    if not basis:
        raise ValueError('the source has no nonzero prepared response')
    phi = s.Matrix.hstack(*basis)
    g = s.simplify(phi.H*phi)
    mu = gram_floor(g)
    c = s.simplify(phi.H*bounded*phi)
    t = s.simplify(g.inv()*c)
    remainder = s.simplify(bounded*phi-phi*t)
    residual_gram = s.simplify(phi.H*bounded**2*phi-c*g.inv()*c)
    assert s.simplify(residual_gram-remainder.H*remainder) == s.zeros(phi.cols)
    coordinates = s.simplify(g.inv()*phi.H*seed)
    assert phi*coordinates == seed
    for power in range(2*depth):
        assert s.simplify(seed.H*bounded**power*seed-coordinates.H*g*t**power*coordinates) == s.zeros(seed.cols)
    return {'S': bounded, 'source_columns': columns, 'basis': phi, 'Gram': g,
            'Gram_floor': mu, 'compression_in_Gram_coordinates': t,
            'seed_coordinates': coordinates, 'residual_Gram': residual_gram,
            'rank': phi.cols, 'input_column_count': columns.cols,
            'moments_matched_through': 2*depth-1,
            'polynomial_preparation_original_success_lower': mu/phi.cols}


def bernstein_coefficients(values):
    vals = tuple(rational(x, 'Bernstein value') for x in values)
    if len(vals) < 2 or any(not 0 <= x <= 1 for x in vals):
        raise ValueError('at least two bounded [0,1] sample values required')
    n = len(vals)-1
    return tuple(F(comb(n, j))*sum(((-1)**(j-k)*comb(j, k)*vals[k]
                                  for k in range(j+1)), F(0)) for j in range(n+1))


def ramp_values(degree, left, right, delta, upper=False):
    if type(degree) is not int or degree < 1 or type(upper) is not bool:
        raise ValueError('positive degree and literal ramp direction required')
    a, b, d = map(rational, (left, right, delta))
    if not 0 <= a <= b <= 1 or d <= 0:
        raise ValueError('ordered unit interval and positive boundary width required')
    values = []
    for j in range(degree+1):
        u = F(j, degree)
        if upper:
            value = max(F(0), 1-max(a-u, u-b, F(0))/d)
        else:
            value = max(F(0), min(F(1), (u-a)/d, (b-u)/d))
        values.append(value)
    return tuple(values)


def window_from_moments(moments, left, right, delta, errors=None):
    """Bernstein brackets retain both boundary resolution and raw score errors."""
    moments = tuple(map(rational, moments))
    if len(moments) < 2 or any(not 0 <= x <= 1 for x in moments):
        raise ValueError('bounded original moment estimates required')
    errors = (F(0),)*len(moments) if errors is None else tuple(map(rational, errors))
    if len(errors) != len(moments) or min(errors) < 0:
        raise ValueError('one nonnegative error for each original moment required')
    n, d = len(moments)-1, rational(delta)
    lower_coeff = bernstein_coefficients(ramp_values(n, left, right, d))
    upper_coeff = bernstein_coefficients(ramp_values(n, left, right, d, upper=True))
    w_upper = min(F(1), moments[0]+errors[0])
    # floor(sqrt(n)) keeps the reported bound rational and conservative.
    analytic = w_upper/(2*d*isqrt(n))
    first = sum(c*m for c, m in zip(lower_coeff, moments, strict=True))
    second = sum(c*m for c, m in zip(upper_coeff, moments, strict=True))
    efirst = sum(abs(c)*e for c, e in zip(lower_coeff, errors, strict=True))
    esecond = sum(abs(c)*e for c, e in zip(upper_coeff, errors, strict=True))
    lower = max(F(0), first-analytic-efirst)
    upper = min(w_upper, second+analytic+esecond)
    if lower > upper:
        raise ValueError('moment/error data give inconsistent interval brackets')
    return {'original_weight_upper': w_upper, 'original_interval_lower': lower,
            'original_interval_upper': upper, 'analytic_Bernstein_error': analytic,
            'lower_score_error': efirst, 'upper_score_error': esecond,
            'boundary_width': d, 'boundary_weight_proved_zero': False,
            'spectral_atom_inferred': False}


def fourier_preparation_budget(alpha, resolution, cutoff, cell_width,
                               commutator_bound, creator_bound=F(1), coefficient_error=F(0)):
    a, z, t, h, ca, norm, ec = map(rational, (alpha, resolution, cutoff, cell_width,
                                            commutator_bound, creator_bound, coefficient_error))
    if a < 0 or min(z, t, h) <= 0 or min(ca, norm, ec) < 0 or norm > 1:
        raise ValueError('valid positive scales and contraction/commutator bounds required')
    # pi>3; positive coefficient weights have total mass at most one.
    tail = min(F(1), 2*a/(3*z*t)) if a else F(0)
    vector_error = norm*(tail+ec)+ca*h/2 if a else F(0)
    return {'Gamma_mean': a/z, 'Fourier_tail_upper': tail,
            'selected_creation_column_error': vector_error,
            'original_moment_score_error': min(F(1), 2*norm*vector_error),
            'arbitrary_complement_channel_error_inferred': False,
            'native_instrument_state_and_sampling_errors_separate': True}


def positive_gap_filter(gap, resolution, power):
    g, z = rational(gap), rational(resolution)
    if min(g, z) <= 0 or type(power) is not int or power < 1:
        raise ValueError('positive gap, resolution and integer filter power required')
    attenuation = (z/(z+g))**power
    return {'ground_projection_operator_error': attenuation,
            'excited_original_weight_retention_factor_lower': (1-attenuation)**2,
            'direct_difference_LCU_normalization': 2,
            'direct_difference_LCU_retention_lower': (1-attenuation)**2/4,
            'unit_normalized_filter_native_dilation_executed': False,
            'continuum_gap_inferred': False}


def invariant_interval(energy_lower, energy_upper, momentum_lower, momentum_upper, speed):
    el, eu, pl, pu, c = map(rational, (energy_lower, energy_upper, momentum_lower, momentum_upper, speed))
    if not 0 <= el <= eu or not 0 <= pl <= pu or c <= 0:
        raise ValueError('ordered nonnegative energy/momentum and independent positive speed required')
    low, high = el*el-c*c*pu*pu, eu*eu-c*c*pl*pl
    return {'energy_momentum_invariant_interval': [low, high],
            'all_points_timelike_or_null': low >= 0,
            'conditional_squared_mass_interval': [low/c**4, high/c**4] if low >= 0 else None,
            'Lorentz_action_supplied_by_interval': False, 'physical_particle_mass_assigned': False}


def dimer_joint_response():
    h, omega, upper, x, z0 = closure._native_dimer()
    ident = s.eye(4)
    k = h+9*ident/8
    bounded = (k+s.eye(4)/8).inv()/8
    psi = (ident-bounded)*(x+z0)*omega/2
    swap = s.zeros(4)
    for mask in range(4):
        target = ((mask & 1) << 1) | ((mask >> 1) & 1)
        swap[target, mask] = 1
    assert swap**2 == ident and swap*h == h*swap and swap*omega == omega
    p0, ptop = omega*omega.H, upper*upper.H
    pmid = ident-p0-ptop
    labels = [(s.Rational(9, 8), 1, pmid), (s.Rational(9, 8), -1, pmid),
              (s.Rational(5, 4), 1, ptop), (s.Rational(5, 4), -1, ptop)]
    rows=[]
    for energy, sign, projector in labels:
        joint = projector*(ident+sign*swap)/2
        assert joint**2 == joint and joint == joint.H and k*joint == energy*joint
        weight = s.simplify((psi.H*joint*psi)[0])
        rows.append({'energy': str(energy), 'site_swap_character': sign, 'original_weight': str(weight)})
    assert [F(row['original_weight']) for row in rows] == [F(81, 4000), F(729, 4000), F(16, 121), F(0)]
    one = krylov(k, psi, 1, F(1, 8))
    two = krylov(k, psi, 2, F(1, 8))
    block = krylov(k, psi.row_join(swap*psi), 2, F(1, 8))
    assert one['residual_Gram'] != s.zeros(1)
    assert two['rank'] == 2 and two['residual_Gram'] == s.zeros(2)
    assert block['rank'] == 3 and block['residual_Gram'] == s.zeros(3)
    assert set(two['compression_in_Gram_coordinates'].eigenvals()) == {s.Rational(1, 10),s.Rational(1, 11)}
    phi = block['basis']
    g = block['Gram']
    representation = s.simplify(g.inv()*phi.H*swap*phi)
    assert phi*representation == swap*phi
    assert representation*block['compression_in_Gram_coordinates'] == block['compression_in_Gram_coordinates']*representation
    ritz = one['compression_in_Gram_coordinates'][0]
    assert ritz not in (s.Rational(1, 10), s.Rational(1, 11))
    return {'source': 'native_ordered_ground.source.matrix(2,((0,1),),J=1,B=3/8)',
            'state': '(1,3,3,1)/sqrt(20)', 'clock': 'K=H_B+9I/8', 'resolution': '1/8',
            'complete_joint_success_rows': rows,
            'preparation_failure_probability': '1/2', 'filter_failure_probability': '7999/48400',
            'total_original_success': str(s.simplify((psi.H*psi)[0])),
            'one_step_Ritz_value': str(ritz), 'one_step_Ritz_is_an_exact_source_level': False,
            'two_step_energy_rank': two['rank'], 'two_step_energy_Gram_floor': str(two['Gram_floor']),
            'joint_block_rank': block['rank'], 'joint_block_input_columns': block['input_column_count'],
            'joint_block_Gram_floor': str(block['Gram_floor']),
            'joint_polynomial_preparation_original_success_lower': str(block['polynomial_preparation_original_success_lower']),
            'joint_translation_compression_commutes_with_clock': True,
            'translation_is_only_a_finite_graph_symmetry': True,
            'native_pulse_tape_executed': False, 'continuum_pole_certified': False,
            'physical_mass_assignment': None}


def certificate():
    return {'actual_source_joint_response': dimer_joint_response(),
            'moment_preparation_example': fourier_preparation_budget(F(3,2),1,10000,F(1,10**6),4,1,F(1,10**6)),
            'positive_field_filter_example': positive_gap_filter(F(1,32),1,64),
            'continuous_score_hierarchy_converges': True,
            'uniform_discontinuous_CDF_convergence_claimed': False,
            'arbitrary_response_projection_implemented': False,
            'uniform_growing_Gram_floor_proved': False,
            'actual_continuum_pole_certified': False}
