"""Exact source homothety controls for the fixed-reference dilation theorem.

No finite control establishes a continuum limit or an original nodal atom.
The finite dimer has no three-dimensional momentum interpretation.
"""
from fractions import Fraction as F

import sympy as s

from model.source_family_compatibility import charge_scaling, local_spectral
from model.source_family_compatibility.defect_budget import exact_source


def _rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational required')
    value = F(value.numerator, value.denominator) if isinstance(value, F) else F(value)
    if positive and value <= 0:
        raise ValueError('positive rational required')
    return value


def _momentum(value):
    value = tuple(value)
    if len(value) != 3:
        raise ValueError('three momentum coordinates required')
    return tuple(_rational(x) for x in value)


def physical_source(spacing, speed=1, ratio=8):
    """Actual homogeneous centered action; its normalized reference is fixed."""
    a, v, u = (_rational(x, positive=True) for x in (spacing, speed, ratio))
    return {'spacing': a, 'speed': v, 'c': v/a, 'U': u*v/a,
            'U_over_c': u, 'literal_hopping_magnitude': v/(2*a),
            'physical_energy_per_normalized_energy': v/a,
            'a_dependent_reference_selected': False}


def scale_coordinates(momentum, energy, factor):
    p = _momentum(momentum)
    e = _rational(energy)
    factor = _rational(factor, positive=True)
    return {'momentum': tuple(factor*x for x in p),
            'energy': factor*e, 'factor': factor}


def dense_power_bracket(spacing, power=1):
    """Exact interpolation bracket for a_n=n**(-power), without float roots."""
    a = _rational(spacing, positive=True)
    if a > 1 or type(power) is not int or power < 1:
        raise ValueError('spacing at most one and positive integer power required')
    lo, hi = 1, 2
    while a * hi**power <= 1:
        lo, hi = hi, 2*hi
    while hi-lo > 1:
        middle = (lo+hi)//2
        if a * middle**power <= 1:
            lo = middle
        else:
            hi = middle
    upper, lower = F(1, lo**power), F(1, (lo+1)**power)
    if not lower < a <= upper:
        raise ArithmeticError('exact refinement bracket failed')
    return {'n': lo, 'power': power, 'a_n': upper, 'a_next': lower,
            'spacing': a, 'lambda': a/upper,
            'lambda_lower': F(lo, lo+1)**power,
            'inverse_test_dilation': upper/a,
            'jacobian': (upper/a)**3}


def geometric_refinement(ratio, index):
    """Discrete dilation of a geometric family; no continuous-scale claim."""
    ratio = _rational(ratio, positive=True)
    if ratio >= 1 or type(index) is not int or index < 0:
        raise ValueError('ratio strictly between zero and one and nonnegative integer index required')
    return {'a_n': ratio**index, 'a_next': ratio**(index+1),
            'ratio': ratio, 'jacobian': ratio**-3,
            'continuous_dilation_covariance_inferred': False,
            'discrete_covariance_requires_whole_sequence_limit': True}


def shell_dilates(momentum, energy, velocity, factors):
    """Geometric control; invariant_mass_squared is rest-energy squared.

    No shell or physical unit conversion is assigned to the source.
    """
    p = _momentum(momentum)
    e, v = _rational(energy), _rational(velocity, positive=True)
    factors = tuple(_rational(x, positive=True) for x in factors)
    if not factors or len(set(factors)) != len(factors):
        raise ValueError('nonempty distinct positive dilations required')
    mass2 = e*e-v*v*sum(x*x for x in p)
    if mass2 <= 0:
        raise ValueError('a positive invariant mass squared is required')
    rows = []
    for factor in factors:
        row = scale_coordinates(p, e, factor)
        actual = row['energy']**2-v*v*sum(x*x for x in row['momentum'])
        if actual != factor*factor*mass2:
            raise ArithmeticError('dilation changed the Lorentz quadratic form')
        rows.append({**row, 'invariant_mass_squared': actual,
                     'joint_measure_weight_factor': factor**3})
    if len({row['invariant_mass_squared'] for row in rows}) != len(rows):
        raise ArithmeticError('positive-mass dilates are not disjoint shells')
    return {'original_invariant_mass_squared': mass2, 'dilates': rows,
            'hypothetical_shell_geometry_only': True,
            'actual_source_shell_proved': False}


def threshold_dilation(momentum, energy, velocity, target):
    """A future-cone support point can be scaled below any positive threshold."""
    p = _momentum(momentum)
    e, v, target = (_rational(x, positive=True) for x in (energy, velocity, target))
    mass2 = e*e-v*v*sum(x*x for x in p)
    if mass2 < 0:
        raise ValueError('future causal-cone support point required')
    factor = min(F(1, 2), target/(2*e))
    row = scale_coordinates(p, e, factor)
    scaled = factor*factor*mass2
    if not scaled < target*target:
        raise ArithmeticError('positive support threshold was not lowered')
    return {**row, 'invariant_mass_squared': scaled, 'target': target,
            'channel_support_dilation_invariance_required': True,
            'actual_source_causal_cone_proved': False}


def _simple(matrix):
    return matrix.applyfunc(s.simplify)


def default_source_scaling(factor=F(3, 2)):
    """All16 actual Fock states, full ground, and every charge-energy record."""
    factor = _rational(factor, positive=True)
    base, rho, a, projectors = local_spectral.dimer_source()
    scaled, _, _ = exact_source(2, ((0, 1, 1),), factor/8, factor)
    q = s.Rational(factor.numerator, factor.denominator)
    if _simple(scaled-q*base) != s.zeros(16):
        raise ArithmeticError('actual centered source fails homogeneity')
    e0 = (4-s.sqrt(17))/8
    if _simple(scaled*rho-q*e0*rho) != s.zeros(16):
        raise ArithmeticError('full source ground was changed by scaling')
    if _simple(rho*rho-rho) != s.zeros(16) or s.trace(rho) != 1:
        raise ArithmeticError('actual reference is not its rank-one projector')
    records = []
    complete = s.zeros(16)
    for label, sign, creator in (('addition', 1, a.H), ('removal', -1, a)):
        effects = s.zeros(16)
        for energy, projector in projectors:
            if _simple(scaled*projector-q*(e0+energy)*projector) != s.zeros(16):
                raise ArithmeticError('full energy projector was not transported')
            effect = _simple(creator.H*projector*creator)
            effects += effect
            weight = s.simplify(s.trace(rho*effect))
            records.append({'charge': label, 'signed_source_energy': str(sign*energy),
                            'signed_scaled_energy': str(s.simplify(sign*q*energy)),
                            'original_probability': str(weight)})
        if _simple(effects-creator.H*creator) != s.zeros(16):
            raise ArithmeticError('a charged energy complement was discarded')
        complete += effects
    if _simple(complete-s.eye(16)) != s.zeros(16):
        raise ArithmeticError('the full two-charge instrument is incomplete')
    weight_sum = s.simplify(sum(s.sympify(row['original_probability']) for row in records))
    if weight_sum != 1:
        raise ArithmeticError('original charge-energy weights fail normalization')
    return {'full_Fock_dimension': 16, 'source_c': F(1, 8), 'source_U': F(1),
            'U_over_c': F(8), 'positive_clock_factor': factor,
            'same_full_ground_projector': True, 'same_energy_projectors': True,
            'complete_charge_energy_records': records, 'total_original_probability': F(1),
            'zero_probability_records_retained': True,
            'conditional_probability_renormalization_applied': False,
            'finite_dimer_is_continuum_momentum_mode': False,
            'native_large_detector_tape_executed': False}


def _detach(value):
    if isinstance(value, F):
        return F(value.numerator, value.denominator)
    if isinstance(value, dict):
        return {key: _detach(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_detach(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_detach(item) for item in value)
    return value


def certificate():
    """Finite source checks and precise conditional scope; no new source premise."""
    return _detach({'actual_default_source_homogeneity': default_source_scaling(),
            'retained_rational_ground_controls': [charge_scaling.homogeneous_dimer(n) for n in (3, 4)],
            'physical_source_example': physical_source(F(1, 16)),
            'dense_refinement_examples': [dense_power_bracket(F(7, 100), p) for p in (1, 2)],
            'geometric_refinement_control': geometric_refinement(F(1, 2), 4),
            'shell_geometry_control': shell_dilates((1, 0, 0), 2, 1, (F(1), F(3, 2), F(2))),
            'future_cone_threshold_control': threshold_dilation((1, 0, 0), 2, 1, F(1, 100)),
            'fixed_reference_and_dimensionless_couplings_required': True,
            'full_joint_vague_limit_required': True,
            'dense_refinement_ratio_tending_to_one_suffices': True,
            'arbitrary_sparse_subsequence_excluded_by_theorem': False,
            'geometric_whole_limit_has_discrete_dilation_covariance': True,
            'geometric_limit_has_continuous_dilation_covariance': False,
            'geometric_positive_shell_requires_entire_mass_orbit': True,
            'positive_original_weight_positive_mass_shell_excluded_if_premises': True,
            'positive_full_channel_mass_threshold_excluded_if_future_cone_and_invariant_support': True,
            'every_unstable_resonance_excluded': False,
            'continuous_timelike_response_excluded': False,
            'actual_full_continuum_limit_proved': False,
            'actual_nodal_atom_positivity_proved': False,
            'running_reference_or_observable_dictionary_excluded': False,
            'physical_species_assigned': 0, 'physical_mass_assigned': None})
