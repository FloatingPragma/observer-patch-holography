"""Actual-clock filtered current bounds and complete finite Pauli records."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import isqrt

import sympy as s

from model.native_ordered_ground import source
from . import bounds


def filtered_current(n, cutoff, p_norm, q_norm):
    row = bounds.family(n)
    m = bounds.rational(cutoff, 'energy cutoff')
    r = bounds.rational(p_norm, 'P-field input norm')
    q = bounds.rational(q_norm, 'Q-field input norm')
    natural = m*r+q/n
    return {**row, 'spectral_multiplier_norm_upper': F(1),
            'spectral_multiplier_support_upper': m,
            'P_vacuum_norm_upper': r, 'Q_vacuum_norm_upper': q,
            'selection_torque_norm_upper': q/n,
            'filtered_current_norm_upper': row['p']**2*natural,
            'filtered_current_squared_norm_upper': row['p']**4*natural**2,
            'p_minus_two_filtered_current_norm_upper': natural,
            'same_H_B_ground_and_filter_required': True,
            'linear_force_residual_vanishes': False,
            'full_history_force_closure_claimed': False}


def score_error(normalization, score_bound, output_error, scale=F(1)):
    c = bounds.rational(normalization, 'current normalization', True)
    score = bounds.rational(score_bound, 'absolute score bound')
    delta = bounds.rational(output_error, 'complete-output half trace error')
    p = bounds.rational(scale, 'scale', True)
    return {'normalization': c, 'absolute_score_bound': score,
            'complete_output_half_trace_error': delta,
            'original_squared_norm_error_upper': 2*c*c*score*delta,
            'p_minus_four_squared_norm_error_upper': 2*c*c*score*delta/p**4,
            'success_conditioning_used': False}


def _simple(value):
    return s.expand(s.radsimp(s.simplify(value)))


def _zero(matrix):
    return all(_simple(value) == 0 for value in matrix)


def _quadratic_interval(value, bits=48):
    """Outward rational interval in Q(sqrt(257)); no floating evidence."""
    value = _simple(value)
    root = s.sqrt(257)
    coefficient = value.coeff(root)
    constant = _simple(value-coefficient*root)
    if not (constant.is_Rational and coefficient.is_Rational):
        raise ArithmeticError('finite source value left its exact quadratic field')
    a, b = F(int(constant.p), int(constant.q)), F(int(coefficient.p), int(coefficient.q))
    denominator = 1 << bits
    numerator = isqrt(257*denominator*denominator)
    low, high = F(numerator, denominator), F(numerator+1, denominator)
    return (a+b*low, a+b*high) if b >= 0 else (a+b*high, a+b*low)


def _nonnegative(value):
    value = _simple(value)
    return value == 0 or _quadratic_interval(value)[0] >= 0


def _tensor(label):
    pauli = {'I': s.eye(2), 'X': s.Matrix([[0, 1], [1, 0]]),
             'Y': s.Matrix([[0, -s.I], [s.I, 0]]), 'Z': s.diag(1, -1)}
    out = s.ones(1, 1)
    for axis in reversed(label):
        out = s.kronecker_product(out, pauli[axis])
    return out


def _pauli_records(effect, ground):
    """Every setting/pointer effect, including identity-axis minus zero labels."""
    labels = [''.join(x) for x in product('IXYZ', repeat=2)]
    coefficients = {label: _simple(s.trace(effect*_tensor(label))/4) for label in labels}
    if any(not value.is_Rational for value in coefficients.values()):
        raise ArithmeticError('finite detector coefficients must be exact rationals')
    reconstructed = sum((coefficients[label]*_tensor(label) for label in labels), s.zeros(4))
    if not _zero(reconstructed-effect):
        raise ArithmeticError('complete Pauli expansion differs from detector effect')
    effects = s.zeros(4)
    probabilities = s.Integer(0)
    mean_score = s.Integer(0)
    records = []
    for label in labels:
        for signs in product((-1, 1), repeat=2):
            projector = s.eye(4)
            parity = 1
            for site, (axis, sign) in enumerate(zip(label, signs)):
                if axis == 'I':
                    if sign == -1:
                        projector = s.zeros(4)
                else:
                    name = axis+'I' if site == 0 else 'I'+axis
                    projector = projector*(s.eye(4)+sign*_tensor(name))/2
                    parity *= sign
            record_effect = projector/16
            probability = _simple(s.trace(ground*record_effect))
            if not _nonnegative(probability):
                raise ArithmeticError('negative exact complete record probability')
            score = 16*coefficients[label]*parity
            effects += record_effect
            probabilities += probability
            mean_score += score*probability
            records.append({'setting': label, 'pointer_signs': list(signs),
                            'original_probability': str(probability), 'score': str(score)})
    expectation = _simple(s.trace(ground*effect))
    if not _zero(effects-s.eye(4)) or _simple(probabilities-1) != 0:
        raise ArithmeticError('complete source Pauli records omit an outcome')
    if _simple(mean_score-expectation) != 0:
        raise ArithmeticError('original Pauli score differs from current response')
    maximum = max(abs(16*value) for value in coefficients.values())
    return {'settings': 16, 'labelled_pointer_records': 64,
            'zero_probability_records_retained': True, 'records': records,
            'Pauli_coefficients': {k: str(v) for k, v in coefficients.items()},
            'score_absolute_bound': str(maximum),
            'mean_score_equals_original_effect_probability': True,
            'original_effect_probability': str(expectation),
            'complete_record_probability_sum': '1',
            'conditional_success_distribution_used': False,
            'native_probe_templates_are_inherited': True,
            'full_detector_pulse_tape_executed_here': False}


@lru_cache(maxsize=1)
def finite_source():
    """Exact selected two-spin joint-family source, full spectrum and all records."""
    n, p, b, j = 2, s.Rational(1, 4), s.Rational(1, 32), s.Integer(1)
    integer, scale = source.matrix(2, ((0, 1),), F(1), F(1, 32))
    marshall = s.diag(1, 1, -1, -1)
    hb = marshall*s.Matrix(integer.tolist())*marshall/int(scale)
    h0 = sum((_tensor(axis*2) for axis in 'XYZ'), s.zeros(4))/4-s.eye(4)/4
    selection = -b*(_tensor('XI')-_tensor('IX'))/2
    if not _zero(hb-h0-selection):
        raise ArithmeticError('finite force source differs from actual H_B')
    minimum = -(16+s.sqrt(257))/32
    maximum = (-16+s.sqrt(257))/32
    spectrum = hb.eigenvals()
    if spectrum != {minimum: 1, s.Integer(0): 2, maximum: 1}:
        raise ArithmeticError('full selected source spectrum differs')
    vector = (hb-minimum*s.eye(4)).nullspace()[0]
    norm = _simple((vector.H*vector)[0])
    ground = (vector*vector.H/norm).applyfunc(_simple)
    if not _zero(ground*ground-ground) or s.trace(ground) != 1:
        raise ArithmeticError('finite selected ground is not a complete rank-one state')
    projector = s.eye(4)-hb*(hb+s.eye(4))/(b*b)
    if not _zero(projector*projector-projector) or s.trace(projector) != 2:
        raise ArithmeticError('the full zero-source-energy filter changed')
    l = (hb-minimum*s.eye(4))/p
    energy = _simple(-minimum/p)
    if not _nonnegative(5-energy) or not _zero(l*projector-energy*projector):
        raise ArithmeticError('source filter lies outside its declared clock band')
    q = p*p*(_tensor('ZI')-2*_tensor('IZ'))/2
    momentum = -p*(_tensor('YI')+2*_tensor('IY'))/2
    current = (_tensor('ZX')-_tensor('XZ'))/4
    z = p*p*current  # G_1-G_0=1, so p^3 D^pG=p^2.
    torque = s.Rational(1, n)
    drift = s.I*(l*momentum-momentum*l)
    if not _zero(drift+z/(p*p)+torque*q):
        raise ArithmeticError('exact actual source current or selection sign changed')
    if not _zero(projector*z*ground+s.I*p*p*l*projector*momentum*ground
                 +p*p*torque*projector*q*ground):
        raise ArithmeticError('same-ground filtered current identity changed')
    if _zero(drift+z/(p*p)):
        raise ArithmeticError('the nonzero finite selection torque was dropped')
    current_norm_squared = _simple(s.trace(ground*z*projector*z))
    bound = filtered_current(n, F(5), F(3, 8), F(3, 32))
    rational_upper = s.Rational(bound['filtered_current_squared_norm_upper'].numerator,
                               bound['filtered_current_squared_norm_upper'].denominator)
    if not _nonnegative(rational_upper-current_norm_squared):
        raise ArithmeticError('actual current exceeds the certified filter bound')
    records = []
    # The second test is an explicit finite comparison observable. It does not
    # assert that this two-site profile supplies a continuum Laplacian.
    q_test = p*p*(_tensor('ZI')+_tensor('IZ'))/2
    for candidate in (s.Integer(0), s.Integer(1)):
        y = z+p*p*candidate*q_test
        c = p*p/2+p**4*abs(candidate)
        effect = (y*projector*y/(c*c)).applyfunc(_simple)
        eigenvalues = effect.eigenvals()
        if any(not _nonnegative(value) or not _nonnegative(1-value)
               for value in eigenvalues):
            raise ArithmeticError('complete two-effect response is not a contraction')
        success = _simple(s.trace(ground*effect))
        failure = _simple(1-success)
        if not _nonnegative(success) or not _nonnegative(failure):
            raise ArithmeticError('complete original effect probabilities are invalid')
        records.append({'candidate_coefficient': str(candidate),
                        'source_norm_normalization': str(c),
                        'full_effect_spectrum': sorted(str(x) for x in eigenvalues),
                        'complete_effect_probabilities': [str(success), str(failure)],
                        'original_squared_response_norm': str(_simple(c*c*success)),
                        'complete_Pauli_score_readout': _pauli_records(effect, ground),
                        'linear_force_equation_claimed': False})
    return {'family_index': n, 'scale': str(p), 'selection_field': str(b),
            'coupling': str(j), 'sites': 2, 'all_source_rows': 4,
            'actual_source_frame_identity_exact': True,
            'complete_source_spectrum': {str(k): v for k, v in spectrum.items()},
            'selected_ground_energy': str(minimum),
            'filter_rank': 2, 'filter_excitation_energy': str(energy),
            'filter_spectral_cutoff': '5', 'current_drift_identity_exact': True,
            'selection_torque_nonzero': True,
            'filtered_ground_current_identity_exact': True,
            'original_filtered_current_squared_norm': str(current_norm_squared),
            'complete_current_and_comparison_records': records,
            'finite_level_identified_as_continuum_mass_shell': False}


def certificate():
    return {'current_bounds': [filtered_current(n, F(5), F(1), F(1))
                               for n in (2, 4, 8, 16)],
            'complete_score_error_example': score_error(F(1, 32), F(4), F(1, 1000), F(1, 4)),
            'finite_source': finite_source(),
            'current_scale_on_fixed_filtered_vacuum_vectors': 'O(p^2)',
            'linear_force_cancellation_proved': False,
            'full_history_force_closure_proved': False,
            'sharp_spectral_branch_or_interaction_inferred': False}
