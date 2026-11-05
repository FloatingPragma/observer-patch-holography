"""Exact finite-amplitude bounds; analytic justification is in PROOF.md."""
from fractions import Fraction as F
import hashlib
from pathlib import Path

from ..certified_reference.packet import SOURCE_PATHS
from ..plaquette_reference.preparation import error_theorem
from ..regional_modular_expansion.expansion import logarithm_tail


ROOT = Path(__file__).resolve().parents[2]
PRIOR = 'model/certified_reference/results.json'
DIMENSION = 16
RADIUS = F(1, 40)


def rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(name+' must be an exact rational')
    return F(value)


def source_bounds(prior):
    """Recompute the regional bounds from the retained literal control proof."""
    if prior.get('schema') != 'oph.certified_reference.v1' or prior.get('all_pass') is not True:
        raise ValueError('a passing literal-source certificate is required')
    if set(prior['source_sha256']) != set(SOURCE_PATHS):
        raise ValueError('literal source certificate custody coverage changed')
    for path, digest in prior['source_sha256'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('literal source certificate custody changed: '+path)
    domain = prior['domain']
    expected = {'geometry': [2, 2, 1], 'field_modes': 8, 'full_Fock_dimension': 256,
                'beta': '1/32', 'speed': '1/8', 'interaction': '1'}
    if any(domain.get(key) != value for key, value in expected.items()):
        raise ValueError('the complete interacting source domain changed')
    if domain.get('all_256_classical_records_and_every_export_abort_retained') is not True:
        raise ValueError('every source record and abort must be retained')
    error = F(prior['preparation_certificate']['complete_preparation_output_error'])
    if not 0 <= error < F(17, 10**13):
        raise ValueError('the proved literal preparation error changed')
    design = error_theorem()
    execution = (2+error)*error
    minimum = design['regional_minimum_weight']-execution
    full_minimum = design['full_minimum_weight']-execution
    bridge = (design['density_trace_error']+execution)/minimum
    tail = logarithm_tail(F(1, 32), F(3))
    total = bridge+tail
    expected_bounds = {
        'actual_regional_success_minimum_weight': minimum,
        'actual_full_success_minimum_weight': full_minimum,
        'regional_modular_operator_error': bridge,
        'regional_P8_remainder': tail,
        'regional_P8_total_operator_and_trace_zero_tangent_bound': total,
    }
    inherited = prior['modular_certificate']
    for key, value in expected_bounds.items():
        if F(inherited[key]) != value:
            raise ValueError('regional bound disagrees with independent recomputation: '+key)
    for key in ('regional_log_is_of_unnormalized_success_block',
                'Gibbs_normalization_scalar_retained_in_absolute_P8_comparison',
                'all_abort_outputs_remain_separate_without_success_renormalization'):
        if inherited.get(key) is not True:
            raise ValueError('success normalization or retained-output domain changed')
    # The full-state lift X tensor I/16 has norm ||X||/16. Its
    # positivity needs the full bound, not only the regional bound.
    lift_minimum = 16*full_minimum
    if min(minimum, lift_minimum) <= F(1, 20) or total >= F(2, 10**9):
        raise ValueError('the claimed common variation radius is not certified')
    return {'preparation_error': error, 'regional_minimum': minimum,
            'full_minimum': full_minimum, 'regional_lift_minimum': lift_minimum,
            'modular_polynomial_error': total, 'abort_probability_upper': error**2}


def remainder_bounds(amplitude, *, minimum, modular_error,
                     operator_norm=1, hs_squared=2, trace_norm=2):
    """Bounds for sigma +/- amplitude*Delta at fixed trace <= 1.

    Input norm values are certified *upper* bounds, except hs_squared in
    the returned lower bound: that lower bound requires its exact value.
    All retained basis directions have the exact norms 1, 2, 2.
    """
    amplitude = rational(amplitude, 'amplitude')
    minimum = rational(minimum, 'minimum')
    modular_error = rational(modular_error, 'modular_error')
    operator_norm = rational(operator_norm, 'operator_norm')
    hs_squared = rational(hs_squared, 'hs_squared')
    trace_norm = rational(trace_norm, 'trace_norm')
    if amplitude < 0 or minimum <= 0 or modular_error < 0:
        raise ValueError('positive spectral floor and nonnegative errors required')
    if min(operator_norm, hs_squared, trace_norm) < 0:
        raise ValueError('norm values must be nonnegative')
    if operator_norm == 0 and (hs_squared != 0 or trace_norm != 0):
        raise ValueError('inconsistent zero direction')
    lower = minimum-amplitude*operator_norm
    if lower <= 0:
        raise ValueError('the entire two-sided path must be strictly positive')
    quadratic = amplitude**2*hs_squared
    upper = quadratic/(2*lower)
    cubic = amplitude**3*operator_norm*hs_squared/(6*lower**2)
    return {'path_minimum': lower,
            'relative_entropy_lower': quadratic/2,
            'relative_entropy_upper': upper,
            'quadratic_Taylor_remainder_upper': cubic,
            'polynomial_finite_entropy_error_upper':
                amplitude*trace_norm*modular_error+upper,
            'symmetric_entropy_derivative_error_upper':
                trace_norm*modular_error+
                amplitude**2*operator_norm*hs_squared/(6*lower**2)}


def direction_certificate():
    """An exact rational basis of every Hermitian trace-zero 16x16 matrix.

    Diagonal directions E_ii-E_15,15 have Gram matrix I+11^T.
    Real and imaginary off-diagonal directions have Gram matrix 2I.
    Their supports prove all cross-sector Gram entries vanish.
    """
    diagonal = [(i, DIMENSION-1) for i in range(DIMENSION-1)]
    off_diagonal = [(i, j) for i in range(DIMENSION) for j in range(i+1, DIMENSION)]
    count = len(diagonal)+2*len(off_diagonal)
    if count != DIMENSION**2-1:
        raise ArithmeticError('the full Hermitian variation space is required')
    even_pairs = [(i, j) for i, j in off_diagonal
                  if i.bit_count() % 2 == j.bit_count() % 2]
    number_pairs = [(i, j) for i, j in off_diagonal
                    if i.bit_count() == j.bit_count()]
    return {'dimension': DIMENSION, 'real_dimension': count,
            'parity_preserving_real_dimension': len(diagonal)+2*len(even_pairs),
            'number_preserving_real_dimension': len(diagonal)+2*len(number_pairs),
            'parity_preserving_off_diagonal_pairs': even_pairs,
            'number_preserving_off_diagonal_pairs': number_pairs,
            'diagonal_pairs': diagonal, 'off_diagonal_pairs': off_diagonal,
            'each_off_diagonal_pair_has_real_and_imaginary_direction': True,
            'diagonal_Gram_determinant': DIMENSION,
            'off_diagonal_Gram_diagonal': 2,
            'full_Gram_determinant': DIMENSION*2**(2*len(off_diagonal)),
            'each_direction_operator_norm': 1,
            'each_direction_Hilbert_Schmidt_squared': 2,
            'each_direction_trace_norm': 2}
