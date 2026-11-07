"""A finite operational test of the identity transfer-clock identification.

The transfer continuation is a proposed map on the working carrier in this
test. The native program does not execute that continuation.
"""
from __future__ import annotations

from fractions import Fraction as F

import sympy as s

from model.native_clock_observables import observables as observable
from model.native_clock_observables import packet as observable_packet
from model.native_record_response import response, transfer
from tools import strict_json


DEGREE = 10
TIME = F(1, 10)
NORM_UPPER = F(103, 200)
DERIVATIVE_LOWER = F(7, 1000)
MEAN_LOWER = F(1, 2000)
CLOCK_ERROR = F(1, 100000)
SQRT_TWO = (F(1414213562373095, 10**15), F(1414213562373096, 10**15))


def rational(value):
    if type(value) is int or type(value) is F:
        return F(value)
    if isinstance(value, s.Rational):
        return F(int(value.p), int(value.q))
    raise ValueError('an exact rational value is required')


def checked_source(data=None, probe=None):
    data = response.data() if data is None else data
    probe = observable.data() if probe is None else probe
    p, pi = data['P'], data['pi']
    if p.shape != (66, 66) or pi.shape != (66, 1):
        raise ValueError('the full native population carrier is required')
    if data['states'] != probe['states'] or s.diag(*pi) != probe['pi']:
        raise ArithmeticError('the transfer and detector source charts differ')
    if (any(x < 0 for x in p) or any(x <= 0 for x in pi) or sum(pi) != 1
            or s.ones(1, 66)*p != s.ones(1, 66) or p*pi != pi):
        raise ArithmeticError('the native stationary population law changed')
    if p*s.diag(*pi) != s.diag(*pi)*p.T:
        raise ArithmeticError('the population transfer is not reversible')
    if min(p[i, i] for i in range(66)) < s.Rational(791, 990):
        raise ArithmeticError('the positive transfer bound does not apply')
    if probe['pairs'] != observable.cap.pairs(0):
        raise ArithmeticError('the native detector pair dictionary changed')
    expected_y = s.zeros(66)
    for i, j, _, _ in probe['pairs']:
        if pi[i]/pi[j] != s.Rational(25, 8):
            raise ArithmeticError('the accessible pair has a different population ratio')
        expected_y[i, j], expected_y[j, i] = -s.I, s.I
    if (probe['Y'].conjugate() != -probe['Y'] or probe['Y'].H != probe['Y']
            or probe['Y']**2 != probe['Q'] or probe['Y'] != expected_y):
        raise ArithmeticError('the imaginary unit-norm native detector changed')
    return data, probe


def derivative_interval(data=None, probe=None, degree=DEGREE):
    """Bound the proposed unitary Y derivative using rational log columns.

    If U(t)=exp(-itH), the derivative is -2 sum_(i,j)(pi_i-pi_j)H_ij.
    The actual probe has pi_i/pi_j=25/8 on every spectator pair.
    """
    response._count(degree, 'logarithm degree')
    if degree == 0:
        raise ValueError('a positive logarithm degree is required')
    data, probe = checked_source(data, probe)
    pairs = probe['pairs']; pi = data['pi']
    columns = s.zeros(66, len(pairs))
    for k, (_, j, _, _) in enumerate(pairs):
        columns[j, k] = 1
    term = s.SparseMatrix(columns); value = s.zeros(66, len(pairs))
    deficit = s.SparseMatrix(s.eye(66)-data['P'])
    for n in range(1, degree+1):
        term = deficit*term
        value += term/s.Integer(n)
    coefficient = rational(-s.Rational(4, 5)*sum(
        (pi[i]-pi[j])*value[i, k] for k, (i, j, _, _) in enumerate(pairs)))
    a, b = SQRT_TWO
    if not (a > 0 and a*a < 2 < b*b):
        raise ArithmeticError('the square-root enclosure is invalid')
    low, high = sorted((coefficient*a, coefficient*b))
    error = 2*transfer.logarithm_tail(degree)
    return {'degree': degree, 'sqrt_two_coefficient': str(coefficient),
            'sqrt_two_interval': [str(a), str(b)],
            'operator_logarithm_tail': str(error/2),
            'derivative_tail_upper': str(error),
            'lower': str(low-error), 'upper': str(high+error),
            'all_spectator_pairs': len(pairs), 'carrier_dimension': 66}


def finite_time_certificate(data=None, probe=None):
    data, probe = checked_source(data, probe)
    derivative = derivative_interval(data, probe)
    log_upper = sum((transfer.DEFICIT**k/k for k in range(1, DEGREE+1)), F())
    log_upper += transfer.logarithm_tail(DEGREE)
    if log_upper >= NORM_UPPER or F(derivative['lower']) <= DERIVATIVE_LOWER:
        raise ArithmeticError('the finite comparison margins failed')
    remainder = F(4, 3)*NORM_UPPER**3*TIME**3
    mean_lower = DERIVATIVE_LOWER*TIME-remainder
    if mean_lower <= MEAN_LOWER:
        raise ArithmeticError('the proposed transfer response has no certified margin')
    return {'derivative': derivative, 'generator_norm_upper': str(NORM_UPPER),
            'scalar_logarithm_upper': str(log_upper),
            'response_is_odd': True, 'time': str(TIME),
            'derivative_lower': str(DERIVATIVE_LOWER),
            'third_order_remainder_upper': str(remainder),
            'mean_lower': str(mean_lower), 'simplified_mean_lower': str(MEAN_LOWER),
            'native_modular_stationary_mean': '0',
            'native_modular_parameter_may_be_any_real': True,
            'identity_identification_of_transfer_and_working_carrier_is_tested': True,
            'transfer_unitary_is_executed_by_native_source': False}


def detector_contract(receipt=None):
    """Reuse the owning executed Y certificate without normalizing its output.

    Its byte-bound proof and per-host acceptance tests own the compiled arrays.
    The comparison independently recomputes every error coefficient it uses.
    """
    receipt = strict_json.load(observable_packet.OUTPUT) if receipt is None else receipt
    if receipt.get('source_bindings') != observable_packet.source_bindings():
        raise ArithmeticError('the executed detector receipt has stale source bindings')
    word = observable.words()['Y']
    admitted = observable.interpreter.source_word(word, 9)
    arithmetic = observable.multiplication_error(len(word))
    _, failure = observable.clock._processor_failure_bound(len(word))
    error = observable.cap.ROTATION_GATE+F(arithmetic['operator_error_upper'])
    delta = (error*(2+error)+failure)/2
    gain = (1+error)**2+failure
    row = receipt['compiler']['words']['Y']
    expected = {'operator_error_upper': str(error),
                'returned_array_operator_error_upper': str(observable.cap.ROTATION_GATE),
                'product_arithmetic': arithmetic,
                'first_processor_failure_mass_upper': str(failure),
                'complete_halfdiamond_upper': str(delta),
                'channel_norm_upper': str(gain),
                'word_digest': observable.interpreter.digest(admitted),
                'resources': admitted['declared_resources']}
    if any(row.get(k) != v for k, v in expected.items()):
        raise ArithmeticError('the executed Y detector allowance differs from its owner')
    rounded_delta, rounded_gain = F(1, 40000000), F(20000001, 20000000)
    if delta >= rounded_delta or gain >= rounded_gain:
        raise ArithmeticError('the Y service exceeds the comparison error budget')
    return {'complete_halfdiamond_upper': str(rounded_delta),
            'channel_norm_upper': str(rounded_gain),
            'rounded_outward_from_exact_owning_allowances': True,
            'raw_native_slots': 25, 'address_export_aborts_retained': True,
            'proper_code_first_failures_retained_in_bound': True,
            'success_conditioning_or_polar_normalization': False,
            'compiler_owner': 'model/native_clock_observables/results.json'}


def operational_certificate(detector=None):
    detector = detector_contract() if detector is None else detector
    epsilon = F(observable.preparation_certificate()['half_diamond_upper'])
    delta = F(detector['complete_halfdiamond_upper'])
    gain = F(detector['channel_norm_upper'])
    if delta < 0 or gain < 1:
        raise ValueError('nonnegative error and gain at least one are required')
    error = gain*(epsilon+CLOCK_ERROR)+delta
    separation = MEAN_LOWER/2-2*error
    if separation <= F(1, 5000):
        raise ArithmeticError('the complete implementation error erases the separation')
    return {'preparation_halftrace_upper': str(epsilon),
            'hypothetical_clock_service_halfdiamond_upper': str(CLOCK_ERROR),
            'each_complete_output_error_upper': str(error),
            'complete_output_halftrace_lower': str(separation),
            'simplified_complete_output_halftrace_lower': '1/5000',
            'score_values': [-1, 0, 1], 'all_abort_scores': 0,
            'all_records_retained': True,
            'preparation_archive_traced_in_replacement_bound': True,
            'adaptation_to_preparation_archive_included': False,
            'clock_service_error_is_a_conditional_acceptance_budget': True,
            'unnormalized_numerical_outputs_are_not_claimed_to_be_probability_laws': True}


def blindness_certificate(data=None):
    data = response.data() if data is None else data
    columns = 0; zero = 0
    for row in data['branches']:
        inputs = [i for i, _ in row.pairs]; outputs = [j for _, j in row.pairs]
        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
            raise ArithmeticError('native complete branch is not a partial injection')
        columns += len(inputs); zero += not inputs
    return {'complete_raw_slots': len(data['branches']), 'zero_slots': zero,
            'partial_injection_columns': columns,
            'interleaved_arbitrary_diagonal_phases_preserve_record_reference_output': True,
            'adaptive_policy_uses_only_retained_records': True,
            'terminal_working_carrier_is_traced': True,
            'coherent_terminal_measurements_included_in_blindness': False,
            'reference_preserving_unitary_intertwining_with_transfer_clock': False,
            'different_reference_or_general_encoding_excluded': False}
