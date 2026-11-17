"""Joint phase/work predictions with complete, unselected output bounds."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from ..clocked_common_observer import construction as work, experiment as budget
from ..common_source import output_certificate as report
from ..native_interface import construction as interface
from . import construction as c
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
NATIVE_GATE = F(1, 10**8)


def source_paths():
    from ..compiled_native_interface import construction as compiled
    return tuple(sorted(set((*budget.SOURCE_PATHS, *interface.SOURCE_PATHS,
        *compiled.SOURCE_PATHS, 'model/clocked_quantum_readout/__init__.py',
        'model/clocked_quantum_readout/__main__.py',
        'model/clocked_quantum_readout/construction.py',
        'model/clocked_quantum_readout/experiment.py'))))


def joint_bound(native, w, bounds):
    native, w = c.validate_native(native), work.validate_work(w)
    if (type(bounds) is not dict or bounds.get('work_binding') != budget.work_binding(w)
            or native.cap > w.horizon):
        raise ValueError('the computed work budget must bind the complete executed family')
    check = c.certify_native(native)
    report._gate(check['error'], NATIVE_GATE, 'actual heralded native density arithmetic')
    total = F(0)
    for (herald, n, current, bit), row in check['per_record'].items():
        error = bounds['preparation'] if herald == 0 else bounds['errors'][n, current]
        total += row['error']+(row['mass']+2*row['error'])*error
    return {'finite_error': budget.up(total), 'native_arithmetic': check['error'],
            'per_record': check['per_record']}


def observations(native, w):
    values = dict.fromkeys(('mass_real', 'mass_imag', 'phase_match_and_success',
        'final_gauge_one_and_success', 'phase_match_and_final_gauge_one',
        'phase_match_and_different_smooth_bits'), F(0))
    count = 0
    for (herald, n, current, bits, phase), nr, wr in c.output_rows(native, w):
        real, imag = budget.probability(nr, wr)
        values['mass_real'] += real
        values['mass_imag'] += imag
        count += 1
        if current:
            if phase == 0:
                values['phase_match_and_success'] += real
                if bits[3] == 1:
                    values['phase_match_and_final_gauge_one'] += real
                if bits[1] != bits[2]:
                    values['phase_match_and_different_smooth_bits'] += real
            if bits[3] == 1:
                values['final_gauge_one_and_success'] += real
    report._gate(abs(values['mass_real']-1)+abs(values['mass_imag']),
                 budget.MASS_GATE, 'complete phase/work binary mass')
    if count != 128*native.cap+3:
        raise ArithmeticError('complete output count changed')
    return values


def exact_native_statistics(law, cap):
    states = c.exact_outputs(law, cap)
    probabilities = {r: c.clock._trace(s) for r, s in states.items()}
    stopped = sum((p for (_, _, j, _), p in probabilities.items() if j), F(0))
    phase = sum((p for (_, _, j, bit), p in probabilities.items() if j and bit == 0), F(0))
    timeout = sum((p for (h, _, j, _), p in probabilities.items() if h and j == 0), F(0))
    if stopped+timeout != c.YIELD or probabilities[0, 0, 0, -1] != 1-c.YIELD:
        raise ArithmeticError('herald failure or timeout population changed')
    return {'probabilities': probabilities, 'success': stopped, 'timeout': timeout,
            'phase_match_and_success': phase}


def difference_interval(left, right, error):
    return [report.outward(left-right-2*error, upper=False),
            report.outward(left-right+2*error)]


def correlation_interval(values, error, success):
    """Unnormalized success-conditioned covariance; no successful outputs removed.

    J*s-P*G equals s^2 Cov(phase-match, gauge-one | successful stop).
    The returned numerator avoids dividing by a possibly rare event.
    """
    def limits(key):
        return max(F(0), values[key]-2*error), min(F(1), values[key]+2*error)
    jl, ju = limits('phase_match_and_final_gauge_one')
    pl, pu = limits('phase_match_and_success')
    gl, gu = limits('final_gauge_one_and_success')
    return [report.outward(jl*success-pu*gu, upper=False),
            report.outward(ju*success-pl*gl)]


def resources(w, cap):
    from ..compiled_native_interface import construction as compiled
    conditional = budget.resources(w, cap)
    native = c.clock.statistics('C', cap)
    prep = dict(w.prefixes[0].program.preparation_resources)['native_pulses']
    pulses = F(conditional['expected_work_source_pulses_upper'])
    expected_mixers = c.YIELD*(2-native['timeout_probability'])
    expected_work = prep+c.YIELD*(pulses-prep)
    mixer_pulses = compiled.resources()['native_pulses']
    return {'native_preparations': 1, 'work_preparations': 1,
        'complete_joint_records': 128*cap+3,
        'maximum_native_attempts_including_herald_and_pointer': cap+2,
        'maximum_balanced_mixer_calls': 2,
        'expected_native_attempts_exact': str(1+c.YIELD*(native['expected_native_attempts']
                                                       +1-native['timeout_probability'])),
        'expected_balanced_mixer_calls_exact': str(expected_mixers),
        'expected_work_source_pulses_upper': str(expected_work),
        'expected_work_and_mixer_source_pulses_upper': str(expected_work+mixer_pulses*expected_mixers),
        'maximum_work_and_mixer_source_pulses': conditional['maximum_work_source_resources']['native_pulses']+2*mixer_pulses,
        'maximum_work_source_resources': conditional['maximum_work_source_resources'],
        'failed_herald_cost': {'native_attempts': 1, 'balanced_mixer_calls': 0,
                              'work_preparations': 1, 'work_waits': 0},
        'timeout_after_n_slots': {'native_attempts': 'n+1', 'balanced_mixer_calls': 1,
                                  'work_preparations': 1, 'work_waits': 'n', 'phase_reads': 0},
        'successful_stop_after_n_slots': {'native_attempts': 'n+2', 'balanced_mixer_calls': 2,
                                          'work_preparations': 1, 'work_waits': 'n', 'phase_reads': 1},
        'physical_service_or_storage_cost_bound': False}


def compiled_service():
    """A separate full-carrier implementation comparison, with leakage retained.

    Intervening exact native instruments extend by projecting onto the native
    code and its complement, with identity on the latter placed in the zero
    current / bit-zero outcome. This is CPTP and fixes native-sector behavior.
    The compiled mixer itself is never compressed or projected in this bound.
    """
    from ..compiled_native_interface import construction as compiled
    actual = compiled.full_matrix()
    local = compiled.local_matrix()
    if (type(actual) is not np.ndarray or actual.shape != (512, 512)
            or not np.array_equal(actual, np.kron(np.eye(64, dtype=int), local))):
        raise ArithmeticError('the actual complete source spectator lift changed')
    certificate = compiled.certify()
    counts = compiled.resources()
    error = certificate['two_call_half_diamond_upper']
    if type(error) is not F or error < 0:
        raise ArithmeticError('an exact complete two-call source bound is required')
    return error, {'maximum_calls': 2,
        'one_call_half_diamond_upper': compiled._outward(certificate['half_diamond_upper']),
        'two_call_half_diamond_upper': compiled._outward(error),
        'source_resources_per_mixer': {key: counts[key] for key in (
            'native_pulses', 'single_code_pulses', 'pair_code_pulses',
            'proper_code_service_events', 'identity_phase_pulses', 'maximum_work_qubits')},
        'full_source_carrier_dimension': 512,
        'source_matrix_compressed_to_native_sector': False,
        'source_padding_leakage_charged_by_full_carrier_bound': True,
        'other_operations_are_exact_CPTP_maps_in_this_implementation_comparison': True,
        'padding_extension': 'dephase native code versus complement; identity on complement in zero/bit-zero outcome',
        'floating_full_512D_history_executed_or_certified': False}


def run():
    w = work.execute_work(max(c.CAPS))
    bounds = budget.work_certificates(w)
    implementation_error, implementation = compiled_service()
    rows, comparisons = [], []
    for cap in c.CAPS:
        retained = {}
        for law in ('C', 'E'):
            native = c.execute_native(law, cap)
            proof = joint_bound(native, w, bounds)
            observed = observations(native, w)
            exact = exact_native_statistics(law, cap)
            finite = proof['finite_error']
            eventual = finite+exact['timeout']
            # Implementation target is the separate full-carrier binary source
            # mixer instrument. Its additional error is not postselected away.
            # Padding can remain unresolved forever under the declared native
            # extension. Do not infer an unlimited binary-source terminal state
            # from the ideal invariant-sector tail theorem.
            source_comparison = finite+implementation_error
            retained[law] = observed, proof, exact
            rows.append({'candidate': law, 'cap': cap, 'complete_joint_records': 128*cap+3,
                'native_density_arithmetic_upper': str(NATIVE_GATE),
                'finite_ideal_service_output_error_upper': report.outward(finite),
                'eventual_ideal_service_output_error_upper': report.outward(eventual),
                'finite_source_mixer_comparison_upper': report.outward(source_comparison),
                'exact_preparation_failure_probability': str(1-c.YIELD),
                'exact_successful_stop_probability': str(exact['success']),
                'exact_unmeasured_timeout_probability': str(exact['timeout']),
                'exact_native_phase_match_and_success': str(exact['phase_match_and_success']),
                'complete_native_record_probabilities': [
                    {'record': list(r), 'probability': str(p)} for r, p in exact['probabilities'].items()],
                'finite_source_mixer_record_weight_intervals': {
                    key: budget.interval(observed[key], source_comparison) for key in observed
                    if key not in ('mass_real', 'mass_imag')},
                'eventual_ideal_service_observation_intervals': {
                    key: budget.interval(observed[key], eventual) for key in observed
                    if key not in ('mass_real', 'mass_imag')},
                'finite_phase_work_covariance_numerator_interval': correlation_interval(
                    observed, finite, exact['success']),
                'resources': resources(w, cap)})
        left, right = retained['C'], retained['E']
        error = left[1]['finite_error']+right[1]['finite_error']
        lp, rp = left[2]['probabilities'], right[2]['probabilities']
        tv = sum((abs(lp[r]-rp[r]) for r in lp), F(0))/2
        # Conditional work laws are the same and normalized for both candidates.
        # Appending their records therefore leaves this complete classical TV unchanged.
        for h in (1, 2):
            for n in range(1, cap+1):
                for j in (-1, 1):
                    if sum(lp[h,n,j,b] for b in (0,1)) != sum(rp[h,n,j,b] for b in (0,1)):
                        raise ArithmeticError('the C/E stopping marginal changed')
        if lp[0,0,0,-1] != rp[0,0,0,-1] or any(lp[h,cap,0,-1] != rp[h,cap,0,-1] for h in (1,2)):
            raise ArithmeticError('failure or timeout selects a different classical law')
        report._gate(abs(left[0]['final_gauge_one_and_success']-right[0]['final_gauge_one_and_success']),
                     NATIVE_GATE, 'C/E work-only marginal')
        phase_difference = left[2]['phase_match_and_success']-right[2]['phase_match_and_success']
        if phase_difference != tv or phase_difference <= 0:
            raise ArithmeticError('the retained complete-record phase discriminator changed')
        comparisons.append({'cap': cap,
            'exact_ideal_complete_joint_classical_total_variation': str(tv),
            'exact_unconditional_phase_match_probability_difference': str(phase_difference),
            'source_binary_complete_record_half_L1_lower': report.outward(max(F(0),tv-2*implementation_error),upper=False),
            'work_only_marginals_equal_in_exact_model': True,
            'joint_phase_and_final_gauge_C_minus_E_interval': difference_interval(
                left[0]['phase_match_and_final_gauge_one'],right[0]['phase_match_and_final_gauge_one'],
                error+2*implementation_error),
            'joint_phase_and_different_smooth_C_minus_E_interval': difference_interval(
                left[0]['phase_match_and_different_smooth_bits'],right[0]['phase_match_and_different_smooth_bits'],
                error+2*implementation_error)})
    return {'schema': 'oph.clocked_quantum_readout.v1', 'all_pass': True,
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths()}),
        'shared_source_mixer': implementation, 'histories': rows, 'candidate_discrimination': comparisons,
        'scope': {'same_stopped_native_and_interacting_work_history': True,
            'both_herald_signs_failure_and_unmeasured_timeouts_retained': True,
            'one_work_preparation_and_no_reset_at_current': True,
            'native_ideal_service_surrogate_executed_in_66_dimensions': True,
            'source_mixer_error_attached_as_separate_full_carrier_comparison': True,
            'unlimited_binary_source_terminal_state_claimed': False,
            'native_phase_preparation_or_target_POVM_independently_supplied': False,
            'stationary_native_input_uniform_scheduler_and_balanced_target_design_supplied': True,
            'coherent_source_services_work_action_and_model_duration_supplied': True,
            'native_clock_records_aligned_with_observer_modular_flow': False,
            'physical_source_service_error_bound': False,
            'infinite_occupation_or_interacting_spacetime_derived': False,
            'new_empirical_validation': False}}
