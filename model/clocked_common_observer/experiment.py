"""Reproducible native-clock/work histories and complete directed budgets."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from ..common_source import exact_output, output_certificate as report, protocol
from ..common_source.infinite_history import canonical
from ..primitive_observer import common, runtime
from . import clock, construction as c, work_bounds as wb
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


SOURCE_PATHS = tuple(sorted(set((*clock.SOURCE_PATHS, *wb.SOURCE_PATHS,
    *report.SOURCE_PATHS, 'model/clocked_common_observer/__init__.py',
    'model/clocked_common_observer/__main__.py',
    'model/clocked_common_observer/construction.py',
    'model/clocked_common_observer/experiment.py', 'tools/strict_json.py'))))
ERROR_GATE = F(1, 10000)
ARITHMETIC_GATE = F(1, 10**8)
MASS_GATE = F(1, 10**8)


def up(value):
    """Widen intermediate bounds to a fixed rational grid, avoiding huge powers."""
    if type(value) is not F or value < 0:
        raise ValueError('a nonnegative exact intermediate bound is required')
    grid = 10**18
    scaled = value*grid
    return F(-(-scaled.numerator//scaled.denominator), grid)


def trace(array):
    """Exact trace of the original entries of a finite base square array."""
    if (type(array) is not np.ndarray or array.ndim != 2
            or array.shape[0] != array.shape[1]):
        raise ValueError('an original base square numeric ndarray is required')
    array = clock._array(array, array.shape)
    parts = [exact_output._component(x) for x in np.diag(array)]
    return sum((r for r, _ in parts), F(0)), sum((i for _, i in parts), F(0))


def probability(native, work):
    a, b = trace(native)
    x, y = trace(work)
    return a*x-b*y, a*y+b*x


def work_binding(work):
    """Bind a trusted computed budget to every executed source/density history."""
    payload = (work.horizon, work.steps, str(work.delta),
               tuple(h.seal for h in work.prefixes),
               tuple((n, j, h.seal) for n, j, h in work.successes))
    return hashlib.sha256(runtime._json(payload).encode()).hexdigest()


def work_certificates(work):
    """Charge source maps, repeated products, preparation and actual densities.

    The only pure reference preparation is the existing finite d2 vacuum.
    The step/suffix operator bounds themselves hold for arbitrary work and
    untouched reference input, not merely for that executed preparation.
    """
    work = c.validate_work(work)
    if (work.prefixes[0].reference_dimension != 1
            or work.prefixes[0].preparation_mode != 'declared'):
        raise ValueError('this output budget executes the declared independent work preparation')
    step = wb.certify_step(work.delta, work.steps)
    suffixes = {j: wb.certify_suffix(j, work.steps) for j in (-1, 1)}
    original = work.prefixes[0].initial_density
    prep = up(wb.certify_preparation(original))
    report._gate(prep, ARITHMETIC_GATE, 'finite work preparation')
    source = step.maps[0][1]
    actual_step = work.prefixes[0].program.events[0].instruments[0].kraus[0][0]
    if not np.array_equal(actual_step, source):
        raise ArithmeticError('waiting execution detached from the certified source step')
    step_error = up(step.stinespring_error_upper)
    prefixes, prefix_errors = [np.eye(16, dtype=complex)], [F(0)]
    errors, arithmetic = {}, {}
    reset = {j: runtime.execute(suffixes[j].compiled_program) for j in (-1, 1)}
    successes = {(n, j): h for n, j, h in work.successes}
    successes.update({(0, j): h for j, h in reset.items()})
    for n in range(work.horizon+1):
        if n:
            matrix = source@prefixes[-1]
            rounded = wb.exact_product_residual(matrix, source, prefixes[-1])
            error = up((1+step_error)*prefix_errors[-1]+step_error+rounded)
            prefixes.append(matrix)
            prefix_errors.append(error)
        prefix, error = prefixes[n], prefix_errors[n]
        if n:
            rounding = wb.exact_density_residual(work.prefixes[n].branches[0].density,
                                                (prefix,), original)
            arithmetic[n, 0] = up(rounding)
            errors[n, 0] = up(error*(2+error)/2+(1+error)**2*prep+rounding)
        for j in (-1, 1):
            suffix = suffixes[j]
            maps = dict(suffix.maps)
            rounding = sum((wb.exact_density_residual(branch.density,
                            (prefix, maps[branch.record]), original)
                            for branch in successes[n, j].branches), F(0))
            delta = up((1+suffix.stinespring_error_upper)*error+suffix.stinespring_error_upper)
            arithmetic[n, j] = up(rounding)
            errors[n, j] = up(delta*(2+delta)/2+(1+delta)**2*prep+rounding)
    report._gate(max(arithmetic.values()), ARITHMETIC_GATE, 'actual work density arithmetic')
    return {'work_binding': work_binding(work),
            'step': step, 'suffixes': suffixes, 'errors': errors, 'arithmetic': arithmetic,
            'preparation': prep, 'reset': reset,
            'prefix_operator_errors': tuple(prefix_errors)}


def joint_certificate(work, bounds, law, native, *, cap, activity=F(1), reset=False):
    """Compose computed bounds; caller-authored bounds are not proof objects."""
    work = c.validate_work(work)
    if type(bounds) is not dict or bounds.get('work_binding') != work_binding(work):
        raise ValueError('the computed work budget must bind this complete executed family')
    if type(reset) is not bool:
        raise ValueError('an explicit reset-control Boolean is required')
    clock._horizon(cap)
    if cap > work.horizon:
        raise ValueError('the requested cap has not been executed')
    native_check = clock.compare_outputs(law, cap, native, activity=activity)
    total = F(0)
    for (n, j), row in native_check['per_record'].items():
        work_error = bounds['errors'][0 if reset and j else n, j]
        # ||R_actual||_1 summed over suffix outcomes is <=1+2*work_error.
        total += row['distance_upper']*(1+2*work_error)+row['ideal_mass']*work_error
    total = up(total)
    report._gate(native_check['distance_upper'], ARITHMETIC_GATE, 'actual native density arithmetic')
    report._gate(total, ERROR_GATE, 'complete native-clock/work finite output')
    stats = clock.statistics(law, cap, activity=activity)
    return {'finite_error': total, 'eventual_error': total+stats['timeout_probability'],
            'native_arithmetic': native_check['distance_upper'], 'clock': stats}


def exact_observations(work, native, *, cap, reset=None):
    """Original binary outputs enter exact trace sums before any rounding."""
    work = c.validate_work(work)
    clock._horizon(cap)
    if (cap > work.horizon or type(native) is not dict
            or any(type(record) is not tuple or len(record) != 2
                   or any(type(v) is not int for v in record) for record in native)
            or set(native) != set(clock.records(cap))):
        raise ValueError('a complete executed native/work cap is required')
    native = {key: clock._array(value, (66, 66)) for key, value in native.items()}
    if reset is not None:
        if (type(reset) is not dict or set(reset) != {-1, 1}
                or any(type(j) is not int for j in reset)):
            raise ValueError('both signed complete fresh-preparation controls are required')
        reset = {j: runtime.validate_history(h) for j, h in reset.items()}
        for j, history in reset.items():
            expected = common.compile_program(protocol.continuation_program(c.DIMENSION, work.steps,
                                                                            native_current=j))
            if (history.program.digest != expected.digest or history.cursor != len(expected.events)
                    or history.reference_dimension != 1 or history.preparation_mode != 'declared'
                    or not np.array_equal(history.initial_density, expected.initial_density)):
                raise ValueError('the reset comparison requires its declared complete signed source suffix')
    values = dict.fromkeys(('mass_real', 'mass_imag', 'final_gauge_one_and_success',
                           'first_slot_gauge_one', 'different_smooth_bits_and_success'), F(0))
    records = []
    histories = work.successful()
    for n, j in clock.records(cap):
        if j:
            history = histories[n, j] if reset is None else reset[j]
            branches = history.branches
        else:
            branches = work.prefixes[cap].branches
        for branch in branches:
            real, imag = probability(native[n, j], branch.density)
            values['mass_real'] += real
            values['mass_imag'] += imag
            records.append((n, j, branch.record))
            if j:
                if branch.record[3] == 1:
                    values['final_gauge_one_and_success'] += real
                    if n == 1:
                        values['first_slot_gauge_one'] += real
                if branch.record[1] != branch.record[2]:
                    values['different_smooth_bits_and_success'] += real
    if len(records) != 32*cap+1 or len(set(records)) != len(records):
        raise ArithmeticError('the complete native-clock/work output menu changed')
    report._gate(abs(values['mass_real']-1)+abs(values['mass_imag']), MASS_GATE,
                 'exact binary complete output mass')
    return values


def interval(value, error):
    # 2*error remains valid for unnormalized/non-Hermitian binary arrays.
    return [report.outward(max(F(0), value-2*error), upper=False),
            report.outward(min(F(1), value+2*error))]


def resources(work, cap, *, activity=F(1)):
    """Actual complete-path maxima, with finite record storage and service scope."""
    work = c.validate_work(work)
    clock._horizon(cap)
    if cap > work.horizon:
        raise ValueError('the requested resource cap has not been executed')
    candidates = [h for n, _, h in work.successes if n <= cap] + [work.prefixes[cap]]
    names = set(key for h in candidates for b in h.branches for key, _ in b.resources)
    maximum = {key: max(dict(b.resources).get(key, 0) for h in candidates for b in h.branches)
               for key in sorted(names)}
    if maximum['work_preparations'] != 1 or maximum['fresh_pointer_blanks'] != 4:
        raise ArithmeticError('repeated preparation or incomplete pointer resource accounting')
    stats = clock.statistics('C', cap, activity=activity)
    step_cost = dict(work.prefixes[0].program.events[0].instruments[0].resources)
    preparation = dict(work.prefixes[0].program.preparation_resources)
    suffix_ceiling = max(dict(b.resources)['native_pulses']-preparation['native_pulses']
                         -n*step_cost['native_pulses']
                         for n, _, h in work.successes if n <= cap for b in h.branches)
    expected_pulses = (preparation['native_pulses']
                       +stats['expected_controller_slots']*step_cost['native_pulses']
                       +(1-stats['timeout_probability'])*suffix_ceiling)
    return {'maximum_work_source_resources': maximum, 'maximum_native_channel_attempts': cap,
            'maximum_controller_slots': cap,
            'expected_native_attempts_exact': str(stats['expected_native_attempts']),
            'expected_controller_slots_exact': str(stats['expected_controller_slots']),
            'expected_work_source_pulses_upper': str(expected_pulses),
            'encoded_clock_record_bits': (2*cap).bit_length(),
            'maximum_work_record_bits': 4,
            'joint_complete_record_alphabet_size': 32*cap+1,
            'native_preparations': 1, 'work_preparations': 1,
            'native_channel_source_pulse_compilation_claimed': False,
            'physical_service_error_upper': None}


def run():
    work = c.validate_work(c.execute_work(max(c.CAPS)))
    bounds = work_certificates(work)
    rows, retained = [], {}
    for law in ('C', 'E'):
        prefixes = clock.execute_prefixes(law, work.horizon)
        for cap in c.CAPS:
            native = prefixes[cap]
            proof = joint_certificate(work, bounds, law, native, cap=cap)
            observed = exact_observations(work, native, cap=cap)
            retained[law, cap] = (proof, observed, native)
            rows.append({'law': law, 'cap': cap, 'complete_outcomes': 32*cap+1,
                'finite_complete_output_error_upper': report.outward(proof['finite_error']),
                'eventual_complete_output_error_upper': report.outward(proof['eventual_error']),
                'exact_timeout_probability': str(proof['clock']['timeout_probability']),
                'uniform_reference_stable_timeout_upper': str(proof['clock']['uniform_tail_bound']),
                'actual_native_density_arithmetic_upper': str(ARITHMETIC_GATE),
                'actual_work_density_arithmetic_upper': str(ARITHMETIC_GATE),
                'complete_mass_defect_upper': str(MASS_GATE),
                'eventual_observation_intervals': {
                    key: interval(observed[key], proof['eventual_error']) for key in
                    ('final_gauge_one_and_success', 'first_slot_gauge_one',
                     'different_smooth_bits_and_success')},
                'resources': resources(work, cap)})
    cap = work.horizon
    baseline, observed, native = retained['C', cap]
    reset_proof = joint_certificate(work, bounds, 'C', native, cap=cap, reset=True)
    reset_observed = exact_observations(work, native, cap=cap, reset=bounds['reset'])
    reset_difference = observed['first_slot_gauge_one']-reset_observed['first_slot_gauge_one']
    reset_margin = abs(reset_difference)-2*(baseline['finite_error']+reset_proof['finite_error'])
    if reset_margin <= 0:
        raise ArithmeticError('the retained-state versus reset witness is unresolved')
    idle_native = clock.execute_native('C', cap, activity=F(1, 2))
    idle = joint_certificate(work, bounds, 'C', idle_native, cap=cap, activity=F(1, 2))
    idle_observed = exact_observations(work, idle_native, cap=cap)
    idle_difference = idle_observed['final_gauge_one_and_success']-observed['final_gauge_one_and_success']
    idle_margin = abs(idle_difference)-2*(baseline['eventual_error']+idle['eventual_error'])
    if idle_margin <= 0:
        raise ArithmeticError('native-only idling versus common work witness is unresolved')
    for cap in c.CAPS:
        for key in ('final_gauge_one_and_success', 'first_slot_gauge_one', 'different_smooth_bits_and_success'):
            # The declared work-only law cannot distinguish C/E; its native
            # coherence remains available in the full factorized output.
            report._gate(abs(retained['C', cap][1][key]-retained['E', cap][1][key]),
                         ARITHMETIC_GATE, 'C/E common work marginal')
    root = Path(__file__).resolve().parents[2]
    return {'schema': 'oph.clocked_common_observer.v1',
        'source_sha256': bound_pins({p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'occupation_dimension': c.DIMENSION, 'source_split_count': c.STEPS,
        'work_model_duration_per_slot': str(c.DELTA), 'caps': list(c.CAPS),
        'finite_work_channel_bounds': {
            'waiting_step_half_diamond_upper': report.outward(bounds['step'].half_diamond_upper),
            'suffix_half_diamond_upper': {str(j): report.outward(b.half_diamond_upper)
                                        for j, b in bounds['suffixes'].items()},
            'arbitrary_work_reference_input': True,
            'actual_density_execution_roundoff_separately_included': True},
        'histories': rows,
        'retained_state_control': {'observable': 'first current at slot one and final gauge bit one',
            'fresh_preparation_comparison_is_not_a_physical_reset_implementation': True,
            'additional_work_preparations_on_a_successful_reset_path': 1,
            'signed_retained_minus_reset_interval':
                [report.outward(reset_difference-2*(baseline['finite_error']+reset_proof['finite_error']), upper=False),
                 report.outward(reset_difference+2*(baseline['finite_error']+reset_proof['finite_error']))],
            'absolute_difference_lower': report.outward(reset_margin, upper=False)},
        'native_only_idling_control': {'activity': '1/2',
            'resources': resources(work, cap, activity=F(1, 2)),
            'observable': 'eventual final gauge bit one, ignoring stop count and sign',
            'finite_error_upper': report.outward(idle['finite_error']),
            'eventual_error_upper': report.outward(idle['eventual_error']),
            'timeout_probability_exact': str(idle['clock']['timeout_probability']),
            'signed_idle_minus_baseline_interval':
                [report.outward(idle_difference-2*(baseline['eventual_error']+idle['eventual_error']), upper=False),
                 report.outward(idle_difference+2*(baseline['eventual_error']+idle['eventual_error']))],
            'absolute_difference_lower': report.outward(idle_margin, upper=False)},
        'scope': {'one_native_clock_and_retained_interacting_work_history': True,
            'every_signed_stopping_outcome_and_timeout_retained': True,
            'work_feedback_uses_committed_records': True,
            'native_channel_semantics_execution': True, 'work_source_words_executed': True,
            'native_source_pulse_words_executed': False,
            'bound_covers_complete_factorized_native_and_work_binary_densities': True,
            'additional_floating_tensor_product_not_certified': True,
            'prepared_native_work_independently': True,
            'native_preparation_is_declared_phase_mixture': True,
            'occupation_refinement_or_infinite_field_comparison': False,
            'physical_clock_calibration': False, 'physical_service_error_upper': None},
        'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    paths = parser.add_mutually_exclusive_group()
    paths.add_argument('--write', type=Path)
    paths.add_argument('--check', type=Path)
    args = parser.parse_args()
    expected = strict_json.load(args.check) if args.check else None
    result = run()
    if args.check:
        report.compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    print(canonical({key: result[key] for key in ('all_pass', 'finite_work_channel_bounds',
                                                'retained_state_control', 'native_only_idling_control')}))
