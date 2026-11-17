"""Certify the actual complete observer outputs against exact history targets.

Every replay executes source words and recomputes the polynomial proof objects.
Only outward rounded guarantees are portable; raw floating arrays, their hashes
and their exact masses are not claimed identical across numerical platforms.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction as F
import hashlib
from pathlib import Path

from tools import strict_json
from . import continuation_bounds, exact_output, infinite_history as ih
from . import open_observer, open_protocol, protocol, smooth_history as sh


STEPS = 512
GRID = 10**6
MASS_GATE = F(1, 10**8)
SOURCE_PATHS = tuple(sorted(set((*open_observer.SOURCE_PATHS,
    'model/common_source/exact_output.py',
    'model/common_source/output_certificate.py',
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json'))))


def outward(value, *, upper=True):
    """Export a directed rational decimal, never a rounded floating bound."""
    if type(value) is not F:
        raise ValueError('exact rational certificate value required')
    scaled = value*GRID
    integer = -(-scaled.numerator//scaled.denominator) if upper else scaled.numerator//scaled.denominator
    return str(F(integer, GRID))


def _gate(value, threshold, label):
    if type(value) is not F or not 0 <= value <= threshold:
        raise ArithmeticError(f'{label} failed its exact acceptance bound')


def compare_receipt(expected, actual):
    """Exact portable guarantees, including original JSON scalar types."""
    if type(expected) is not type(actual):
        raise ArithmeticError('receipt value type differs')
    if isinstance(actual, dict):
        if expected.keys() != actual.keys():
            raise ArithmeticError('receipt fields differ')
        for key in actual:
            compare_receipt(expected[key], actual[key])
    elif isinstance(actual, list):
        if len(expected) != len(actual):
            raise ArithmeticError('receipt list length differs')
        for left, right in zip(expected, actual, strict=True):
            compare_receipt(left, right)
    elif expected != actual:
        raise ArithmeticError('recomputed outward guarantees or source bindings differ')


def _summary(result, gate):
    _gate(result['distance_upper'], gate, 'complete output')
    _gate(abs(result['source_mass_real']-1)+abs(result['source_mass_imag']), MASS_GATE,
          'complete output mass')
    return {
        'half_trace_norm_distance_interval': [outward(result['distance_lower'], upper=False),
                                               outward(result['distance_upper'])],
        'exact_arithmetic_mass_defect_upper': str(MASS_GATE),
        'acceptance_half_trace_norm_upper': str(gate),
        'retained_records': len(result['records']),
    }


def feedback_control(program):
    """Locked intervention: set both feedback angles to zero, retain all reads."""
    program = protocol.validate_program(program)
    return replace(program, operations=tuple(
        replace(op, enabled=False) if op.kind == 'feedback' else op
        for op in program.operations))


def _event_probability(vectors, slot=3, bit=1):
    # exact_output validates the arrays before this helper is used. Convert
    # each original binary component before any norm or probability reduction.
    return sum((F(float(x.real))**2+F(float(x.imag))**2
                for record, vector in vectors.items() if record[slot] == bit
                for x in vector), F(0))


def _event_comparison(vectors, enclosure):
    observed = _event_probability(vectors)
    lower, upper = sh.probability_interval(enclosure, 3, 1)
    return {'event': 'last gauge record r3=1, all earlier records retained',
            'computed_event_mass_interval': [outward(observed, upper=False), outward(observed)],
            'ideal_event_probability_interval': [outward(lower, upper=False), outward(upper)],
            'computed_event_error_upper': outward(max(abs(observed-lower), abs(observed-upper)))}


def run():
    rows, controls = [], []
    environment = open_observer.canonical_environment()
    for current in (-1, 0, 1):
        base = protocol.continuation_program(2, STEPS, native_current=current)
        infinite = sh.enclose(base)
        for dimension in (2, 4):
            program = replace(base, dimension=dimension)
            finite = sh.enclose(program, dimension=dimension)
            source = protocol.execute(program)
            vectors = {b.record: b.state[:, 0] for b in source.branches}
            to_finite = exact_output.compare_vectors(program, vectors, finite)
            to_infinite = exact_output.compare_vectors(program, vectors, infinite)
            noisy = open_protocol.execute(program, environment)
            densities = {b.record: b.density for b in noisy.branches}
            mixed_finite = exact_output.compare_densities(program, densities, finite)
            _, occupation = sh.distance_interval(finite, infinite)
            noisy_infinite = mixed_finite['distance_upper']+occupation
            _gate(noisy_infinite, F(3, 4), 'noisy output to noiseless infinite target')
            rows.append({
                'dimension': dimension, 'native_current': current,
                'canonical_program': program.specification(),
                'computed_pure_to_exact_finite_continuous': _summary(to_finite, F(1, 1000)),
                'computed_pure_to_exact_infinite': _summary(to_infinite, F(3, 10)),
                'computed_mixed_to_exact_NOISELESS_finite_continuous': _summary(mixed_finite, F(1, 4)),
                'computed_mixed_to_exact_NOISELESS_infinite_upper': outward(noisy_infinite),
                'clock_event_to_infinite': _event_comparison(vectors, infinite),
                'resources': {
                    'pure_path_component_maxima': {key: max(dict(b.resources)[key] for b in source.branches)
                                                   for key in protocol.RESOURCE_KEYS},
                    'pure_expanded_tree': dict(source.expanded_tree_resources),
                    'mixed_path_component_maxima': {key: max(dict(b.resources)[key] for b in noisy.branches)
                                                    for key in open_protocol.RESOURCE_KEYS},
                    'mixed_expanded_tree': dict(noisy.expanded_tree_resources)},
            })
            if dimension == 4 and current == 0:
                control = feedback_control(program)
                control_source = protocol.execute(control)
                control_vectors = {b.record: b.state[:, 0] for b in control_source.branches}
                control_infinite = sh.enclose(control)
                control_certificate = exact_output.compare_vectors(control, control_vectors, control_infinite)
                a, b = sh.probability_interval(infinite, 3, 1)
                c, d = sh.probability_interval(control_infinite, 3, 1)
                low, high = a-d, b-c
                if low <= 0 <= high:
                    raise ArithmeticError('locked feedback intervention has no certified clock contrast')
                controls.append({
                    'dimension': dimension, 'native_current': current,
                    'intervention': 'set both feedback angles to zero; keep preparation, reads and elapsed model time',
                    'computed_control_to_exact_infinite': _summary(control_certificate, F(3, 10)),
                    'control_clock_event_to_infinite': _event_comparison(control_vectors, control_infinite),
                    'exact_infinite_feedback_ON_minus_OFF_clock_probability_interval':
                        [outward(low, upper=False), outward(high)],
                    'nonzero_contrast_certified': True,
                    'control_path_component_maxima': {key: max(dict(branch.resources)[key]
                                                              for branch in control_source.branches)
                                                     for key in protocol.RESOURCE_KEYS},
                    'control_expanded_tree': dict(control_source.expanded_tree_resources),
                })
    root = Path(__file__).resolve().parents[2]
    return {'schema': 'oph.common_source.exact_output_certificate.v1',
            'source_sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                              for name in SOURCE_PATHS},
            'steps': STEPS, 'outward_report_grid': str(F(1, GRID)),
            'environment': environment.specification(),
            'remaining_work_service_stage_keys': list(continuation_bounds.required_service_keys(base.operations)),
            'remaining_environment_service_stage_keys': [
                f'environment_{index}:site_{event.site}@{event.after}'
                for index, event in enumerate(environment.events)],
            'executions': rows, 'feedback_interventions': controls,
            'scope': {'all_complete_work_outputs_compared': True,
                      'numerical_preparation_evolution_read_and_feedback_errors_included': True,
                      'input_arrays_or_proof_bounds_loaded_from_receipt': False,
                      'raw_arrays_renormalized': False,
                      'arbitrary_input_channel_certificate': False,
                      'mixed_comparison_target_contains_noise': False,
                      'floating_native_attachment_certified': False,
                      'physical_service_error_upper': None,
                      'physical_clock_calibrated': False,
                      'temporal_error_added_again_to_direct_comparison': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    paths = parser.add_mutually_exclusive_group()
    paths.add_argument('--write', type=Path)
    paths.add_argument('--check', type=Path)
    args = parser.parse_args()
    # Parse first; never let a malformed retained claim trigger costly work.
    expected = strict_json.load(args.check) if args.check else None
    result = run()
    if args.check:
        compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(ih.canonical(result))
    print(ih.canonical({'all_pass': result['all_pass'],
                        'feedback_interventions': result['feedback_interventions']}))


if __name__ == '__main__':
    main()
