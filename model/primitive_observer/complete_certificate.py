"""Exact output bounds for the executed native/work observer channel.

The certified object is the factorized map represented by the original binary
native maps and work density arrays. All logical/reference input states are
covered, with independently prepared work. No numerical positivity, trace
normalization or unrounded source gate is assumed for these arrays.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from ..common_source import continuation_bounds, exact_output, open_observer, protocol, smooth_history
from ..common_source import output_certificate as report
from ..common_source.open_environment import Environment
from . import exact_native, native, runtime


STEPS = 512
WEIGHTS = (F(1, 4), F(1, 2), F(1, 4))
SOURCE_PATHS = tuple(sorted(set((*report.SOURCE_PATHS,
    'model/primitive_observer/runtime.py', 'model/primitive_observer/common.py',
    'model/primitive_observer/native.py', 'model/primitive_observer/exact_native.py',
    'model/primitive_observer/complete_certificate.py'))))


def nuclear_upper(array):
    """Directed sqrt(D)*Frobenius upper bound on the original square array."""
    # Check the concrete type before reading metadata or copying entries;
    # ndarray subclasses may replace either operation with arbitrary code.
    if type(array) is not np.ndarray:
        raise ValueError('an original base ndarray is required')
    if (len(array.shape) != 2 or array.shape[0] < 1
            or array.shape[0] != array.shape[1]):
        raise ValueError('a square original numerical array is required')
    array = exact_output._array(array, array.shape)
    squared = sum((a*a+b*b for a, b in map(exact_output._component, array.flat)), F(0))
    return exact_output._sqrt_interval(F(array.shape[0])*squared)[1]


def native_error(history):
    """All-input/reference bound for replacing only the numerical native maps.

    The actual mixed work blocks remain fixed in both outputs. A complete
    runtime history is required; no postselection or positive-matrix repair is
    performed. Integrity checks are not authentication of hostile proof objects.
    """
    history = native.validate_history(history)
    if history.native_mode not in ('protected_qubit', 'passive_qubit'):
        raise ValueError('the two-column encoded native contract is required')
    exact_native.completeness(history.law, history.native_mode)
    maps, norms = {}, {}
    total = maximum = abort = F(0)
    for branch in history.branches:
        ideal = exact_native.ideal_map(history.law, history.native_mode, branch.record)
        key = id(branch.native_map), tuple(sorted(ideal.items()))
        if key not in maps:
            maps[key] = exact_native.compare_map(branch.native_map, ideal)
        difference = maps[key]['distance_upper']
        maximum = max(maximum, maps[key]['frobenius_error_upper'])
        key = id(branch.work_density)
        if key not in norms:
            norms[key] = nuclear_upper(branch.work_density)
        term = difference*norms[key]
        total += term
        if branch.abort:
            abort += term
    return {'distance_upper': total, 'maximum_map_frobenius_error_upper': maximum,
            'abort_distance_upper': abort, 'complete_raw_outcomes': len(history.branches)}


def work_error(worker, program, finite, infinite):
    """Compare the actual shared-runtime matrices, not the older interpreter."""
    program = protocol.validate_program(program)
    smooth_history._validate(finite)
    smooth_history._validate(infinite)
    if finite.dimension != program.dimension or infinite.dimension is not None:
        raise ValueError('distinct matched finite and infinite target enclosures required')
    worker = runtime.validate_history(worker)
    if (worker.cursor != len(worker.program.events) or worker.reference_dimension != 1
            or worker.preparation_mode != 'declared'
            or strict_json.loads(worker.program.binding).get('canonical_program') != program.specification()):
        raise ValueError('complete declared work preparation and canonical program required')
    densities = {branch.record: branch.density for branch in worker.branches}
    comparison = exact_output.compare_densities(program, densities, finite)
    _, occupation = smooth_history.distance_interval(finite, infinite)
    mass_defect = abs(comparison['source_mass_real']-1)+abs(comparison['source_mass_imag'])
    report._gate(mass_defect, report.MASS_GATE, 'shared runtime complete mass')
    return {'finite_upper': comparison['distance_upper'],
            'infinite_upper': comparison['distance_upper']+occupation,
            'mass_defect_upper': mass_defect,
            'finite_to_infinite_upper': occupation}


def compose(native_bound, work_bounds):
    """The encoded first-current law is (1/4,1/2,1/4) for every input.

    These are trusted computed bounds, not a verifier for caller-authored proof
    claims. run() obtains all of them from actual arrays and fresh enclosures.
    """
    if (type(native_bound) is not F or native_bound < 0 or type(work_bounds) is not tuple
            or len(work_bounds) != 3 or any(type(x) is not F or x < 0 for x in work_bounds)):
        raise ValueError('one exact native bound and three nonnegative exact work bounds required')
    return native_bound+sum((weight*error for weight, error in zip(WEIGHTS, work_bounds, strict=True)), F(0))


def run():
    """Reexecute both native laws and memory modes, all currents and records."""
    enclosures, comparisons = {}, {}
    work_rows, rows, witnesses = [], [], []
    for dimension in (2, 4):
        for noisy in (False, True):
            env = open_observer.canonical_environment() if noisy else Environment(())
            for law in ('C', 'E'):
                pair = {}
                for mode in ('protected_qubit', 'passive_qubit'):
                    history = native.execute(env, law, dimension, STEPS, native_mode=mode)
                    history = native.validate_history(history)
                    errors = []
                    for current, worker in zip((-1, 0, 1), history.work_histories, strict=True):
                        program = protocol.continuation_program(dimension, STEPS, native_current=current)
                        for cutoff in (None, dimension):
                            key = current, cutoff
                            if key not in enclosures:
                                enclosures[key] = smooth_history.enclose(program, dimension=cutoff)
                        # Cache only validated, byte-bound histories. Different
                        # runtime outputs must obtain their own comparisons.
                        if worker.seal not in comparisons:
                            comparison = work_error(worker, program, enclosures[current, dimension],
                                                    enclosures[current, None])
                            comparisons[worker.seal] = comparison
                            work_rows.append({'dimension': dimension, 'noise': noisy, 'native_current': current,
                                'canonical_program': program.specification(),
                                'finite_NOISELESS_distance_upper': report.outward(comparison['finite_upper']),
                                'infinite_NOISELESS_distance_upper': report.outward(comparison['infinite_upper']),
                                'mass_defect_upper': str(report.MASS_GATE),
                                'all_work_records': len(worker.records)})
                        errors.append(comparisons[worker.seal])
                    attachment = native_error(history)
                    report._gate(attachment['distance_upper'], F(1, 10**8), 'complete native attachment')
                    finite = compose(attachment['distance_upper'], tuple(e['finite_upper'] for e in errors))
                    infinite = compose(attachment['distance_upper'], tuple(e['infinite_upper'] for e in errors))
                    report._gate(finite, F(1, 4) if noisy else F(1, 1000), 'complete finite observer')
                    report._gate(infinite, F(3, 4), 'complete infinite observer')
                    mass = 2*attachment['distance_upper']+compose(F(0), tuple(e['mass_defect_upper'] for e in errors))
                    report._gate(mass, F(1, 10**7), 'all-input complete observer mass')
                    pair[mode] = (finite, infinite)
                    rows.append({'dimension': dimension, 'noise': noisy, 'law': law, 'native_mode': mode,
                        'complete_raw_outcomes': attachment['complete_raw_outcomes'],
                        'native_map_replacement_upper': '1/100000000',
                        'all_input_abort_output_nuclear_norm_upper': '1/50000000',
                        'complete_finite_NOISELESS_distance_upper': report.outward(finite),
                        'complete_infinite_NOISELESS_distance_upper': report.outward(infinite),
                        'all_input_absolute_complex_mass_defect_upper': '1/10000000',
                        'source_resources': native.resources(history)})
                # Existing exact memory theorem: D(active +/-)=1,
                # D(passive +/-)=1-q, q=sum_current p_j P_j(r2=0).
                for cutoff in (dimension, None):
                    target = 'finite' if cutoff is not None else 'infinite'
                    index = int(cutoff is None)
                    intervals = [smooth_history.probability_interval(enclosures[c, cutoff], 2, 0)
                                 for c in (-1, 0, 1)]
                    qlo = sum((p*x[0] for p, x in zip(WEIGHTS, intervals, strict=True)), F(0))
                    qhi = sum((p*x[1] for p, x in zip(WEIGHTS, intervals, strict=True)), F(0))
                    active, passive = pair['protected_qubit'][index], pair['passive_qubit'][index]
                    advantage = qlo-2*(active+passive)
                    # Both noisy/noiseless finite comparisons resolve memory;
                    # the infinite comparison resolves d=4 without noise.
                    if (cutoff is not None or (dimension == 4 and not noisy)) and advantage <= 0:
                        raise ArithmeticError('complete computed memory advantage was not certified')
                    witnesses.append({'dimension': dimension, 'noise': noisy, 'law': law, 'target': target,
                        'ideal_passive_loss_interval': [report.outward(qlo, upper=False), report.outward(qhi)],
                        'computed_active_phase_distance_lower': report.outward(max(F(0), 1-2*active), upper=False),
                        'computed_passive_phase_distance_upper': report.outward(1-qlo+2*passive),
                        'computed_active_minus_passive_distance_lower': report.outward(advantage, upper=False),
                        'positive_advantage_certified': advantage > 0})
    root = Path(__file__).resolve().parents[2]
    return {'schema': 'oph.primitive_observer.complete_certificate.v1',
        'source_sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
        'steps': STEPS, 'outward_report_grid': str(F(1, report.GRID)),
        'noise_schedule': open_observer.canonical_environment().specification(),
        'unbounded_physical_service_obligations': {
            'work_stage_keys': list(continuation_bounds.required_service_keys(program.operations)),
            'environment_stage_keys': [f'environment_{i}:site_{event.site}@{event.after}'
                                       for i, event in enumerate(open_observer.canonical_environment().events)],
            'native_and_storage': ['logical_input_preparation', 'native_encoding',
                'first_native_repair_and_record', 'first_conditional_recovery_if_active',
                'logical_storage_during_work', 'final_record_selected_repair',
                'final_conditional_recovery_if_active', 'record_storage_and_feedforward'],
            'total_error_upper': None},
        'exact_first_current_weights': [str(x) for x in WEIGHTS],
        'work_comparisons': work_rows, 'complete_observers': rows, 'memory_witnesses': witnesses,
        'scope': {'arbitrary_logical_reference_input': True, 'work_prepared_independently': True,
            'all_raw_records_abort_and_padding_retained': True,
            'actual_shared_runtime_work_and_native_maps_certified': True,
            'factorized_binary_array_map_is_the_certified_output': True,
            'subsequent_floating_tensor_multiplication_certified': False,
            'matched_noisy_infinite_target': False, 'arbitrary_work_input': False,
            'physical_service_error_upper': None, 'physical_clock_calibrated': False,
            'full_78_state_native_input_certified': False,
            'terminal_destructive_read_certified': False},
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
            stream.write(smooth_history.ih.canonical(result))
    print(smooth_history.ih.canonical({'all_pass': result['all_pass'],
        'complete_observers': len(result['complete_observers']), 'memory_witnesses': result['memory_witnesses']}))


if __name__ == '__main__':
    main()
