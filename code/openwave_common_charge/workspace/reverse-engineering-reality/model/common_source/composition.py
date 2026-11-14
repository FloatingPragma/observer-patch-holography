"""Replay a continued observer, its native attachment and complete error ledger.

One source preparation supports four committed quantum reads. The two scalar
reads and subsequent clock probe act on the disturbed retained state. Complete
branch maps, independent Hamiltonian comparisons and reset controls are rebuilt
on every replay; this receipt stores summaries, not a large trajectory dump.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from tools import strict_json
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..observer.source_scalar import field
from . import continuation_bounds as bounds
from . import experiment as old
from . import native_join, protocol, scalar_coupling as scalar, thermal_clock


TOLERANCE = 2e-9
STEPS = 512


@lru_cache(maxsize=64)
def _executed(program, backend='source', stop_after=None):
    return protocol.execute(program, backend=backend, stop_after=stop_after)


def _states(history):
    return {branch.record: branch.state for branch in history.branches}


def complete_distance(left, right):
    """Trace distance of the full CQ outputs, including a retained reference."""
    if not left or set(left) != set(right):
        raise ValueError('both complete outcome menus are required')
    values = []
    for record in sorted(left):
        x, y = np.asarray(left[record]), np.asarray(right[record])
        if (x.shape != y.shape or x.size == 0 or x.ndim not in (1, 2)
                or not np.isfinite(x).all() or not np.isfinite(y).all()):
            raise ValueError('matching finite work/reference branches required')
        x, y = x.ravel(), y.ravel()
        difference = np.outer(x, x.conj())-np.outer(y, y.conj())
        values.append(float(np.sum(abs(np.linalg.eigvalsh(difference))))/2)
    return math.fsum(values)


def _direct_continuous(program, *, reset_after_first_scalar=False):
    """Independent complete experiment from the supplied action and effects.

The negative control explicitly discards the disturbed work state just after
the first scalar read, preserving its old record probabilities. The actual
execution vocabulary contains no such reset instruction.
"""
    d, size = program.dimension, program.work_dimension
    h = scalar.ideal_hamiltonian(d, float(program.coupling))
    q = np.kron(np.kron(np.eye(d), field.quadratures(d)[0]), np.eye(4))
    x = np.zeros((size, size))
    for index in range(size):
        x[index ^ 1, index] = 1
    initial = np.kron(field.vacuum_vector(d)[0], [1., 0., 0., 0.])
    initial = expm(-1j*program.native_current*math.atan2(4, 3)*x)@initial
    gauge = tuple(np.diag([int(index % 2 == bit) for index in range(size)])
                  for bit in (0, 1))
    minus, plus = expm(-.5j*q), expm(.5j*q)
    smooth = ((minus-1j*plus)/2, (minus+1j*plus)/2)
    current = {(): initial}
    scalar_reads = 0
    for operation in program.operations:
        if operation.kind == 'evolve':
            u = expm(-1j*float(operation.time)*h)
            current = {record: u@state for record, state in current.items()}
        elif operation.kind == 'feedback':
            current = {record: expm(-1j*(math.atan2(4, 3) if operation.enabled
                        and record[operation.record] else 0.)*x)@state
                       for record, state in current.items()}
        elif operation.kind in ('gauge_read', 'smooth_read'):
            maps = gauge if operation.kind == 'gauge_read' else smooth
            current = {(*record, bit): k@state for record, state in current.items()
                       for bit, k in enumerate(maps)}
            if operation.kind == 'smooth_read':
                scalar_reads += 1
                if reset_after_first_scalar and scalar_reads == 1:
                    current = {record: math.sqrt(float(np.vdot(state, state).real))*initial
                               for record, state in current.items()}
        else:
            raise ValueError('the continued observer uses gauge and smooth reads only')
    return {record: state[:, None] for record, state in current.items()}


def _marginal(states, slot):
    return [math.fsum(float(np.vdot(state, state).real)
                     for record, state in states.items() if record[slot] == bit)
            for bit in (0, 1)]


def _prefix_error(before, after):
    length = len(before.branches[0].record)
    effects = protocol.prefix_effects(after, length)
    return max(float(np.linalg.norm(effects[b.record]-b.kraus.conj().T@b.kraus, 2))
               for b in before.branches)


def _clock(program):
    """Read the common Hamiltonian after the same two sequential scalar reads."""
    if (tuple(op.kind for op in program.operations if op.kind in protocol.READS)
            != ('gauge_read', 'smooth_read', 'smooth_read', 'gauge_read')
            or tuple(op.kind for op in program.operations[-2:]) != ('evolve', 'gauge_read')):
        raise ValueError('the smooth continuation with terminal evolve/gauge read is required')
    # This is the prefix just before the last evolution and gauge measurement.
    prefix = _executed(program, stop_after=len(program.operations)-2)
    ref = thermal_clock.reference(program.dimension, beta=1., coupling=float(program.coupling))
    effect = np.diag([index % 2 for index in range(program.work_dimension)])
    rows = []
    for time in (F(0), F(1, 2), F(1)):
        setting = replace(program, operations=(*program.operations[:-2],
                          protocol.Operation('evolve', time), program.operations[-1]))
        result = _executed(setting)
        modular_effect = thermal_clock.modular_flow(ref, effect, -float(time))
        independent_effect = expm(1j*float(time)*ref.hamiltonian)@effect@expm(-1j*float(time)*ref.hamiltonian)
        predicted = math.fsum(float(np.vdot(b.state, modular_effect@b.state).real)
                              for b in prefix.branches)
        observed = _marginal(_states(result), 3)[1]
        time_bound = F(scalar.time_error_bound(program.dimension, program.steps,
                       time=time, coupling=program.coupling)
                       ['ideal_channel_half_diamond_error_upper'])
        defect = float(np.linalg.norm(modular_effect-independent_effect, 2))
        if defect > TOLERANCE or abs(observed-predicted) > float(time_bound)+TOLERANCE:
            raise AssertionError('clock continuation is detached from the same Hamiltonian')
        rows.append({'model_time': str(time), 'modular_parameter_beta_one': str(-time),
                     'source_clock_bit_one_probability': observed,
                     'continuous_prediction_on_same_retained_prefix': predicted,
                     'gibbs_identity_operator_error': defect,
                     'ideal_appended_time_error_upper': str(time_bound)})
    variation = max(row['source_clock_bit_one_probability'] for row in rows)-min(
        row['source_clock_bit_one_probability'] for row in rows)
    if variation <= 1e-5:
        raise AssertionError('the declared continued clock readout is constant')
    return {'settings': rows, 'readout_probability_range': variation,
            'shared_prefix_includes_both_disturbing_scalar_reads': True,
            'settings_are_separate_runs_of_that_same_prefix': True,
            'physical_duration_or_temperature_calibrated': False,
            'gibbs_reference_prepared_as_actual_state': False}


def _one(program):
    result = _executed(program)
    product = _executed(program, 'product')
    continuous = _direct_continuous(program)
    reset = _direct_continuous(program, reset_after_first_scalar=True)
    prefix = _executed(program, stop_after=5)
    checkpoint = protocol.checkpoint_from_json(program, protocol.checkpoint_to_json(prefix))
    resumed = protocol.execute(program, checkpoint=checkpoint)
    comparison = protocol.compare_histories(result, product)
    resumed_comparison = protocol.compare_histories(result, resumed)
    retained_error = _prefix_error(prefix, result)
    direct_distance = complete_distance(_states(result), continuous)
    reset_distance = complete_distance(continuous, reset)
    budget = bounds.program_bound(program)
    temporal = F(budget['ideal_temporal_instrument_error_upper'])
    if (comparison['maximum_kraus_operator_error'] > TOLERANCE
            or comparison['complete_flagged_vector_error'] > TOLERANCE
            or resumed_comparison['maximum_kraus_operator_error'] > TOLERANCE
            or resumed_comparison['complete_flagged_vector_error'] > TOLERANCE
            or retained_error > TOLERANCE
            or direct_distance > float(temporal)+TOLERANCE
            or reset_distance <= .01):
        raise AssertionError('continued-source composition or its retained-state control failed')
    h = scalar.ideal_hamiltonian(program.dimension, float(program.coupling))
    return {
        'program': program.specification(),
        'branches': [{'record': list(b.record), 'probability': b.probability,
                      'unnormalized_work_energy': float(np.vdot(b.state, h@b.state).real)}
                     for b in result.branches],
        'source_versus_independent_product': comparison,
        'source_versus_independent_continuous_complete_trace_distance': direct_distance,
        'reset_control_complete_trace_distance': reset_distance,
        'reset_control_last_scalar_probability_difference': _marginal(continuous, 2)[1]-_marginal(reset, 2)[1],
        'old_record_effect_operator_error': retained_error,
        'serialized_checkpoint_resumption': resumed_comparison,
        'resources_per_realized_path': dict(result.branches[0].resources),
        'maximum_resources_per_realized_path': {key: max(dict(b.resources)[key] for b in result.branches)
                                               for key, _ in result.branches[0].resources},
        'expanded_tree_resources_not_physical_path_cost': dict(result.expanded_tree_resources),
        'number_of_retained_quantum_outcomes': len(result.branches),
        'working_qubits': scalar.qubits(program.dimension),
        'ideal_finite_temporal_error_upper': str(temporal),
        'all_pass': True,
    }


def native_experiment(programs):
    """Close the classical native -> quantum -> native loop on every outcome.

The smooth scalar record has the opposite sign convention to the old sharp
record: bit zero is the plus effect. Only that bit is decoded for actuation.
The native load is carried unchanged during the quantum program. This is a
classical native input ensemble, not a compiled arbitrary M78 quantum input.
"""
    if (type(programs) is not dict or set(programs) != {-1, 0, 1}
            or any(type(j) is not int or type(p) is not protocol.Program for j, p in programs.items())):
        raise ValueError('all three native-current programs are required')
    signature = {(p.dimension, p.steps, p.coupling, p.operations) for p in programs.values()}
    if len(signature) != 1 or any(p.native_current != current for current, p in programs.items()):
        raise ValueError('native currents must select the preparation of one common program')
    for program in programs.values():
        if (tuple(op.kind for op in program.operations if op.kind in protocol.READS)
                != ('gauge_read', 'smooth_read', 'smooth_read', 'gauge_read')
                or tuple(op.kind for op in program.operations[-2:]) != ('evolve', 'gauge_read')):
            raise ValueError('the native attachment needs the declared smooth/clock bit meanings')
    histories = {current: protocol.validate_checkpoint(program, _executed(program), 'source')
                 for current, program in programs.items()}
    if any(len(h.records[0]) != 4 for h in histories.values()):
        raise ValueError('the canonical four-read continuation is required')
    totals = {'probability_mass': 0., 'mean_final_port0_load': 0.,
              'mean_final_port0_load_disconnected_actuator': 0.,
              'last_smooth_plus_probability': 0., 'clock_bit_one_probability': 0.}
    paths = 0
    bits = set()
    for loads in native.states(binary=True):
        for first_coin in (0, 1):
            retained_loads, first_record = native.native(loads, native.E, first_coin)
            history = histories[first_record[1]]
            for branch in history.branches:
                # Smooth bit 0 is plus; native_continuation expects plus as 1.
                actuator = 1-branch.record[2]
                for second_coin in (0, 1):
                    final, final_record = native_join.native_continuation(retained_loads, actuator, second_coin)
                    disconnected, _ = native_join.native_continuation(retained_loads, 0, second_coin)
                    weight = branch.probability/264
                    totals['probability_mass'] += weight
                    totals['mean_final_port0_load'] += weight*final[0]
                    totals['mean_final_port0_load_disconnected_actuator'] += weight*disconnected[0]
                    totals['last_smooth_plus_probability'] += weight*actuator
                    totals['clock_bit_one_probability'] += weight*branch.record[3]
                    bits.add((first_record, branch.record, final_record))
                    paths += 1
    if abs(totals['probability_mass']-1) > TOLERANCE or paths != 4224:
        raise AssertionError('the common observer omitted a native/quantum outcome')
    return {'complete_joint_outcome_count': paths,
            'declared_distinct_public_record_words': len(bits),
            'fixed_width_public_record_storage_bits': 10,
            'native_events_per_realized_path': 2,
            'quantum_work_preparations_per_realized_path': 1,
            'classical_native_input_states': 66,
            'actuator_input': 'one minus the second smooth record bit; no other record or work state',
            'quantum_terminal_states_are_retained': True,
            'arbitrary_quantum_native_input_claimed': False,
            'scalar_and_clock_bits_are_joint_observations': True,
            'summary': totals}


def _compact_budget(budget):
    # The operation sequence determines each intermediate recurrence; the
    # receipt retains its complete final bounds once per cutoff, not once
    # per current or a duplicate of the readily recomputed stage ledger.
    return {**budget, 'occupation_certificate': {
        key: value for key, value in budget['occupation_certificate'].items() if key != 'stages'}}


def run(steps=STEPS):
    if type(steps) is not int or steps < 1:
        raise ValueError('a positive integer split count required')
    levels = []
    for d in (2, 4):
        programs = {j: protocol.continuation_program(d, steps, native_current=j) for j in (-1, 0, 1)}
        required = (*bounds.required_service_keys(programs[0].operations),
                    'native_preparation_and_first_record', 'native_actuation_and_last_record')
        histories = [_one(programs[j]) for j in (-1, 0, 1)]
        budget = _compact_budget(bounds.program_bound(programs[0]))
        levels.append({'dimension': d, 'histories': histories,
                       'ideal_error_budget_uniform_over_native_currents': budget,
                       'native_observer': native_experiment(programs),
                       'physical_service_coverage': bounds.service_error_budget(
                           F(budget['complete_ideal_product_to_infinite_instrument_error_upper']),
                           required_stage_keys=required),
                       'clock_after_retained_smooth_reads': _clock(programs[0])})
    canonical_program = protocol.continuation_program(2, steps)
    limit = [_compact_budget(bounds.complete_error_budget(d, canonical_program.operations, steps=d*d))
             for d in (256, 1024, 4096)]
    return {'schema': 'oph.common_source.composed_observer.v1',
            'levels': levels,
            'analytic_refinement_not_large_cutoff_source_execution': limit,
            'scope': {'common_source_preparation_count': 1, 'quantum_reads': 4,
                      'retained_sequential_scalar_reads': 2,
                      'exact_tower_or_full_Lean_OperationalObserver_instantiated': False,
                      'sharp_detector_moment_bound_transferred': False,
                      'model_time_is_physical_duration': False,
                      'independent_empirical_evidence': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    print(json.dumps({'all_pass': result['all_pass'],
                      'retained_joint_histories_per_cutoff': [level['native_observer']['complete_joint_outcome_count']
                                                           for level in result['levels']],
                      'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()
