"""Separate source-executed Ramsey variant of the common interacting history.

The old sharp-sign protocol is unchanged. This detector retains both outcomes
and both conditional quantum states. Target matrices below are independent
checks of source execution, never replacements for emitted operations.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c, field
from ..observer.source_scalar.experiment import encoded
from . import experiment as e
from . import scalar_coupling as scalar


TOLERANCE = 2e-10


def pointer_word(d):
    """H; exp(-i Z_pointer Q0/2); S^dag; H, using the existing source compiler."""
    qubits, _ = e._register(d)
    bits = d.bit_length()-1
    values, basis = field.canonical_eigenbasis(field.quadratures(d)[0])
    basis_word = c.compile_real_orthogonal(basis)
    hadamard = [c.pulse('H', math.pi/2, [bits]), c.pulse('I', -math.pi/2, [bits])]
    local = (hadamard + c.inverse(basis_word)
             + c.compile_diagonal(np.concatenate((-.5*values, .5*values)))
             + basis_word
             + [c.pulse('Z', math.pi/4, [bits]), c.pulse('I', -math.pi/4, [bits])]
             + hadamard)
    return c.remap(local, [*range(2, 2+bits), qubits])


def coordinate(d):
    e._register(d)
    return np.kron(np.kron(np.eye(d), field.quadratures(d)[0]), np.eye(4))


def ideal_kraus(d):
    """Independent analytic maps, including each physical circuit phase."""
    q = coordinate(d)
    minus, plus = expm(-.5j*q), expm(.5j*q)
    return ((minus-1j*plus)/2, (minus+1j*plus)/2)


@lru_cache(maxsize=2)
def _pointer(d):
    qubits, size = e._register(d)
    word = pointer_word(d)
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, qubits+1))
    branches = tuple(unitary[b*size:(b+1)*size, :size].copy() for b in (0, 1))
    expected = ideal_kraus(d)
    discrepancy = max(float(np.linalg.norm(actual-target, 2))
                      for actual, target in zip(branches, expected, strict=True))
    complete = float(np.linalg.norm(sum(k.conj().T@k for k in branches)-np.eye(size), 2))
    if max(discrepancy, complete) > TOLERANCE:
        raise ArithmeticError('emitted Ramsey instrument differs from its complete analytic map')
    for branch in branches:
        branch.flags.writeable = False
    return branches, {'guard': guard, 'resources': e._word_resources(word),
                      'emitted_word': word, 'blank_pointer_bit': qubits,
                      'maximum_full_kraus_operator_error': discrepancy,
                      'full_kraus_completeness_error': complete,
                      'circuit': 'H; exp(-i Z_pointer Q0/2); Sdag; H; binary read',
                      'effects': ['(I+sin(Q0))/2', '(I-sin(Q0))/2'],
                      'kraus': ['exp(-i*pi/4) cos(Q0/2-pi/4)',
                                'exp(i*pi/4) cos(Q0/2+pi/4)'],
                      'pointer_auxiliary_parameter_interval': [0., 1.],
                      'pointer_parameter_is_physical_elapsed_time': False}


def pointer(d):
    branches, metadata = _pointer(d)
    return tuple(branch.copy() for branch in branches), deepcopy(metadata)


def _energy(states, h):
    return math.fsum(float(np.vdot(state, h@state).real) for _, state in states)


def _stage(name, states, h):
    return {'stage': name, 'energy': _energy(states, h),
            'probability_mass': math.fsum(float(np.vdot(v, v).real) for _, v in states),
            'branches': [{'record': list(record),
                          'probability': float(np.vdot(v, v).real),
                          'unnormalized_energy': float(np.vdot(v, h@v).real)}
                         for record, v in states]}


def _evolution_drift(before, after, independent, h):
    ideal_after = [(record, independent@state) for record, state in before]
    source_delta = _energy(after, h)-_energy(before, h)
    ideal_delta = _energy(ideal_after, h)-_energy(before, h)
    return {'source_executed_energy_change': source_delta,
            'independent_ideal_split_energy_change': ideal_delta,
            'source_minus_ideal_split_energy_defect': source_delta-ideal_delta,
            'exact_continuous_same_H_energy_change': 0.,
            'physical_work_or_energy_unit_calibration_claimed': False}


def history(d, steps=512, coupling=.25, feedback=True, native_current=0):
    """Complete four-branch source history with the smooth final detector."""
    qubits, size = e._register(d)
    if type(steps) is not int or steps < 1 or type(feedback) is not bool:
        raise ValueError('positive integer split count and explicit feedback switch required')
    if type(native_current) is not int or native_current not in (-1, 0, 1):
        raise ValueError('native current -1, 0 or 1 required')
    coupling = scalar._real(coupling, 'coupling')
    from . import protocol
    _, prepare_meta = scalar.preparation(d)
    kick_word = [c.pulse('X', -native_current*e.FEEDBACK_ANGLE, [0])]
    _, kick_guard = c.stabilize_executed_unitary(c.execute(kick_word, qubits))
    kick_meta = {'guard': kick_guard, 'resources': e._word_resources(kick_word),
                 'native_current': native_current}
    _, first_meta = scalar.source_unitary(d, steps, time=1., coupling=coupling)
    _, second_meta = scalar.source_unitary(d, steps, time=.5, coupling=coupling)
    _, a_meta = e.pointer(d, 0)
    read_b, b_meta = pointer(d)
    ideal_b = ideal_kraus(d)
    independent_first = scalar.independent_split_unitary(d, steps, time=1., coupling=coupling)
    independent_second = scalar.independent_split_unitary(d, steps, time=.5, coupling=coupling)
    # One canonical execution supplies every disturbed prefix and terminal
    # state. The retained energy ledger and independent oracle use those
    # states; they do not run a second production history.
    program = protocol.smooth_program(d, steps, native_current=native_current,
                                      feedback=feedback, coupling=F(coupling))
    prefixes = protocol.execution_prefixes(program)
    prepared = prefixes[0].initial_state[:, 0]
    initial = prefixes[0].branches[0].state[:, 0]
    ensemble = lambda cursor: [(b.record, b.state[:, 0]) for b in prefixes[cursor].branches]
    terminal = {branch.record: branch for branch in prefixes[-1].branches}
    target = np.kron(field.vacuum_vector(d)[0], np.array([1., 0., 0., 0.]))
    x0 = np.zeros((size, size))
    for index in range(size):
        x0[index ^ 1, index] = 1
    target = expm(-1j*native_current*e.FEEDBACK_ANGLE*x0)@target
    h = scalar.ideal_hamiltonian(d, coupling)
    first_ensemble, read_ensemble = ensemble(1), ensemble(2)
    feedback_ensemble, second_ensemble, terminal_ensemble = ensemble(3), ensemble(4), ensemble(5)
    first_margin = [float(np.vdot(v, v).real) for _, v in read_ensemble]
    q = coordinate(d)
    sinq = (expm(1j*q)-expm(-1j*q))/(2j)
    complete = np.zeros((size, size), complex)
    record_errors, independent_errors, law_errors = [], [], []
    branches, conditional, control_meta = [], [], []
    for a in (0, 1):
        _, metadata = e.feedback(d, a, feedback)
        control_meta.append(metadata)
        oracle_control = e._independent_feedback(size, a, feedback)
        before_b = second_ensemble[a][1]
        prefix_map = prefixes[2].branches[a].kraus
        effects = np.zeros((size, size), complex)
        marginal = float(np.vdot(before_b, before_b).real)
        signal = float(np.vdot(before_b, sinq@before_b).real)
        actual_probabilities = []
        for b in (0, 1):
            k = terminal[a, b].kraus
            state = terminal[a, b].state[:, 0]
            oracle = (ideal_b[b]@independent_second@oracle_control
                      @e.projector(size, 0, a)@independent_first@target)
            independent_errors.append(float(np.linalg.norm(state-oracle)))
            probability = float(np.vdot(state, state).real)
            expected_probability = (marginal+(1 if b == 0 else -1)*signal)/2
            law_errors.append(abs(probability-expected_probability))
            actual_probabilities.append(probability)
            effects += k.conj().T@k
            complete += k.conj().T@k
            branches.append({'record': [a, b], 'probability': probability,
                             'unnormalized_terminal_state': encoded(state),
                             'unnormalized_energy': float(np.vdot(state, h@state).real)})
        record_errors.append(float(np.linalg.norm(effects-prefix_map.conj().T@prefix_map, 2)))
        conditional.append({'first_record': a, 'first_record_probability': marginal,
                            'unnormalized_sine_expectation': signal,
                            'joint_probabilities': actual_probabilities,
                            'conditional_probabilities': [p/marginal for p in actual_probabilities]
                                if marginal > 0 else None,
                            'law': 'p(a,b)=[p(a)+(-1)^b <sin Q0>_a]/2'})
    stages = [_stage('prepared_scalar_and_gauge_blank', [((), prepared)], h),
              _stage('after_native_current_kick', [((), initial)], h),
              _stage('after_first_evolution', first_ensemble, h),
              _stage('after_gauge_readout_nonselective', read_ensemble, h),
              _stage('after_record_conditioned_feedback', feedback_ensemble, h),
              _stage('after_second_evolution', second_ensemble, h),
              _stage('after_smooth_readout_nonselective', terminal_ensemble, h)]
    changes = [{'from': a['stage'], 'to': b['stage'], 'energy_change': b['energy']-a['energy']}
               for a, b in zip(stages, stages[1:])]
    dual_actual = sum(k.conj().T@h@k for k in read_b)-h
    dual_ideal = sum(k.conj().T@h@k for k in ideal_b)-h
    injection = stages[-1]['energy']-stages[-2]['energy']
    predicted_injection = _energy(second_ensemble, dual_ideal)
    boundary_correction = _energy(second_ensemble, dual_ideal-np.eye(size)/8)
    actual_dual_energy = _energy(second_ensemble, dual_actual)
    mass_before = stages[-2]['probability_mass']
    energy = {'units': 'declared finite Hamiltonian model units; no physical energy calibration',
              'hamiltonian': 'existing scalar.ideal_hamiltonian with projected Q-square boundary terms',
              'stages': stages, 'stage_changes': changes,
              'energy_telescope_residual': stages[-1]['energy']-stages[0]['energy']
                  -math.fsum(row['energy_change'] for row in changes),
              'evolution_drift': {
                  'first': _evolution_drift([((), initial)], first_ensemble, independent_first, h),
                  'second': _evolution_drift(feedback_ensemble, second_ensemble, independent_second, h)},
              'smooth_readout': {
                  'source_executed_unconditional_energy_change': injection,
                  'source_dual_operator_expectation': actual_dual_energy,
                  'source_dual_energy_identity_residual': injection-actual_dual_energy,
                  'independent_analytic_finite_energy_change': predicted_injection,
                  'source_minus_analytic_finite_energy_residual': injection-predicted_injection,
                  'full_dual_energy_operator_error': float(np.linalg.norm(dual_actual-dual_ideal, 2)),
                  'infinite_normalized_unconditional_injection': 1/8,
                  'ideal_infinite_injection_times_actual_mass': mass_before/8,
                  'analytic_finite_boundary_correction': boundary_correction,
                  'finite_boundary_identity_residual': predicted_injection-mass_before/8-boundary_correction,
                  'source_injection_minus_infinite_one_eighth': injection-1/8,
                  'finite_or_conditional_injection_assumed_to_equal_one_eighth': False}}
    marginals_a = [sum(row['probability'] for row in branches if row['record'][0] == a) for a in (0, 1)]
    marginals_b = [sum(row['probability'] for row in branches if row['record'][1] == b) for b in (0, 1)]
    norm_error = abs(sum(row['probability'] for row in branches)-1)
    completeness_error = float(np.linalg.norm(complete-np.eye(size), 2))
    checks = {'normalized_complete_history': norm_error < TOLERANCE,
              'full_input_channel_complete': completeness_error < TOLERANCE,
              'first_record_effect_preserved_on_all_inputs': max(record_errors) < TOLERANCE,
              'independent_full_terminal_vectors': max(independent_errors) < TOLERANCE,
              'conditional_sine_law': max(law_errors) < TOLERANCE,
              'old_record_marginal_preserved': max(abs(x-y) for x, y in zip(first_margin, marginals_a)) < TOLERANCE,
              'finite_source_energy_identity': abs(energy['smooth_readout']['source_dual_energy_identity_residual']) < TOLERANCE,
              'independent_analytic_finite_energy': abs(injection-predicted_injection) < TOLERANCE}
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'schema': 'oph.common_source.smooth_history.v1',
            'dimension': d, 'steps_per_segment': steps, 'coupling': coupling,
            'native_current': native_current, 'feedback_enabled': feedback,
            'model_time_segments': [1., .5],
            'meaning': {'first_bit': 'species-0 right occupation after first evolution',
                        'second_bit': 'Ramsey scalar detector: bit zero has effect (I+sin Q0)/2'},
            'branches': branches, 'first_record_marginal': marginals_a, 'second_record_marginal': marginals_b,
            'conditional_readout_law': conditional,
            'normalization_error': norm_error, 'full_kraus_completeness_error': completeness_error,
            'maximum_old_record_effect_error': max(record_errors),
            'maximum_independent_terminal_vector_error': max(independent_errors),
            'maximum_conditional_law_error': max(law_errors),
            'resources_per_realized_history': e._resources(
                [prepare_meta, kick_meta, first_meta, a_meta, control_meta[0], second_meta, b_meta], d),
            'execution': {'preparation': prepare_meta, 'native_current_preparation': kick_meta,
                          'first_evolution': first_meta, 'first_pointer': a_meta,
                          'feedback': control_meta, 'second_evolution': second_meta, 'smooth_pointer': b_meta},
            'energy_ledger': energy,
            'scope': {'old_sharp_sign_protocol_changed': False,
                      'physical_clock_or_energy_calibration_claimed': False,
                      'physical_service_error_certified': False,
                      'ideal_compression_bound_includes_floating_execution': False},
            'checks': checks, 'all_pass': True}


def run(steps=512):
    levels = []
    for d in (2, 4):
        on = history(d, steps)
        off = history(d, steps, feedback=False)
        zero_on = history(d, steps, coupling=0.)
        zero_off = history(d, steps, coupling=0., feedback=False)
        response = on['second_record_marginal'][1]-off['second_record_marginal'][1]
        disconnected = zero_on['second_record_marginal'][1]-zero_off['second_record_marginal'][1]
        if abs(disconnected) > TOLERANCE:
            raise AssertionError('disconnected scalar marginal responds to gauge feedback')
        levels.append({'dimension': d, 'coupled_with_feedback': on, 'coupled_without_feedback': off,
                       'disconnected_with_feedback': zero_on, 'disconnected_without_feedback': zero_off,
                       'smooth_scalar_feedback_response': response,
                       'disconnected_scalar_feedback_response': disconnected})
    return {'schema': 'oph.common_source.smooth_protocol.v1', 'levels': levels,
            'old_sharp_sign_protocol_preserved': True, 'all_pass': True}
