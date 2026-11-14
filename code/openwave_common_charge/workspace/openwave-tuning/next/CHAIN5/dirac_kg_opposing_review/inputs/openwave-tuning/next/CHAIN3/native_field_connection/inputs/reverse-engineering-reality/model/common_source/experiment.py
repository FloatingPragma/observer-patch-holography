"""Joint observations and feedback on one interacting quantum register.

Gauge matter occupies bits 0,1 and the two neutral scalar sites occupy the
higher bits. Both recorded outcomes come from native source pointer words.
No component supplies a probability to a separate simulator. The Hamiltonian,
preparation and model-time schedule are specified construction inputs.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from ..observer.source_scalar.experiment import encoded, sources
from ..observer.source_scalar.field import quadratures, canonical_eigenbasis
from . import gauge_matter as gauge
from . import scalar_coupling as scalar
from . import records

DIMENSIONS = (2, 4)
STEPS = (128, 512)
COUPLING = .25
FEEDBACK_ANGLE = math.atan2(4, 3)
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')


def _word_resources(word):
    # Exact-zero classification of computed angles varies with LAPACK.
    # Every request, including zero-angle requests, stays in the charged word.
    return {key: value for key, value in c.resources(word).items()
            if key != 'zero_angle_requests'}


def _register(d):
    if type(d) is not int or d not in DIMENSIONS:
        raise ValueError('declared scalar dimension 2 or 4 required')
    return 2+2*(d.bit_length()-1), 4*d*d


def projector(size, bit, outcome):
    return np.diag([int((index >> bit & 1) == outcome) for index in range(size)])


@lru_cache(maxsize=8)
def pointer(d, bit):
    """Full source CNOT, blank pointer and both unnormalized readout maps."""
    qubits, size = _register(d)
    if type(bit) is not int or not 0 <= bit < qubits:
        raise ValueError('a system qubit required')
    word = c.cnot(bit, qubits)
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, qubits+1))
    branches = tuple(unitary[outcome*size:(outcome+1)*size, :size]
                     for outcome in (0, 1))
    discrepancy = max(float(np.linalg.norm(k-projector(size, bit, outcome), 2))
                      for outcome, k in enumerate(branches))
    if discrepancy > 2e-10:
        raise ArithmeticError('source pointer differs from its physical occupation effect')
    return branches, {'guard': guard, 'resources': _word_resources(word),
                      'projector_operator_error': discrepancy}


def scalar_projectors(d):
    values, basis = np.linalg.eigh(quadratures(d)[0])
    return tuple(np.kron(np.kron(np.eye(d),
                     (basis*((values > 0) == outcome))@basis.T), np.eye(4))
                 for outcome in (0, 1))


@lru_cache(maxsize=2)
def scalar_pointer(d):
    """Native QND measurement of the sign of the first scalar coordinate."""
    qubits, size = _register(d)
    bits = d.bit_length()-1
    values, basis = canonical_eigenbasis(quadratures(d)[0])
    if not (np.all(values[:d//2] < 0) and np.all(values[d//2:] > 0)):
        raise ArithmeticError('declared even scalar code has no zero coordinate eigenvalue')
    basis_word = c.remap(c.compile_real_orthogonal(basis), list(range(2, 2+bits)))
    word = c.inverse(basis_word)+c.cnot(2+bits-1, qubits)+basis_word
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, qubits+1))
    branches = tuple(unitary[outcome*size:(outcome+1)*size, :size]
                     for outcome in (0, 1))
    discrepancy = max(float(np.linalg.norm(k-target, 2))
                      for k, target in zip(branches, scalar_projectors(d)))
    if discrepancy > 2e-10:
        raise ArithmeticError('source scalar pointer differs from coordinate-sign projector')
    return branches, {'guard': guard, 'resources': _word_resources(word),
                      'projector_operator_error': discrepancy}


@lru_cache(maxsize=8)
def feedback(d, outcome, enabled):
    """The committed bit alone selects a dressed first-species hopping pulse."""
    qubits, _ = _register(d)
    if type(outcome) is not int or outcome not in (0, 1) or type(enabled) is not bool:
        raise ValueError('a committed bit and explicit feedback switch required')
    angle = -FEEDBACK_ANGLE if enabled and outcome == 1 else 0.
    word = [c.pulse('X', angle, [0])]
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, qubits))
    return unitary, {'guard': guard, 'resources': _word_resources(word),
                     'angle': angle, 'controller_input': outcome}


def _independent_feedback(size, outcome, enabled):
    x = np.zeros((size, size))
    for index in range(size):
        x[index ^ 1, index] = 1
    return expm(-1j*(FEEDBACK_ANGLE if enabled and outcome else 0.)*x)


def _resources(metadata, d):
    rows = [item['resources'] for item in metadata]
    return {**{key: sum(row[key] for row in rows) for key in COUNT_KEYS},
            'system_qubits': _register(d)[0],
            'pointer_qubits_per_read': 1, 'committed_classical_bits': 2,
            'fresh_pointer_blanks': 2,
            'valid_code_failure_probability_in_exact_arithmetic': 0,
            'physical_initialization_storage_readout_costs_included': False}


def history(d, steps=512, coupling=COUPLING, feedback_enabled=True, current=0):
    """Complete four-branch law of two observers on the same quantum history."""
    _, size = _register(d)
    if type(steps) is not int or steps < 1 or type(feedback_enabled) is not bool:
        raise ValueError('positive integer split count and explicit feedback switch required')
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('native current -1, 0 or 1 required')
    from . import protocol
    _, prepare_meta = scalar.preparation(d)
    kick_word = [c.pulse('X', -current*FEEDBACK_ANGLE, [0])]
    _, kick_guard = c.stabilize_executed_unitary(c.execute(kick_word, _register(d)[0]))
    kick_meta = {'guard': kick_guard, 'resources': _word_resources(kick_word),
                 'native_current': current}
    _, first_meta = scalar.source_unitary(d, steps, time=1., coupling=coupling)
    _, second_meta = scalar.source_unitary(d, steps, time=.5, coupling=coupling)
    _, read_a_meta = pointer(d, 0)
    _, read_b_meta = scalar_pointer(d)
    independent_first = scalar.independent_split_unitary(d, steps, time=1., coupling=coupling)
    independent_second = scalar.independent_split_unitary(d, steps, time=.5, coupling=coupling)
    # Production states and prefix maps come from the canonical engine. The
    # original metadata and independent target calculations remain this
    # compatibility receipt's diagnostic surface.
    program = protocol.sharp_program(d, steps, native_current=current,
                                     feedback=feedback_enabled, coupling=F(coupling))
    prefixes = protocol.execution_prefixes(program)
    terminal = {branch.record: branch for branch in prefixes[-1].branches}
    # The independent preparation is the declared scalar target with gauge |00>.
    from ..observer.source_scalar.field import vacuum_vector
    target = np.kron(vacuum_vector(d)[0], np.array([1., 0., 0., 0.]))
    x0 = np.zeros((size, size))
    for index in range(size):
        x0[index ^ 1, index] = 1
    target = expm(-1j*current*FEEDBACK_ANGLE*x0)@target
    completeness = np.zeros((size, size), complex)
    unconditional = np.zeros((size, size), complex)
    branches = []
    errors = []
    bit_a = np.array([index & 1 for index in range(size)])
    flux = np.array([-(index & 1)+(index >> 1 & 1) for index in range(size)])
    q0 = np.kron(np.kron(np.eye(d), quadratures(d)[0]), np.eye(4))
    for a in (0, 1):
        _, control_meta = feedback(d, a, feedback_enabled)
        oracle_control = _independent_feedback(size, a, feedback_enabled)
        for b in (0, 1):
            kraus = terminal[a, b].kraus
            state = terminal[a, b].state[:, 0]
            oracle = (scalar_projectors(d)[b]@independent_second@oracle_control
                      @projector(size, 0, a)@independent_first@target)
            probability = float(np.vdot(state, state).real)
            density = np.outer(state, state.conj())
            errors.append(float(np.linalg.norm(density-np.outer(oracle, oracle.conj()), 1)))
            completeness += kraus.conj().T@kraus
            unconditional += density
            branches.append({'record': [a, b], 'probability': probability,
                             'unnormalized_terminal_state': encoded(state),
                             'source_right_occupation_weight': float(np.dot(abs(state)**2, bit_a)),
                             'electric_flux_weight': float(np.dot(abs(state)**2, flux)),
                             'scalar_coordinate_weight': float(np.vdot(state, q0@state).real)})
    margin_a = [sum(row['probability'] for row in branches if row['record'][0] == a)
                for a in (0, 1)]
    margin_b = [sum(row['probability'] for row in branches if row['record'][1] == b)
                for b in (0, 1)]
    prefix_marginal = [branch.probability for branch in prefixes[2].branches]
    normalization_error = abs(sum(row['probability'] for row in branches)-1)
    kraus_error = float(np.linalg.norm(completeness-np.eye(size), 2))
    checks = {'normalized_complete_history': normalization_error < 2e-10,
              'full_input_channel_complete': kraus_error < 2e-10,
              'independent_complete_branch_densities': max(errors) < 2e-10,
              'later_operations_preserve_first_record_marginal':
                  max(abs(a-b) for a, b in zip(margin_a, prefix_marginal)) < 2e-10,
              'both_observers_read_one_joint_law':
                  abs(sum(margin_a)-sum(margin_b)) < 2e-10}
    if not all(checks.values()):
        raise AssertionError(checks)
    errors_time = [scalar.time_error_bound(d, steps, time=t, coupling=F(coupling))
                   for t in (F(1), F(1, 2))]
    ideal_error = min(F(1), sum((F(row['ideal_channel_half_diamond_error_upper'])
                               for row in errors_time), F(0)))
    return {'dimension': d, 'steps_per_segment': steps, 'coupling': coupling,
            'native_current': current,
            'feedback_enabled': feedback_enabled,
            'model_time_segments': [1., .5],
            'meaning': {'first_bit': 'species-0 right occupation after first evolution',
                        'second_bit': 'positive site-0 scalar coordinate after feedback and second evolution'},
            'branches': branches, 'first_record_marginal': margin_a,
            'second_record_marginal': margin_b,
            'normalization_error': normalization_error,
            'full_kraus_completeness_error': kraus_error,
            'maximum_branch_density_matrix_one_norm_error': max(errors),
            'terminal_gauge_scalar_covariance': float(
                np.trace(unconditional@(np.diag(bit_a)@q0)).real
                -np.trace(unconditional@np.diag(bit_a)).real*np.trace(unconditional@q0).real),
            'resources_per_realized_history': _resources(
                [prepare_meta, kick_meta, first_meta, read_a_meta, control_meta, second_meta, read_b_meta], d),
            'execution': {'preparation': prepare_meta, 'first_evolution': first_meta,
                          'native_current_preparation': kick_meta,
                          'second_evolution': second_meta,
                          'feedback': [feedback(d, a, feedback_enabled)[1] for a in (0, 1)],
                          'first_pointer': read_a_meta, 'second_pointer': read_b_meta},
            'ideal_finite_continuous_comparison': {
                'segment_certificates': errors_time,
                'complete_instrument_half_diamond_error_upper': str(ideal_error),
                'classical_history_tv_error_upper': str(ideal_error),
                'scope': 'ideal split versus continuous evolution of this same finite interacting Hamiltonian with the same interventions; numerical, physical service, occupation and spatial errors excluded'},
            'checks': checks, 'all_pass': True}


def run():
    from . import finite_signal, native_join
    source = sources()
    levels = []
    for d in DIMENSIONS:
        refinements = [history(d, steps) for steps in STEPS]
        coupled_off = history(d, STEPS[-1], feedback_enabled=False)
        disconnected_on = history(d, STEPS[-1], coupling=0.)
        disconnected_off = history(d, STEPS[-1], coupling=0., feedback_enabled=False)
        response = refinements[-1]['second_record_marginal'][1]-coupled_off['second_record_marginal'][1]
        disconnected = disconnected_on['second_record_marginal'][1]-disconnected_off['second_record_marginal'][1]
        if abs(disconnected) > 2e-10:
            raise AssertionError('gauge-only feedback changes a disconnected scalar marginal')
        if abs(response) < .02:
            raise AssertionError('declared interacting protocol lost its finite scalar response')
        levels.append({'dimension': d, 'refinements': refinements,
                       'coupled_without_feedback': coupled_off,
                       'disconnected_with_feedback': disconnected_on,
                       'disconnected_without_feedback': disconnected_off,
                       'scalar_feedback_response': response,
                       'disconnected_scalar_feedback_response': disconnected})
    exact_controlled = tuple(tuple(gauge.matmul(gauge.occupation_projector(0, bit), u)
                                  for bit in (0, 1))
                            for u in (gauge.ideal_electric(), gauge.ideal_echo()))
    return {'schema': 'oph.common_source.interacting_observer.v1',
            'source': source, 'gauge_matter_exact_audit': gauge.audit(),
            'checkpoint_adapter': records.native_checkpoint_certificate(),
            'exact_controlled_record_instrument': records.controlled_certificate(exact_controlled),
            'levels': levels, 'finite_response_certificate': finite_signal.certificate(steps=STEPS[-1]),
            'native_quantum_join': native_join.run(d=4, steps=STEPS[-1]), 'all_pass': True,
            'scope': {'carrier': 'one flux-dressed four-state CAR code times two neutral scalar codes and two committed classical bits',
                      'operations': 'quantum preparation, field/matter/electric evolution, pointers and feedback are emitted twelve-port source words; native integer repair remains the existing separately declared primitive in the joint instrument',
                      'interaction': 'g sqrt(2) Q0 (nL0-nR0), with g=1/4 a declared model coefficient',
                      'time': 'declared Hamiltonian model coordinate with intervals 1 and 1/2; no physical clock calibration',
                      'reference': 'faithful tracial reference for block-algebra repair; the separately prepared pure state is not selected by that unconstrained reference',
                      'errors': 'rational all-input finite temporal channel bounds; floating source execution checked against independent finite split matrices; separate interacting occupation and infinite coordinate-sign certificates are in infinite_results.json; numerical execution, physical services and spatial continuum errors remain separate',
                      'completion': 'finite common-process construction; not a complete OPH physical model or empirical validation'}}
