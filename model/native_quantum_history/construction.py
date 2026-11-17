"""Reconstruct native operator transfer and execute its coherent readout.

The uniform scheduler, binary carrier, coherent preparation and readout are
declared. The logarithm acts on operator space, not on carrier kets. Native
attempt counts are discrete; no physical clock or fractional CP flow is used.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..observer import native, quantum_memory, quantum_mixing as mixing
from ..primitive_observer import runtime


ROOT = Path(__file__).resolve().parents[2]
STATES = tuple(state for state in mixing.STATES if max(state) == 1)
INDEX = {state: i for i, state in enumerate(STATES)}
FULL_TO_BINARY = {mixing.INDEX[state]: i for i, state in enumerate(STATES)}
RATES = quantum_memory.COHERENCE_RATES
ATTEMPTS = (0, 1, 2, 4, 8)
TOLERANCE = 2e-10
SOURCE_PATHS = (
    'model/native_quantum_history/__init__.py',
    'model/native_quantum_history/__main__.py',
    'model/native_quantum_history/construction.py',
    'model/observer/native.py', 'model/observer/quantum_mixing.py',
    'model/observer/quantum_memory.py', 'model/observer/QUANTUM_MIXING.md',
    'model/primitive_observer/runtime.py', 'model/common_source/open_protocol.py',
    'model/common_source/open_environment.py', 'model/common_source/protocol.py',
    'model/common_source/slow_geometry.py',
    'model/observer/scalar_response/experiment.py', 'model/requirements.txt',
    'tools/strict_json.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/SourceLaw/ReflectionPositivity.lean',
)


def candidate(value):
    if type(value) is not str or value not in RATES:
        raise ValueError('native candidate C or E required')
    return value


def attempts(value):
    if type(value) is not int or not 0 <= value <= 64:
        raise ValueError('an integer native attempt count from zero to sixty-four required')
    return value


def bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= TOLERANCE:
        raise ArithmeticError(f'{name} exceeds the declared numerical tolerance')
    return value


def immutable(array):
    array = np.asarray(array)
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


@dataclass(frozen=True)
class Block:
    pairs: tuple[tuple[int, int], ...]
    matrix: np.ndarray
    values: np.ndarray
    vectors: np.ndarray


@lru_cache(maxsize=2)
def reconstruction(law):
    """Check every exact native matrix unit, then evaluate its spectral log.

    Each block is independently identified with the killed exclusion formula
    before its spectrum is evaluated. Rational diagonal-dominance bounds prove
    strict positivity without relying on the floating diagonalization.
    """
    law = candidate(law)
    blocks, lower, digest = [], F(1), hashlib.sha256()
    spectral_min, nonconstant_max, fixed, maximum_residual = 1., 0., 0, 0.
    for sector in mixing.sectors(binary=True):
        exact = mixing.native_block(law, sector.difference, binary=True)
        expected = mixing.binary_reduction(law, sector.difference)
        if exact != expected or exact != tuple(zip(*exact, strict=True)):
            raise ArithmeticError('native quantum block violates the exact exclusion attachment')
        if sector.difference == mixing.ZERO:
            reached, pending = {0}, [0]
            while pending:
                i = pending.pop()
                for j, value in enumerate(exact[i]):
                    if value and j not in reached:
                        reached.add(j)
                        pending.append(j)
            if len(reached) != 66:
                raise ArithmeticError('the diagonal native block is disconnected')
        elif mixing.killing_rate(law, sector.difference) <= 0:
            raise ArithmeticError('a nonconstant coherence sector has no strict loss')
        for i, row in enumerate(exact):
            lower = min(lower, row[i]-sum(abs(x) for j, x in enumerate(row) if j != i))
            if sum(row, F(0)) != 1-mixing.killing_rate(law, sector.difference):
                raise ArithmeticError('native block does not have its exact coherence loss')
        digest.update(json.dumps([sector.difference, [[str(x) for x in row] for row in exact]],
                                 separators=(',', ':')).encode())
        matrix = np.asarray(exact, float)
        values, vectors = np.linalg.eigh(matrix)
        if values[0] <= 0 or values[-1] > 1+TOLERANCE:
            raise ArithmeticError('positive native contraction has an inconsistent spectrum')
        maximum_residual = max(maximum_residual, float(np.linalg.norm(
            (vectors*values)@vectors.T-matrix, 'fro')))
        spectral_min = min(spectral_min, float(values[0]))
        fixed += int(np.count_nonzero(np.abs(values-1) < TOLERANCE))
        below = values[values < 1-TOLERANCE]
        if len(below):
            nonconstant_max = max(nonconstant_max, float(below[-1]))
        pairs = tuple((FULL_TO_BINARY[i], FULL_TO_BINARY[j]) for i, j in sector.pairs)
        blocks.append(Block(pairs, immutable(matrix), immutable(values), immutable(vectors)))
    wanted = F(2, 3) if law == 'C' else F(11, 30)
    if lower != wanted or fixed != 1 or sum(len(b.pairs) for b in blocks) != 66**2:
        raise ArithmeticError('full native operator-space positivity or ground multiplicity changed')
    expected_gap = -math.log((55+math.sqrt(5))/60)
    return tuple(blocks), {
        'all_binary_matrix_units_checked': 66**2,
        'block_dimensions': {str(k): v for k, v in sorted(Counter(len(b.pairs) for b in blocks).items())},
        'exact_native_block_sha256': digest.hexdigest(),
        'exact_lower_operator_bound': str(lower), 'upper_operator_bound': '1',
        'unique_ground_operator': 'I66/sqrt(66)',
        'ground_uniqueness_proof': 'connected diagonal block; every other sector has strictly subunit row sum',
        'operator_space_dimension': 66**2,
        'minimum_transfer_eigenvalue_diagnostic': spectral_min,
        'gap_per_Euclidean_native_attempt_diagnostic': -math.log(nonconstant_max),
        'gap_vs_native_exclusion_formula_residual': bounded(
            abs(-math.log(nonconstant_max)-expected_gap), 'operator-space gap'),
        'block_spectral_reconstruction_residual': bounded(maximum_residual, 'spectral block'),
        'Hamiltonian': 'H_L=-log(Phi_L), acting on Hilbert-Schmidt operator space',
        'all_complex_operator_inputs_and_reference_extensions_covered_by_matrix_unit_identity': True,
    }


def fractional_obstruction():
    """Exact diagonal input excluding positive quantum half-step evolution.

    Both laws preserve diagonals with the same K. The functional-calculus
    square root restricts to sqrt(K), so one negative population excludes
    positivity of the entire fractional quantum map, without a Choi matrix.
    """
    k = mixing.native_block('C', mixing.ZERO, binary=True)
    if k != mixing.native_block('E', mixing.ZERO, binary=True):
        raise ArithmeticError('the two laws no longer share their diagonal channel')
    d = tuple(tuple(F(int(i == j))-k[i][j] for j in range(66)) for i in range(66))
    if (min(k[i][i]-sum(abs(k[i][j]) for j in range(66) if j != i)
            for i in range(66)) != F(2, 3)):
        raise ArithmeticError('the half-step series norm bound changed')
    def address(ports):
        return INDEX[tuple(int(i in ports) for i in range(12))]
    row, column = address((3, 11)), address((2, 4))
    vector = tuple(F(int(i == column)) for i in range(66))
    partial = F(0)
    for n in range(1, 9):
        vector = tuple(sum((value*vector[j] for j, value in enumerate(values) if value), F(0))
                       for values in d)
        coefficient = F(math.comb(2*n, n), 4**n*(2*n-1))
        partial -= coefficient*vector[row]
    remainder = F(1, 3)**9/(1-F(1, 3))
    upper = partial+remainder
    if upper >= 0:
        raise ArithmeticError('the finite Taylor certificate did not exclude a positive half-step')
    return {'positive_input': '|ports 2,4><ports 2,4|',
        'negative_output_population': 'ports 3,11',
        'series_order': 8, 'exact_partial_sum': str(partial),
        'absolute_tail_bound': str(remainder), 'strictly_negative_upper_bound': str(upper),
        'both_candidates_excluded': ['C', 'E'],
        'fractional_parameter': 'one half Euclidean native attempt',
        'fractional_map_fails_positivity_already_on_a_diagonal_density': True,
        'integer_native_channels_affected': False}


def operator_action_summary():
    """Evaluate the repair-derived positive quadratic form on a phase mode.

    The form is built from H=-log(Phi) on Hilbert-Schmidt operator space.
    This summary deliberately distinguishes that transfer-space form from a
    carrier-ket Hamiltonian or a space-time action.
    """
    _, probe, _, _ = phase_apparatus()
    norm_squared = float(np.vdot(probe, probe).real)
    if abs(norm_squared-20.) > TOLERANCE:
        raise ArithmeticError('declared phase probe no longer has Hilbert-Schmidt norm squared 20')
    rows = {}
    for law in ('C', 'E'):
        rate = RATES[law]
        actual = apply_transfer(probe, law, 1)
        residual = bounded(float(np.linalg.norm(actual-float(rate)*probe, 'fro')),
                           'repair phase-mode eigenidentity')
        if residual > TOLERANCE:
            raise ArithmeticError('native transfer no longer preserves the declared phase mode')
        log_rate = {'C': 'log(6/5)', 'E': 'log(10/7)'}[law]
        rows[law] = {
            'probe': 'S=A+A* from the declared first-seam phase readout',
            'hilbert_schmidt_norm_squared': 20,
            'exact_one_attempt_transfer_eigenvalue': str(rate),
            'operator_log_eigenvalue_per_native_attempt': log_rate,
            'positive_quadratic_form_half_H_on_probe': f'10*{log_rate}',
            'norm_squared_after_n_native_attempts_over_initial': f'({rate})^(2*n)',
            'unitary_conjugation_excluded_for_forgotten_record_map': True,
            'repair_mode_residual_diagnostic': residual,
        }
    return {
        'space': 'Hilbert-Schmidt operator space of the 66-state binary carrier',
        'form': 'Q_L[O]=<O,(-log(Phi_L))O>_HS/2 >= 0',
        'parameter': 'integer Euclidean native repair-attempt count',
        'phase_mode_results': rows,
        'physical_space_time_action_derived': False,
        'carrier_unitary_dynamics_derived': False,
        'principal_log_generates_a_continuous_CP_semigroup': False,
        'other_continuous_generator_or_environmental_dilation_excluded': False,
    }


def apply_transfer(matrix, law, count, *, reference_dimension=1, spectral=False):
    """Independent block action on arbitrary operators, including references."""
    law, count = candidate(law), attempts(count)
    if (type(reference_dimension) is not int or reference_dimension < 1
            or reference_dimension > 3 or type(spectral) is not bool):
        raise ValueError('reference dimension one to three and Boolean spectral switch required')
    matrix = np.asarray(matrix, complex)
    size = 66*reference_dimension
    if matrix.shape != (size, size) or not np.isfinite(matrix).all():
        raise ValueError('finite matrix of the declared work/reference dimension required')
    tensor = matrix.reshape(66, reference_dimension, 66, reference_dimension)
    output = np.zeros_like(tensor)
    for block in reconstruction(law)[0]:
        left, right = np.array(block.pairs).T
        transfer = ((block.vectors*np.exp(count*np.log(block.values)))@block.vectors.T
                    if spectral else np.linalg.matrix_power(block.matrix, count))
        output[left, :, right, :] = np.einsum('ab,brs->ars', transfer, tensor[left, :, right, :])
    return output.reshape(matrix.shape)


def native_instrument(law):
    """Compile the actual native Kraus maps, summing every committed record."""
    law = candidate(law)
    groups = {}
    for number, seam in enumerate(mixing.seams()):
        for coin in (0, 1):
            for i, state in enumerate(STATES):
                output, record = native.native(state, seam, coin)
                if record != (coin, state[seam[0]]-output[seam[0]]) or output not in INDEX:
                    raise ArithmeticError('native coin/current or binary carrier changed')
                label = record[1] if law == 'C' else state[seam[0]]-state[seam[1]]
                key = number, coin, label
                groups.setdefault(key, np.zeros((66, 66), complex))[INDEX[output], i] = 1/math.sqrt(60)
    return runtime.Instrument(kraus=(tuple(groups[key] for key in sorted(groups)),),
        resources=(('native_attempts', 1),),
        binding={'candidate': law, 'scheduler': 'uniform thirty seams and fair placement coin',
                 'native_record': 'current' if law == 'C' else 'input endpoint imbalance',
                 'all_native_records_forgotten': True, 'native_pulse_words_emitted': False})


def phase_apparatus(sign=1):
    """The declared common-spectator preparation and complete binary pointer."""
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('preparation phase sign must be minus or plus one')
    u, v = mixing.seams()[0]
    a = np.zeros((66, 66), complex)
    rho = np.zeros_like(a)
    pairs = []
    for spectator in range(12):
        if spectator in (u, v):
            continue
        left, right = [0]*12, [0]*12
        left[u] = left[spectator] = right[v] = right[spectator] = 1
        i, j = INDEX[tuple(left)], INDEX[tuple(right)]
        pairs.append((i, j))
        a[i, j] = 1
        vector = np.zeros(66, complex)
        vector[i], vector[j] = 1/math.sqrt(2), sign/math.sqrt(2)
        rho += np.outer(vector, vector.conj())/10
    s, identity = a+a.conj().T, np.eye(66)
    q = s@s
    plus, minus = (q+s)/2+(identity-q)/math.sqrt(2), (q-s)/2+(identity-q)/math.sqrt(2)
    read = runtime.Instrument(kraus=((plus,), (minus,)), visible=True,
        resources=(('phase_readouts', 1),),
        binding={'outcomes': ['plus', 'minus'], 'effects': '(I +/- A +/- A*)/2',
                 'postmeasurement_maps': 'positive square roots of the declared effects',
                 'coherent_preparation_and_readout_supplied': True})
    return rho, s, read, tuple(pairs)


def program(law, count, *, sign=1, compiled=None):
    law, count = candidate(law), attempts(count)
    rho, _, read, _ = phase_apparatus(sign)
    hidden = native_instrument(law) if compiled is None else compiled
    events = [runtime.Event(f'unread_native_attempt_{i+1}', instruments=(hidden,)) for i in range(count)]
    events.append(runtime.Event('complete_phase_readout', instruments=(read,)))
    return runtime.Program(f'native_quantum_history_{law}_{count}_{sign}', 66, rho, tuple(events),
        binding={'candidate': law, 'native_attempts': count,
                 'preparation': 'uniform ten-spectator phase-pair mixture', 'phase_sign': sign,
                 'complete_record_dictionary': ['plus', 'minus'],
                 'parameter': 'integer native repair attempts; no physical calibration'})


def reference_preparation():
    """Entangled mixture with a complex relative phase and untouched qubit."""
    _, _, _, pairs = phase_apparatus()
    rho = np.zeros((132, 132), complex)
    for i, j in pairs:
        vector = np.zeros(132, complex)
        vector[2*i], vector[2*j+1] = 1/math.sqrt(2), 1j/math.sqrt(2)
        rho += np.outer(vector, vector.conj())/10
    return rho


def terminal_oracle(initial, law, count, read, reference_dimension=1):
    evolved = apply_transfer(initial, law, count, reference_dimension=reference_dimension)
    outputs = []
    for group in read.kraus:
        k = np.kron(group[0], np.eye(reference_dimension))
        outputs.append(k@evolved@k.conj().T)
    return tuple(outputs)


def run():
    quantum_memory.phase_pair_witness()
    rows, candidates, preparations, channel_outputs = [], {}, {}, {}
    for law in ('C', 'E'):
        _, candidates[law] = reconstruction(law)
        hidden = native_instrument(law)
        for count in ATTEMPTS:
            for sign in (1, -1):
                p = program(law, count, sign=sign, compiled=hidden)
                prefix = runtime.execute(p, stop_after=count)
                history = runtime.execute(p, checkpoint=prefix)
                if history.records != ((0,), (1,)) or prefix.records != ((),):
                    raise ArithmeticError('the complete phase pointer dictionary changed')
                expected_plus = (1+sign*RATES[law]**count)/2
                targets = (expected_plus, 1-expected_plus)
                read = p.events[-1].instruments[0]
                oracle = terminal_oracle(p.initial_density, law, count, read)
                residual = sum(np.linalg.norm(branch.density-wanted, 'fro')
                               for branch, wanted in zip(history.branches, oracle, strict=True))
                spectral = apply_transfer(p.initial_density, law, count, spectral=True)
                spectral_error = np.linalg.norm(prefix.branches[0].density-spectral, 'fro')
                error = max(abs(branch.probability-float(wanted)) for branch, wanted in
                            zip(history.branches, targets, strict=True))
                preparations[law, count, sign] = prefix.branches[0].density
                channel_outputs[law, count, sign] = history
                rows.append({'candidate': law, 'native_attempts': count, 'preparation_sign': sign,
                    'complete_records': ['plus', 'minus'],
                    'exact_probabilities': list(map(str, targets)),
                    'runtime_probabilities': [branch.probability for branch in history.branches],
                    'quantum_output_dimension_per_record': 66,
                    'complete_terminal_density_residual_diagnostic': bounded(residual, 'complete phase output'),
                    'probability_residual_diagnostic': bounded(error, 'exact phase probability'),
                    'exp_minus_n_H_vs_runtime_density_residual_diagnostic': bounded(spectral_error, 'quantum transfer'),
                    'expanded_tree_resources': dict(history.expanded_tree_resources)})
    discrimination = []
    for count in ATTEMPTS:
        delta = quantum_memory.phase_pair_signal(count)
        state_difference = preparations['C', count, 1]-preparations['E', count, 1]
        trace_distance = float(np.linalg.norm(state_difference, 'nuc')/2)
        cq_distance = runtime.complete_distance(channel_outputs['C', count, 1],
                                               channel_outputs['E', count, 1])
        memory = {}
        for law in ('C', 'E'):
            distance = float(np.linalg.norm(preparations[law, count, 1]-preparations[law, count, -1], 'nuc')/2)
            memory[law] = {'exact_phase_memory_distance': str(RATES[law]**count),
                          'runtime_distance': distance,
                          'distance_residual_diagnostic': bounded(abs(distance-float(RATES[law]**count)), 'phase memory')}
        discrimination.append({'native_attempts': count,
            'exact_C_minus_E_plus_probability': str(delta),
            'exact_equal_prior_success_probability': str((1+delta)/2),
            'complete_pre_readout_trace_distance': trace_distance,
            'complete_post_readout_cq_distance': cq_distance,
            'optimal_binary_readout_residual_diagnostic': bounded(max(abs(trace_distance-float(delta)),
                    abs(cq_distance-float(delta))), 'complete C/E discrimination'), 'phase_memory': memory})
    reference_rows = []
    for law in ('C', 'E'):
        p, initial = program(law, 4), reference_preparation()
        prefix = runtime.execute(p, initial=initial, stop_after=2)
        restored = runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix))
        continued = runtime.execute(p, checkpoint=restored)
        direct = runtime.execute(p, initial=initial)
        oracle = terminal_oracle(initial, law, 4, p.events[-1].instruments[0], 2)
        residual = sum(np.linalg.norm(branch.density-wanted, 'fro')
                       for branch, wanted in zip(direct.branches, oracle, strict=True))
        total = sum(branch.density for branch in direct.branches).reshape(66, 2, 66, 2)
        marginal = np.trace(total, axis1=0, axis2=2)
        reference_rows.append({'candidate': law, 'native_attempts': 4, 'reference_dimension': 2,
            'quantum_output_dimension_per_record': 132,
            'complete_record_probabilities': [branch.probability for branch in direct.branches],
            'exact_complete_record_probabilities': ['1/2', '1/2'],
            'complete_record_probability_residual_diagnostic': bounded(
                max(abs(branch.probability-.5) for branch in direct.branches), 'reference pointer probability'),
            'complete_terminal_reference_density_residual_diagnostic': bounded(residual, 'reference output'),
            'untouched_reference_marginal_residual_diagnostic': bounded(np.linalg.norm(marginal-np.eye(2)/2), 'reference marginal'),
            'checkpoint_after_native_attempts': 2, 'checkpoint_independently_replayed': True,
            'checkpoint_complete_cq_distance_diagnostic': bounded(runtime.complete_distance(continued, direct), 'checkpoint')})
    return {'schema': 'oph.native_quantum_history.v1', 'all_pass': True,
        'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'operator_transfer': candidates, 'executed_phase_measurements': rows,
        'operator_space_quadratic_form': operator_action_summary(),
        'fractional_channel_obstruction': fractional_obstruction(),
        'discrimination': discrimination,
        'optimal_count_for_declared_phase_pair_experiment': 4,
        'maximum_C_minus_E_plus_probability': str(quantum_memory.phase_pair_signal(4)),
        'reference_and_checkpoint': reference_rows,
        'scope': {'carrier_and_uniform_scheduler_declared': True,
            'coherent_preparation_and_binary_effects_declared': True,
            'native_channel_semantics_executed': True,
            'native_pulse_words_emitted': False,
            'all_native_records_retained': False,
            'complete_final_measurement_records_and_quantum_outputs_retained': True,
            'spectral_log_is_unique_on_full_binary_operator_space': True,
            'log_is_a_carrier_ket_Hamiltonian_or_Lindblad_generator': False,
            'fractional_transfer_powers_are_physical_channels': False,
            'all_quantum_history_reflection_axioms_proved': False,
            'physical_clock_or_continuum_QFT_derived': False,
            'empirical_dataset_used_or_empirical_agreement_claimed': False,
            'roundoff_or_physical_service_certificate': False}}
