"""Reduce the coherent apparatus premise and check its native boundary.

Actual native records supply the herald and phase effects once one balanced
real mixer is admitted. Native partial injections cannot supply that mixer
from diagonal preparations by classical conditioning or feedback alone.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from ..observer import native
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/native_interface/__init__.py', 'model/native_interface/__main__.py',
    'model/native_interface/construction.py', *quantum.SOURCE_PATHS,
    'Lean/Geometry/SourceRepairInstrument.lean',
    'Lean/QFT/SourcePhaseLiftBridge.lean', 'Lean/QFT/ConjugationGauge.lean',
    'Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean',
)))
TOLERANCE = 2e-10
SUCCESS = F(5, 66)
COUNTS = (1, 2, 4)


def bounded(value, label):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= TOLERANCE:
        raise ArithmeticError(f'{label} exceeds the native interface tolerance')
    return value


@dataclass(frozen=True)
class RecordedMap:
    coin: int
    label: int
    current: int
    mapping: tuple[tuple[int, int], ...]
    matrix: np.ndarray


def fixed_maps(law, seam=None):
    """Actual native fair-coin branches on one fixed seam, including zeros."""
    law = quantum.candidate(law)
    seam = quantum.mixing.seams()[0] if seam is None else seam
    if seam not in quantum.mixing.seams():
        raise ValueError('one of the declared oriented native seams required')
    groups, currents = {}, {}
    for coin in (0, 1):
        for i, state in enumerate(quantum.STATES):
            output, record = native.native(state, seam, coin)
            if output not in quantum.INDEX or record != (coin, state[seam[0]]-output[seam[0]]):
                raise ArithmeticError('native record grammar or binary domain changed')
            current = record[1]
            label = current if law == 'C' else state[seam[0]]-state[seam[1]]
            key = coin, label
            if current not in (-1, 0, 1) or key in currents and currents[key] != current:
                raise ArithmeticError('the committed record fails to determine the current')
            currents[key] = current
            groups.setdefault(key, {})[i] = quantum.INDEX[output]
    result = []
    for (coin, label), mapping in sorted(groups.items()):
        if len(set(mapping.values())) != len(mapping):
            raise ArithmeticError('a recorded native map is not a partial injection')
        matrix = np.zeros((66, 66), complex)
        for i, j in mapping.items():
            matrix[j, i] = 1/math.sqrt(2)
        result.append(RecordedMap(coin, label, currents[coin, label], tuple(sorted(mapping.items())), matrix))
    return tuple(result)


@lru_cache(maxsize=1)
def phase_pairs():
    """Pairing read directly from the actual fixed-seam moving map."""
    branch, = (branch for branch in fixed_maps('C') if branch.current == -1)
    pairs = tuple((output, source) for source, output in branch.mapping)
    if len(pairs) != 10 or len({i for pair in pairs for i in pair}) != 20:
        raise ArithmeticError('the native current does not define ten disjoint orientation pairs')
    return pairs


def mixer():
    """One declared balanced mixer; native branch support fixes its pairing.

    Membership in the generated matrix algebra does not make this an
    available operation of the committed native channel grammar.
    """
    result = np.eye(66, dtype=complex)
    for left, right in phase_pairs():
        result[np.ix_((left, right), (left, right))] = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    return result


def herald_and_pointer(law):
    """Both instruments use only committed native coin/current processing."""
    branches = fixed_maps(law)
    herald = runtime.Instrument(
        kraus=tuple(tuple(b.matrix for b in branches if b.current == current)
                    for current in (0, -1, 1)), visible=True,
        resources=(('native_attempts', 1), ('native_record_readouts', 1)),
        binding={'candidate': law, 'outcomes': ['failed_herald', 'plus_phase_herald', 'minus_phase_herald'],
                 'success_rule': 'nonzero native signed current; retain its sign',
                 'all_failure_branches_retained': True})
    pointer = runtime.Instrument(
        kraus=tuple(tuple(b.matrix for b in branches
                          if int(2*b.current+2*b.coin-1 < 0) == outcome) for outcome in (0, 1)),
        visible=True, resources=(('native_attempts', 1), ('native_record_readouts', 1)),
        binding={'candidate': law, 'outcomes': ['plus', 'minus'],
                 'plus_rule': '2*native_current + 2*native_coin - 1 > 0',
                 'hidden_load_access': False, 'native_pulse_words_emitted': False})
    return herald, pointer


def exact_native_audit(law):
    """Check monomial closure and exact herald/effect counts from all sources."""
    branch_count, addresses = 0, 0
    for seam in quantum.mixing.seams():
        for branch in fixed_maps(law, seam):
            branch_count += 1
            addresses += len(branch.mapping)
            # One entry in each nonzero column and row proves the diagonal
            # commutation and diagonal Heisenberg-effect identities exactly.
            inputs, outputs = zip(*branch.mapping, strict=True)
            if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
                raise ArithmeticError('native coherence obstruction lost its monomial premise')
    maps = fixed_maps(law)
    effects = [[0]*66 for _ in (0, 1)]
    success_population = {sign: [F(0)]*66 for sign in (-1, 1)}
    for branch in maps:
        outcome = int(2*branch.current+2*branch.coin-1 < 0)
        for i, j in branch.mapping:
            effects[outcome][i] += 1
            if branch.current:
                success_population[branch.current][j] += F(1, 132)
    pairs = phase_pairs()
    lefts, rights = {i for i, _ in pairs}, {j for _, j in pairs}
    expected_plus = tuple(F(1) if i in lefts else F(0) if i in rights else F(1, 2) for i in range(66))
    if tuple(F(value, 2) for value in effects[0]) != expected_plus:
        raise ArithmeticError('native coin/current readout has the wrong orientation effect')
    if any(a+b != 2 for a, b in zip(*effects, strict=True)):
        raise ArithmeticError('the complete native orientation pointer is not normalized')
    for current, support in ((-1, lefts), (1, rights)):
        population = success_population[current]
        if sum(population) != SUCCESS or any(
                value != (SUCCESS/10 if i in support else 0) for i, value in enumerate(population)):
            raise ArithmeticError('the source herald does not prepare the uniform orientation mixture')
    return {'candidate': law, 'native_branch_addresses_checked': addresses,
            'recorded_partial_injection_maps_checked': branch_count,
            'exact_herald_probability_per_sign': str(SUCCESS),
            'exact_total_useful_herald_probability': str(2*SUCCESS),
            'exact_failure_probability': str(1-2*SUCCESS),
            'conditional_herald_population_per_selected_orientation_state': '1/10',
            'orientation_pointer_plus_effect_diagonal': list(map(str, expected_plus)),
            'native_outcomes_commute_with_work_dephasing': True,
            'native_classical_feedback_cannot_create_coherence_from_diagonal_inputs': True,
            'native_record_effects_are_diagonal_on_all_finite_adaptive_paths': True,
            'joint_quantum_outputs_are_claimed_equal_after_dephasing': False}


def program(law, count):
    law, count = quantum.candidate(law), quantum.attempts(count)
    if count == 0:
        raise ValueError('at least one unread native attempt is required in this experiment')
    herald, pointer = herald_and_pointer(law)
    identity = runtime.Instrument(kraus=((np.eye(66),),))
    balanced = runtime.Instrument(kraus=((mixer(),),), resources=(('supplied_balanced_mixers', 1),),
        binding={'one_additional_coherent_service': True, 'source_selected_or_native_reachable': False})
    hidden = quantum.native_instrument(law)
    failure_read = runtime.Instrument(kraus=((np.eye(66),), (np.zeros((66, 66)),)), visible=True,
        binding={'failed_preparation_record': 0, 'phase_measurement_executed': False})
    opposite_pointer = runtime.Instrument(kraus=tuple(reversed(pointer.kraus)), visible=True,
        resources=pointer.resources,
        binding={'candidate': law, 'outcomes': ['matches_prepared_phase', 'opposes_prepared_phase'],
                 'classical_postprocessing': 'reverse native pointer bit for minus preparation'})
    events = [runtime.Event('native_preparation_herald', instruments=(herald,)),
              runtime.Event('prepare_phase_if_heralded', instruments=(identity, balanced, balanced), control=0)]
    for n in range(count):
        events.append(runtime.Event(f'unread_native_attempt_{n+1}', instruments=(identity, hidden, hidden), control=0))
    events.extend((runtime.Event('analyze_phase_if_heralded', instruments=(identity, balanced, balanced), control=0),
                   runtime.Event('native_coin_current_pointer', instruments=(failure_read, pointer, opposite_pointer), control=0)))
    return runtime.Program(f'native_interface_{law}_{count}', 66, np.eye(66)/66, tuple(events),
        binding={'stationary_diagonal_input_supplied': True,
                 'complete_records': ['failed,0', 'failed,1 impossible', 'plus_herald,match',
                                      'plus_herald,opposite', 'minus_herald,match', 'minus_herald,opposite'],
                 'herald_conditioning_is_reported_with_its_probability': True,
                 'balanced_mixer_is_additional_interface_service': True,
                 'parameter': 'native attempts; physical duration not supplied'})


def oracle(program, law, count, initial=None):
    """Independent block evolution with explicit preparation/analyzer outputs."""
    initial = program.initial_density if initial is None else initial
    reference = initial.shape[0]//66
    def apply(group, density):
        return sum((np.kron(k, np.eye(reference))@density@np.kron(k.conj().T, np.eye(reference))
                    for k in group), np.zeros_like(density))
    herald, pointer = herald_and_pointer(law)
    failure = apply(herald.kraus[0], initial)
    unitary = np.kron(mixer(), np.eye(reference))
    output = {(0, 0): failure, (0, 1): np.zeros_like(failure)}
    for flag in (1, 2):
        prepared = apply(herald.kraus[flag], initial)
        prepared = unitary@prepared@unitary.conj().T
        evolved = quantum.apply_transfer(prepared, law, count, reference_dimension=reference)
        analyzed = unitary@evolved@unitary.conj().T
        for bit in (0, 1):
            raw = bit if flag == 1 else 1-bit
            output[flag, bit] = apply(pointer.kraus[raw], analyzed)
    return output


def pointer_audit(law):
    """Same phase effects, explicitly different postmeasurement operation."""
    rho, s, old_pointer, _ = quantum.phase_apparatus()
    _, pointer = herald_and_pointer(law)
    b = mixer()
    error = 0.
    for outcome, group in enumerate(pointer.kraus):
        effect = sum((b.conj().T@k.conj().T@k@b for k in group), np.zeros((66, 66), complex))
        target = (np.eye(66)+(1 if outcome == 0 else -1)*s)/2
        error = max(error, float(np.linalg.norm(effect-target, 'fro')))
    actual = sum(k@b@rho@b.conj().T@k.conj().T for k in pointer.kraus[0])
    old = old_pointer.kraus[0][0]@rho@old_pointer.kraus[0][0].conj().T
    distance = float(np.linalg.norm(actual-old, 'nuc')/2)
    if abs(distance-.5) > TOLERANCE:
        raise ArithmeticError('same-effect native detector no longer has the expected disturbance control')
    return {'candidate': law, 'complete_effect_matrix_residual_diagnostic': bounded(error, 'native phase effects'),
            'certain_plus_state_native_vs_old_pointer_trace_distance': distance,
            'same_phase_outcome_laws_for_every_work_reference_input': True,
            'same_complete_postmeasurement_instrument_as_old_square_root_pointer': False,
            'whole_carrier_effects_are_projective': False,
            'projective_selection_theorem_domain': 'twenty-dimensional active orientation corner',
            'certain_state_invariance_required_for_Lueders_selection_satisfied': False}


def real_phase_boundary(law):
    """The added real mixer still cannot supply imaginary-phase tomography."""
    _, pointer = herald_and_pointer(law)
    b = mixer()
    states = []
    for sign in (1, -1):
        density = np.zeros((66, 66), complex)
        for i, j in phase_pairs():
            vector = np.zeros(66, complex)
            vector[i], vector[j] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
            density += np.outer(vector, vector.conj())/10
        states.append(density)
    outputs = []
    for state in states:
        read = runtime.Instrument(kraus=tuple(tuple(k@b for k in group) for group in pointer.kraus), visible=True)
        p = runtime.Program('real_native_interface_imaginary_phase_control', 66, state,
                            (runtime.Event('actual_native_readout', instruments=(read,)),))
        outputs.append(runtime.execute(p))
    probabilities = [[branch.probability for branch in history.branches] for history in outputs]
    residual = max(abs(value-.5) for row in probabilities for value in row)
    return {'candidate': law, 'opposite_imaginary_phase_state_trace_distance': float(np.linalg.norm(states[0]-states[1], 'nuc')/2),
            'complete_native_readout_probabilities': probabilities,
            'equal_half_probability_residual_diagnostic': bounded(residual, 'real-effect blindness'),
            'full_complex_tomography_supplied_by_one_real_mixer': False,
            'conjugation_fixed_effect_obstruction_applies': True}


def reference_input():
    """Coherence within each herald support survives its recorded branch."""
    pairs = phase_pairs()
    vector = np.zeros(132, complex)
    vector[2*pairs[0][1]], vector[2*pairs[1][1]+1] = .5, .5j
    vector[2*pairs[2][0]], vector[2*pairs[3][0]+1] = .5, .5
    return np.outer(vector, vector.conj())


def run():
    rows, references = [], []
    audits = {law: exact_native_audit(law) for law in ('C', 'E')}
    for law in ('C', 'E'):
        for count in COUNTS:
            p = program(law, count)
            history = runtime.execute(p)
            expected = oracle(p, law, count)
            plus = (1+quantum.RATES[law]**count)/2
            probabilities = (1-2*SUCCESS, F(0), SUCCESS*plus, SUCCESS*(1-plus), SUCCESS*plus, SUCCESS*(1-plus))
            if history.records != ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)):
                raise ArithmeticError('complete herald and native pointer dictionary changed')
            residual = max(np.linalg.norm(branch.density-expected[branch.record], 'fro') for branch in history.branches)
            probability_error = max(abs(branch.probability-float(value)) for branch, value in zip(history.branches, probabilities, strict=True))
            for branch in history.branches:
                if branch.probability > TOLERANCE:
                    costs = dict(branch.resources)
                    native_cost = count+2 if branch.record[0] else 1
                    if costs.get('native_attempts') != native_cost or costs.get('supplied_balanced_mixers', 0) != 2*int(bool(branch.record[0])):
                        raise ArithmeticError('failed herald or successful apparatus has wrong native cost')
            prefix = runtime.execute(p, stop_after=2)
            preparation_error = 0.
            for flag, sign in ((1, 1), (2, -1)):
                prepared = next(branch.density for branch in prefix.branches if branch.record == (flag,))/float(SUCCESS)
                preparation_error = max(preparation_error, np.linalg.norm(prepared-quantum.phase_apparatus(sign)[0], 'fro'))
            resumed = runtime.execute(p, checkpoint=runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix)))
            rows.append({'candidate': law, 'unread_native_attempts_after_herald': count,
                'complete_records': [list(record) for record in history.records],
                'exact_joint_probabilities': list(map(str, probabilities)),
                'runtime_joint_probabilities': [branch.probability for branch in history.branches],
                'conditional_prepared_phase_match_probabilities': [str(plus), str(1-plus)],
                'herald_probability_per_sign': str(SUCCESS),
                'total_useful_herald_probability': str(2*SUCCESS), 'failure_probability': str(1-2*SUCCESS),
                'full_terminal_density_dimension': 66,
                'conditional_native_preparation_vs_previous_phase_service_residual': bounded(preparation_error, 'native heralded preparation'),
                'complete_terminal_density_residual_diagnostic': bounded(residual, 'native interface terminal state'),
                'complete_probability_residual_diagnostic': bounded(probability_error, 'native interface law'),
                'checkpoint_complete_cq_distance_diagnostic': bounded(runtime.complete_distance(resumed, history), 'herald checkpoint'),
                'expanded_tree_resources': dict(history.expanded_tree_resources)})
        p = program(law, 1)
        # The arbitrary entangled input exercises the complete instrument;
        # it is not claimed to obey the stationary-input herald law.
        initial = reference_input()
        actual = runtime.execute(p, initial=initial)
        expected = oracle(p, law, 1, initial)
        error = max(np.linalg.norm(branch.density-expected[branch.record], 'fro') for branch in actual.branches)
        tensor = sum(branch.density for branch in actual.branches).reshape(66, 2, 66, 2)
        references.append({'candidate': law, 'reference_dimension': 2,
            'terminal_density_dimension_per_record': 132,
            'complete_reference_output_residual_diagnostic': bounded(error, 'reference terminal states'),
            'untouched_reference_marginal_residual_diagnostic': bounded(np.linalg.norm(np.trace(tensor, axis1=0, axis2=2)-np.eye(2)/2), 'reference marginal'),
            'stationary_herald_probability_asserted_for_entangled_override': False})
    delta = quantum.quantum_memory.phase_pair_signal(4)
    return {'schema': 'oph.native_interface.v1', 'all_pass': True,
        'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'native_reachability_and_record_audit': audits,
        'complete_heralded_measurements': rows,
        'native_phase_pointer': [pointer_audit(law) for law in ('C', 'E')],
        'reference_checks': references,
        'four_attempt_discrimination': {
            'conditional_prepared_phase_match_probability_difference': str(delta),
            'unconditional_complete_record_total_variation': str(2*SUCCESS*delta),
            'useful_preparation_yield': str(2*SUCCESS),
            'both_phase_heralds_used_with_classical_bit_interpretation': True},
        'remaining_complex_phase_boundary': [real_phase_boundary(law) for law in ('C', 'E')],
        'scope': {'native_preparation_herald_and_binary_effects_constructed': True,
            'one_shared_balanced_mixer_service_supplied': True,
            'mixer_derived_as_an_executable_native_operation': False,
            'stationary_binary_input_and_uniform_evolution_scheduler_supplied': True,
            'all_failed_heralds_and_complete_outputs_retained': True,
            'only_successful_preparations_used_for_conditional_phase_law': True,
            'independent_target_density_or_POVM_received_by_runtime': False,
            'native_source_pulse_words_emitted': False,
            'full_complex_observer_interface_or_physical_clock_derived': False,
            'new_empirical_measurement_or_validation': False,
            'new_Lean_proof_or_outward_roundoff_certificate': False}}
