"""Program the active-corner modular flow of the native Palm reference.

The reference is diag(d/600), while the probe preparation is still the
existing uniform native herald. Modular parameter, programmed phase and
native event count remain distinct; no identification with physical duration
or the natural dissipative native evolution is made.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..common_source import exact_output as exact
from ..common_source import output_certificate as report
from ..compiled_native_interface import construction as compiled
from ..event_circulation import construction as events
from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..primitive_observer import runtime


ROOT = Path(__file__).resolve().parents[2]
QUARTERS = (0, 1, 2, 4)
PHASE_ERROR_GATE = F(1, 10**10)
SOURCE_PATHS = tuple(sorted(set((*compiled.SOURCE_PATHS, *events.SOURCE_PATHS,
    'model/native_modular_reference/__init__.py', 'model/native_modular_reference/__main__.py',
    'model/native_modular_reference/construction.py',
    'Lean/Dynamics/PrivateInner.lean', 'docs/INTERPRETATION_AND_TERMS.md'))))


def quarter(value):
    if type(value) is not int or value not in QUARTERS:
        raise ValueError('a declared phase setting 0,1,2,4 times pi/2 is required')
    return value


def reference():
    degree, palm, _, _ = events.embedded()
    pairs = interface.phase_pairs()
    census = Counter((degree[a], degree[b]) for a, b in pairs)
    if census != {(8, 10): 2, (10, 8): 2, (8, 8): 2, (10, 10): 4}:
        raise ArithmeticError('the native modular pair census changed')
    active = {i for pair in pairs for i in pair}
    mask = tuple(sorted(compiled.circuit.encode(quantum.STATES[i]) for i in active if degree[i] == 8))
    if len(active) != 20 or len(mask) != 8 or any(code >= 128 for code in mask):
        raise ArithmeticError('active modular corner no longer has the declared seven-bit embedding')
    if min(palm) <= 0 or sum(palm) != 1 or sum(palm[i] for i in active) != F(23, 75):
        raise ArithmeticError('native Palm reference or active weight changed')
    return degree, palm, pairs, mask


def word(setting):
    setting = quarter(setting)
    phases = np.zeros(128)
    phases[list(reference()[3])] = setting*math.pi/2
    return compiled.compiler.compile_diagonal(phases)


@lru_cache(maxsize=4)
def _executed(setting):
    matrix = compiled.compiler.execute(word(setting), 7)
    return matrix.dtype.str, matrix.tobytes()


def local_matrix(setting):
    dtype, payload = _executed(quarter(setting))
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(128, 128)


def full_matrix(setting):
    matrix = np.kron(np.eye(4, dtype=int), local_matrix(setting))
    return np.frombuffer(matrix.tobytes(), dtype=matrix.dtype).reshape(512, 512)


def ideal_local(setting):
    result = np.eye(128, dtype=complex)
    indices = list(reference()[3])
    result[indices, indices] = 1j**quarter(setting)
    return result


def certify(setting, matrix=None):
    """Exact original-binary comparison on all 128 local/padding coordinates.

    These four diagonal targets have rational complex entries. Tensoring the
    map difference with any spectator/reference identity leaves operator
    norms unchanged, so no Frobenius dimension multiplier is introduced.
    """
    setting = quarter(setting)
    matrix = local_matrix(setting) if matrix is None else matrix
    if type(matrix) is not np.ndarray:
        raise ValueError('an original base ndarray is required before dispatch')
    matrix = exact._array(matrix, (128, 128))
    target = ideal_local(setting)
    squared = F(0)
    for actual, wanted in zip(matrix.flat, target.flat, strict=True):
        ar, ai = exact._component(actual)
        wr, wi = exact._component(wanted)
        squared += (ar-wr)**2+(ai-wi)**2
    epsilon = exact._sqrt_interval(squared)[1]
    return {'operator_error_upper': epsilon, 'frobenius_error_squared': squared,
            'half_diamond_upper': epsilon*(2+epsilon)/2, 'channel_norm_upper': (1+epsilon)**2,
            'full_spectator_dimension': 512, 'arbitrary_untouched_references': True}


def probabilities(setting):
    cosine = {0: F(1), 1: F(0), 2: F(-1), 4: F(1)}[quarter(setting)]
    return F(4, 5)+cosine/5


def composite_bound(modular, mixer):
    """Amplification-aware half-diamond bound for raw B D B execution.

    Every ideal map is CPTP. For actual norm gains g and errors delta,
    E_next <= g_next E_previous + delta_next, including untouched references.
    Exact adaptive failed branches and ideal native instruments contract it.
    """
    delta_b, g_b = mixer['half_diamond_upper'], mixer['channel_norm_upper']
    delta_d, g_d = modular['half_diamond_upper'], modular['channel_norm_upper']
    return delta_b*(1+g_d*g_b)+delta_d*g_b


def _lift_native(instrument, *, herald=False):
    """Retain every native branch, with explicit complete padding outcomes.

    Padding is failed preparation for the herald and an unbiased pointer bit
    for the read. The ideal prepared history stays inside the active corner.
    No source leakage is discarded by compressing back to 66 dimensions.
    """
    embedding = compiled.embedding()
    outside = np.eye(512)-embedding@embedding.T
    groups = [tuple(embedding@k@embedding.T for k in group) for group in instrument.kraus]
    if herald:
        groups[0] += (outside,)
    else:
        groups = [group+(outside/math.sqrt(2),) for group in groups]
    return runtime.Instrument(kraus=tuple(groups), visible=True, resources=instrument.resources,
        binding={'native_binding': instrument.binding, 'padding_outcomes_retained': True,
                 'padding_rule': 'failed herald' if herald else 'unbiased pointer',
                 'native_maps_are_direct_finite_repair_instruments': True})


def program(law, setting, *, uniform_reference=False):
    law, setting = quantum.candidate(law), quarter(setting)
    if type(uniform_reference) is not bool:
        raise ValueError('Boolean reference control required')
    herald, pointer = interface.herald_and_pointer(law)
    herald, pointer = _lift_native(herald, herald=True), _lift_native(pointer)
    identity = runtime.Instrument(kraus=((np.eye(512),),))
    mixed = runtime.Instrument(kraus=((compiled.full_matrix(),),),
        resources=(('compiled_mixer_pulses', compiled.resources()['native_pulses']),),
        binding={'source_word': 'compiled_native_interface.mixer', 'polar_projection_used': False})
    selected = 0 if uniform_reference else setting
    modular = runtime.Instrument(kraus=((full_matrix(selected),),),
        resources=(('compiled_modular_pulses', compiled.compiler.resources(word(selected))['native_pulses']),),
        binding={'quarter_turns': selected, 'programmed_active_corner_modular_probe': True,
                 'natural_native_waiting_evolution': False, 'polar_projection_used': False})
    failure_read = runtime.Instrument(kraus=((np.eye(512),), (np.zeros((512, 512)),)), visible=True)
    reverse = runtime.Instrument(kraus=tuple(reversed(pointer.kraus)), visible=True,
                                 resources=pointer.resources, binding={'minus_herald_decoding': 'reverse pointer bit'})
    embedding = compiled.embedding()
    return runtime.Program(f'native_modular_reference_{law}_{setting}', 512, embedding@embedding.T/66,
        (runtime.Event('native_signed_current_herald', instruments=(herald,)),
         runtime.Event('prepare_phase', instruments=(identity, mixed, mixed), control=0),
         runtime.Event('program_modular_phase', instruments=(identity, modular, modular), control=0),
         runtime.Event('analyze_phase', instruments=(identity, mixed, mixed), control=0),
         runtime.Event('native_coin_current_pointer', instruments=(failure_read, pointer, reverse), control=0)),
        binding={'preparation_state': 'I66/66, not Palm',
                 'reference_state': 'I66/66' if uniform_reference else 'diag(degree/600)',
                 'parameter': 'modular parameter s; programmed phase theta=s log(5/4)',
                 'source_reference_has_no_supplied_Hamiltonian_or_Gibbs_fit': True})


def oracle(law, setting):
    """Independent ideal active-pair formula, including terminal native maps."""
    degree, _, pairs, _ = reference()
    embedding = compiled.embedding()
    outputs = {(0, 0): np.zeros((512, 512), complex), (0, 1): np.zeros((512, 512), complex)}
    # Integer native paths give the failed-herald diagonal directly.
    counts = [0]*66
    for branch in interface.fixed_maps(law):
        if branch.current == 0:
            for _, output in branch.mapping:
                counts[output] += 1
    outputs[0, 0] = embedding@np.diag([float(F(n, 132)) for n in counts])@embedding.T
    _, pointer = interface.herald_and_pointer(law)
    for flag, sign in ((1, 1), (2, -1)):
        analyzed = np.zeros((66, 66), complex)
        for left, right in pairs:
            u = 1j**setting if degree[left] == 8 else 1
            v = 1j**setting if degree[right] == 8 else 1
            vector = np.zeros(66, complex)
            vector[left], vector[right] = (u+sign*v)/2, (u-sign*v)/2
            analyzed += np.outer(vector, vector.conj())/132
        for bit in (0, 1):
            raw = bit if flag == 1 else 1-bit
            density = sum((k@analyzed@k.conj().T for k in pointer.kraus[raw]), np.zeros((66, 66), complex))
            outputs[flag, bit] = embedding@density@embedding.T
    return outputs


def factor_execute(p):
    """Execute the same complete Program as Kraus factors, without projection.

    rho=w sum a a* starts with all 66 columns of E and the original binary
    input weight w=1/66 (no rounded square root of that weight). Every hidden
    Kraus alternative becomes another column; no branch or small amplitude
    is dropped. This avoids subtractive cancellation in nominally zero read
    outcomes of dense K rho K* propagation. Arithmetic remains diagnostic.
    """
    p = runtime.validate_program(p)
    if p.work_dimension != 512:
        raise ValueError('factor execution requires the declared 512-dimensional carrier')
    embedding = compiled.embedding()
    first = compiled.circuit.encode(quantum.STATES[0])
    weight = float(p.initial_density[first, first].real)
    if not np.array_equal(p.initial_density, weight*(embedding@embedding.T)):
        raise ValueError('factor execution requires exactly the declared binary uniform native input')
    branches = [((), embedding, p.preparation_resources)]
    for event in p.events:
        if event.sequence or event.repeat != 1:
            raise ValueError('this finite factor executor supports the declared single-step events')
        following = []
        for record, amplitudes, costs in branches:
            instrument = event.instruments[0 if event.control is None else record[event.control]]
            if instrument.schur is not None:
                raise ValueError('factor execution requires explicit Kraus instruments')
            for outcome, group in enumerate(instrument.kraus):
                factors = np.concatenate([k@amplitudes for k in group], axis=1)
                next_record = (*record, outcome) if event.visible else record
                following.append((next_record, factors, runtime.add_resources(costs, instrument.resources)))
        branches = following
    return tuple(runtime.Branch(record, weight*(amplitudes@amplitudes.conj().T), costs)
                 for record, amplitudes, costs in branches)


def execute_setting(law, setting, *, uniform_reference=False, checkpoint=False):
    p = program(law, setting, uniform_reference=uniform_reference)
    branches = factor_execute(p)
    selected = 0 if uniform_reference else setting
    match = probabilities(selected)
    success = interface.SUCCESS
    exact_probabilities = (1-2*success, F(0), success*match, success*(1-match), success*match, success*(1-match))
    if tuple(branch.record for branch in branches) != ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)):
        raise ArithmeticError('a complete raw herald/read outcome was lost')
    probability_error = max(abs(branch.probability-float(expected))
                            for branch, expected in zip(branches, exact_probabilities, strict=True))
    expected = oracle(law, selected)
    state_error = max(float(np.linalg.norm(branch.density-expected[branch.record], 'fro')) for branch in branches)
    for branch in branches:
        costs = dict(branch.resources)
        if branch.record[0] == 0 and any(costs.get(key, 0) for key in ('compiled_mixer_pulses', 'compiled_modular_pulses')):
            raise ArithmeticError('failed preparation acquired a phase-control cost')
    checkpoint_error = None
    generic_error = None
    if checkpoint:
        history = runtime.execute(p)
        generic_error = interface.bounded(max(float(np.linalg.norm(a.density-b.density, 'fro'))
            for a, b in zip(branches, history.branches, strict=True)), 'factor versus generic execution')
        prefix = runtime.execute(p, stop_after=2)
        resumed = runtime.execute(p, checkpoint=runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix)))
        checkpoint_error = interface.bounded(runtime.complete_distance(history, resumed), 'modular checkpoint')
    return {'law': law, 'quarter_turns': setting, 'uniform_reference_control': uniform_reference,
        'exact_conditional_match_probability': str(match), 'exact_joint_probabilities': list(map(str, exact_probabilities)),
        'executed_joint_probabilities': [branch.probability for branch in branches],
        'complete_probability_residual_diagnostic': interface.bounded(probability_error, 'modular readout'),
        'full_512_terminal_density_residual_diagnostic': interface.bounded(state_error, 'modular terminal state'),
        'checkpoint_complete_distance_diagnostic': checkpoint_error,
        'generic_runtime_factor_distance_diagnostic': generic_error,
        'execution_backend': 'complete Kraus-factor propagation through the declared Program',
        'realized_branch_resources': [{'record': list(branch.record), 'resources': dict(branch.resources)}
                                      for branch in branches]}


def run():
    degree, palm, pairs, mask = reference()
    actual_mixer = compiled.certify()
    report._gate(actual_mixer['operator_error_upper'], compiled.ERROR_GATE, 'raw mixer source word')
    mixer = {'half_diamond_upper': compiled.ERROR_GATE*(2+compiled.ERROR_GATE)/2,
             'channel_norm_upper': (1+compiled.ERROR_GATE)**2}
    sources = []
    for setting in QUARTERS:
        bounds = certify(setting)
        report._gate(bounds['operator_error_upper'], PHASE_ERROR_GATE, 'raw modular source word')
        # Fixed validated envelopes are portable across BLAS last-bit changes.
        portable = {'operator_error_upper': PHASE_ERROR_GATE,
                    'half_diamond_upper': PHASE_ERROR_GATE*(2+PHASE_ERROR_GATE)/2,
                    'channel_norm_upper': (1+PHASE_ERROR_GATE)**2}
        sources.append({'quarter_turns': setting, 'phase': f'{setting} pi/2',
            'operator_error_upper': compiled._outward(portable['operator_error_upper']),
            'half_diamond_upper': compiled._outward(portable['half_diamond_upper']),
            'channel_norm_upper': compiled._outward(portable['channel_norm_upper']),
            'mixer_phase_mixer_half_diamond_upper': compiled._outward(composite_bound(portable, mixer)),
            'exact_target_entries': 'rational complex, with diagonal values 1 or i^quarter_turns',
            'source_word_sha256': hashlib.sha256(canonical(word(setting)).encode()).hexdigest(),
            'resources': compiled.compiler.resources(word(setting))})
    rows = [execute_setting(law, setting, checkpoint=(law == 'C' and setting == 1))
            for law in ('C', 'E') for setting in QUARTERS]
    control = execute_setting('C', 2, uniform_reference=True)
    # The native unread channel has a different contraction law on this probe.
    native_control = {law: str((1+quantum.RATES[law])/2) for law in ('C', 'E')}
    return {'schema': 'oph.native_modular_reference.v1', 'all_pass': True,
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'reference': {'Palm_configuration_probabilities': list(map(str, palm)),
            'faithful': True, 'eigenvalue_ratio': '5/4', 'state_derived_modular_frequency': 'log(5/4)',
            'frequency_float_diagnostic': math.log(5/4), 'active_corner_reference_weight': '23/75',
            'ordered_pair_degree_census': {'8,10': 2, '10,8': 2, '8,8': 2, '10,10': 4},
            'eight_degree_8_active_code_addresses': list(mask), 'reference_is_probe_preparation': False},
        'readout': {'conditional_match_law': '4/5 + cos(s log(5/4))/5',
            'exact_law_parameter': 'modular parameter s', 'programmed_phase': 'theta=s log(5/4)',
            'probability_inversion_window': '0 <= s log(5/4) <= pi',
            'period': '2 pi/log(5/4)', 'phase_sign_and_period_aliases': True,
            'single_outcome_is_not_an_exact_clock_reading': True,
            'native_herald_useful_probability': '5/33', 'failed_preparation_probability': '28/33'},
        'compiled_modular_words': sources,
        'compiled_mixer': {'half_diamond_upper': compiled._outward(mixer['half_diamond_upper']),
            'channel_norm_upper': compiled._outward(mixer['channel_norm_upper']),
            'uses_per_successful_branch': 2, 'half_diamond_composition': 'delta_B*(1+g_D*g_B)+delta_D*g_B'},
        'complete_executions': rows,
        'uniform_reference_control': control, 'one_unread_native_attempt_match_probabilities': native_control,
        'scope': {'state_reference_constructed_from_native_event_statistics': True,
            'native_reference_inferred_by_fitting_a_supplied_Hamiltonian': False,
            'full_128_phase_map_and_512_spectator_lift_certified': True,
            'modular_flow_identification_domain': 'known twenty-dimensional active corner only',
            'full_66_state_modular_flow_implemented': False,
            'raw_pulse_products_polar_normalized': False,
            'native_instruments_and_floating_state_propagation_numerically_certified': False,
            'programmed_modular_probe_is_natural_native_waiting_evolution': False,
            'physical_service_error_bound': None,
            'continuing_clock_and_effective_field_dynamics_identified': False,
            'new_empirical_observations': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native-reference modular readout replayed with all preparation failures.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result)+'\n')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
