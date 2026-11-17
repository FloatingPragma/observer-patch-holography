"""Join native event-time continuation to the Palm modular protocol algebra.

The actual stopped CP channel is not replaced by modular automorphisms.
The reference, invariant algebra, full finite KMS identity and continuing
record predictions are attached separately, with an operational obstruction
to identifying their two dynamics. Native finite instruments and declared
coherent/projection apparatus execute; exact first-event limits do not.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..clocked_quantum_readout import construction as prepared
from ..compiled_native_interface import construction as compiled
from ..event_circulation import construction as events
from ..native_interface import construction as interface
from ..native_modular_reference import construction as modular
from ..native_relational_readout import construction as relational
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical, compare_receipt
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*modular.SOURCE_PATHS, *relational.SOURCE_PATHS,
    *prepared.clock.SOURCE_PATHS,
    'model/clocked_quantum_readout/construction.py',
    'model/relative_modular_dynamics/__init__.py',
    'model/relative_modular_dynamics/__main__.py',
    'model/relative_modular_dynamics/construction.py',
    'Lean/QFT/StructuralInheritance.lean', 'Lean/Dynamics/PrivateInner.lean',
    'Lean/Time/TimeOrderLedger.lean'))))
PAIRS = interface.phase_pairs()
ALGEBRA = frozenset((*((i, i) for i in range(66)), *PAIRS, *((j, i) for i, j in PAIRS)))
PROBE = (64, 54)  # source pair with equal Palm weights, checked below
ZERO = (F(0), F(0))
CAP = 2


def z(real=0, imag=0):
    return F(real), F(imag)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conjugate(a):
    return a[0], -a[1]


def scale(a, value):
    return a[0]*value, a[1]*value


def sparse(state):
    if type(state) is not dict or any(type(key) is not tuple or len(key) != 2
            or any(type(i) is not int or not 0 <= i < 66 for i in key)
            or type(value) is not tuple or len(value) != 2
            or any(type(c) is not F for c in value) for key, value in state.items()):
        raise ValueError('native matrix units with exact rational-complex coefficients required')
    return {key: value for key, value in state.items() if value != ZERO}


def accumulate(state, key, value):
    state[key] = add(state.get(key, ZERO), value)


def trace(state):
    value = z()
    for (i, j), entry in sparse(state).items():
        if i == j:
            value = add(value, entry)
    if value[1]:
        raise ArithmeticError('a purported density has nonreal trace')
    return value[0]


def dense(state):
    result = np.zeros((66, 66), complex)
    for key, value in sparse(state).items():
        result[key] = complex(float(value[0]), float(value[1]))
    return result


def marked(law):
    return _marked(relational.quantum.candidate(law))


@lru_cache(maxsize=2)
def _marked(law):
    return relational.source(law)


def event_columns(law):
    """Every matrix-unit column of the existing convergent stopped CP map.

    T0 is a Schur multiplier W=count/60. Summing T1 T0^n gives each
    moving partial-injection contribution 1/(60-count), entrywise.
    Immutable tuples keep callers from altering the source-derived map.
    """
    return _event_columns(relational.quantum.candidate(law))


@lru_cache(maxsize=2)
def _event_columns(law):
    source = marked(law)
    result = [defaultdict(F) for _ in range(66*66)]
    for mapping in source.moving_maps:
        for i, a in mapping:
            for j, b in mapping:
                result[66*i+j][a, b] += F(1, 60-source.no_current_counts[i][j])
    return tuple(tuple((a, b, coefficient) for (a, b), coefficient in sorted(column.items()))
                 for column in result)


def apply(law, state, *, first_at=None, through=None, no_current=None):
    """Exact linear native event map, or a named finite prefix/tail."""
    if sum(value is not None for value in (first_at, through, no_current)) > 1:
        raise ValueError('choose at most one finite event or no-current prefix')
    for value in (first_at, through, no_current):
        if value is not None:
            relational.horizon(value)
    source = marked(relational.quantum.candidate(law))
    columns = event_columns(law)
    result = {}
    for (i, j), value in sparse(state).items():
        w = F(source.no_current_counts[i][j], 60)
        if no_current is not None:
            accumulate(result, (i, j), scale(value, w**no_current))
            continue
        factor = ((1-w)*w**(first_at-1) if first_at is not None
                  else 1-w**through if through is not None else F(1))
        for a, b, coefficient in columns[66*i+j]:
            accumulate(result, (a, b), scale(value, coefficient*factor))
    return sparse(result)


def native_attempt_apply(law, state):
    """Exact sum of both outcomes of one actual native repair attempt.

    The no-current branches are identity partial maps whose total coefficient
    on E_ij is count[i][j]/60. Each moving partial injection contributes
    1/60 to its mapped matrix unit. This is the native-attempt channel T0+T1,
    not the first-current stopped map and not modular evolution.
    """
    source = marked(law)
    result = {}
    for (i, j), value in sparse(state).items():
        accumulate(result, (i, j), scale(value, F(source.no_current_counts[i][j], 60)))
        for mapping in source.moving_maps:
            mapping = dict(mapping)
            if i in mapping and j in mapping:
                accumulate(result, (mapping[i], mapping[j]), scale(value, F(1, 60)))
    return sparse(result)


def native_attempt_clock_audit(law, attempts=8):
    """Test whether native repair generates the existing modular clock hand.

    The declared coherent plus preparation and all native branch matrices are
    real. Consequently every iterate of their actual CP channel is real, so
    its Y readout stays unbiased for every integer native-attempt count. The
    finite execution checks this boundary through ``attempts``; real-linearity
    of each native Kraus branch proves the same statement for every n.
    """
    relational.horizon(attempts, 256)
    source = marked(law)
    state = prepared_states(law)[1]
    if trace(state) != interface.SUCCESS or any(value[1] for value in state.values()):
        raise ArithmeticError('the existing successful coherent preparation is no longer real')
    # The source instrument consists exclusively of real partial-injection
    # matrices. Check the actual arrays used by the bounded runtime as well as
    # the exact integer-address map used for the all-n argument.
    if any(np.max(np.abs(matrix.imag)) != 0 for group in source.instrument.kraus for matrix in group):
        raise ArithmeticError('the committed native repair instrument gained complex amplitudes')
    actual = dense(state)
    rows = []
    maximum_residual = 0.0
    for n in range(attempts+1):
        plus = trace(y_read(state, 0))
        minus = trace(y_read(state, 1))
        total = trace(state)
        if plus+minus != total or plus != total/2 or minus != total/2:
            raise ArithmeticError('native repair generated a Y phase from a real preparation')
        rows.append({'native_attempts': n, 'successful_preparation_mass': str(total),
                     'plus_probability_mass': str(plus), 'minus_probability_mass': str(minus),
                     'conditional_plus_probability': str(plus/total)})
        if n == attempts:
            break
        actual_outcomes = runtime._apply(source.instrument, actual, 1)
        actual = sum(actual_outcomes, np.zeros((66, 66), dtype=complex))
        state = native_attempt_apply(law, state)
        maximum_residual = max(maximum_residual,
            float(np.linalg.norm(actual-dense(state), 'fro')))
    quarter_turn = math.pi/(2*math.log(5/4))
    modular_plus = modular_clock_response(law, quarter_turn)
    expected_modular_plus = float(F(1, 2)-(F(1, 140) if law == 'C' else F(1, 220)))
    if abs(modular_plus-expected_modular_plus) > 2e-16 or modular_plus == 0.5:
        raise ArithmeticError('the ideal modular phase boundary witness changed')
    return {'law': law, 'parameter': 'integer native repair-attempt count n',
        'initial_state': 'existing successful coherent plus branch, mass 5/66',
        'source_step': 'exact T0+T1 from native zero-current identity branches and moving partial injections',
        'native_kraus_matrices_are_real': True,
        'initial_density_is_exactly_real': True,
        'conditional_plus_probabilities': rows,
        'all_integer_attempt_counts_unbiased_by_real-linearity': True,
        'runtime_vs_exact_native_attempt_maximum_frobenius_residual': maximum_residual,
        'quarter_modular_turn_conditional_plus_probability': str(F(1, 2)-(
            F(1, 140) if law == 'C' else F(1, 220))),
        'native_repair_generates_modular_phase': False,
        'physical_duration_or_lab_clock_identified': False}


def port_factor_audit():
    """Identify the source-derived local support of each nonzero frequency.

    The twelve ports label fixed-total-two load configurations, not twelve
    independent qubit factors in the current 66-state observer. Restricting
    the Palm state to a triple carrying one of its existing cross-degree
    coherent pairs gives a faithful three-state conditional reference. This
    makes the modular frequency local to three port labels without asserting
    that the triple is an autonomous subsystem of the repair dynamics.
    """
    degree, palm, _, _ = events.embedded()
    states = relational.quantum.STATES
    rows = []
    for i, j in PAIRS:
        if degree[i] == degree[j]:
            continue
        left, right = states[i], states[j]
        ports = tuple(k for k in range(12) if left[k] or right[k])
        if len(ports) != 3:
            raise ArithmeticError('a coherent cross-degree mode escaped its three-port support')
        local = tuple(k for k, state in enumerate(states)
            if max(state) == 1 and all(state[p] == 0 for p in range(12) if p not in ports))
        if len(local) != 3 or i not in local or j not in local:
            raise ArithmeticError('the three-port binary two-load subspace changed')
        mass = sum((palm[k] for k in local), F(0))
        weights = tuple(palm[k]/mass for k in local)
        if (mass != F(13, 300) or min(weights) <= 0
                or max(weights)/min(weights) != F(5, 4)):
            raise ArithmeticError('the local Palm restriction lost faithfulness or modular ratio')
        mode = {(i, j): z(1)}
        if trace(mode) != 0:
            raise ArithmeticError('a coherence mode acquired a diagonal component')
        phased = modular_action(mode)
        actuals = [native_attempt_apply(law, mode) for law in ('C', 'E')]
        if any(any(value[1] for value in actual.values()) or any(key not in ALGEBRA for key in actual)
                for actual in actuals):
            raise ArithmeticError('a native attempt no longer preserves real protocol coherences')
        phase_value = phased[i, j]
        if abs(phase_value[0]) != 0 or abs(phase_value[1]) != 1:
            raise ArithmeticError('the local source modular quarter turn changed')
        rows.append({'coherent_pair_indices': [i, j], 'ports': list(ports),
            'local_binary_two-load_indices': list(local),
            'conditional_Palm_weights': [str(weight) for weight in weights],
            'Palm_mass_of_local_two-load_condition': str(mass),
            'faithful_on_three-state_restricted_algebra': True,
            'cross-pair_eigenvalue_ratio': '5/4',
            'source_modular_frequency': 'log(5/4)',
            'native_attempt_image_is_real_and_stays_in_protocol_algebra': True})
    if len(rows) != 4:
        raise ArithmeticError('the native phase preparation no longer has four cross-degree pairs')
    return {'carrier_representation': 'fixed-total-two load configurations; 66 binary basis states, not 12 independent qubits',
        'factor_scope': 'three-port binary subspace with both loads inside the selected triple; conditional reference only',
        'local_coherent_modes': rows,
        'both_native_laws_keep_every_tested_mode_real': True,
        'autonomous_three-port_factor_dynamics_derived': False,
        'native_repair_generates_any_local_modular_phase': False}


def diagonal(values):
    return {(i, i): z(value) for i, value in enumerate(values) if value}


def modular_action(state, quarter=1, *, full=True):
    """Dual state action sigma^(-is) rho sigma^(is), theta=s log(5/4).

    The KMS observable convention is the inverse action
    alpha_s(A)=sigma^(is) A sigma^(-is).
    """
    if type(quarter) is not int or not -4 <= quarter <= 4 or type(full) is not bool:
        raise ValueError('bounded exact quarter-turn and Boolean domain required')
    degree, _, _, _ = events.embedded()
    active = {i for pair in PAIRS for i in pair}
    powers = (z(1), z(0, 1), z(-1), z(0, -1))
    result = {}
    for (i, j), value in sparse(state).items():
        a = int(degree[i] == 8 and (full or i in active))
        b = int(degree[j] == 8 and (full or j in active))
        result[i, j] = mul(powers[(quarter*(a-b)) % 4], value)
    return sparse(result)


def modular_clock_response(law, modular_parameter):
    """Analytic Y-readout probability for the existing coherent preparation.

    This evaluates the exact one-frequency response derived in
    ``modular_clock_audit``.  It is a prediction conditional on applying the
    reference modular action as the existing ideal target service; it does
    not claim that native repair itself implements that action.
    """
    if (type(modular_parameter) not in (int, float)
            or not math.isfinite(modular_parameter)):
        raise ValueError('a finite real modular parameter is required')
    law = relational.quantum.candidate(law)
    amplitude = 1/140 if law == 'C' else 1/220
    return 0.5-amplitude*math.sin(float(modular_parameter)*math.log(5/4))


def modular_clock_audit(law):
    """Derive a local, periodic modular-phase readout from existing source data.

    The initial coherent native preparation has only modular frequencies
    0,+/-log(5/4). The stopped native channel and existing Y read are linear,
    so the response has one sine/cosine harmonic. Four exact quarter-turn
    evaluations fix its coefficients; the remaining claims are checked
    against that harmonic and the source support.
    """
    law = relational.quantum.candidate(law)
    degree, _, _, _ = events.embedded()
    start = prepared_states(law)[1]
    if any(key not in ALGEBRA for key in start):
        raise ArithmeticError('the coherent source preparation left the protocol algebra')
    frequencies = sorted({int(degree[i] == 8)-int(degree[j] == 8)
                          for (i, j), value in start.items() if value != ZERO})
    if frequencies != [-1, 0, 1]:
        raise ArithmeticError('the source coherent preparation has new modular frequencies')
    probabilities = {}
    for quarter in (-1, 0, 1, 2):
        prepared_phase = modular_action(start, quarter, full=False)
        if modular_action(start, quarter) != prepared_phase:
            raise ArithmeticError('the active target is not the full modular action on this state')
        stopped = apply(law, prepared_phase)
        if any(key not in ALGEBRA for key in stopped):
            raise ArithmeticError('the stopped source channel left the protocol algebra')
        probabilities[quarter] = trace(y_read(stopped, 0))/interface.SUCCESS
    center = (probabilities[0]+probabilities[2])/2
    cosine = (probabilities[0]-probabilities[2])/2
    sine = (probabilities[1]-probabilities[-1])/2
    expected_amplitude = F(1, 140) if law == 'C' else F(1, 220)
    if (center, cosine, sine) != (F(1, 2), F(0), -expected_amplitude):
        raise ArithmeticError('the exact modular phase-response coefficients changed')
    for quarter, actual in probabilities.items():
        expected = center+cosine*(0 if quarter % 2 else (-1 if quarter % 4 == 2 else 1))
        # The sine quadrature is exact on quarter turns; even turns have zero sine.
        if quarter % 4 == 1:
            expected += sine
        elif quarter % 4 == 3:
            expected -= sine
        if actual != expected:
            raise ArithmeticError('quarter-turn phase prediction disagrees with its derived harmonic')
    return {
        'law': law,
        'initial_state': 'existing fixed-seam coherent native preparation, successful plus branch',
        'phase_service': 'ideal target for the native Palm reference modular action on the protocol algebra',
        'readout': 'existing Y readout after the actual first-moving-current stopped channel, conditioned on its 5/33 herald',
        'input_modular_frequency_support_in_units_of_log_5_over_4': frequencies,
        'exact_quarter_turn_conditional_plus_probabilities': {
            str(q): str(probabilities[q]) for q in (-1, 0, 1, 2)},
        'derived_exact_response': (
            '1/2 - sin(s*log(5/4))/140' if law == 'C'
            else '1/2 - sin(s*log(5/4))/220'),
        'exact_sine_amplitude': str(expected_amplitude),
        'exact_cosine_coefficient': '0',
        'modular_frequency': 'log(5/4)',
        'modular_phase_period': '2*pi/log(5/4)',
        'local_probability_slope_at_s_zero': (
            '-log(5/4)/140' if law == 'C' else '-log(5/4)/220'),
        'locally_phase_sensitive_near_zero': True,
        'single_probability_globally_identifies_phase': False,
        'native_repair_generates_modular_parameter': False,
        'modular_parameter_to_native_attempt_conversion_derived': False,
        'continuous_modular_gate_and_detector_error_bounded': False,
    }


def modular_lab_scale_audit(law):
    """Show that the dimensionless modular parameter has no lab-time scale.

    For any positive candidate rate kappa, s=kappa*t gives the same internal
    response as a function of s, but changes its period and slope in seconds.
    The candidates are scale witnesses, not fitted clock rates.
    """
    law = relational.quantum.candidate(law)
    amplitude = F(1, 140) if law == 'C' else F(1, 220)
    base_quarter_time = math.pi/(2*math.log(5/4))
    candidates = (F(1), F(2), F(5, 3))
    rows = []
    expected = float(F(1, 2)-amplitude)
    for ds_per_second in candidates:
        lab_time = base_quarter_time/float(ds_per_second)
        recovered_parameter = float(ds_per_second)*lab_time
        actual = modular_clock_response(law, recovered_parameter)
        if abs(actual-expected) > 2e-16:
            raise ArithmeticError('time-rescaled modular quarter-turn response changed')
        rows.append({
            'illustrative_modular_parameter_rate_per_second': str(ds_per_second),
            'implied_lab_period_ratio_to_unit_rate': str(1/ds_per_second),
            'quarter_turn_lab_time_ratio_to_unit_rate': str(1/ds_per_second),
            'exact_conditional_plus_probability_at_corresponding_quarter_turn': str(F(1, 2)-amplitude),
            'floating_response_residual_diagnostic': abs(actual-expected),
        })
    if len({row['implied_lab_period_ratio_to_unit_rate'] for row in rows}) != len(rows):
        raise ArithmeticError('the modular lab-time rescaling witness lost distinct periods')
    return {
        'source_parameter': 'dimensionless modular parameter s',
        'source_frequency_per_unit_s': 'log(5/4)',
        'candidate_map_if_calibrated': 's=kappa*t with kappa in modular-parameter units per second',
        'same_internal_response_at_matched_s_for_every_positive_kappa': True,
        'candidate_rates_are_estimates': False,
        'illustrative_rescalings': rows,
        'modular_parameter_to_seconds_map_identified': False,
        'native_attempts_generate_s': False,
        'matched_public_lab_dataset_identified': False,
        'reason_no_direct_dataset_test_is_defined': (
            'the current source provides an ideal modular phase service but no physical process mapped to s'),
        'minimal_calibration_and_operational_link': (
            'identify a physical implementation of the modular phase operation and pair a nonzero delta-s '
            'with an independent lab-time interval; justify constant ds/dt over the tested range. '
            'Calibration alone does not supply a native generator for s'),
    }


def balanced(state):
    # Reuse the existing exact native block support check, including its
    # rejection of unaccounted irrational cross-corner coefficients.
    real = prepared._balanced_exact({key: value[0] for key, value in sparse(state).items()})
    imag = prepared._balanced_exact({key: value[1] for key, value in sparse(state).items()})
    return sparse({key: z(real.get(key, 0), imag.get(key, 0)) for key in real.keys() | imag.keys()})


def prepared_states(law):
    groups = ({}, {}, {})
    for branch in interface.fixed_maps(law):
        flag = {0: 0, -1: 1, 1: 2}[branch.current]
        for _, output in branch.mapping:
            accumulate(groups[flag], (output, output), z(F(1, 132)))
    return groups[0], balanced(groups[1]), balanced(groups[2])


def degree_read(state, outcome):
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('degree outcome 0 or 1 required')
    degree = events.embedded()[0]
    wanted = (8, 10)[outcome]
    return {key: value for key, value in sparse(state).items()
            if degree[key[0]] == degree[key[1]] == wanted}


def y_read(state, outcome):
    """Exact Lüders read on one pair within the invariant protocol algebra."""
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError('phase outcome 0 or 1 required')
    state = sparse(state)
    if any(key not in ALGEBRA for key in state):
        raise ValueError('the rational read formula requires the protocol algebra')
    a, b = PROBE
    sign = 1-2*outcome
    projector = {(a, a): z(F(1, 2)), (b, b): z(F(1, 2)),
                 (a, b): z(0, F(-sign, 2)), (b, a): z(0, F(sign, 2))}
    result = {key: scale(value, F(1, 2)) for key, value in state.items() if key[0] not in PROBE}
    for (i, j), value in state.items():
        if i in PROBE:
            for r in PROBE:
                for s in PROBE:
                    accumulate(result, (r, s), mul(mul(projector[r, i], value), conjugate(projector[s, j])))
    return sparse(result)


def native_targets():
    degree = events.embedded()[0]
    active = {i for pair in PAIRS for i in pair}
    d = np.diag([1j if degree[i] == 8 and i in active else 1 for i in range(66)])
    phase = runtime.Instrument(kraus=((d,),), resources=(('declared_modular_phase_calls', 1),),
        binding={'target_surrogate': True, 'source_comparison': 'native_modular_reference',
                 'ideal_target_equals_full_modular_state_action_on_protocol_algebra': True})
    projectors = tuple(np.diag([int(value == wanted) for value in degree]) for wanted in (8, 10))
    degree_instrument = runtime.Instrument(kraus=tuple((p,) for p in projectors), visible=True,
        resources=(('declared_degree_reads', 1),), binding={'supplied_projection_apparatus': True})
    a, b = PROBE
    projection = np.zeros((66, 66), complex); projection[a, a] = projection[b, b] = 1
    y = np.zeros((66, 66), complex); y[a, b], y[b, a] = -1j, 1j
    outside = np.eye(66)-projection
    reads = tuple((projection+sign*y)/2+outside/math.sqrt(2) for sign in (1, -1))
    phase_instrument = runtime.Instrument(kraus=tuple((k,) for k in reads), visible=True,
        resources=(('declared_phase_reads', 1),), binding={'supplied_phase_sensitive_apparatus': True})
    return phase, degree_instrument, phase_instrument


def finite_preparation(law, *, coherent):
    if type(coherent) is not bool:
        raise ValueError('Boolean preparation selection required')
    law = relational.quantum.candidate(law)
    if coherent:
        return runtime.execute(prepared.preparation_program(law))
    source = marked(law)
    program = runtime.Program('actual_native_Palm_preparation', 66, np.eye(66)/66,
        (runtime.Event('uniform_native_attempt_current', instruments=(source.instrument,)),),
        binding={'success': 'movement on this specified attempt; not waiting for first movement'})
    return runtime.execute(program)


def degree_records(cap):
    relational.horizon(cap)
    return ((0, 0, -1, 0, -1), (1, cap+1, -1, 0, -1))+tuple(
        (1, n, a, m, b) for n in range(1, cap+1) for a in (0, 1)
        for m, b in (*((m, b) for m in range(1, cap+1) for b in (0, 1)), (cap+1, -1)))


def degree_experiment(law, cap=CAP):
    """Two continuing first-event blocks; retain preparation and both tails."""
    relational.horizon(cap, 16)
    source = marked(law)
    history = finite_preparation(law, coherent=False)
    states = {branch.record[0]: branch.density for branch in history.branches}
    _, detector, _ = native_targets()
    outputs = {(0, 0, -1, 0, -1): states[0]}
    waiting = states[1]
    for n in range(1, cap+1):
        waiting, arrived = runtime._apply(source.instrument, waiting, 1)
        for a, conditional in enumerate(runtime._apply(detector, arrived, 1)):
            second = conditional  # actual retained postmeasurement state; no reset
            for m in range(1, cap+1):
                second, arrived_again = runtime._apply(source.instrument, second, 1)
                for b, density in enumerate(runtime._apply(detector, arrived_again, 1)):
                    outputs[1, n, a, m, b] = density
            outputs[1, n, a, cap+1, -1] = second
    outputs[1, cap+1, -1, 0, -1] = waiting
    degree, palm, _, _ = events.embedded()
    good = diagonal([F(5, 33)*p for p in palm])
    failed = diagonal([F(60-d, 3960) for d in degree])
    exact = {(0, 0, -1, 0, -1): failed,
             (1, cap+1, -1, 0, -1): apply(law, good, no_current=cap)}
    for n in range(1, cap+1):
        arrived = apply(law, good, first_at=n)
        for a in (0, 1):
            conditional = degree_read(arrived, a)
            for m in range(1, cap+1):
                again = apply(law, conditional, first_at=m)
                for b in (0, 1):
                    exact[1, n, a, m, b] = degree_read(again, b)
            exact[1, n, a, cap+1, -1] = apply(law, conditional, no_current=cap)
    if outputs.keys() != exact.keys():
        raise ArithmeticError('two-event experiment dropped a raw record')
    if set(outputs) != set(degree_records(cap)) or sum(map(trace, exact.values())) != 1:
        raise ArithmeticError('two-event experiment lost complete probability')
    residual = max(float(np.linalg.norm(outputs[key]-dense(exact[key]), 'fro')) for key in outputs)
    rows = []
    for key in degree_records(cap):
        flag, n, a, m, b = key
        attempts = 1+(min(n, cap) if flag else 0)+(min(m, cap) if a >= 0 else 0)
        rows.append({'record': list(key), 'probability': str(trace(exact[key])),
                     'native_attempts': attempts, 'degree_reads': int(a >= 0)+int(b >= 0)})
    return {'law': law, 'attempt_cap_per_event': cap, 'complete_records': rows,
        'full_native_terminal_state_dimension': 66,
        'full_state_residual_diagnostic': relational.bounded(residual, 'continuing degree output'),
        'execution': 'existing primitive runtime CP instruments on retained postmeasurement matrices',
        'reset_between_events': False}


def covariance_records(cap):
    relational.horizon(cap)
    return ((0, 0, -1),)+tuple((flag, n, bit) for flag in (1, 2)
        for n in range(1, cap+1) for bit in (0, 1))+tuple((flag, cap+1, -1) for flag in (1, 2))


def covariance_experiment(law, order, cap=4):
    if type(order) is not str or order not in ('phase_before_event', 'phase_after_event'):
        raise ValueError('declared intervention order required')
    relational.horizon(cap, 32)
    source = marked(law)
    phase, _, detector = native_targets()
    history = finite_preparation(law, coherent=True)
    initial = {b.record[0]: b.density for b in history.branches}
    exact_initial = prepared_states(law)
    actual, exact = {(0, 0, -1): initial[0]}, {(0, 0, -1): exact_initial[0]}
    for flag in (1, 2):
        waiting, start = initial[flag], exact_initial[flag]
        if order == 'phase_before_event':
            waiting = runtime._apply(phase, waiting, 1)[0]
            start = modular_action(start, full=False)
        for n in range(1, cap+1):
            waiting, arrived = runtime._apply(source.instrument, waiting, 1)
            exact_arrived = apply(law, start, first_at=n)
            if order == 'phase_after_event':
                arrived = runtime._apply(phase, arrived, 1)[0]
                exact_arrived = modular_action(exact_arrived, full=False)
            read = runtime._apply(detector, arrived, 1)
            for bit in (0, 1):
                raw = bit if flag == 1 else 1-bit
                actual[flag, n, bit] = read[raw]
                exact[flag, n, bit] = y_read(exact_arrived, raw)
        actual[flag, cap+1, -1] = waiting
        exact[flag, cap+1, -1] = apply(law, start, no_current=cap)
    keys = covariance_records(cap)
    if set(actual) != set(keys) or set(exact) != set(keys) or sum(map(trace, exact.values())) != 1:
        raise ArithmeticError('covariance probe dropped a preparation failure or timeout')
    residual = max(float(np.linalg.norm(actual[k]-dense(exact[k]), 'fro')) for k in keys)
    return {'law': law, 'intervention_order': order, 'attempt_cap': cap,
        'complete_records': [{'record': list(k), 'probability': str(trace(exact[k])),
                              'native_attempts': 1+(min(k[1], cap) if k[0] else 0),
                              'balanced_mixer_calls': int(k[0] != 0),
                              'phase_calls': int(k[0] != 0 and (order == 'phase_before_event' or k[2] >= 0)),
                              'phase_reads': int(k[2] >= 0)} for k in keys],
        'full_state_residual_diagnostic': relational.bounded(residual, 'ordered covariance probe'),
        'ideal_native_phase_service_with_separate_full_carrier_source_comparison': True}


def algebra_audit(law):
    degree, palm, _, _ = events.embedded()
    columns = event_columns(law)
    if len(ALGEBRA) != 86 or PROBE not in PAIRS or degree[PROBE[0]] != degree[PROBE[1]]:
        raise ArithmeticError('protocol algebra or modular-fixed phase probe changed')
    covariance_defects = 0
    for i in range(66):
        for j in range(66):
            for a, b, coefficient in columns[66*i+j]:
                if ((i, j) in ALGEBRA) != ((a, b) in ALGEBRA):
                    raise ArithmeticError('native stopped map or its adjoint escapes the protocol algebra')
                covariance_defects += int(F(degree[i], degree[j]) != F(degree[a], degree[b]))
    if not covariance_defects:
        raise ArithmeticError('native/modular covariance obstruction vanished')
    sigma = diagonal(palm)
    if apply(law, sigma) != sigma:
        raise ArithmeticError('native first-event channel does not preserve the Palm reference')
    attempt = apply(law, sigma, no_current=1)
    for key, value in apply(law, sigma, first_at=1).items():
        accumulate(attempt, key, value)
    attempt_tv = sum(abs(attempt.get((i, i), ZERO)[0]-palm[i]) for i in range(66))/2
    if attempt_tv != F(1, 150):
        raise ArithmeticError('attempt-time and event-time stationary references were conflated')
    # The fullM66 finite KMS boundary identity is bilinear. All potentially
    # nonzero matrix-unit pairs Eij,Eji reduce to this exact ratio identity;
    # all other pairs have zero trace on both sides.
    for i in range(66):
        for j in range(66):
            if palm[i]*(palm[j]/palm[i]) != palm[j]:
                raise ArithmeticError('finite modular KMS ratio identity failed')
            unit = {(i, j): z(1)}
            if modular_action(apply(law, unit, no_current=1)) != apply(law, modular_action(unit), no_current=1):
                raise ArithmeticError('native no-current map lost full modular covariance')
    for key in ALGEBRA:
        unit = {key: z(1)}
        if modular_action(unit) != modular_action(unit, full=False):
            raise ArithmeticError('compiled corner phase fails the full protocol-algebra action')
    active = {i for pair in PAIRS for i in pair}
    omitted = next(i for i in range(66) if i not in active and degree[i] == 8)
    other = next(i for i in range(66) if i not in active and degree[i] == 10)
    unit = {(omitted, other): z(1)}
    if modular_action(unit) == modular_action(unit, full=False):
        raise ArithmeticError('fullM66 omission control no longer distinguishes the two actions')
    # Instrument-level support checks: every native partial injection maps
    # this matrix-unit algebra into itself; all declared apparatus is block
    # diagonal in its ten2x2 and46scalar blocks.
    mappings = (*marked(law).moving_maps, *(branch.mapping for branch in interface.fixed_maps(law)))
    for mapping in mappings:
        mapping = dict(mapping)
        inverse = {a: i for i, a in mapping.items()}
        for direction in (mapping, inverse):
            for i, j in ALGEBRA:
                if i in direction and j in direction and (direction[i], direction[j]) not in ALGEBRA:
                    raise ArithmeticError('a native Kraus branch or its adjoint leaves the protocol algebra')
    for state in prepared_states(law):
        if any(key not in ALGEBRA for key in state):
            raise ArithmeticError('native preparation leaves the protocol algebra')
    for key in ALGEBRA:
        unit = {key: z(1)}
        for output in (balanced(unit), degree_read(unit, 0), degree_read(unit, 1), y_read(unit, 0), y_read(unit, 1)):
            if any(k not in ALGEBRA for k in output):
                raise ArithmeticError('a declared intervention leaves the protocol algebra')
    return {'law': law, 'matrix_unit_columns_checked': 4356,
        'protocol_algebra_dimension': 86, 'protocol_algebra': 'direct_sum(10 copies M2, 46 copies C)',
        'source_map_and_adjoint_preserve_protocol_algebra': True,
        'actual_native_and_declared_apparatus_preserve_protocol_algebra': True,
        'Palm_reference_stationary_for_first_event_channel': True,
        'Palm_reference_one_native_attempt_TV_change': str(attempt_tv),
        'nonzero_full_channel_modular_frequency_mixing_entries': covariance_defects,
        'no_current_Schur_map_commutes_with_full_modular_flow': True,
        'full_M66_KMS_nonzero_basis_pairs_checked': 4356,
        'KMS_sign_convention': 'alpha_s(A)=sigma^(is) A sigma^(-is); omega(A alpha_(-i)(B))=omega(B A)',
        'ideal_active_phase_target_equals_full_modular_state_action_on_protocol_algebra': True,
        'outside_protocol_algebra_counterexample_matrix_unit': [omitted, other],
        'source_modular_channel_covariance': False}


def predictions(law):
    _, palm, _, _ = events.embedded()
    sigma = diagonal(palm)
    first = apply(law, sigma)
    joint = [[trace(degree_read(apply(law, degree_read(first, a)), b)) for b in (0, 1)] for a in (0, 1)]
    if joint != [[F(1, 5), F(1, 5)], [F(1, 5), F(2, 5)]]:
        raise ArithmeticError('continuing stationary degree record law changed')
    start = prepared_states(law)[1]
    before = apply(law, modular_action(start))
    after = modular_action(apply(law, start))
    before_p, after_p = trace(y_read(before, 0))/interface.SUCCESS, trace(y_read(after, 0))/interface.SUCCESS
    expected = F(69, 140) if law == 'C' else F(109, 220)
    if before_p != expected or after_p != F(1, 2):
        raise ArithmeticError('heralded operational covariance witness changed')
    limits = []
    for cap in (2, 8, 32, 64, 128):
        finite_first = apply(law, sigma, through=cap)
        finite_joint = [[trace(degree_read(apply(law, degree_read(finite_first, a), through=cap), b))
                        for b in (0, 1)] for a in (0, 1)]
        unresolved = 1-sum(map(sum, finite_joint))
        bound = min(F(1), 2*F(13, 15)**cap)
        if not 0 <= unresolved <= bound or any(finite_joint[a][b] > joint[a][b] for a in (0, 1) for b in (0, 1)):
            raise ArithmeticError('complete finite two-event truncation failed its uniform tail bound')
        limits.append({'attempt_cap_per_event': cap, 'conditional_unresolved_mass': str(unresolved),
                       'unconditional_unresolved_mass': str(F(5, 33)*unresolved),
                       'uniform_two_event_record_TV_bound': str(bound)})
    return {'law': law, 'two_event_degree_joint_probabilities': [[str(v) for v in row] for row in joint],
        'degree8_marginal': '2/5', 'degree8_two_event_connected_correlation': '1/25',
        'reset_control_joint_probabilities': [['4/25', '6/25'], ['6/25', '9/25']],
        'phase_before_event_conditional_plus': str(before_p), 'phase_after_event_conditional_plus': str(after_p),
        'pooled_herald_joint_probability_defect': str(F(5, 33)*(after_p-before_p)),
        'herald_success_probability': '5/33', 'preparation_failure_probability': '28/33',
        'two_event_finite_cap_enclosures': limits}


def run():
    phase_bound = modular.certify(1)
    mixer_bound = compiled.certify()
    if phase_bound['operator_error_upper'] > modular.PHASE_ERROR_GATE or mixer_bound['operator_error_upper'] > compiled.ERROR_GATE:
        raise ArithmeticError('fresh full-carrier source comparison failed')
    results = {'schema': 'oph.relative_modular_dynamics.v1', 'all_pass': True,
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'algebra_audits': [algebra_audit(law) for law in ('C', 'E')],
        'modular_clock_readouts': [modular_clock_audit(law) for law in ('C', 'E')],
        'modular_lab_time_scale_identifiability': [modular_lab_scale_audit(law) for law in ('C', 'E')],
        'native_attempt_clock_audits': [native_attempt_clock_audit(law) for law in ('C', 'E')],
        'port_factor_audit': port_factor_audit(),
        'exact_predictions': [predictions(law) for law in ('C', 'E')],
        'continuing_degree_executions': [degree_experiment(law) for law in ('C', 'E')],
        'ordered_phase_executions': [covariance_experiment(law, order) for law in ('C', 'E')
                                     for order in ('phase_before_event', 'phase_after_event')],
        'fresh_source_gate_scope': {'phase_operator_error_upper': str(modular.PHASE_ERROR_GATE),
            'preparation_mixer_operator_error_upper': str(compiled.ERROR_GATE),
            'full_512_source_and_padding_comparisons_recomputed': True,
            'runtime_uses_native_target_surrogates': True,
            'detectors_native_matrix_production_and_floating_state_arithmetic_certified': False},
        'scope': {'reference_algebra': 'native M66 with stationary Palm sigma=diag(degree/600)',
            'continued_source_channel': 'actual first-nonzero-current stopped CP map',
            'parameter_of_continuation': 'recorded nonzero native events, retaining attempt waits',
            'modular_parameter': 's, theta=s log(5/4)',
            'transfer_Hilbert_space_reference_is_the_same_operational_state': False,
            'full_M66_modular_implementation_claimed': False,
            'dissipative_event_channel_identified_with_modular_automorphisms': False,
            'infinite_first_event_limit_executed_in_finite_budget': False,
            'physical_apparatus_errors_bounded': False, 'new_empirical_measurement': False}}
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native continuing records and modular protocol-algebra receipt replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
