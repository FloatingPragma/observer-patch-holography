"""Preparation, native-generated seam choices and continuing complete trials.

Coin/current scratch records are discarded at their stated interfaces. Every
public herald, pointer and scheduler-timeout record remains. The auxiliary
is traced out by the exact coin-instrument identity, not reset or replaced.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from model.closed_prediction_protocol import construction as sequential
from model.compiled_native_interface import construction as compiled
from model.native_equilibration import sweeps
from model.native_interface import construction as interface
from model.native_quantum_history import construction as quantum
from model.observer import native
from model.primitive_observer import runtime
from model.readout_controls import prospective
from . import coins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
PROBE_ATTEMPTS = 3
TRIALS = 3927
ROUND_CAP = 7
PREPARATION_ALLOWANCE = F(1, 1570800)
CERTIFIED_PREPARATION_BOUND = F(259, 10**9)
SOURCE_PATHS = tuple(sorted(set((*sequential.SOURCE_PATHS,
    'model/native_equilibration/sweeps.py', 'model/repair_scheduler/__init__.py',
    'model/repair_scheduler/__main__.py', 'model/repair_scheduler/coins.py',
    'model/repair_scheduler/construction.py'))))


def public_records(attempts=PROBE_ATTEMPTS):
    """(herald, stage, bit); stage 1..n is timeout, n+1 completed pointer."""
    if type(attempts) is not int or not 1 <= attempts <= 8:
        raise ValueError('one to eight declared probe attempts required')
    return ((0, 0, 0), (0, 0, 1),
            *((h, stage, 0) for h in (1, 2) for stage in range(1, attempts+1)),
            *((h, attempts+1, bit) for h in (1, 2) for bit in (0, 1)))


def _density(density, dimension):
    if (type(density) is not np.ndarray or density.dtype.kind not in 'fc'
            or density.ndim != 2 or density.shape[0] != density.shape[1]
            or density.shape[0] not in (dimension, 2*dimension, 3*dimension)
            or not np.isfinite(density).all()):
        raise ValueError('plain finite floating/complex source density with up to three reference levels required')
    trace = np.trace(density)
    if (np.max(abs(density-density.conj().T)) > 1e-10 or abs(trace.imag) > 1e-10
            or not -1e-10 <= trace.real <= 1+1e-10
            or np.linalg.eigvalsh(density).min() < -1e-10):
        raise ValueError('positive Hermitian subnormalized density required')
    return density


def uniform_native_step(density, law):
    """Reduced accepted selector/work channel, enumerating all source seams."""
    result = np.zeros_like(density, dtype=complex)
    for seam in quantum.mixing.seams():
        for coin in (0, 1):
            result += coins.apply_coin(density, law, coin, seam=seam)/30
    return result


def trial_outputs(density, law, *, attempts=PROBE_ATTEMPTS, rounds=ROUND_CAP):
    """Complete ideal-target experiment with native-generated seam choices.

    Auxiliary reduction is exact on arbitrary correlated inputs at every
    call. Numerical native work channels are enumerated, never substituted
    by a spectral propagator. The supplied balanced mixer remains explicit.
    """
    density = _density(density, 66)
    records = public_records(attempts)
    selector = coins.selector_certificate(rounds=rounds)
    timeout = float(F(selector['timeout_probability']))
    reference = density.shape[0]//66
    herald, pointer = interface.herald_and_pointer(law)
    branches = runtime._apply(herald, density, reference)
    result = {record: np.zeros_like(density, dtype=complex) for record in records}
    result[0, 0, 0] = branches[0]
    mixer = runtime.Instrument(kraus=((interface.mixer(),),))
    for h in (1, 2):
        rho, = runtime._apply(mixer, branches[h], reference)
        for attempt in range(1, attempts+1):
            # No work attempt runs after a selector timeout. This state is
            # the actual disturbed prefix, not a normalized/reset state.
            result[h, attempt, 0] = timeout*rho
            rho = (1-timeout)*uniform_native_step(rho, law)
        rho, = runtime._apply(mixer, rho, reference)
        measured = runtime._apply(pointer, rho, reference)
        for bit in (0, 1):
            result[h, attempts+1, bit] = measured[bit if h == 1 else 1-bit]
    if any(not np.isfinite(rho).all() for rho in result.values()):
        raise ArithmeticError('complete scheduled trial produced a nonfinite output')
    total = sum(np.trace(rho).real for rho in result.values())
    if not np.isfinite(total) or abs(total-np.trace(density).real) > 2e-10:
        raise ArithmeticError('complete scheduled trial lost probability mass')
    return result


def _distance(a, b):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ArithmeticError('nonfinite output in distance diagnostic')
    result = float(np.linalg.norm(a-b, 'nuc')/2)
    if not np.isfinite(result):
        raise ArithmeticError('nonfinite distance diagnostic')
    return result


def _binary_indices():
    full = native.states()
    return np.array([full.index(x) for x in quantum.STATES], dtype=int)


def _embed_binary(density):
    reference = density.shape[0]//66
    ix = (_binary_indices()[:, None]*reference+np.arange(reference)).reshape(-1)
    result = np.zeros((78*reference, 78*reference), complex)
    result[np.ix_(ix, ix)] = density
    return result


def _restrict_binary(density):
    reference = density.shape[0]//78
    ix = (_binary_indices()[:, None]*reference+np.arange(reference)).reshape(-1)
    result = density[np.ix_(ix, ix)]
    restored = _embed_binary(result)
    if np.any(restored != density):
        raise ArithmeticError('native absorption failed; restriction would discard actual output')
    return result


def prepare_full(density, law, count):
    """Apply actual sweeps, then reindex their exactly binary-supported output."""
    density = _density(density, 78)
    output = sweeps.execute(density, law=law, sweeps=count,
                            reference_dimension=density.shape[0]//78, full=True)
    return _restrict_binary(output)


def _ideal_record_law(law):
    old = prospective.exact_law(law, PROBE_ATTEMPTS)
    out = {r: F(0) for r in public_records()}
    out[0, 0, 0], out[0, 0, 1] = old[:2]
    for h in (1, 2):
        for bit in (0, 1):
            out[h, PROBE_ATTEMPTS+1, bit] = old[2*h+bit]
    return out


def classify(counts):
    """No timeout or impossible record can be silently selected away."""
    if (type(counts) is not dict
            or any(type(key) is not tuple or len(key) != 3
                   or any(type(item) is not int for item in key) for key in counts)
            or set(counts) != set(public_records())
            or any(type(n) is not int or n < 0 for n in counts.values())
            or sum(counts.values()) != TRIALS):
        raise ValueError('all twelve nonnegative integer record counts at the fixed sample size required')
    if counts[0, 0, 1]:
        return 'outside_both_declared_models'
    if any(counts[h, stage, 0] for h in (1, 2) for stage in range(1, PROBE_ATTEMPTS+1)):
        return 'scheduler_timeout_inconclusive'
    coarse = tuple(counts[0, 0, b] for b in (0, 1))+tuple(
        counts[h, PROBE_ATTEMPTS+1, b] for h in (1, 2) for b in (0, 1))
    return prospective.classify(coarse, prospective.exact_law('C', PROBE_ATTEMPTS),
                                 prospective.exact_law('E', PROBE_ATTEMPTS), TRIALS)


def _json_exact(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    return value


def execute_history(law, count):
    """Two trials on one continuing carrier, beginning outside the binary shell."""
    vector = np.zeros(78, complex)
    full = native.states()
    collision = next(i for i, x in enumerate(full) if max(x) == 2)
    binary = _binary_indices()[0]
    vector[collision], vector[binary] = .6, .8j
    initial = np.outer(vector, vector.conj())
    prepared = prepare_full(initial, law, count)
    prepare_error = _distance(prepared, np.eye(66)/66)
    first = trial_outputs(prepared, law)
    ideal = _ideal_record_law(law)
    rows = []
    first_tv = sum(abs(float(np.trace(rho).real)-float(ideal[record]))
                   for record, rho in first.items())/2
    total_tv, conditional_max = 0., 0.
    for record, rho in first.items():
        mass = float(np.trace(rho).real)
        second = trial_outputs(prepare_full(_embed_binary(rho), law, count), law)
        if mass > 0:
            conditional_max = max(conditional_max, sum(abs(float(np.trace(output).real)/mass-float(ideal[r]))
                                                       for r, output in second.items())/2)
        for r, output in second.items():
            probability = float(np.trace(output).real)
            expected = ideal[record]*ideal[r]
            total_tv += abs(probability-float(expected))/2
            rows.append({'record': [*record, *r], 'probability': probability,
                         'ideal_product_probability': str(expected)})
    delta = PREPARATION_ALLOWANCE+PROBE_ATTEMPTS*F(16)**(-ROUND_CAP)
    if (prepare_error > float(PREPARATION_ALLOWANCE)+2e-10
            or first_tv > float(delta)+2e-10 or total_tv > 2*float(delta)+2e-10
            or conditional_max > float(delta)+2e-10):
        raise ArithmeticError('continuing complete history violates its declared bound')
    return {'law': law, 'executed_trials': 2, 'initial_full_sector_collision_mass': '9/25',
            'preparation_distance_diagnostic': prepare_error,
            'first_trial_TV_diagnostic': first_tv, 'two_trial_TV_diagnostic': total_tv,
            'maximum_second_conditional_TV_diagnostic': conditional_max,
            'complete_two_trial_law': rows, 'work_replaced_between_trials': False,
            'auxiliary_handling': 'exact reduced instrument; no auxiliary reset',
            'source_balanced_mixers': 'ideal target; fresh full-code comparison retained separately'}


def reference_history(law, count):
    """Full-sector input entangled with an untouched qubit and all output states."""
    vector = np.zeros(156, complex)
    full = native.states()
    collision = next(i for i, x in enumerate(full) if max(x) == 2)
    vector[2*collision], vector[2*_binary_indices()[0]+1] = .6, .8j
    initial = np.outer(vector, vector.conj())
    marginal = np.trace(initial.reshape(78, 2, 78, 2), axis1=0, axis2=2)
    prepared = prepare_full(initial, law, count)
    target = np.kron(np.eye(66)/66, marginal)
    actual, reference = trial_outputs(prepared, law), trial_outputs(target, law)
    error = sum(_distance(actual[key], reference[key]) for key in actual)
    marginal_error = _distance(sum(np.trace(value.reshape(66, 2, 66, 2), axis1=0, axis2=2)
                                   for value in actual.values()), marginal)
    oracle_error = _distance(uniform_native_step(prepared, law), sequential.prepare(prepared, law, 1))
    if error > float(PREPARATION_ALLOWANCE)+2e-10 or max(marginal_error, oracle_error) > 2e-10:
        raise ArithmeticError('reference-complete source/history comparison failed')
    return {'law': law, 'reference_dimension': 2,
            'complete_output_distance_diagnostic': error,
            'untouched_reference_marginal_error_diagnostic': marginal_error,
            'independent_uniform_source_channel_error_diagnostic': oracle_error}


def _compact_history(result):
    """Store every probability once, with one shared row/column record order."""
    result = dict(result)
    rows = result.pop('complete_two_trial_law')
    order = public_records()
    lookup = {tuple(row['record']): row['probability'] for row in rows}
    if len(rows) != len(order)**2 or len(lookup) != len(order)**2:
        raise ArithmeticError('complete history lost or duplicated a joint record')
    result['public_record_order'] = [list(record) for record in order]
    result['joint_probabilities'] = [[lookup[(*first, *second)] for second in order] for first in order]
    ideal = _ideal_record_law(result['law'])
    result['ideal_single_trial_probabilities'] = [str(ideal[key]) for key in order]
    return result


def _compact_auxiliary(result):
    result = dict(result)
    rows = result.pop('all_accepted_round_address_records')
    lookup = {(row['round'], row['address']): row['probability'] for row in rows}
    if len(rows) != 30*ROUND_CAP or len(lookup) != 30*ROUND_CAP:
        raise ArithmeticError('auxiliary execution lost or duplicated a selector record')
    result['probabilities_by_round_then_address'] = [[lookup[r, e] for e in range(30)]
                                                  for r in range(1, ROUND_CAP+1)]
    return result


def run():
    preparation = sweeps.budget(PREPARATION_ALLOWANCE, law='all', full=True)
    if preparation['trace_distance_upper'] > CERTIFIED_PREPARATION_BOUND:
        raise ArithmeticError('native sweep certificate exceeds its fixed outward budget')
    # The fixed existing readout menu is retained for comparison, not
    # re-advertised as optimal after accounting for auxiliary repairs.
    count = preparation['sweeps']
    selector = coins.selector_certificate(rounds=ROUND_CAP)
    decision = prospective.trial_design(prospective.exact_law('C', PROBE_ATTEMPTS),
                                        prospective.exact_law('E', PROBE_ATTEMPTS), F(1, 200))
    if decision['fixed_total_trials'] != TRIALS:
        raise ArithmeticError('existing likelihood design changed')
    scheduler_error = TRIALS*PROBE_ATTEMPTS*F(selector['timeout_probability'])
    error = F(1, 200)+TRIALS*CERTIFIED_PREPARATION_BOUND+scheduler_error
    if error >= F(1, 100):
        raise ArithmeticError('finite source decision budget exceeds one percent')
    source = compiled.certify()
    if source['two_call_half_diamond_upper'] > F(1, 10**10):
        raise ArithmeticError('source mixer comparison exceeds existing envelope')
    resources = compiled.resources()
    return {'schema': 'oph.repair_scheduler.v1',
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'native_coin_effects': coins.effect_certificate(), 'selector': selector,
        'sweep_structure': _json_exact(sweeps.structural_audit()), 'preparation': _json_exact(preparation),
        'auxiliary_executions': [_compact_auxiliary(coins.auxiliary_execution(law)) for law in ('C', 'E')],
        'continuing_histories': [_compact_history(execute_history(law, count)) for law in ('C', 'E')],
        'reference_histories': [reference_history(law, count) for law in ('C', 'E')],
        'decision': {'fixed_trials': TRIALS, 'probe_attempts': PROBE_ATTEMPTS,
            'ideal_likelihood_certificate': decision,
            'preparation_error_upper': str(TRIALS*CERTIFIED_PREPARATION_BOUND),
            'scheduler_timeout_union_upper': str(scheduler_error),
            'wrong_candidate_probability_upper': str(error),
            'entire_sampled_campaign_executed': False,
            'new_optimality_or_new_frozen_prediction_claim': False},
        'resources': {'work_native_attempts_per_trial_upper': 30*count+PROBE_ATTEMPTS+2,
            'auxiliary_native_attempts_per_trial_upper': PROBE_ATTEMPTS*5*ROUND_CAP,
            'total_native_attempts_campaign_upper': TRIALS*(30*count+PROBE_ATTEMPTS+2+PROBE_ATTEMPTS*5*ROUND_CAP),
            'source_pulses_campaign_upper': 2*TRIALS*resources['native_pulses'],
            'proper_code_service_events_campaign_upper': 2*TRIALS*resources['proper_code_service_events'],
            'mixer_raw_array_two_call_half_diamond_upper': '1/10000000000',
            'physical_service_error_authenticated': False},
        'scope': {'external_IID_seam_sampler_required': False,
            'exact_stationary_or_binary_initial_preparation_required': False,
            'auxiliary_initial_independence_or_reset_required': False,
            'working_total_two_and_native_C_or_E_instrument_supplied': True,
            'additional_auxiliary_carrier_and_classical_controller_supplied': True,
            'primitive_fair_coin_and_coherent_services_derived_from_endpoint_rule': False,
            'all_paper_premises_or_physical_data_agreement_established': False,
            'source_clock_counts_are_physical_durations': False}}
