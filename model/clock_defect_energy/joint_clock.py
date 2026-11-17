"""A finite complete joint native-event clock and protected-field instrument.

The same committed incident count drives repeated transported parameter-four
source words and the existing coherent native arc walk.  The clock and field
generators need not agree.  Source arrays are executed; the whole arc walk,
its preparation, controlled branch and Ramsey read remain declared services.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np

from model.native_modular_reference import history_clock as clock
from model.observer import coherent_continuation as wave
from .joins_contract import binding, check_receipt, source_bindings
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "joint_clock_results.json"
FIELD = tuple(6*state[0]-1 for state in wave.STATES)
ZERO = (F(), F(), F())
COEFFICIENT_BASIS = ("1", "cos(4*log(8/5))", "cos(8*log(8/5))")


def _law(law):
    return clock.records.candidate(law)


def tick_count(law, history):
    """Zero-inclusive incidences in the predecessor's two-attempt policy."""
    if type(history) is not tuple or len(history) != 2:
        raise ValueError("two complete native records required")
    law = _law(law)
    clock.transport._history(law, history)
    if history[0][:2] != clock.FIRST_EDGE or history[1][:2] != clock.next_edge(law, history[0]):
        raise ValueError("history does not follow the committed predecessor policy")
    return sum(0 in record[:2] for record in history)


def field_probabilities(ticks, phase=-1):
    if type(ticks) is not int or ticks not in (1, 2):
        raise ValueError("one or two admitted incident ticks required")
    if type(phase) is not int or phase not in (-1, 1):
        raise ValueError("the declared integer continuation phase required")
    return _field_probabilities(ticks, phase)


@lru_cache(maxsize=4)
def _field_probabilities(ticks, phase):
    outputs = wave.ramsey_outputs(FIELD, ticks, phase)
    answer = tuple(wave.norm_squared(x)/330 for x in outputs)
    if sum(answer) != 1:
        raise ArithmeticError("a field Ramsey branch was lost")
    return answer


def _plus(a, b):
    return tuple(x+y for x, y in zip(a, b, strict=True))


def _scale(a, value):
    return tuple(x*value for x in a)


def exact_clock_rows(law, *, event_clock=True):
    """Complete attained native paths, with exact trigonometric coefficients.

    Input is the ideal faithful target, followed by the existing two-address
    probe rotation.  This does not replace the finite preparation channel.
    Unattained raw-register slots have exactly zero maps on every input.
    """
    law = _law(law)
    if type(event_clock) is not bool:
        raise ValueError("an explicit event-clock control is required")
    weights = clock.transport.reference()
    a, b = clock.probe_pair()
    diagonal = list(weights)
    diagonal[a] = diagonal[b] = (weights[a]+weights[b])/2
    coherence = (weights[b]-weights[a])/2
    rows = {}
    for branch in clock.branches(law):
        ticks = tick_count(law, branch.history)
        populations = {out: (diagonal[i]*branch.weight, F(), F())
                       for i, out in zip(branch.inputs, branch.outputs, strict=True)}
        if a in branch.inputs and b in branch.inputs:
            left, right = (branch.outputs[branch.inputs.index(i)] for i in (a, b))
            average = _scale(_plus(populations[left], populations[right]), F(1, 2))
            term = [F(), F(), F()]
            term[ticks if event_clock else 1] = coherence*branch.weight
            populations[left] = _plus(average, _scale(term, -1))
            populations[right] = _plus(average, term)
        for read in clock.interface.fixed_maps(law, clock.READ_EDGE):
            value = ZERO
            for i, _ in read.mapping:
                value = _plus(value, _scale(populations.get(i, ZERO), F(1, 2)))
            key = (branch.history, read.coin, read.label, read.current, ticks)
            rows[key] = value
    total = ZERO
    for value in rows.values():
        total = _plus(total, value)
    if total != (F(1), F(), F()):
        raise ArithmeticError("the exact complete clock law is not normalized")
    return rows


@lru_cache(maxsize=1)
def cosine_intervals():
    result = []
    for count in (1, 2):
        re, _, radius = clock.common.phase_enclosure(F(4*count), F(8, 5))
        result.append((max(F(-1), re-radius), min(F(1), re+radius)))
    return tuple(result)


def evaluate_interval(coefficients):
    lo = hi = coefficients[0]
    for coefficient, interval in zip(coefficients[1:], cosine_intervals(), strict=True):
        a, b = sorted(coefficient*x for x in interval)
        lo, hi = lo+a, hi+b
    return lo, hi


def _outward(interval, places=12):
    scale = 10**places
    a, b = interval
    return [str(F(math.floor(a*scale), scale)), str(F(math.ceil(b*scale), scale))]


def exact_joint(law):
    """A changed coupling with exactly the same two separate marginals."""
    rows = exact_clock_rows(law)
    ticks_mass = defaultdict(F)
    grouped = defaultdict(lambda: ZERO)
    for (history, coin, label, current, ticks), coefficients in rows.items():
        grouped[ticks, current] = _plus(grouped[ticks, current], coefficients)
        ticks_mass[ticks] += coefficients[0]
    # The clock phases do not change any pre-detector history probability.
    for ticks in (1, 2):
        total = ZERO
        for (n, _), coefficients in grouped.items():
            if n == ticks:
                total = _plus(total, coefficients)
        if total != (ticks_mass[ticks], F(), F()):
            raise ArithmeticError("the phase word changed the incident count")
    field_plus = sum((mass*field_probabilities(n)[0] for n, mass in ticks_mass.items()), F())
    independent_tv = sum((mass*abs(field_probabilities(n)[0]-field_plus)
                          for n, mass in ticks_mass.items()), F())
    negative_read = joint_negative_plus = ZERO
    retained = []
    for (ticks, current), coefficients in sorted(grouped.items()):
        plus = _scale(coefficients, field_probabilities(ticks)[0])
        minus = _scale(coefficients, field_probabilities(ticks)[1])
        if current == -1:
            negative_read = _plus(negative_read, coefficients)
            joint_negative_plus = _plus(joint_negative_plus, plus)
        retained.append({"incident_ticks": ticks, "clock_detector_current": current,
                         "clock_coefficients": list(map(str, coefficients)),
                         "joint_plus_coefficients": list(map(str, plus)),
                         "joint_minus_coefficients": list(map(str, minus)),
                         "joint_plus_interval": _outward(evaluate_interval(plus))})
    covariance = _plus(joint_negative_plus, _scale(negative_read, -field_plus))
    # Run only the first phase word while retaining the identical native/wave
    # schedule: this changes the clock detector and tests the event coupling.
    one_word = exact_clock_rows(law, event_clock=False)
    tv_lo = tv_hi = F()
    for key, value in rows.items():
        lo, hi = evaluate_interval(_plus(value, _scale(one_word[key], -1)))
        tv_lo += min(abs(lo), abs(hi)) if lo*hi > 0 else 0
        tv_hi += max(abs(lo), abs(hi))
    if independent_tv <= 0 or tv_lo <= 0 or evaluate_interval(covariance)[1] >= 0:
        raise ArithmeticError("a joint or triggered-clock control disappeared")
    return {"coefficient_basis": list(COEFFICIENT_BASIS),
            "conditional_target_input": "diag(transported local reference), probe rotation, no reset",
            "complete_attained_clock_read_records": len(rows),
            "incident_count_probabilities": {str(n): str(p) for n, p in sorted(ticks_mass.items())},
            "field_Ramsey_plus_marginal": str(field_plus),
            "grouped_joint_law": retained,
            "clock_negative_current_field_plus_covariance_coefficients": list(map(str, covariance)),
            "clock_negative_current_field_plus_covariance_interval": _outward(evaluate_interval(covariance)),
            "independent_field_resampling_control": {
                "complete_clock_record_marginal_preserved": True,
                "field_Ramsey_marginal_preserved": True,
                "joint_total_variation_exact": str(independent_tv),
                "complete_ideal_preparation_joint_TV_exact": str(F(5, 33)*independent_tv),
                "definition": "Replace P(Ramsey|complete clock record) by its unconditional field marginal."},
            "one_phase_word_control_complete_clock_TV_interval": _outward((tv_lo/2, tv_hi/2)),
            "full_record_coefficient_digest": _digest([
                [key, list(map(str, value))] for key, value in sorted(rows.items())])}


@dataclass(frozen=True)
class Frame:
    """An abort's padding output is traced, retaining every external reference."""
    key: tuple
    ticks: int
    abort: bool
    matrix: np.ndarray


def _matrix(branch, raw):
    if raw:
        return branch.matrix
    value = np.zeros((66, 66), complex)
    for i, j in branch.mapping:
        value[j, i] = 1/math.sqrt(2)
    return value


def clock_frames(law, *, raw=True, event_clock=True):
    """Execute complete source operators, including each off-code export.

    Successful read frame J_h is paired with (I +/- W^N)/2.  An abort
    frame is paired with the accumulated W^N and no Ramsey measurement.
    Tensoring these operators with I_R defines every untouched-reference
    extension, including inputs entangled between clock, field and R.
    """
    law = _law(law)
    if type(raw) is not bool or type(event_clock) is not bool:
        raise ValueError("explicit source and trigger Booleans required")
    codes = np.array(clock.full.support()[0])
    outside = np.array([i for i in range(512) if i not in set(codes)])
    result = []

    def export(array, incoming, key, ticks):
        result.append(Frame(key, ticks, True, array[np.ix_(outside, codes)]@incoming))
        return array[np.ix_(codes, codes)]@incoming

    probe = export(clock.rotation(clock.probe_pair(), 1, raw), np.eye(66, dtype=complex),
                   ("abort", "probe"), 0)
    first_branches = {b.history[0]: b for b in
                      clock.transport.word_instrument(law, (clock.FIRST_EDGE,))}
    admitted = {b.history: b for b in clock.branches(law)}
    for first in clock.interface.fixed_maps(law, clock.FIRST_EDGE):
        h1 = (*clock.FIRST_EDGE, first.coin, first.label)
        middle = _matrix(first, raw)@probe
        middle = export(clock.clock(F(4), clock.ratios(first_branches[h1], law=law), raw),
                        middle, ("abort", "clock1", h1), 0)
        edge = clock.next_edge(law, h1)
        for second in clock.interface.fixed_maps(law, edge):
            history = (h1, (*edge, second.coin, second.label))
            ticks = tick_count(law, history)
            branch = admitted.get(history)
            final = _matrix(second, raw)@middle
            if branch is not None and ticks == 2 and event_clock:
                final = export(clock.clock(F(4), clock.ratios(branch, law=law), raw),
                               final, ("abort", "clock2", history), 1)
            else:
                result.append(Frame(("abort", "clock2", history), 1, True,
                                    np.zeros((446, 66), complex)))
            if branch is not None and all(i in branch.inputs for i in clock.probe_pair()):
                pair = tuple(branch.outputs[branch.inputs.index(i)] for i in clock.probe_pair())
                final = export(clock.rotation(pair, -1, raw), final,
                               ("abort", "analyzer", history), ticks)
            else:
                result.append(Frame(("abort", "analyzer", history), ticks, True,
                                    np.zeros((446, 66), complex)))
            for read in clock.interface.fixed_maps(law, clock.READ_EDGE):
                result.append(Frame(("read", history, read.coin, read.label, read.current),
                                    ticks, False, _matrix(read, raw)@final))
    if len({frame.key for frame in result}) != len(result):
        raise ArithmeticError("complete record collision")
    return tuple(result)


def factor_probabilities(frames, density):
    """All joint probabilities for a clock input and the declared field input."""
    if (type(density) is not np.ndarray or density.shape != (66, 66)
            or not np.isfinite(density).all()
            or not np.allclose(density, density.conj().T, atol=1e-13, rtol=0)
            or np.linalg.eigvalsh(density).min() < -1e-13):
        raise ValueError("a positive original finite native input required")
    result = {}
    for frame in frames:
        probability = float(np.einsum('ai,ai->', frame.matrix@density, frame.matrix.conj()).real)
        if frame.abort:
            result[frame.key] = probability
        else:
            for sign, weight in zip((1, -1), field_probabilities(frame.ticks), strict=True):
                result[(*frame.key, sign)] = probability*float(weight)
    return result


@lru_cache(maxsize=32)
def _field_gram(vertices, ticks, sign):
    if sign is None:
        vectors = tuple(wave.embed(v) for v in vertices)
    else:
        vectors = tuple(wave.ramsey_outputs(v, ticks)[int(sign == -1)] for v in vertices)
    return tuple(tuple(sum((a*b for a, b in zip(v, w, strict=True)), F())/60
                       for w in vectors) for v in vectors)


def schmidt_reference(frame, source_columns, vertices, sign=None):
    """Actual reduced untouched-reference output for a finite Schmidt input.

    Input is sum_k C[:, :, k] tensor A|v_k>, in source, reference, field
    order. The coherent source/field output is defined by the same factors;
    only this requested reference marginal is traced here. No branch is
    normalized. An abort keeps the walk's live arm and its reference arm;
    their common Gram matrix is unchanged by unitary accumulation.
    """
    if (type(frame) is not Frame or type(frame.ticks) is not int or not 0 <= frame.ticks <= 2
            or type(frame.abort) is not bool or type(frame.matrix) is not np.ndarray
            or frame.matrix.ndim != 2 or frame.matrix.shape[1] != 66
            or not np.isfinite(frame.matrix).all()
            or type(source_columns) is not np.ndarray or source_columns.ndim != 3
            or source_columns.shape[0] != 66 or not np.isfinite(source_columns).all()
            or min(source_columns.shape[1:]) < 1
            or type(vertices) is not tuple or len(vertices) != source_columns.shape[2]
            or (frame.abort and sign is not None)
            or (not frame.abort and (type(sign) is not int or sign not in (-1, 1)))):
        raise ValueError("finite Schmidt columns and the declared complete output required")
    vertices = tuple(wave._vector(v, 66) for v in vertices)
    gram = np.array(_field_gram(vertices, frame.ticks, sign), dtype=float)
    transformed = np.einsum("ab,bri->ari", frame.matrix, source_columns)
    return np.einsum("ari,asj,ji->rs", transformed, transformed.conj(), gram)


def _digest(value):
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def source_budget():
    predecessor = clock.source_budget()  # executes and certifies the actual words
    error, gain = F(), F(1)
    for eps in (clock.ROTATION_GATE, clock.NATIVE_GATE, clock.WORD_GATE,
                clock.NATIVE_GATE, clock.WORD_GATE, clock.ROTATION_GATE, clock.NATIVE_GATE):
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    word_size = len(clock.common.diagonal_word(F(4), clock.ratios()))
    pulses = predecessor['maximum_source_pulses_per_successful_history']+word_size
    prep = clock.preparation.budget()
    joined = gain*F(prep['complete_flagged_reference_halftrace_upper'])+F(5, 33)*error
    return {"source_halfdiamond_upper": str(error), "source_channel_norm_upper": str(gain),
            "complete_preparation_and_joint_halftrace_upper": str(joined),
            "initial_prepared_reference_parameter8_stability_only": clock.transport.preparation_certificate(F(8)),
            "unknown_prepared_reference_history_covariance_claimed": False,
            "maximum_clock_word_calls": 2, "clock_parameter_per_incident": "4",
            "maximum_source_pulses": pulses, "maximum_coherent_load_export_calls": 2*pulses,
            "native_attempts_after_preparation_including_detector": 3,
            "preparation_native_attempts_upper": prep['total_native_attempts_upper'],
            "clock_address_qubits": 9, "protected_wave_arc_coordinates": wave.DIMENSION,
            "wave_address_qubits_upper": wave.DIMENSION.bit_length(), "Ramsey_control_qubits": 1,
            "maximum_controlled_wave_steps": 2, "field_preparation_calls": 1, "field_readout_calls": 1,
            "complete_native_record_bits_upper": 27, "transported_core_bits_upper": 8,
            "incident_counter_bits": 2, "terminal_kind_and_location_bits_upper": 4,
            "source_program_bytes": predecessor['source_program_utf8_bytes'],
            "controller_policy_bytes": predecessor['controller_policy_utf8_bytes'],
            "extra_event_controller_rule": "after successful port0-incident clock export, advance protected live wave arm once",
            "whole_wave_source_compiler_executed": False,
            "joint_halfdiamond_with_wave_service_errors":
                "source_bound + source_gain*(epsilon_field_prepare + 2*epsilon_controlled_W + epsilon_Ramsey_read)",
            "wave_errors_are_declared_CPTP_halfdiamond_allowances": True,
            "physical_apparatus_error_bound": None,
            "host_workspace_bounded": False}


def run():
    budget = source_budget()
    rows = []
    reference = np.diag(list(map(float, clock.transport.reference()))).astype(complex)
    for law in ('C', 'E'):
        ideal = clock_frames(law, raw=False)
        actual = clock_frames(law, raw=True)
        ideal_probs = factor_probabilities(ideal, reference)
        actual_probs = factor_probabilities(actual, reference)
        variation = sum(abs(actual_probs[key]-value) for key, value in ideal_probs.items())/2
        effect = sum((frame.matrix.conj().T@frame.matrix for frame in ideal), np.zeros((66, 66), complex))
        if (np.linalg.norm(effect-np.eye(66), ord=2) > 2e-12
                or variation > 2e-10 or abs(sum(actual_probs.values())-1) > 2e-10):
            raise ArithmeticError("complete executed joint source instrument failed")
        exact = exact_clock_rows(law)
        for frame in ideal:
            if not frame.abort:
                key = (*frame.key[1:], frame.ticks)
                value = exact.get(key, ZERO)
                lo, hi = evaluate_interval(value)
                probability = float(np.trace(frame.matrix@reference@frame.matrix.conj().T).real)
                if not float(lo)-2e-12 <= probability <= float(hi)+2e-12:
                    raise ArithmeticError("executed clock law differs from exact coefficients")
        preparation = clock.preparation.population_maps(law=law)
        finite = factor_probabilities(actual, np.diag(preparation['success'][:, 0]).astype(complex))
        failures = {flag: float(np.sum(matrix[:, 0])) for flag, matrix in preparation.items() if flag != 'success'}
        if abs(sum(finite.values())+sum(failures.values())-1) > 2e-9:
            raise ArithmeticError("finite preparation flags were dropped")
        prepared_tv = sum(abs(finite[key]-float(F(5, 33))*value)
                          for key, value in ideal_probs.items())/2
        prepared_tv += sum(abs(value-(float(F(28, 33)) if flag == 'herald_failure' else 0))
                           for flag, value in failures.items())/2
        if prepared_tv > float(F(budget['complete_preparation_and_joint_halftrace_upper']))+2e-10:
            raise ArithmeticError("executed preparation plus joint record law exceeds its source bound")
        rows.append({"law": law, **exact_joint(law),
                     "complete_source_export_and_read_frames": len(actual),
                     "complete_joint_read_and_abort_records": len(actual_probs),
                     "raw_source_against_ideal_joint_TV_interval": ["0", "1/5000000000"],
                     "finite_preparation_flags": sorted(preparation),
                     "finite_preparation_and_joint_normalization_error_upper": "1/500000000",
                     "finite_preparation_and_joint_TV_within_complete_bound": True,
                     "raw_register_zero_extension": {
                         "per_attempt_coin_and_raw_slots": len(clock.raw_menu(law)),
                         "three_attempt_raw_words": len(clock.raw_menu(law))**3,
                         "rule": "history_clock.raw_read_key; unattained paths remain zero on every reference"}})
    paths = set(clock.transport.source_bindings()) | set(wave.SOURCE_PATHS)
    paths.update({"model/clock_defect_energy/joins_contract.py", "model/clock_defect_energy/contract.py",
                  "model/clock_defect_energy/joint_clock.py", "model/clock_defect_energy/test_joint_clock.py",
                  "model/clock_defect_energy/JOINT_CLOCK.md", "model/native_modular_reference/history_clock.py",
                  "model/native_modular_reference/history_clock_results.json.gz", "tools/strict_json.py",
                  "model/relative_observer_geometry/temporal_records.py", "model/common_source/THERMAL_CLOCK.md",
                  "Lean/Time/ClockReadout.lean", "Lean/Time/ClockComparison.lean"})
    return {"schema": "oph.clock_defect_energy.joint_clock.v1", **binding(),
            "positive_result": "Complete C/E records jointly drive repeated transported clock words and a protected field walk.",
            "field_Ramsey_by_incident_count": {str(n): list(map(str, field_probabilities(n))) for n in (1, 2)},
            "source_and_service_budget": budget, "laws": rows,
            "all_input_statement": "Each successful map is J_h tensor (I +/- W^N)/2 tensor I_R; both signs and all aborts retained.",
            "preparation_scope": "Native clock preparation uses all four existing flags; field preparation is a declared coherent service.",
            "clock_readout_scope": "Known finite N in {1,2}; detector distributions differ, no one-shot or global phase unwrapping.",
            "physical_seconds_calibrated": False, "equal_generators_required": False,
            "single_native_clock_field_primitive_derived": False, "new_empirical_predictions": 0,
            "source_sha256": bound_pins(source_bindings(paths))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    parser.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        check_receipt(args.check, result)
    if not args.write and not args.check:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print('Complete joint clock/field instrument and matched-marginal control verified.')


if __name__ == '__main__':
    main()
