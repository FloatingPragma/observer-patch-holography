"""Continue the source-prepared edge reference with a clock and feedback.

The complete evidence keeps the sampler request, every first herald failure,
both native seam/coin/current records, all sixteen detector outcomes and a
twelve-dimensional untouched reference. Modular identification concerns the
successful work marginal, not the larger classical-record/reference algebra.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from tools import strict_json
from ..observer import native, quantum_mixing
from ..observer.source_scalar import compiler
from ..physical_space_time import local_continuum as spatial_limit
from ..physical_space_time import source_green as spatial_source
from . import inverse_clock as preparation
from . import prepared_clock_source as source

ROOT = Path(__file__).resolve().parents[2]
INITIAL = (0,)*10+(1, 1)
CLOCK_PORT = 10
BAND = 2  # One vector of the lambda=6 sector, not its multiplicity.
FILTER_GUARD = F(1, 10**9)
ARITHMETIC_ALLOWANCE = F(1, 10**10)
WEYL_FACTORS = ((3, "X", F(1, 10)), (3, "Y", F(1, 10)),
                (1, "X", F(1, 50)), (2, "Y", F(1, 40)))
WEYL_COPIES = 4096
OPTICAL_DATA = "model/predictions/data/modular_optical_ratios.json"


def _array(value, shape):
    if (type(value) is not np.ndarray or value.shape != shape
            or value.dtype != np.dtype("complex128") or not np.isfinite(value).all()):
        raise ValueError(f"original finite complex128 array of shape {shape} required")
    return np.frombuffer(value.tobytes(), dtype=np.complex128).reshape(shape)


def _choice(axis, intervention):
    if type(axis) is not str or axis not in ("X", "Y"):
        raise ValueError("the complete X or Y detector setting is required")
    if type(intervention) is not int or intervention not in (0, 1):
        raise ValueError("an original binary intervention setting is required")


def _native_record(state, edge, coin):
    seam = quantum_mixing.seams()[edge]
    after, (_, current) = native.native(state, seam, coin)
    return after, (edge, coin, current, int(CLOCK_PORT in seam))


@lru_cache(maxsize=1)
def native_histories():
    """All 3600 native tapes, grouped only by identical work operations.

    Sampling the request first uses 64 fair placement coins on the empty
    seam (0,1). It leaves INITIAL fixed for every word, so it introduces no
    unrecorded dependence between that request and the later native carrier.
    """
    if quantum_mixing.seams()[0] != (0, 1):
        raise ArithmeticError("the fixed empty sampler seam changed")
    for coin in (0, 1):
        after, record = native.native(INITIAL, (0, 1), coin)
        if after != INITIAL or record != (coin, 0):
            raise ArithmeticError("sampling coins disturb the subsequent native carrier")
    groups = defaultdict(list)
    for e1 in range(30):
        for c1 in (0, 1):
            middle, first = _native_record(INITIAL, e1, c1)
            # The first committed local current is the entire feedback input.
            left, right = quantum_mixing.seams()[e1]
            local_current = first[2]*(int(right == CLOCK_PORT)-int(left == CLOCK_PORT))
            feedback = int(local_current < 0)
            for e2 in range(30):
                for c2 in (0, 1):
                    final, second = _native_record(middle, e2, c2)
                    control = first[3], feedback, second[3]
                    groups[control].append((first, second, final))
    return tuple((key, tuple(rows)) for key, rows in sorted(groups.items()))


def _ideal_filter():
    """Independent endpoint dilation in the declared event chart.

    Only the event basis completion is reused from the existing preparation.
    The shift is evaluated directly from seam endpoints, without its circuit,
    permutation executor, filter compression or native transfer matrix.
    """
    _, event_basis = preparation._lazy_filter_program(F(1, 2))
    output = np.zeros((1 << 15, 12), dtype=complex)
    for event in range(60):
        seam = quantum_mixing.seams()[event//2]
        for site in range(12):
            destination = seam[0]+seam[1]-site if event & 1 and site in seam else site
            for final_event in range(64):
                output[destination+16*final_event, site] += (
                    event_basis[event, final_event]*event_basis[event, 0])
    return output


@dataclass(frozen=True)
class Preparation:
    weights: tuple[F, ...]
    rejection: F
    successes: tuple[np.ndarray, ...]
    failures: tuple[np.ndarray, ...]


def prepare(*, ideal=False):
    """Unnormalized amplitudes on all work/environment/reference outputs."""
    if type(ideal) is not bool:
        raise ValueError("an explicit ideal comparison flag is required")
    _, data = preparation._prepared_edge_gibbs_data()
    raw = _ideal_filter() if ideal else data["raw_filter"]
    raw = _array(raw, (1 << 15, 12))
    weights = preparation.subordinated_filter_weights()
    successful = np.eye(12, dtype=complex)/np.sqrt(12.)
    successes, failures = [], []
    for n in range(len(weights)):
        embedded = np.zeros((16, 12), dtype=complex)
        embedded[:12] = successful
        successes.append(_array(embedded, (16, 12)))
        if n+1 < len(weights):
            full = raw@successful
            failures.append(_array(full[12:].copy(), ((1 << 15)-12, 12)))
            successful = full[:12]
    return Preparation(weights, 1-sum(weights), tuple(successes), tuple(failures))


def _maps(ideal):
    compiled = source.compiled_source()
    return compiled["targets" if ideal else "raw"]


def _prefix(amplitude, control, intervention, maps):
    tick, feedback = control[:2]
    result = maps["clock1"]@amplitude if tick else amplitude
    if feedback and intervention:
        result = maps[f"kick_plus_{BAND}"]@result
    return result


def live_prefix(*, intervention=1, ideal=False):
    """Prepare, execute the first native tick and its committed feedback."""
    _choice("X", intervention)
    prepared, maps = prepare(ideal=ideal), _maps(ideal)
    controls = sorted({key[:2] for key, _ in native_histories()})
    return {(n, *control): _array(_prefix(amplitude, control, intervention, maps), (16, 12))
            for n, amplitude in enumerate(prepared.successes) for control in controls}


def checkpoint(*, intervention=1):
    """Serialize a replay-bound live prefix without normalizing its arrays.

    Failed/rejected outputs are the immutable completed preparation log; the
    checkpoint names that log by its raw binary digest. Replay checks its
    custody before continuing the serialized live arrays.
    """
    _choice("X", intervention)
    prefix, prepared = live_prefix(intervention=intervention), prepare()
    rows = []
    for key, value in sorted(prefix.items()):
        rows.append([list(key), value.astype("<c16", copy=False).tobytes().hex()])
    payload = {
        "schema": "oph.prepared_clock.prefix.v1", "intervention": intervention,
        "source_sha256": bindings(),
        "completed_preparation_sha256": _completed_digest(prepared),
        "live_arrays_complex128_little_endian": rows,
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _completed_digest(prepared):
    digest = hashlib.sha256()
    digest.update(str((prepared.weights, prepared.rejection)).encode())
    # Rejection leaves the supplied joint work/reference input untouched.
    digest.update(prepared.successes[0].tobytes())
    for value in prepared.failures:
        digest.update(value.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class History:
    axis: str
    intervention: int
    preparation: Preparation
    # Each row retains every native tape with the same CP work operation.
    # The terminal address read projects each of the sixteen rows separately.
    successful: dict
    resumed: bool = False


def execute(*, axis="Y", intervention=1, ideal=False, saved=None):
    _choice(axis, intervention)
    if type(ideal) is not bool:
        raise ValueError("an explicit ideal comparison flag is required")
    prepared, maps = prepare(ideal=ideal), _maps(ideal)
    prefix = live_prefix(intervention=intervention, ideal=ideal)
    if saved is not None:
        if ideal or type(saved) is not str or saved != checkpoint(intervention=intervention):
            raise ValueError("checkpoint differs from the complete original source-prefix replay")
        payload = json.loads(saved)
        prefix = {tuple(key): _array(np.frombuffer(bytes.fromhex(raw), dtype="<c16")
                                    .reshape(16, 12).copy(), (16, 12))
                  for key, raw in payload["live_arrays_complex128_little_endian"]}
    result = {}
    for n in range(len(prepared.weights)):
        for control, _ in native_histories():
            amplitude = prefix[n, *control[:2]]
            if control[2]:
                amplitude = maps["clock1"]@amplitude
            amplitude = maps[f"detector_{axis.lower()}_{BAND}"]@amplitude
            result[n, *control] = _array(amplitude, (16, 12))
    return History(axis, intervention, prepared, result, saved is not None)


def _pure_distance(a, b):
    """Half trace norm of |a><a|-|b><b|, using a stable upper bound.

    Vector phases are irrelevant to a classical quantum-output block.
    The independent algebraic rank-two norm formula is checked in tests.
    """
    if np.array_equal(a, b):
        return 0.
    # Finite input entries do not imply finite dot products or norms. A NaN
    # must never bypass the caller's error > allowance rejection branch.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        overlap = np.vdot(b, a)
        if not (np.isfinite(na) and np.isfinite(nb) and np.isfinite(overlap)):
            raise ArithmeticError("quantum-output norm calculation is not finite")
        phase = overlap/abs(overlap) if overlap else 1.
        result = float((na+nb)*np.linalg.norm(a-phase*b)/2)
    if not math.isfinite(result):
        raise ArithmeticError("quantum-output distance calculation is not finite")
    return result


def _validate_preparation(p):
    weights = preparation.subordinated_filter_weights()
    if (type(p) is not Preparation or type(p.weights) is not tuple
            or any(type(w) is not F or type(w.numerator) is not int
                   or type(w.denominator) is not int or w.denominator <= 0
                   or math.gcd(w.numerator, w.denominator) != 1 for w in p.weights)
            or p.weights != weights
            or type(p.rejection) is not F or p.rejection != 1-sum(weights)
            or type(p.successes) is not tuple or len(p.successes) != len(weights)
            or type(p.failures) is not tuple or len(p.failures) != len(weights)-1):
        raise ValueError("the complete unnormalized preparation menu is required")
    for state in p.successes:
        _array(state, (16, 12))
    for state in p.failures:
        _array(state, ((1 << 15)-12, 12))


def _validate(actual):
    if type(actual) is not History:
        raise ValueError("a complete typed clock history is required")
    if type(actual.resumed) is not bool:
        raise ValueError("the continuation marker must be an original Boolean")
    _choice(actual.axis, actual.intervention)
    _validate_preparation(actual.preparation)
    weights = actual.preparation.weights
    expected = {(n, *control) for n in range(len(weights)) for control, _ in native_histories()}
    if (type(actual.successful) is not dict
            or any(type(key) is not tuple or any(type(x) is not int for x in key)
                   for key in actual.successful)
            or set(actual.successful) != expected):
        raise ValueError("every request and native-control outcome must be retained")
    for state in actual.successful.values():
        _array(state, (16, 12))


def complete_distance(actual, expected):
    """Compare every failure and every terminal quantum-reference block."""
    _validate(actual); _validate(expected)
    if (actual.axis, actual.intervention) != (expected.axis, expected.intervention):
        raise ValueError("identical observation menus are required for comparison")
    p, q = actual.preparation, expected.preparation
    result = float(p.rejection)*_pure_distance(p.successes[0], q.successes[0])
    for n, weight in enumerate(p.weights):
        for k in range(n):
            result += float(weight)*_pure_distance(p.failures[k], q.failures[k])
        for control, tapes in native_histories():
            a, b = actual.successful[n, *control], expected.successful[n, *control]
            result += float(weight*F(len(tapes), 3600))*sum(
                _pure_distance(a[j], b[j]) for j in range(16))
    return result


def probabilities(actual):
    """Complete observer law with native tapes factored by identical outputs."""
    _validate(actual)
    result = {("sampler_reject",): float(actual.preparation.rejection)*float(
        np.linalg.norm(actual.preparation.successes[0])**2)}
    for n, weight in enumerate(actual.preparation.weights):
        for k in range(n):
            result["herald_failure", n, k+1] = float(weight)*float(
                np.linalg.norm(actual.preparation.failures[k])**2)
        for control, tapes in native_histories():
            output = actual.successful[n, *control]
            for address in range(16):
                result["read", n, *control, address] = float(weight*F(len(tapes), 3600))*float(
                    np.vdot(output[address], output[address]).real)
    return result


def bindings():
    paths = {*preparation.SOURCE_PATHS, *spatial_limit.SOURCE_PATHS,
             "model/predictions/prepared_clock.py", "model/predictions/prepared_clock_source.py",
             "tools/strict_json.py", OPTICAL_DATA}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def error_budget():
    """Complete CP composition bound under the checked numerical guards.

    These are finite numerical checks, not directed-rounding certificates.
    All first failures are outcomes of each filter. Thus n filters have gain
    at most g**n and half-diamond error (g**n-1)/2, g=(1+delta)**2.
    A triangle through ideal preparation charges continuation error only on
    its successful mass; no successful branch is renormalized.
    """
    _, data = preparation._prepared_edge_gibbs_data()
    residual = float(np.linalg.norm(data["raw_filter"]-_ideal_filter(), "fro"))
    if not math.isfinite(residual) or residual > float(FILTER_GUARD):
        raise ArithmeticError("the complete raw filter exceeds its numerical guard")
    compiled = source.compiled_source()
    guards = compiled["guards"]
    # Both axes and every branch use at most two clocks, one kick and one read.
    delta = ((1+guards["clock1"])**2*(1+guards[f"kick_plus_{BAND}"])
             *(1+max(guards[f"detector_{a}_{BAND}"] for a in ("x", "y")))-1)
    gain = (1+delta)**2
    continuation_error = (gain-1)/2
    weights = preparation.subordinated_filter_weights()
    g = (1+FILTER_GUARD)**2
    prep_error = sum((w*(g**n-1)/2 for n, w in enumerate(weights)), F())
    success = F(source.prepared_reference()["normalization"])/12
    combined = gain*prep_error+success*continuation_error+ARITHMETIC_ALLOWANCE
    # Outward rational reporting keeps the receipt reviewable without changing
    # the exact product calculation that establishes each upper allowance.
    upper = lambda value: F(math.ceil(value*10**12), 10**12)
    return {"filter_operator_guard": FILTER_GUARD,
        "filter_complete_Frobenius_residual": residual,
        "stopped_preparation_half_diamond_upper_under_guards": upper(prep_error),
        "continuation_half_diamond_upper_under_guards": upper(continuation_error),
        "continuation_channel_norm_upper_under_guards": upper(gain),
        "complete_half_trace_upper_under_guards": upper(combined),
        "complete_event_weight_error_upper_under_guards": upper(2*combined),
        "numerical_arithmetic_allowance": ARITHMETIC_ALLOWANCE,
        "ideal_success_probability": success,
        "conditioned_success_half_trace_upper_under_guards": min(F(1), upper(2*combined/success)),
        "directed_rounding_or_physical_fidelity_certificate": False}


def audit(actual, target):
    """The acceptance gate checks complete quantum outputs, not just masses."""
    error = complete_distance(actual, target)
    budget = error_budget()
    allowance = float(budget["complete_half_trace_upper_under_guards"])
    probabilities_actual = probabilities(actual)
    mass = sum(probabilities_actual.values())
    if (not math.isfinite(error) or not math.isfinite(mass)
            or error > allowance or abs(mass-1) > 2*allowance
            or any(not math.isfinite(p) or p < 0 for p in probabilities_actual.values())):
        raise ArithmeticError("complete source clock/history quantum outputs exceed the budget")
    # Returned preparation states are part of the exposed continuation data.
    for a, b in zip(actual.preparation.successes, target.preparation.successes, strict=True):
        if _pure_distance(a, b) > 4*float(FILTER_GUARD)+float(ARITHMETIC_ALLOWANCE):
            raise ArithmeticError("a source preparation checkpoint exceeds its allowance")
    return {"axis": actual.axis, "intervention": actual.intervention,
        "complete_quantum_half_trace_upper_diagnostic": error,
        "complete_raw_mass_diagnostic": mass,
        "retained_coarse_records": len(probabilities_actual),
        "native_tapes_per_successful_request": 3600,
        "complete_uncompressed_records": 1+10+len(actual.preparation.weights)*3600*16,
        "all_terminal_reference_entries_compared": True,
        "all_first_failure_environment_rows_compared": True,
        "resumed_from_serialized_live_prefix": actual.resumed}


def resources():
    edge, _ = preparation._prepared_edge_gibbs_data()
    compiled = source.compiled_source()
    words = compiled["words"]
    continuation_word = (words["clock1"]+words[f"kick_plus_{BAND}"]
                         +words["clock1"]+words[f"detector_y_{BAND}"])
    continuation_pulses = max(
        2*len(words["clock1"])+len(words[f"kick_plus_{BAND}"])+len(words[f"detector_{a}_{BAND}"])
        for a in ("x", "y"))
    weights = preparation.subordinated_filter_weights()
    reference = source.prepared_reference()
    lazy = reference["lazy_filter"]
    expected_filters = sum((w*sum((F((lazy**(2*k)).trace())/12 for k in range(n)), F())
                            for n, w in enumerate(weights)), F())
    success = F(reference["normalization"])/12
    expected_continuation = success*(
        F(1, 3)*len(words["clock1"])+F(1, 15)*len(words[f"kick_plus_{BAND}"])
        +max(len(words[f"detector_{a}_{BAND}"]) for a in ("x", "y")))
    filter_pulses = edge["one_filter_resources"]["native_pulses"]
    maximum = 4*filter_pulses+continuation_pulses
    info = compiler.resources(continuation_word)
    filter_program, _ = preparation._lazy_filter_program(F(1, 2))
    program_bytes = len(json.dumps(
        {"filter": filter_program, "continuation_words": words},
        sort_keys=True, separators=(",", ":"), allow_nan=False).encode())
    thresholds = [row["past_last_word"] for row in sampler_contract()
                  if row["request"] is not None]
    if any(not 0 < value < 1 << 64 for value in thresholds):
        raise ArithmeticError("the sampler thresholds do not fit their declared storage")
    return {"maximum_filter_requests": 4,
        "one_filter_source_pulses": filter_pulses,
        "maximum_successful_continuation_source_pulses": continuation_pulses,
        "maximum_total_source_pulses": maximum,
        "maximum_proper_code_load_events": maximum,
        "maximum_proper_code_export_events": maximum,
        "mean_executed_filters_with_first_failure_stopping": expected_filters,
        "ideal_mean_total_source_pulses_upper": expected_filters*filter_pulses+expected_continuation,
        "fair_native_sampler_coin_events": 64,
        "native_carrier_placement_coin_events_maximum": 66,
        "native_carrier_placement_coin_events_ideal_mean": 64+2*success,
        "uniform_thirty_seam_selection_service_calls_maximum": 2,
        "uniform_thirty_seam_selection_service_calls_ideal_mean": 2*success,
        "exact_uniform_seam_selector_is_a_supplied_service": True,
        "finite_worst_case_fair_bit_budget_for_seam_selection_claimed": False,
        "subsequent_native_attempts_maximum": 2,
        "ideal_subsequent_native_attempts_mean": 2*success,
        "clock_ticks_per_successful_trial_mean": F(1, 3),
        "feedback_kicks_per_successful_trial_mean": F(1, 15),
        "native_carriers": 1, "native_ports": 12,
        "protected_work_qubits": 4, "reference_storage_qubits": 4,
        "reference_populated_dimension": 12,
        "event_and_scratch_qubits_during_filter": 11,
        "peak_encoded_qubits_including_reference": 19,
        "source_processor_dimension": 6,
        "classical_record_storage_bits": 64+2*8+4+3+3+1,
        "record_encoding": "64 sampler bits; two seam(5)/coin(1)/current(2) records; "
                           "address(4), request(3), failure index(3), status(1); schema distinguishes branches",
        "fixed_attempt_horizon_not_wait_for_two_ticks": True,
        "zero_tick_trials_retained": True,
        "mode_basis_and_detector_costs_included": True,
        "source_word_packet_utf8_bytes": program_bytes,
        "source_word_packet_byte_cursor_bits": (program_bytes-1).bit_length(),
        "expanded_source_pulse_ordinal_bits": maximum.bit_length(),
        "sampler_cumulative_threshold_storage_bits": 64*len(thresholds),
        "observer_setting_storage_bits": 2,
        "controller_interpreter_working_memory_included": False,
        "filter_request_counter_bits": 3,
        "attempt_counter_bits": 2,
        "herald_measurements_maximum": 4,
        "herald_measurements_ideal_mean": expected_filters,
        "terminal_address_measurements_maximum": 1,
        "terminal_address_measurements_ideal_mean": success,
        "initial_native_state_work_I12_and_purification_supplied": True,
        "herald_and_terminal_readout_services_supplied": True,
        "source_filter_routine_reuse_removes_no_service_invocations": True,
        "example_continuation_load_export_events": info["proper_code_service_events"],
        "physical_duration_or_error_of_supplied_services_bounded": False}


def sampler_contract():
    weights = preparation.subordinated_filter_weights()
    rows, boundary = [], 0
    for label, weight in enumerate((*weights, 1-sum(weights))):
        width = weight*(1 << 64)
        if width.denominator != 1:
            raise ArithmeticError("the request sampler is not on its declared coin grid")
        width = int(width)
        if not width:
            continue
        selected = label if label < len(weights) else None
        for integer in (boundary, boundary+width-1):
            state, records = INITIAL, []
            for bit in map(int, f"{integer:064b}"):
                state, record = native.native(state, (0, 1), bit)
                records.append(record)
            decoded = preparation.select_subordinated_request(
                tuple(coin for coin, _ in records), weights)
            if decoded != selected or state != INITIAL:
                raise ArithmeticError("request sampler does not continue the same native carrier")
        rows.append({"request": selected, "first_word": boundary,
                     "past_last_word": boundary+width, "probability": weight})
        boundary += width
    if boundary != 1 << 64:
        raise ArithmeticError("request and rejection intervals omit native coin words")
    return rows


def shared_fluctuation_modes():
    """Eleven fluctuation pairs from ONE successful twelve-state marginal.

    No independent state is assigned to a mode. Matrices here are expressed
    in the declared spectral basis; only three band representatives have
    the source-controlled pair readout used below. The other modes are an
    algebraic statement, not a claim of executed access to eleven modes.
    """
    reference = source.prepared_reference()
    populations = np.array([float(value/reference["normalization"])
                            for value in reference["masses"]])
    density = np.diag(np.r_[populations, np.zeros(4)]).astype(complex)
    modes, observables = [], []
    for address in range(1, 12):
        delta = populations[0]-populations[address]
        if not delta > 0:
            raise ArithmeticError("the prepared fluctuation pair is not polarized")
        x = np.zeros((16, 16), complex)
        x[0, address] = x[address, 0] = 1
        y = np.zeros((16, 16), complex)
        y[0, address], y[address, 0] = -1j, 1j
        q, p = x/math.sqrt(2*delta), y/math.sqrt(2*delta)
        gap = reference["modular_gaps"][address]
        variance = (populations[0]+populations[address])/(2*delta)
        modes.append({"address": address, "delta": float(delta), "gap": gap,
            "variance": float(variance), "Q": q, "P": p,
            "equal_mode_commutator_L2_variance_times_copy_count":
                float((populations[0]+populations[address])/delta**2-1)})
        observables.extend((q, p))
    two_point = np.array([[np.trace(density@a@b) for b in observables]
                          for a in observables])
    covariance, symplectic = two_point.real, 2*two_point.imag
    wanted_symplectic = np.kron(np.eye(11), [[0., 1.], [-1., 0.]])
    wanted_covariance = np.diag([mode["variance"] for mode in modes for _ in range(2)])
    if (np.max(abs(symplectic-wanted_symplectic)) > 1e-12
            or np.max(abs(covariance-wanted_covariance)) > 1e-12
            or any(abs(mode["variance"]-.5/math.tanh(mode["gap"]/2)) > 2e-11
                   for mode in modes)):
        raise ArithmeticError("the shared-state fluctuation dictionary failed")
    return {"density": density, "modes": tuple(modes),
        "covariance": covariance, "symplectic": symplectic}


def _weyl_choice(factors, copies):
    if type(copies) is not int or not 1 <= copies <= 1 << 40:
        raise ValueError("an original positive copy count at most 2^40 is required")
    if (type(factors) is not tuple or not 1 <= len(factors) <= 8
            or any(type(row) is not tuple or len(row) != 3
                   or type(row[0]) is not int or row[0] not in (1, 2, 3)
                   or type(row[1]) is not str or row[1] not in ("X", "Y")
                   or type(row[2]) is not F or abs(row[2]) > 1 for row in factors)):
        raise ValueError("one to eight declared band/axis/rational-coefficient factors required")


def ordered_weyl_law(factors=WEYL_FACTORS, copies=WEYL_COPIES):
    """Direct finite-copy law and a quantitative Gaussian limit.

    Source factors are chronological: U=exp(i A_m/sqrt(N))...exp(i A_1/sqrt(N)).
    Independent identical copies give z_N**N exactly. The parity commuting
    with rho and anticommuting with each A eliminates all odd Taylor terms.
    Taylor's integral remainder and a contraction telescoping bound give
    |z_N**N-exp(-C)| <= (D**4/24+|C|**2/2)/N. This is a finite bound, not
    an assertion that any finite register has exact canonical commutators.
    """
    _weyl_choice(factors, copies)
    data, reference = shared_fluctuation_modes(), source.prepared_reference()
    addresses = source.PAIR_ADDRESSES
    matrices, source_factors = [], []
    unitary = np.eye(16, dtype=complex)
    norm_sum = 0.
    for band, axis, coefficient in factors:
        address = addresses[band-1]
        mode = data["modes"][address-1]
        operator = mode["Q" if axis == "X" else "P"]
        matrices.append(float(coefficient)*operator)
        norm_sum += abs(float(coefficient))/math.sqrt(2*mode["delta"])
        exact_delta = ((reference["masses"][0]-reference["masses"][address])
                       /reference["normalization"])
        angle = sp.Rational(coefficient)/sp.sqrt(2*exact_delta*copies)
        source_factors.append((band, axis, source._dyadic(angle)))
        theta = float(sp.N(angle, 60))
        pair = operator*math.sqrt(2*mode["delta"])
        projector = pair@pair
        unitary = (np.eye(16)+(math.cos(theta)-1)*projector
                   +1j*math.sin(theta)*pair)@unitary
    density = data["density"]
    c = sum((np.trace(density@a@a)/2 for a in matrices), 0j)
    c += sum((np.trace(density@matrices[later]@matrices[earlier])
              for later in range(len(matrices)) for earlier in range(later)), 0j)
    if c.real < -1e-13 or not (math.isfinite(norm_sum) and np.isfinite(c)):
        raise ArithmeticError("the ordered fluctuation exponent is not finite and positive")
    single = complex(np.trace(density@unitary))
    limit = complex(np.exp(-c))
    bound = (norm_sum**4/24+abs(c)**2/2)/copies
    return {"source_factors": tuple(source_factors), "matrices": tuple(matrices),
        "one_copy_unitary_spectral": unitary, "one_copy_characteristic": single,
        "independent_copy_characteristic": single**copies,
        "gaussian_characteristic": limit, "gaussian_exponent": complex(c),
        "operator_norm_sum": norm_sum, "finite_copy_error_upper": min(2., bound),
        "uncapped_finite_copy_error_upper": bound,
        "angle_rounding_operator_error_upper": F(len(factors), 1 << (source.SYNTHESIS_BITS+1))}


@dataclass(frozen=True)
class PointerHistory:
    factors: tuple
    copies: int
    axis: str
    preparation: Preparation
    successful: tuple[np.ndarray, ...]


def execute_weyl(*, factors=WEYL_FACTORS, copies=WEYL_COPIES, axis="X", ideal=False):
    """One source-prepared register and one actual controlled pointer read.

    ``copies`` sets the pulse angle, not the number of executed registers.
    All first failures and 32 joint pointer/work addresses retain the same
    twelve reference columns. Product-copy expectations are derived later.
    """
    _choice(axis, 0)
    if type(ideal) is not bool:
        raise ValueError("an explicit ideal comparison flag is required")
    law = ordered_weyl_law(factors, copies)
    routine = source.pair_pointer_source(law["source_factors"], axis=axis)
    matrix = routine["target" if ideal else "raw"]
    prepared, outputs = prepare(ideal=ideal), []
    for amplitude in prepared.successes:
        embedded = np.zeros((32, 12), complex)
        embedded[:16] = amplitude
        outputs.append(_array(matrix@embedded, (32, 12)))
    return PointerHistory(factors, copies, axis, prepared, tuple(outputs))


def _validate_pointer(history):
    if type(history) is not PointerHistory:
        raise ValueError("a complete typed pointer history is required")
    _choice(history.axis, 0)
    _weyl_choice(history.factors, history.copies)
    _validate_preparation(history.preparation)
    if (type(history.successful) is not tuple
            or len(history.successful) != len(history.preparation.weights)):
        raise ValueError("every pointer request and terminal address must be retained")
    for amplitude in history.successful:
        _array(amplitude, (32, 12))


def pointer_probabilities(history):
    _validate_pointer(history)
    p = history.preparation
    result = {("sampler_reject",): float(p.rejection)*float(np.linalg.norm(p.successes[0])**2)}
    for n, weight in enumerate(p.weights):
        for k in range(n):
            result["herald_failure", n, k+1] = float(weight)*float(np.linalg.norm(p.failures[k])**2)
        for address, row in enumerate(history.successful[n]):
            result["read", n, address] = float(weight)*float(np.vdot(row, row).real)
    if any(not math.isfinite(value) or value < 0 for value in result.values()):
        raise ArithmeticError("the complete pointer law is not finite and nonnegative")
    return result


def pointer_audit(actual, target):
    _validate_pointer(actual); _validate_pointer(target)
    if (actual.factors, actual.copies, actual.axis) != (target.factors, target.copies, target.axis):
        raise ValueError("identical pointer observation menus are required")
    law = ordered_weyl_law(actual.factors, actual.copies)
    routine = source.pair_pointer_source(law["source_factors"], axis=actual.axis)
    budget = error_budget()
    p, q = actual.preparation, target.preparation
    error = float(p.rejection)*_pure_distance(p.successes[0], q.successes[0])
    for n, weight in enumerate(p.weights):
        for k in range(n):
            error += float(weight)*_pure_distance(p.failures[k], q.failures[k])
        error += float(weight)*sum(_pure_distance(a, b) for a, b in
                                   zip(actual.successful[n], target.successful[n], strict=True))
        if _pure_distance(p.successes[n], q.successes[n]) > 4*float(FILTER_GUARD)+float(ARITHMETIC_ALLOWANCE):
            raise ArithmeticError("the retained pointer preparation violates its numerical guard")
    guard = routine["guard"]
    gain = (1+guard)**2
    success = budget["ideal_success_probability"]
    allowance = (gain*budget["stopped_preparation_half_diamond_upper_under_guards"]
                 +success*(gain-1)/2+ARITHMETIC_ALLOWANCE)
    probabilities_ = pointer_probabilities(actual)
    total = sum(probabilities_.values())
    if (not math.isfinite(error) or error > float(allowance)
            or abs(total-1) > 2*float(allowance)):
        raise ArithmeticError("complete pointer quantum outputs violate the numerical guards")
    contrast = sum(value*(1 if key[2] < 16 else -1)
                   for key, value in probabilities_.items() if key[0] == "read")
    return {"axis": actual.axis, "complete_raw_mass_diagnostic": total,
        "complete_quantum_half_trace_upper_diagnostic": error,
        "complete_half_trace_upper_under_guards": allowance,
        "unconditional_pointer_contrast": contrast,
        "retained_coarse_records": len(probabilities_),
        "joint_pointer_work_addresses_per_success": 32,
        "all_terminal_reference_and_first_failure_environment_entries_compared": True,
        "source_word_sha256": routine["word_sha256"],
        "source_pulses_after_success": routine["pulse_count"]}


def shared_state_field():
    """Same-state noise/response/clock law, with a source pointer witness."""
    data, law = shared_fluctuation_modes(), ordered_weyl_law()
    rows, histories = [], {}
    for axis in ("X", "Y"):
        actual, target = execute_weyl(axis=axis), execute_weyl(axis=axis, ideal=True)
        rows.append(pointer_audit(actual, target))
        histories[axis] = actual
    success = error_budget()["ideal_success_probability"]
    estimator = complex(rows[0]["unconditional_pointer_contrast"],
                        rows[1]["unconditional_pointer_contrast"])/float(success)
    one_error = (2*math.sqrt(2)*max(float(row["complete_half_trace_upper_under_guards"])
                                  for row in rows)/float(success)
                 +float(law["angle_rounding_operator_error_upper"]))
    # The raw numerical estimator need not be a contraction. Its possible
    # amplification is charged, not removed by normalization or clipping.
    product_error = math.expm1(WEYL_COPIES*math.log1p(one_error))
    if abs(estimator-law["one_copy_characteristic"]) > one_error:
        raise ArithmeticError("the source pointer fails the ordered one-copy law")
    combined = law["finite_copy_error_upper"]+product_error
    if abs(estimator**WEYL_COPIES-law["gaussian_characteristic"]) > combined:
        raise ArithmeticError("the source estimator fails the finite-copy bound")
    # The request-zero branch is already retained in the same experiment.
    zero = []
    for axis in ("X", "Y"):
        amplitude = histories[axis].successful[0]
        zero.append(float(np.linalg.norm(amplitude[:16])**2-np.linalg.norm(amplitude[16:])**2))
    modes = [{key: value for key, value in mode.items() if key not in ("Q", "P")}
             for mode in data["modes"]]
    pair = lambda value: [float(value.real), float(value.imag)]
    filter_pulses = resources()["one_filter_source_pulses"]
    maximum_pointer = max(row["source_pulses_after_success"] for row in rows)
    return {"reference": "the ideal successful work marginal M/Tr(M); raw source error is separately bounded",
        "mode_specific_preparations_supplied": False,
        "formal_fluctuation_pairs": modes,
        "executed_representative_addresses": list(source.PAIR_ADDRESSES),
        "covariance_rule": "Cov(Qj,Qk)=Cov(Pj,Pk)=delta_jk*coth(gj/2)/2; Cov(Qj,Pk)=0",
        "expected_commutator_rule": "< [Qj,Pk] >= i*delta_jk; cross Q/Q and P/P vanish",
        "clock_rule": "Qj(t)=cos(gj*t)Qj+sin(gj*t)Pj under the same H=-log(M/m0)",
        "ordered_time_correlation": "<Qj(t)Qk>=delta_jk*(coth(gj/2)*cos(gj*t)-i*sin(gj*t))/2",
        "finite_copy_rule": "z_N^N; |z_N^N-exp(-C)| <= (D^4/24+|C|^2/2)/N",
        "factors_chronological": [[band, axis, str(coefficient)] for band, axis, coefficient in WEYL_FACTORS],
        "angle_copy_parameter": WEYL_COPIES, "executed_work_registers_per_trial": 1,
        "one_copy_characteristic": pair(law["one_copy_characteristic"]),
        "independent_copy_characteristic": pair(law["independent_copy_characteristic"]),
        "gaussian_characteristic": pair(law["gaussian_characteristic"]),
        "raw_pointer_estimator": pair(estimator),
        "finite_copy_error_upper": law["finite_copy_error_upper"],
        "one_copy_source_estimator_error_upper_under_guards": one_error,
        "source_estimator_power_error_upper_under_guards": product_error,
        "combined_gaussian_error_upper_under_guards": combined,
        "resolved_ordered_phase_lower_under_guards": law["gaussian_characteristic"].imag-combined,
        "request_zero_pointer_characteristic_control": zero,
        "request_zero_gaussian_symplectic_form_vanishes": True,
        "request_zero_finite_ordering_phase_required_to_vanish": False,
        "executions": rows,
        "resources": {"maximum_filters_per_attempt": 4,
            "maximum_source_pulses_per_executed_trial": 4*filter_pulses+maximum_pointer,
            "one_load_and_export_per_source_pulse": True,
            "sampler_native_coin_events_per_attempt": 64,
            "additional_native_attempts_in_pointer_protocol": 0,
            "protected_pointer_qubits_after_success": 1,
            "peak_encoded_qubits_including_reference": 19,
            "pointer_prepared_after_filter_environment_is_released": True,
            "ideal_expected_preparation_attempts_for_N_successful_copies": F(WEYL_COPIES)/success,
            "independent_copy_collection_executed": False,
            "unbounded_retry_has_no_finite_worst_case_cost": True},
        "independent_identical_successful_copies_are_a_preparation_premise": True,
        "multicopy_pointer_circuit_or_growing_spatial_encoder_executed": False,
        "all_eleven_mode_controls_executed": False,
        "central_group_image_algebra_is_the_fluctuation_algebra": False,
        "source_numerical_guards_are_physical_fidelities": False,
        "empirical_seconds_interacting_QFT_or_spacetime_identified": False}


def native_full_port_eigenvalue_lower(level, index):
    """Finite positive eigenvalue bound for D_l=4^l L_l, counted from one.

    K1 and Schur minimization control the coarse cell residual by
    624*4^-k E(v); the native internal gap controls its orthogonal port
    residual by E(v)/(2*4^l). The centered coarse range has rank20*4^k-1,
    so min-max bounds eigenvalue20*4^k. At k=l the cell residual is zero.
    The whole-domain gap handles the first nineteen positive eigenvalues.
    """
    level = preparation._integer(level, 0, 100)
    scale = 4**level
    index = preparation._integer(index, 1, 240*scale-1)
    if index < 20:
        return F(scale, 4096*scale+1)
    if index >= 20*scale:
        return F(2*scale)
    coarse_scale = 1
    while 80*coarse_scale <= index:
        coarse_scale *= 4
    return F(coarse_scale)/(spatial_limit.COARSE_ENERGY_TAIL_CONSTANT+F(1, 2))


def native_prepared_trace_bound(level, inverse_temperature, heat_step, retained_mass_lower):
    """Trace, anchor gap and success bounds from the finite native source.

    Every positive eigenvalue of D_l is >=n/81940. Thus its square-root
    heat trace is <=1+163880/b², uniformly BEFORE taking a continuum limit.
    For the native lazy filter F=(1-p)I+pK=I-hD_l/2, p=630h16^l<=1.
    Its full-domain norm obeys h||D_l||<=12h4^l<2. Positive retained
    subordination gives M<=exp[-b sqrt(D_l-hD_l²/4)]<=exp[-b sqrt(D_l/2)].
    Z=Tr(M)/m0 and the successful probability is Tr(M)/(240*4^l).
    The sampler coefficients must be retained below that positive law;
    this function does not implement a growing preparation or compiler.
    """
    level = preparation._integer(level, 0, 100)
    beta, step, mass = (preparation._fraction(value, 0)
                        for value in (inverse_temperature, heat_step, retained_mass_lower))
    if not beta or not step or not 0 < mass <= 1:
        raise ValueError("positive temperature, step and retained mass at most one are required")
    scale, count_constant = 4**level, F(81940)
    dimension = 240*scale
    laziness = 630*step*scale**2
    if laziness > 1:
        raise ValueError("the native lazy filter requires630*heat_step*16^level<=1")
    unnormalized_upper = 1+4*count_constant/beta**2
    logarithmic_mass_loss = (1-mass)/mass
    gap = beta/182 if logarithmic_mass_loss <= beta/182 else None
    return {"level": level, "full_port_dimension": dimension,
        "generator": "D_l=4^l L_l on the complete native port space",
        "positive_eigenvalue_index_denominator": count_constant,
        "positive_generator_gap_lower": F(scale, 4096*scale+1),
        "square_root_Gibbs_trace_upper": min(F(dimension), 1+2*count_constant/beta**2),
        "inverse_temperature": beta, "heat_step": step,
        "native_lazy_event_probability": laziness,
        "full_domain_heat_step_times_generator_norm_upper": 12*step*scale,
        "retained_sampler_mass_lower": mass,
        "unnormalized_filter_trace_upper": min(F(dimension), unnormalized_upper),
        "normalized_partition_upper": min(F(dimension), unnormalized_upper/mass),
        "uniform_normalized_partition_upper": unnormalized_upper/mass,
        "ideal_herald_success_lower": mass/dimension,
        "ideal_herald_success_upper": min(F(1), unnormalized_upper/dimension),
        "logarithmic_zero_mass_loss_upper": logarithmic_mass_loss,
        "prepared_and_target_positive_gap_floor": gap,
        "gap_floor_includes_no_raw_compiler_error": True,
        "finite_native_counting_proved_before_continuum_limit": True,
        "positive_sampler_and_full_native_domain_required": True,
        "growing_state_preparation_executed": False}


@lru_cache(maxsize=1)
def native_trace_controls():
    """Actual level-zero/one spectra and the forbidden endpoint alias."""
    rows = []
    for level in (0, 1):
        graph, laplacian, *_ = spatial_source.matrices(level)
        values = np.linalg.eigvalsh((4**level*laplacian).toarray().astype(float))
        positive = values[1:]
        lower = np.array([float(native_full_port_eigenvalue_lower(level, index))
                          for index in range(1, graph.ports)])
        if np.min(positive-lower) < -1e-10 or abs(values[0]) > 1e-10:
            raise ArithmeticError("the actual full-port spectrum violates its finite count bound")
        # The exact native constant contributes one. Applying sqrt to its
        # floating residual would amplify backend noise at the zero mode.
        trace = 1+float(np.exp(-np.sqrt(positive)).sum())
        bound = native_prepared_trace_bound(level, F(1), F(1, 1024*16**level), F(1023, 1024))
        if trace > float(bound["square_root_Gibbs_trace_upper"])+1e-10:
            raise ArithmeticError("the finite source trace exceeds its proved upper bound")
        rows.append({"level": level, "dimension": graph.ports,
            "minimum_positive_eigenvalue_over_index_diagnostic":
                float(np.min(positive/np.arange(1, graph.ports))),
            "square_root_Gibbs_trace_b1_diagnostic": trace,
            "bound": bound})
    # A contractive filter alone can preserve arbitrarily high energy:
    # D=diag(0,4/h), F=diag(1,-1), hence every F^(2n)=I.
    return {"finite_source_spectra": rows,
        "contractive_endpoint_alias_control": {
            "heat_step_times_high_energy": 4, "high_energy_filter_eigenvalue": -1,
            "every_even_filter_power_eigenvalue": 1,
            "prepared_high_to_anchor_ratio": 1,
            "target_high_to_anchor_ratio_tends_to_zero_as_heat_step_vanishes": True,
            "native_lazy_realizability_admits_this_endpoint": False},
        "uniform_partition_bound_removes_no_herald_cost": True}


def growing_fluctuation_bound(dimension, copies, factor_radii, *, gap_floor,
        inverse_temperature, clock_gap_error=F(0), maximum_time=F(0),
        complete_preparation_error=F(0), success_lower=F(1),
        one_copy_word_error=F(0), spectral_pair_error=F(0),
        cutoff_tail_radius=F(0), partition_upper=None):
    """A conditional finite bound for one growing shared-state field.

    Each star factor has phase-space radius R_k, where the field variables
    are q_j=Q_j/sqrt(omega_j), pi_j=sqrt(omega_j)P_j. The same reference
    rho has p0=max p_j, delta_j=p0(1-exp(-g_j)); Z=1/p0<=dimension.
    Both g_j and b*omega_j are >=gap_floor on the certified band. The
    spectral error bounds the normalized two-point forms under the actual
    common embedding; it is not inferred from a supplied eigenvalue list.
    The cutoff tail belongs to the target field, with b*omega>=gap_floor
    also on that tail. Finite-source observables are band-filtered; no
    thermal bound on uncontrolled high source modes is silently assumed.

    Products and their derivatives preserve anchor+span(v_1,...,v_m), of
    dimension at most m+1. Hence |Tr rho V''''| <= (m+1)*p0*D^4. This
    improves the crude norm-only d^2/N allowance to Z/N. Supplied source
    and conditioning errors are charged before the N-copy power. No growing
    source preparation, spatial encoder or independent copies are inferred.
    """
    dimension = preparation._integer(dimension, 2)
    copies = preparation._integer(copies, 1)
    if type(factor_radii) is not tuple or not 1 <= len(factor_radii) <= 64:
        raise ValueError("a fixed nonempty tuple of phase-space radii is required")
    radii = tuple(preparation._fraction(value, 0) for value in factor_radii)
    gamma, beta = (preparation._fraction(value, 0)
                   for value in (gap_floor, inverse_temperature))
    if not gamma or not beta:
        raise ValueError("positive gap and inverse-temperature bounds are required")
    error, duration, prep_error, word_error, spatial_error, tail = (
        preparation._fraction(value, 0) for value in
        (clock_gap_error, maximum_time, complete_preparation_error,
         one_copy_word_error, spectral_pair_error, cutoff_tail_radius))
    success = preparation._fraction(success_lower, 0, 1)
    if not success:
        raise ValueError("a positive herald probability lower bound is required")
    partition = (F(dimension) if partition_upper is None else
                 preparation._fraction(partition_upper, 1, dimension))
    radius, factors = sum(radii), len(radii)
    # exp(gamma)>=1+gamma gives 1-exp(-gamma)>=gamma/(1+gamma).
    theta = gamma/(1+gamma)
    gaussian_error = radius**4*((factors+1)*partition/F(96)+F(1, 8))/(theta**2*copies)
    normalized_state_error = 4*prep_error/success if 2*prep_error < success else None
    single_error = (normalized_state_error+word_error
                    if normalized_state_error is not None else None)
    amplified = copies*single_error if single_error is not None else None
    source_error = (amplified/(1-amplified)
                    if amplified is not None and amplified < 1 else None)
    # |d[coth(g/2)/2]/dg|<=1/gamma^2, and clock coefficient
    # rotations differ by at most T*|g-b*omega|/b on the same band.
    covariance_error = radius**2*(error/gamma**2+2*duration*error/(beta*theta)
                                  +spatial_error)/2
    cutoff_error = (radius*tail+tail**2/2)/theta
    total = (gaussian_error+source_error+covariance_error+cutoff_error
             if source_error is not None else None)
    return {"dimension": dimension, "copies": copies, "factor_count": factors,
        "phase_space_radius_sum": radius, "inverse_anchor_weight_upper": partition,
        "positive_gap_floor": gamma, "inverse_temperature": beta,
        "finite_copy_gaussian_error_upper": gaussian_error,
        "conditional_state_trace_norm_error_upper": normalized_state_error,
        "one_copy_characteristic_error_upper": single_error,
        "N_copy_raw_source_error_upper": source_error,
        "clock_and_common_spectral_error_upper": covariance_error,
        "spectral_cutoff_error_upper": cutoff_error,
        "total_ordered_characteristic_error_upper": total,
        "herald_success_lower": success,
        "expected_preparation_attempts_upper": F(copies)/success,
        "clock_parameter": "model time t for H/b; modular parameter s=-t/b",
        "native_event_count_to_model_time_error_requires_its_separate_bound": True,
        "partition_bound_is_a_supplied_or_separately_proved_input": partition_upper is not None,
        "common_source_state_encoder_and_independent_copies_executed": False,
        "finite_band_convergence_alone_implies_local_smear_convergence": False}


def native_growing_field_contract(level, inverse_temperature=F(1), band=F(8)):
    """Join finite-source trace, preparation and same-state detector budgets.

    This is a two-factor fixed-band contract, with phase-space radii1/4,
    |t|<=1 and N=2^level independent successful copies. It uses the same
    finite D_l on both sides, so spatial/cutoff errors are zero here, not
    asserted zero in the continuum. The growing words and source accuracy
    are constructive resource obligations, not an executed apparatus.
    """
    resources = source.native_growing_source_resources(level, inverse_temperature, band)
    mass = 1-resources["missing_sampler_mass_upper"]
    result = {"source_resources": resources,
        "factor_count": 2, "copies": 2**level,
        "spatial_and_cutoff_errors_zero_only_for_the_same_finite_generator": True,
        "independent_successful_copies_and_precise_services_are_required": True,
        "growing_words_or_continuum_identification_executed": False}
    if mass <= 0:
        return {**result, "native_trace_bound": None, "joint_field_bound": None}
    traced = native_prepared_trace_bound(level, inverse_temperature, resources["heat_step"], mass)
    result["native_trace_bound"] = traced
    gap, ideal_error = (traced["prepared_and_target_positive_gap_floor"],
                        resources["ideal_prepared_generator_error_upper"])
    if gap is None or ideal_error is None:
        return {**result, "joint_field_bound": None}
    words = resources["resources"]
    pulses = (words["ordered_pointer_common_encoder_and_readout_pulses_upper"]
        +2*words["controlled_complex_star_factor_pulses_upper_in_spectral_basis"]
        +4*words["one_same_state_spectral_clock_word_pulses_upper"])
    # Clock sandwiches are included. Reserve factors for nonunitary raw
    # gain and both X/Y detector contrasts, not just a nominal unitary word.
    bits = 2*resources["sampler_index"]+(8*pulses-1).bit_length()
    linear_error = F(pulses, 1 << bits)
    complex_error = 4*linear_error/(1-linear_error)
    bound = growing_fluctuation_bound(traced["full_port_dimension"], 2**level,
        (F(1, 4), F(1, 4)), gap_floor=gap, inverse_temperature=inverse_temperature,
        clock_gap_error=ideal_error, maximum_time=F(1),
        complete_preparation_error=resources["complete_preparation_half_trace_error_upper"],
        success_lower=traced["ideal_herald_success_lower"], one_copy_word_error=complex_error,
        partition_upper=traced["normalized_partition_upper"])
    result.update(ordered_pointer_and_clock_sandwich_pulses_upper=pulses,
        required_detector_with_services_error_base2_negative_exponent=bits,
        one_copy_complex_estimator_word_error_upper=complex_error,
        joint_field_bound=bound,
        partition_gap_and_success_inputs_derived_from_the_same_native_source=True,
        nominal_clock_vs_actual_state_log_error_is_separately_retained=True)
    return result


def growing_field_controls():
    """Exact finite-family countercontrols and actual clock normalization.

    The equal-gap d-state family is a mathematical sharpness control, not
    a source-prepared growing Gibbs state. The separate twelve-state mode
    checks use the existing source preparation and clock on their actual
    register, with their complete preparation/error contract unchanged.
    """
    gaussian, rows = math.exp(-.75), []
    for dimension in (16, 64, 256, 1024):
        pair = []
        # rho=diag(2,1,...,1)/(d+1), delta=1/(d+1), q=X/sqrt(2delta).
        # The exact support trace is 3/(d+1), giving (d+10)/(32N).
        for copies in (dimension, dimension**2):
            angle = math.sqrt((dimension+1)/(2*copies))
            one_minus = 6*math.sin(angle/2)**2/(dimension+1)
            characteristic = math.exp(copies*math.log1p(-one_minus))
            bound = F(dimension+10, 32*copies)
            if abs(characteristic-gaussian) > float(bound)+1e-13:
                raise ArithmeticError("the rank-supported finite-copy bound failed")
            pair.append({"copies": copies, "characteristic": characteristic,
                "gaussian_error": abs(characteristic-gaussian),
                "exact_rank_supported_error_upper": bound})
        rows.append({"dimension": dimension, "copy_controls": pair})
    data = shared_fluctuation_modes()
    reference, compiled, actual = source.prepared_reference(), source.compiled_source(), prepare()
    raw_state = sum(float(weight)*(amplitude@amplitude.conj().T)
                    for weight, amplitude in zip(actual.weights, actual.successes, strict=True))
    raw_state /= np.trace(raw_state).real  # Explicit successful work marginal.
    basis, raw_clock = compiled["spectrum"]["basis"], compiled["raw"]["clock1"]
    guard = float(compiled["guards"]["clock1"])
    budget = error_budget()
    state_error = 4*float(budget["stopped_preparation_half_diamond_upper_under_guards"]
                         +ARITHMETIC_ALLOWANCE)/float(
        budget["ideal_success_probability"])
    normalization = []
    for address in source.PAIR_ADDRESSES:
        mode = data["modes"][address-1]
        omega = mode["gap"]  # The finite H clock, with declared b=1.
        q = basis@mode["Q"]@basis.T/math.sqrt(omega)
        raw_two_point = complex(np.trace(raw_state@raw_clock.conj().T@q@raw_clock@q))
        target = complex(mode["variance"]*math.cos(omega)/omega, -math.sin(omega)/(2*omega))
        norm_squared = 1/(2*mode["delta"]*omega)
        error = norm_squared*((1+guard)**2-1+state_error*(1+guard)**2
                              +float(ARITHMETIC_ALLOWANCE))
        if abs(raw_two_point-target) > error:
            raise ArithmeticError("the actual prepared clock violates its normalized field dictionary")
        normalization.append({"address": address, "omega_equals_prepared_gap_at_b1": omega,
            "raw_operator_norm_scaled_two_point":
                [raw_two_point.real/norm_squared, raw_two_point.imag/norm_squared],
            "operator_norm_squared": norm_squared,
            "ideal_ordered_two_point": [target.real, target.imag],
            "source_numerical_error_upper": error,
            "normalized_retarded_derivative_at_zero": 1,
            "unscaled_Q_retarded_derivative_at_zero_control": omega})
    schedules = []
    for level in (0, 1, 2):
        dimension = 240*4**level  # Candidate full-port family; no Gibbs encoder asserted.
        copies = dimension**2
        schedules.append(growing_fluctuation_bound(dimension, copies, (F(1, 4), F(1, 4)),
            gap_floor=F(1, 2), inverse_temperature=F(1),
            clock_gap_error=F(1, dimension), maximum_time=F(1),
            complete_preparation_error=F(1, 64*dimension**4),
            success_lower=F(1, 2*dimension), one_copy_word_error=F(1, 8*dimension**3),
            spectral_pair_error=F(1, dimension), cutoff_tail_radius=F(1, dimension)))
    return {"equal_gap_family_is_a_mathematical_countercontrol_not_a_native_encoder": True,
        "finite_family": rows, "gaussian_characteristic": gaussian,
        "N_equals_dimension_limit": math.exp(3*(math.cos(1/math.sqrt(2))-1)),
        "fixed_number_of_copies_has_gaussian_growing_dimension_limit": False,
        "finite_prepared_clock_normalization": normalization,
        "finite_clock_frequencies_identified_with_sqrt_native_L": False,
        "native_full_port_trace": native_trace_controls(),
        "native_joint_preparation_clock_field_contract": [
            native_growing_field_contract(level) for level in (0, 1, 2)],
        "sufficient_full_port_resource_schedule": schedules,
        "schedule_error_inputs_are_obligations_not_measured_source_rates": True,
        "angular_A_uses_native_constant_anchor_and_positive_modes": True,
        "positive_Dirichlet_product_B_requires_a_separate_zero_energy_anchor": True,
        "physical_zero_frequency_stationary_thermal_mode_constructed": False,
        "interacting_field_or_observed_spacetime_identified": False}


def prepared_gap_intervals():
    """Exact rational enclosures, independently of observations and libm.

    Positive coefficients make m(lambda)=sum w_n(1-lambda/120)^(2n)
    decreasing on [0,8].  For x=(m0-m)/(m0+m), log(m0/m)=2*atanh(x);
    the positive series has remainder at most 2*x^17/(17*(1-x*x)).
    """
    scale = 10**24
    lower = F(math.isqrt(5*scale*scale), scale)
    upper = lower+F(1, scale)
    eigenvalues = ((5-upper, 5-lower), (F(6), F(6)), (5+lower, 5+upper))
    weights, result = preparation.subordinated_filter_weights(), []
    mass0 = sum(weights)
    polynomial = lambda value: sum((w*(1-value/120)**(2*n)
                                   for n, w in enumerate(weights)), F())
    series = lambda value: 2*sum((value**(2*k+1)/F(2*k+1) for k in range(8)), F())
    for lo, hi in eigenvalues:
        mlo, mhi = polynomial(hi), polynomial(lo)
        xlo, xhi = (mass0-mhi)/(mass0+mhi), (mass0-mlo)/(mass0+mlo)
        if not 0 < xlo <= xhi < 1:
            raise ArithmeticError("the prepared positive-log enclosure lost its domain")
        result.append((series(xlo), series(xhi)+2*xhi**17/(17*(1-xhi*xhi))))
    return tuple(result)


def prepared_optical_comparison():
    """Retrospective test of nine direct assignments, with no retuning.

    The total ratio error and per-gap attachment error are explicit hard
    hypotheses. Published standard uncertainties and numerical compiler
    guards do not by themselves establish physical error bounds.
    """
    gaps = prepared_gap_intervals()  # Deliberately evaluated before data access.
    gap_error, ratio_error = F(1, 10**6), F(1, 10**12)
    extended = tuple((lo-gap_error, hi+gap_error) for lo, hi in gaps)
    if any(lo <= 0 for lo, _ in extended):
        raise ArithmeticError("the attachment error exceeds a positive clock gap")
    data = strict_json.load(ROOT/OPTICAL_DATA)
    expected = [(40, "1.05776938758748094", "0.00000000000000011"),
                (36, "1.05776646273518748", "0.00000000000000013")]
    if (data["schema"] != "oph.optical_absolute_ratios.v1"
            or [(row["isotope"], row["ratio"], row["standard_uncertainty"])
                for row in data["rows"]] != expected):
        raise ValueError("the pinned primary-table optical observations changed")
    assignments = []
    # Even same-gap assignments receive independent allowed attachment errors;
    # this enlarges the tested family instead of assuming their cancellation.
    for i, (lo1, hi1) in enumerate(extended, 1):
        for j, (lo2, hi2) in enumerate(extended, 1):
            assignments.append((i, j, lo1/hi2, hi1/lo2))
    comparisons = []
    for row in data["rows"]:
        measured = F(row["ratio"])
        separations = [max(F(), lo-measured-ratio_error, measured-ratio_error-hi)
                       for _, _, lo, hi in assignments]
        comparisons.append({"isotope": row["isotope"], "published_ratio": row["ratio"],
            "quoted_standard_uncertainty": row["standard_uncertainty"],
            "all_nine_direct_assignments_incompatible": all(value > 0 for value in separations),
            "minimum_absolute_ratio_separation_lower":
                F((min(separations)*10**12).__floor__(), 10**12)})
    return {"data": OPTICAL_DATA, "primary_url": data["primary_url"],
        "source_table": data["source_table"], "fixed_generator": "H=-log(M/m0)",
        "gap_enclosures": [[F((lo*10**18).__floor__(), 10**18),
                            F((hi*10**18).__ceil__(), 10**18)] for lo, hi in gaps],
        "declared_absolute_per_gap_attachment_error": gap_error,
        "declared_absolute_ratio_error": ratio_error,
        "numerically_checked_clock_operator_guard": source.compiled_source()["guards"]["clock1"],
        "assignments": [{"numerator_band": i, "denominator_band": j,
            "ratio_enclosure_with_attachment_error": [F((lo*10**12).__floor__(), 10**12),
                                                       F((hi*10**12).__ceil__(), 10**12)]}
                        for i, j, lo, hi in assignments],
        "comparisons": comparisons, "one_common_positive_frequency_scale": True,
        "sampler_gaps_or_species_scales_fitted_to_observations": False,
        "offsets_harmonics_phase_aliases_or_excited_to_excited_transitions_included": False,
        "standard_uncertainties_or_numerical_guards_are_physical_hard_bounds": False,
        "status": "retrospective conditional direct-gap attachment exclusion",
        "prospective_prediction_or_OPH_exclusion": False,
        "atomic_transition_dictionary_derived": False}


def run():
    budget = error_budget()
    compiled = source.compiled_source()
    rows, executions = [], {}
    for axis in ("X", "Y"):
        for intervention in (0, 1):
            saved = checkpoint(intervention=1) if (axis, intervention) == ("Y", 1) else None
            actual = execute(axis=axis, intervention=intervention, saved=saved)
            target = execute(axis=axis, intervention=intervention, ideal=True)
            rows.append(audit(actual, target))
            executions[axis, intervention] = actual
    # A source-specific nonzero phase signal on one *unconditional* native
    # record event: first tick loses the port load, second attempt also ticks.
    def contrast(history):
        law = probabilities(history)
        address = compiled["pair_addresses"][BAND]
        return sum(law["read", n, 1, 1, 1, 0]-law["read", n, 1, 1, 1, address]
                   for n in range(len(history.preparation.weights)))
    reference = compiled["reference"]
    mode = compiled["pair_addresses"][BAND]
    m0, mj = float(reference["m0"]), float(reference["masses"][mode])
    gap = reference["modular_gaps"][mode]
    analytic = (m0-mj)*math.sin(gap)/(12*90)
    signal = contrast(executions["X", 1])-contrast(executions["X", 0])
    allowance = 4*float(budget["complete_half_trace_upper_under_guards"])
    if analytic <= allowance or abs(signal-analytic) > allowance:
        raise ArithmeticError("the retained source clock/feedback signal is not resolved")
    # Execute the clock-off control through the same source preparation and
    # noncommuting feedback; do not infer its absence from a clock marginal.
    null_histories = []
    for ideal in (False, True):
        prepared, maps = prepare(ideal=ideal), _maps(ideal)
        outputs = {}
        for n, state in enumerate(prepared.successes):
            for control, _ in native_histories():
                kicked = maps[f"kick_plus_{BAND}"]@state if control[1] else state
                outputs[n, *control] = _array(maps[f"detector_x_{BAND}"]@kicked, (16, 12))
        null_histories.append(History("X", 1, prepared, outputs))
    audit(null_histories[0], null_histories[1])
    null_signal = contrast(null_histories[0])-contrast(executions["X", 0])
    if abs(null_signal) > allowance:
        raise ArithmeticError("source clock-off control falsely retains a phase signal")
    # Common H is not modular time of the retained request/work algebra:
    # request zero has a tracial success state, while H is nonconstant.
    request_zero_gap = 0.
    if gap <= .01 or request_zero_gap == gap:
        raise ArithmeticError("the retained-request modular countercontrol disappeared")
    groups = native_histories()
    native_json = json.dumps(groups, separators=(",", ":"))
    return preparation._json({
        "schema": "oph.prepared_clock.v1", "source_sha256": bindings(),
        "complete_sampler_intervals": sampler_contract(),
        "native_record_groups": [{"control": list(key), "tape_count": len(tapes)}
                                  for key, tapes in groups],
        "complete_native_tape_sha256": hashlib.sha256(native_json.encode()).hexdigest(),
        "executions": rows, "error_budget": budget, "resources": resources(),
        "shared_state_field": shared_state_field(),
        "growing_shared_state_field": growing_field_controls(),
        "prepared_positive_log_limit": source.prepared_log_limit_bound(filter_operator_error=FILTER_GUARD),
        "positive_log_refinement_schedule": [
            {key: value for key, value in source.positive_log_refinement_schedule(
                level, source_preserves_band=True).items()
             if key in ("level", "rate", "heat_step", "jump_cap", "count_cap", "sampling_bits",
                        "maximum_filter_requests", "one_filter_operator_error",
                        "raw_generator_to_b_sqrt_D_upper", "required_phase_word_error_upper",
                        "growing_sampler_evaluated", "bounded_band_then_spectral_cutoff_limit",
                        "raw_implementation_preserves_band_or_entire_encoded_domain_is_in_band")}
            for level in (4, 6, 8)],
        "prepared_gap_optical_comparison": prepared_optical_comparison(),
        "source_words": {"sha256": compiled["word_sha256"],
                         "pulse_counts": compiled["pulse_counts"],
                         "full_register_guards": compiled["guards"]},
        "clock_field_response": {
            "event": "first incident attempt loses port10 load, second attempt is incident",
            "unconditional_native_event_probability": F(1, 90),
            "observable": "X detector address0 minus selected lambda6 address, feedback on minus off",
            "analytic_signal": analytic, "executed_signal": signal,
            "resolved_signal_lower_under_guards": analytic-allowance,
            "executed_clock_off_signal": null_signal,
            "zero_signal_when_clock_or_noncommuting_feedback_is_absent": True,
            "formula": "(m0-m6)*sin(log(m0/m6))/(12*90)"},
        "state_and_clock": {
            "prepared_state": "M/Tr(M), M=sum_n w_n(I-L/120)^(2n)",
            "same_register_clock": "exp(-iH), H=-log(M/m0)",
            "modular_identification": "sigma_s=alpha_(-s) on successful work B(C^12)",
            "reference_and_sampler_request_are_traced_for_this_algebra": True,
            "sampler_request_retained_in_complete_evidence": True,
            "request_zero_modular_gap": request_zero_gap,
            "work_marginal_lambda6_modular_gap": gap,
            "one_clock_is_modular_for_every_retained_request_block": False,
            "faithful_modular_state_on_padding_or_reference_algebra": False,
            "sqrt_L_comparison_clock_compiled": True,
            "finite_prepared_reference_is_exact_sqrt_L_Gibbs_state": False},
        "scope": {
            "one_source_prepared_register_through_clock_feedback_and_readout": True,
            "complete_failed_preparation_and_native_histories_retained": True,
            "same_native_carrier_continues_after_sampler": True,
            "raw_arrays_normalized_or_projected_to_ideal_unitaries": False,
            "source_preparation_filter_and_clock_are_supplied_observer_services": True,
            "growing_product_field_encoder_or_original_scalar_benchmark_completed": False,
            "physical_time_or_empirical_field_dictionary_identified": False}})


def check_receipt(recorded, actual):
    """Strict replay additionally enforces the producer's scientific domains."""
    from ..observer.scalar_response.experiment import compare_receipt
    compare_receipt(recorded, actual)
    budget = actual["error_budget"]
    allowance = float(F(budget["complete_half_trace_upper_under_guards"]))
    filter_guard = float(F(budget["filter_operator_guard"]))
    residual = recorded["error_budget"]["filter_complete_Frobenius_residual"]
    if not 0 <= residual <= filter_guard:
        raise ValueError("stored filter norm violates its scientific guard")
    for row in recorded["executions"]:
        if (not 0 <= row["complete_quantum_half_trace_upper_diagnostic"] <= allowance
                or abs(row["complete_raw_mass_diagnostic"]-1) > 2*allowance):
            raise ValueError("stored complete quantum output violates its scientific guard")
    signal = recorded["clock_field_response"]
    if (signal["resolved_signal_lower_under_guards"] <= 0
            or abs(signal["executed_signal"]-signal["analytic_signal"]) > 4*allowance
            or abs(signal["executed_clock_off_signal"]) > 4*allowance):
        raise ValueError("stored clock signal violates its scientific guard")
    field = recorded["shared_state_field"]
    if (field["finite_copy_error_upper"] < 0
            or field["one_copy_source_estimator_error_upper_under_guards"] < 0
            or field["resolved_ordered_phase_lower_under_guards"] <= 0):
        raise ValueError("stored fluctuation bound violates its scientific guard")
    # The Nth power is a derived diagnostic: tiny permissible backend drift
    # in one-copy amplitudes is amplified by N. Retain its inputs and check
    # the scientifically propagated bound, without assigning it an unrelated
    # fixed absolute replay tolerance or weakening any other field's check.
    estimator = complex(*field["raw_pointer_estimator"])
    one_copy = complex(*field["one_copy_characteristic"])
    gaussian = complex(*field["gaussian_characteristic"])
    if (abs(estimator-one_copy) > field["one_copy_source_estimator_error_upper_under_guards"]
            or abs(estimator**field["angle_copy_parameter"]-gaussian)
            > field["combined_gaussian_error_upper_under_guards"]):
        raise ValueError("stored pointer estimator violates its propagated scientific guard")
    for row in field["executions"]:
        bound = float(F(row["complete_half_trace_upper_under_guards"]))
        if (not 0 <= row["complete_quantum_half_trace_upper_diagnostic"] <= bound
                or abs(row["complete_raw_mass_diagnostic"]-1) > 2*bound):
            raise ValueError("stored pointer outputs violate their scientific guard")
    growing = recorded["growing_shared_state_field"]
    for row in growing["finite_prepared_clock_normalization"]:
        raw = complex(*row["raw_operator_norm_scaled_two_point"])*row["operator_norm_squared"]
        target = complex(*row["ideal_ordered_two_point"])
        if (row["source_numerical_error_upper"] < 0
                or abs(raw-target) > row["source_numerical_error_upper"]):
            raise ValueError("stored normalized field violates its scientific guard")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path); group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        check_receipt(strict_json.load(args.check), result)
        print("Prepared reference, native clock, feedback and complete outputs replayed.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
