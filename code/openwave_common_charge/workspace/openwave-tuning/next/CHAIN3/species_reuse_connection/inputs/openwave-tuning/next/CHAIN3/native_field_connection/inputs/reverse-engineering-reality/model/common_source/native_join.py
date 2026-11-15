"""Native repair records drive and read the same interacting quantum history.

The first native current selects an emitted source pulse before quantum
evolution. Both quantum records and each unnormalized terminal work state are
retained. Only the final scalar bit selects the second native seam. Native
integer repair is the existing declared primitive; this module does not claim
to compile the entire M78 channel into the M6 proper-code processor.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
import math
from numbers import Real

from model.observer import native
from model.observer.p_response.experiment import decoder_seam, controlled_native_step
from . import experiment

CURRENT_VALUES = (-1, 0, 1)
ALL_STATES = native.states()
BINARY_STATES = native.states(binary=True)
STATE_INDEX = {state: i for i, state in enumerate(ALL_STATES)}
TOLERANCE = 2e-10


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def _parameters(d, steps, coupling, feedback_enabled):
    if type(d) is not int or d not in (2, 4):
        raise ValueError("Executed scalar dimension must be 2 or 4")
    if type(steps) is not int or steps < 1:
        raise ValueError("A positive integer split count is required")
    if type(feedback_enabled) is not bool:
        raise ValueError("Quantum feedback switch must be Boolean")
    return d, steps, _number(coupling, "coupling"), feedback_enabled


def actuator_seam(final_bit):
    """The native actuator's only decision input is the committed scalar bit."""
    if type(final_bit) is not int or final_bit not in (0, 1):
        raise ValueError("A committed binary scalar record is required")
    return decoder_seam(2*final_bit-1)


def native_continuation(loads, final_bit, coin):
    """Apply the existing bit-only controller to the actual native work state."""
    actuator_seam(final_bit)
    if (type(loads) is not tuple or len(loads) != 12
            or any(type(value) is not int for value in loads) or loads not in STATE_INDEX):
        raise ValueError("The actual native total-two source state is required")
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("A binary native placement coin is required")
    output, retained, record = controlled_native_step(native, loads, 2*final_bit-1, coin)
    if retained != 2*final_bit-1:
        raise AssertionError("The native continuation altered its committed bit")
    return output, record


@lru_cache(maxsize=48)
def _quantum_history(d, steps, coupling, feedback_enabled, current):
    return experiment.history(d, steps, coupling=coupling,
                              feedback_enabled=feedback_enabled, current=current)


def _validated_histories(histories):
    """Bind every probability to a retained vector on one declared work space."""
    if type(histories) is not dict or set(histories) != set(CURRENT_VALUES):
        raise ValueError("One complete quantum history per native current is required")
    first = histories[0]
    signature = (first["dimension"], first["steps_per_segment"],
                 first["coupling"], first["feedback_enabled"])
    d, _, _, _ = _parameters(*signature)
    indexed = {}
    max_error = 0.
    for current in CURRENT_VALUES:
        history = histories[current]
        if (type(history.get("native_current")) is not int or history["native_current"] != current
                or (history["dimension"], history["steps_per_segment"],
                    history["coupling"], history["feedback_enabled"]) != signature):
            raise ValueError("Quantum histories must share one model and current attachment")
        if history.get("all_pass") is not True:
            raise ValueError("Quantum history execution checks did not pass")
        if type(history["branches"]) is not list or len(history["branches"]) != 4:
            raise ValueError("Every quantum history must retain four outcome branches")
        seen = {}
        for index, branch in enumerate(history["branches"]):
            record = branch["record"]
            if (type(record) is not list or len(record) != 2
                    or any(type(bit) is not int or bit not in (0, 1) for bit in record)):
                raise ValueError("Quantum records must be two binary symbols")
            key = tuple(record)
            if key in seen:
                raise ValueError("Duplicate quantum branch")
            state = branch["unnormalized_terminal_state"]
            if type(state) is not list or len(state) != 4*d*d:
                raise ValueError("Terminal state does not belong to the declared common carrier")
            squared = []
            for amplitude in state:
                if type(amplitude) is not list or len(amplitude) != 2:
                    raise ValueError("Terminal amplitudes need real and imaginary components")
                real = _number(amplitude[0], "state real part")
                imaginary = _number(amplitude[1], "state imaginary part")
                squared.append(real*real+imaginary*imaginary)
            probability = math.fsum(squared)
            reported = _number(branch["probability"], "branch probability")
            if reported < 0. or not math.isfinite(probability):
                raise ValueError("Invalid quantum branch probability")
            discrepancy = abs(probability-reported)
            max_error = max(max_error, discrepancy)
            if discrepancy > TOLERANCE:
                raise ValueError("Quantum probability is detached from its retained terminal state")
            seen[key] = (index, probability)
        if set(seen) != {(0, 0), (0, 1), (1, 0), (1, 1)}:
            raise ValueError("Incomplete quantum outcome menu")
        if abs(math.fsum(probability for _, probability in seen.values())-1) > TOLERANCE:
            raise ValueError("Quantum branches are not a complete normalized history")
        indexed[current] = seen
    return signature, indexed, max_error


def _summary(rows):
    return {
        "history_probability_sum": math.fsum(row["history_probability"] for row in rows),
        "mean_first_native_current": math.fsum(
            row["history_probability"]*row["first_native_record"][1] for row in rows),
        "mean_second_native_current": math.fsum(
            row["history_probability"]*row["second_native_record"][1] for row in rows),
        "mean_final_native_port0_load": math.fsum(
            row["history_probability"]*row["final_native_loads"][0] for row in rows),
        "final_scalar_bit_one_probability": math.fsum(
            row["history_probability"] for row in rows if row["quantum_record"][1] == 1),
        "first_gauge_bit_one_probability": math.fsum(
            row["history_probability"] for row in rows if row["quantum_record"][0] == 1),
        "signed_scalar_bit_times_first_native_current": math.fsum(
            row["history_probability"]*(2*row["quantum_record"][1]-1)
            *row["first_native_record"][1] for row in rows),
    }


def _control(rows, plus_probability):
    """Independent actuator, retaining every original quantum state/record."""
    means = []
    loads = []
    masses = []
    choices = ((0, 1-plus_probability), (1, plus_probability))
    for row in rows:
        source = ALL_STATES[row["after_first_source"]]
        coin = row["second_native_record"][0]
        for bit, probability in choices:
            final, second = native_continuation(source, bit, coin)
            weight = row["history_probability"]*probability
            masses.append(weight)
            means.append(weight*second[1])
            loads.append(weight*final[0])
    return {"actuator_bit_one_probability": plus_probability,
            "history_probability_sum": math.fsum(masses),
            "mean_second_native_current": math.fsum(means),
            "mean_final_native_port0_load": math.fsum(loads),
            "original_quantum_records_and_terminal_states_retained": True,
            "actuator_independent_of_native_state_and_quantum_records": True}


def _temporal_transfer(histories):
    key = "ideal_finite_continuous_comparison"
    if any(key not in histories[current] for current in CURRENT_VALUES):
        return {"status": "not_supplied", "joint_error_bound_claimed": False}
    bounds = {current: Fraction(histories[current][key][
        "complete_instrument_half_diamond_error_upper"]) for current in CURRENT_VALUES}
    if any(not 0 <= error <= 1 for error in bounds.values()):
        raise ValueError("Invalid inherited finite instrument error bound")
    error = max(bounds.values())
    return {"status": "conditional_bound_transferred", "joint_error_bound_claimed": True,
            "complete_joint_instrument_half_diamond_error_upper": str(error),
            "complete_classical_history_tv_error_upper": str(error),
            "mean_second_native_current_error_upper": str(2*error),
            "per_current_quantum_bounds": {str(current): str(bound)
                                           for current, bound in bounds.items()},
            "scope": "same finite occupation code and interventions; exact native instruments and complete record-dependent composition contract the imported ideal temporal bound; numerical, physical service, occupation and spatial errors remain separate"}


def compose_histories(histories):
    """Join complete quantum operators to all actual native branches.

    The terminal state for each row is its referenced quantum vector tensored
    with the listed final native basis state and complete classical record.
    Its density has exact prefactor 1/264. Arrays are retained once per
    (current,a,b), rather than copied for 1056 histories.
    """
    signature, indexed, state_error = _validated_histories(histories)
    rows = []
    first_mass = defaultdict(list)
    for initial in BINARY_STATES:
        for first_coin in (0, 1):
            after_first, first = native.native(initial, native.E, first_coin)
            current = first[1]
            for record in ((0, 0), (0, 1), (1, 0), (1, 1)):
                branch_index, probability = indexed[current][record]
                bit = record[1]
                for second_coin in (0, 1):
                    final, second = native_continuation(after_first, bit, second_coin)
                    weight = probability/264
                    first_mass[first].append(weight)
                    rows.append({
                        "initial_source": STATE_INDEX[initial],
                        "after_first_source": STATE_INDEX[after_first],
                        "first_native_record": list(first),
                        "quantum_record": list(record),
                        "actuator_input": bit,
                        "second_seam": list(actuator_seam(bit)),
                        "second_native_record": list(second),
                        "final_source": STATE_INDEX[final],
                        "final_native_loads": list(final),
                        "quantum_branch_reference": {"current": current, "branch_index": branch_index,
                                                     "record": list(record)},
                        "terminal_density_prefactor": "1/264",
                        "history_probability": weight,
                    })
    summary = _summary(rows)
    fixed = _control(rows, 0.)
    independent = _control(rows, summary["final_scalar_bit_one_probability"])
    summary["feedback_current_contrast_against_independent_same_marginal"] = (
        summary["mean_second_native_current"]-independent["mean_second_native_current"])
    summary["feedback_port0_contrast_against_independent_same_marginal"] = (
        summary["mean_final_native_port0_load"]-independent["mean_final_native_port0_load"])
    prefix = [{"first_native_record": list(record), "probability": math.fsum(weights)}
              for record, weights in sorted(first_mass.items())]
    expected = {(0, 0): 14/33, (1, 0): 14/33, (0, 1): 5/66, (1, -1): 5/66}
    checks = {
        "all_1056_histories_retained": len(rows) == 1056,
        "complete_joint_history_normalized": abs(summary["history_probability_sum"]-1) < TOLERANCE,
        "old_native_record_marginal_preserved": set(first_mass) == set(expected)
            and all(abs(math.fsum(first_mass[record])-probability) < TOLERANCE
                    for record, probability in expected.items()),
        "native_load_conserved_in_every_branch": all(sum(row["final_native_loads"]) == 2 for row in rows),
        "controller_uses_only_retained_final_bit": all(
            tuple(row["second_seam"]) == actuator_seam(row["quantum_record"][1]) for row in rows),
        "probabilities_bound_to_unnormalized_terminal_states": state_error < TOLERANCE,
        "independent_control_has_same_actuator_marginal": (
            independent["actuator_bit_one_probability"] == summary["final_scalar_bit_one_probability"]),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    d, steps, coupling, enabled = signature
    return {
        "schema": "oph.common_source.native_quantum_join.v1",
        "dimension": d, "steps_per_segment": steps, "coupling": coupling,
        "quantum_feedback_enabled": enabled,
        "quantum_histories": [histories[current] for current in CURRENT_VALUES],
        "paths": rows, "summary": summary, "first_native_record_marginal": prefix,
        "controls": {"fixed_native_seam_E": fixed, "independent_bit_same_marginal": independent},
        "finite_temporal_transfer": _temporal_transfer(histories),
        "maximum_probability_state_norm_discrepancy": state_error,
        "checks": checks, "all_pass": True,
        "scope": {
            "joint_instrument": "native record -> emitted current-controlled gauge pulse -> common interacting quantum history -> final-bit-only native continuation",
            "quantum_state": "every terminal gauge/scalar state is retained once and referenced by every complete native history, with density prefactor 1/264",
            "native_primitive": "the existing full M78 integer repair and record channel is a declared primitive; no M6 compilation of that entire channel is claimed",
            "reference": "one prepared joint experiment, not a proof that its state is the unconstrained A3 tracial optimizer",
            "time": "native event order and the supplied Hamiltonian model coordinate; no physical duration calibration",
            "errors": "the supplied ideal finite temporal bound transfers through native composition; floating source execution and branch-state checks are diagnostics, with no inherited free-scalar occupation or spatial continuum certificate",
            "prediction": "conditional finite model law; no empirical validation or derived physical coupling",
        },
    }


def run(d=4, steps=512, coupling=.25, feedback_enabled=True):
    d, steps, coupling, feedback_enabled = _parameters(d, steps, coupling, feedback_enabled)
    histories = {current: deepcopy(_quantum_history(d, steps, coupling, feedback_enabled, current))
                 for current in CURRENT_VALUES}
    return compose_histories(histories)
