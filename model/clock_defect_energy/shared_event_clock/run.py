"""Continuing C/E source clocks with four common event checkpoints.

One admissible preferred-seam feedback path keeps a two-state coherence
through three FIFO exports. Every alternative native record stops as a
retained outcome. A selected record requests a freshly certified parameter-1
modular source word and one protected density-wave step. Independent copies
terminate at four declared checkpoints for source analyzer/Ramsey reads.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from model.completion import observer
from model.completion.family import FamilyParameters
from model.native_modular_reference import common_caps, history_clock, full_flow
from model.clock_defect_energy import clock_services, continuation, ramsey_service, joint_clock
from model.clock_defect_energy.joins_contract import check_receipt, source_bindings
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE / "campaign.json"
OUTPUT = HERE / "results.json"
SEAMS = ((0, 7), (6, 7), (2, 6), (2, 3))
COINS_LABELS = ((0, 1), (1, -1), (1, -1), (0, 1))
INITIAL_PAIR = (55, 65)
PARAMETER = F(1)
SOURCE_GATE = F(1, 10**10)
CODES = np.array(full_flow.support()[0])
OUTSIDE = np.array([i for i in range(512) if i not in set(CODES)])


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _predecessor(name):
    receipt = strict_json.load(HERE.parent / name)
    if source_bindings(bound_pins(receipt["source_sha256"])) != bound_pins(receipt["source_sha256"]):
        raise ValueError(f"{name} source closure changed")
    return receipt


def contract():
    value = strict_json.load(CAMPAIGN)
    expected = {
        "schema": "oph.clock_defect_energy.shared_event_clock_campaign.v1",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "frozen_before_execution": True,
        "whole_program_frozen": False,
        "native_family_changed": False,
        "integrated_gate_accepted": False,
        "physical_comparison_count_added": 0,
        "new_empirical_comparisons": 0,
        "native_laws": ["C", "E"],
        "preferred_seam_intervention_by_attempt": [list(x) for x in SEAMS],
        "selected_feedback_coin_label_by_attempt": [list(x) for x in COINS_LABELS],
        "selected_initial_binary_indices": list(INITIAL_PAIR),
        "clock_parameter_per_selected_event": "1",
        "terminal_horizons": [1, 2, 3, 4],
        "affine_anchor_events": [1, 2],
        "third_event_check": 4,
        "extra_unfitted_event_check": 3,
        "individual_shot_clock_readout_claimed": False,
        "source_selected_geometry_claimed": False,
        "physical_seconds_calibrated": False,
    }
    if any(value.get(k) != v for k, v in expected.items()):
        raise ValueError("shared-event frozen contract changed")
    for key, name in (("clock_count_result_sha256", "clock_count_inference_results.json"),
                      ("clock_service_result_sha256", "clock_services_results.json"),
                      ("ramsey_result_sha256", "ramsey_service_results.json"),
                      ("continuation_result_sha256", "continuation_results.json")):
        if value[key] != hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest():
            raise ValueError("predecessor result bytes changed")
    return value


def _outward(interval, places=12):
    scale = 10**places
    return [str(F(math.floor(interval[0]*scale), scale)),
            str(F(math.ceil(interval[1]*scale), scale))]


def _float_box(value):
    if not math.isfinite(value) or value < -1e-9:
        raise ArithmeticError("invalid source probability diagnostic")
    scale = 10**8
    center = round(value*scale)
    return [str(F(max(0, center-1), scale)), str(F(center+1, scale))]


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _sub(a, b):
    return a[0]-b[1], a[1]-b[0]


def _mul(a, b):
    values = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(values), max(values)


def _div(a, b):
    if b[0] <= 0:
        raise ArithmeticError("nonpositive clock anchor separation")
    return _mul(a, (1/b[1], 1/b[0]))


def _q_clock(n, parameter=PARAMETER):
    re, _, radius = common_caps.phase_enclosure(parameter*n, F(8, 5))
    return (1-re-radius)/2, (1-re+radius)/2


def _q_field(n):
    output = joint_clock.wave.ramsey_outputs(joint_clock.FIELD, n, -1)
    return joint_clock.wave.norm_squared(output[1])/330


def _cross(clock, field, third):
    first, second = 1, 2
    return _sub(_mul(_sub(field[third], field[first]), _sub(clock[second], clock[first])),
                _mul(_sub(field[second], field[first]), _sub(clock[third], clock[first])))


def exact_comparison(bounds):
    clock = {n: _q_clock(n) for n in range(1, 5)}
    field = {n: (v, v) for n in range(1, 5) for v in (_q_field(n),)}
    if any(clock[n][1] >= clock[n+1][0] or field[n][1] >= field[n+1][0]
           for n in range(1, 4)):
        raise ArithmeticError("observed mean readouts do not give strict clocks")
    anchor_gap = _sub(clock[2], clock[1])
    rate = _div(_sub(field[2], field[1]), anchor_gap)
    offset = _sub(field[1], _mul(rate, clock[1]))
    third = _cross(clock, field, 4)
    extra = _cross(clock, field, 3)
    if third[0] <= 0 or extra[0] <= 0:
        raise ArithmeticError("unfitted third-event affine test did not fail")
    residual = _div(third, anchor_gap)
    clock_s4 = {n: _q_clock(n, F(4)) for n in range(1, 5)}
    if clock_s4[2][0] <= clock_s4[3][1]:
        raise ArithmeticError("the inherited parameter-four alias control disappeared")
    robust_clock, robust_field = {}, {}
    for n in range(1, 5):
        eta = bounds[n]
        robust_clock[n] = (max(F(), clock[n][0]-eta), min(F(1), clock[n][1]+eta))
        robust_field[n] = (max(F(), field[n][0]-eta), min(F(1), field[n][1]+eta))
    robust = _cross(robust_clock, robust_field, 4)
    robust_residual = _div(robust, _sub(robust_clock[2], robust_clock[1]))
    if (robust[0] <= 0 or any(robust_clock[n][1] >= robust_clock[n+1][0]
                             or robust_field[n][1] >= robust_field[n+1][0]
                             for n in range(1, 4))):
        raise ArithmeticError("source bounds cannot certify strict clocks and nonaffinity")
    return {
        "record_order": "retained selected native prefix 1<2<3<4",
        "shared_event_map": "clock and field terminal checkpoint n both map to the same complete native prefix and FIFO export history",
        "readout_scope": "conditional ensemble terminal probabilities, not individual stochastic shot results",
        "clock_first_code_probabilities": {str(n): _outward(clock[n]) for n in clock},
        "field_Ramsey_minus_probabilities": {str(n): str(field[n][0]) for n in field},
        "both_readouts_strictly_increasing_on_four_events": True,
        "positive_affine_two_event_rate": _outward(rate),
        "positive_affine_two_event_offset": _outward(offset),
        "third_event": 4,
        "third_event_cross_product_nonzero": _outward(third),
        "third_event_field_probability_above_affine_prediction": _outward(residual),
        "extra_unfitted_event_3_cross_product_nonzero": _outward(extra),
        "source_bound_robust_third_cross_product": _outward(robust),
        "source_bound_robust_field_probability_above_affine_prediction": _outward(robust_residual),
        "source_bound_robust_strict_readouts": True,
        "parameter_four_control_clock_probabilities": {str(n): _outward(clock_s4[n]) for n in clock_s4},
        "parameter_four_control_loses_monotonicity_between_events_2_and_3": True,
        "response_inverse_scope": "disjoint source-certified mean bands identify n in {1,2,3,4}; the lookup is model-conditioned and is not a physical clock calibration",
        "raw_affine_test_passed": False,
        "calibrated_count_lookup_affine_by_construction": True,
        "formal_source_specific_Lean_instance_claimed": False,
        "physical_affine_timelike_atlas_supplied": False,
    }


@lru_cache(maxsize=4)
def _word(core):
    core = continuation._core(core)
    ratios = clock_services.ratios(core)
    tape = common_caps.diagonal_word(PARAMETER, ratios)
    array = full_flow.compiler.execute(tape, 9)
    bound = common_caps.certify_diagonal(PARAMETER, ratios, array)
    if bound > SOURCE_GATE:
        raise ArithmeticError("parameter-one source word exceeds gate")
    return array, tape, bound


@lru_cache(maxsize=6)
def _rotation(pair, sign):
    tape = history_clock.rotation_word(pair, sign)
    array = history_clock.rotation(pair, sign, True)
    bound = history_clock.certify_rotation(pair, sign, array)
    if bound > SOURCE_GATE:
        raise ArithmeticError("pair source rotation exceeds gate")
    return array, tape, bound


def _native(vector, branch):
    result = np.zeros_like(vector)
    if branch.weight:
        amplitude = math.sqrt(float(branch.weight))
        for i, j in branch.pairs:
            if i in continuation.FULL_TO_BINARY:
                result[continuation.FULL_TO_BINARY[j]] += amplitude*vector[continuation.FULL_TO_BINARY[i]]
    return result


def _selected(rows, seam, coin, label):
    hits = [b for b in rows if b.record.kind == "repair"
            and b.record.selector == "feedback" and b.record.seam == seam
            and b.record.coin == coin and b.record.label == label
            and b.weight > 0]
    if len(hits) != 1 or hits[0].weight != F(1, 4):
        raise ArithmeticError("the selected native feedback branch changed")
    return hits[0]


def _budget(n):
    clock_delta = (1+SOURCE_GATE)**(n+2)-1  # one prep, n words, one analyzer
    clock_halftrace = clock_delta*(2+clock_delta)/2
    field = F(ramsey_service.budget(n)["complete_halftrace_error_upper"])
    full = clock_halftrace+(1+clock_delta)**2*field
    selected = F(1, 4**n)
    if full >= selected:
        raise ArithmeticError("conditional source service bound is vacuous")
    conditional = 2*full/(selected-full)
    return {"clock_operator_product_upper": str(clock_delta),
            "clock_complete_halftrace_upper": str(clock_halftrace),
            "field_complete_halftrace_upper": str(field),
            "full_flagged_halftrace_upper": str(full),
            "selected_ideal_path_mass": str(selected),
            "conditional_joint_read_TV_upper": str(conditional),
            "physical_pulse_or_transfer_errors_included": False}


def _source_field():
    reads, waits = {}, {}
    for n in range(5):
        waiting = ramsey_service.field_columns(n, read=False)
        waits[n] = float(np.vdot(waiting, waiting).real)
        if n:
            plus, minus = ramsey_service.field_columns(n)
            reads[n] = (float(np.vdot(plus, plus).real),
                        float(np.vdot(minus, minus).real))
    return reads, waits


def source_law(law, field_reads, field_waits):
    params = FamilyParameters(c_weight=F(int(law == "C")))
    prep, prep_tape, prep_bound = _rotation(INITIAL_PAIR, +1)
    input_address = np.zeros(512, complex)
    input_address[CODES[INITIAL_PAIR[1]]] = 1
    output = prep@input_address
    prep_abort = float(np.vdot(output[OUTSIDE], output[OUTSIDE]).real)
    active = output[CODES]
    initial_reference_mass = float(np.vdot(active, active).real)
    memory, core, pair = (), (0, 1), INITIAL_PAIR
    exports, records = [], []
    stopped = prep_abort*field_waits[0]
    semantic = positive = clock_aborts = 0
    stopped_keys = [(0, "source_preparation_export_abort")]
    rows_out = []
    for position, (seam, (coin, label)) in enumerate(zip(SEAMS, COINS_LABELS, strict=True), 1):
        menu = observer.branches(memory, params, preferred=seam)
        selected = _selected(menu, seam, coin, label)
        native_total = 0.
        selected_native = None
        for branch in menu:
            value = _native(active, branch)
            mass = float(np.vdot(value, value).real)
            native_total += mass
            semantic += 1
            if branch is selected:
                selected_native = value
            else:
                stopped += mass*field_waits[position-1]
                positive += int(mass > 1e-20)
                stopped_keys.append((position, branch.record.__dict__,
                                     [r.__dict__ for r in branch.memory],
                                     branch.exported.__dict__ if branch.exported else None))
        if abs(native_total-float(np.vdot(active, active).real)) > 1e-10:
            raise ArithmeticError("a fixed preferred native menu lost work mass")
        if selected_native is None:
            raise ArithmeticError("selected native record absent")
        mapping = {continuation.FULL_TO_BINARY[i]: continuation.FULL_TO_BINARY[j]
                   for i, j in selected.pairs if i in continuation.FULL_TO_BINARY}
        if pair[0] not in mapping or pair[1] not in mapping:
            raise ArithmeticError("continuing coherent pair left selected branch")
        pair = mapping[pair[0]], mapping[pair[1]]
        core = continuation.update_core(core, selected.record)
        memory = selected.memory
        if selected.exported is not None:
            exports.append(selected.exported)
        records.append(selected.record.__dict__)
        if clock_services.ratios(core)[CODES[pair[0]]] / clock_services.ratios(core)[CODES[pair[1]]] != F(5, 8):
            raise ArithmeticError("transported coherence changed reference ratio")
        word, tape, word_bound = _word(core)
        advanced = word[:, CODES]@selected_native
        clock_abort = float(np.vdot(advanced[OUTSIDE], advanced[OUTSIDE]).real)
        clock_aborts += 1
        stopped += clock_abort*field_waits[position-1]
        active = advanced[CODES]
        analyzer, analyzer_tape, analyzer_bound = _rotation(pair, -1)
        analyzed = analyzer[:, CODES]@active
        analyzer_abort = float(np.vdot(analyzed[OUTSIDE], analyzed[OUTSIDE]).real)
        first = float(abs(analyzed[CODES[pair[0]]])**2)
        other = float(np.vdot(analyzed[CODES], analyzed[CODES]).real)-first
        if min(first, other) < -1e-10:
            raise ArithmeticError("terminal analyzer has negative read mass")
        field_plus, field_minus = field_reads[position]
        terminal = (first+other)*(field_plus+field_minus)
        complete = stopped + analyzer_abort*field_waits[position] + terminal
        if abs(complete-1) > 2e-7:
            raise ArithmeticError("complete prefix source instrument lost mass")
        source_clock = first/(first+other)
        source_field = field_minus/(field_plus+field_minus)
        qlo, qhi = _q_clock(position)
        ideal_field = _q_field(position)
        conditional = F(_budget(position)["conditional_joint_read_TV_upper"])
        if (not float(qlo-conditional)-1e-10 <= source_clock <= float(qhi+conditional)+1e-10
                or abs(source_field-float(ideal_field)) > float(conditional)+1e-10):
            raise ArithmeticError("source checkpoint readouts exceed complete bounds")
        event = {"law": law, "native_prefix": records.copy(),
                 "exported_records": [r.__dict__ for r in exports],
                 "retained_fifo": [r.__dict__ for r in memory]}
        rows_out.append({
            "checkpoint": position, "shared_event_sha256": _digest(event),
            "transported_core": list(core), "source_pair_indices": list(pair),
            "native_prefix": records.copy(), "exported_record_count": len(exports),
            "retained_fifo_size": len(memory),
            "semantic_native_slots_through_checkpoint": semantic,
            "positive_stopped_native_outcomes_through_checkpoint": positive,
            "all_stopped_native_record_keys_sha256": _digest(stopped_keys),
            "clock_export_abort_slots_through_checkpoint": clock_aborts,
            "source_preparation_export_abort_slot": True,
            "terminal_analyzer_export_abort_slot": True,
            "terminal_joint_clock_field_read_slots": 4,
            "source_selected_path_mass_diagnostic": _float_box(float(np.vdot(active, active).real)),
            "source_complete_mass_diagnostic": _float_box(complete),
            "source_clock_first_code_conditional_diagnostic": _float_box(source_clock),
            "source_field_Ramsey_minus_conditional_diagnostic": _float_box(source_field),
            "source_clock_word_sha256": hashlib.sha256(full_flow.corner.canonical(tape).encode()).hexdigest(),
            "source_clock_word_pulses": len(tape),
            "source_analyzer_word_pulses": len(analyzer_tape),
            "source_parameter_one_word_operator_bound_below": str(SOURCE_GATE),
            "source_analyzer_operator_bound_below": str(SOURCE_GATE),
            "clock_source_pulses_this_horizon": len(prep_tape)+position*len(tape)+len(analyzer_tape),
            "field_source_pulses_this_horizon": ramsey_service.budget(position)["source_pulses"],
            "bounds": _budget(position),
        })
    return {"law": law, "fixed_preferred_schedule": [list(x) for x in SEAMS],
            "selected_feedback_records": records,
            "source_preparation_word_pulses": len(prep_tape),
            "source_preparation_operator_bound_below": str(SOURCE_GATE),
            "pure_source_probe_is_faithful_reference_preparation": False,
            "initial_prepared_native_mass_diagnostic": _float_box(initial_reference_mass),
            "checkpoints": rows_out}


def run():
    campaign = contract()
    for name in ("clock_count_inference_results.json", "clock_services_results.json",
                 "ramsey_service_results.json", "continuation_results.json"):
        _predecessor(name)
    bounds = {n: F(_budget(n)["conditional_joint_read_TV_upper"]) for n in range(1, 5)}
    exact = exact_comparison(bounds)
    field_reads, field_waits = _source_field()
    laws = [source_law(law, field_reads, field_waits) for law in ("C", "E")]
    if ([row["shared_event_sha256"] for row in laws[0]["checkpoints"]] ==
            [row["shared_event_sha256"] for row in laws[1]["checkpoints"]]):
        raise ArithmeticError("C/E complete-record event maps were conflated")
    paths = {
        "model/clock_defect_energy/shared_event_clock/__init__.py",
        "model/clock_defect_energy/shared_event_clock/run.py",
        "model/clock_defect_energy/shared_event_clock/test_shared_event_clock.py",
        "model/clock_defect_energy/shared_event_clock/README.md",
        "model/clock_defect_energy/shared_event_clock/campaign.json",
        "model/clock_defect_energy/clock_count_inference_results.json",
        "model/clock_defect_energy/clock_services_results.json",
        "model/clock_defect_energy/ramsey_service_results.json",
        "model/clock_defect_energy/continuation_results.json",
        "Lean/Time/ClockReadout.lean", "Lean/Time/ClockComparison.lean",
        "Lean/Dynamics/PrivateInner.lean", "Lean/Dynamics/ChoiCPTP.lean",
    }
    return {"schema": "oph.clock_defect_energy.shared_event_clock.v1",
            "campaign_id": campaign["campaign_id"],
            "campaign_sha256": hashlib.sha256(CAMPAIGN.read_bytes()).hexdigest(),
            "native_laws": ["C", "E"], "exact_comparison": exact,
            "source_complete_instruments": laws,
            "parameter_one_new_clock_words": 4,
            "parameter_one_clock_word_pulses_each": 2042,
            "clock_and_field_share_selected_native_event_prefix": True,
            "preparation_finite_mixed_reference_flags_reused": False,
            "source_processor_proper_code_failure_contract": "all complement outputs retained; ideal probability zero on admitted proper-code inputs",
            "physical_seconds_calibrated": False,
            "integrated_gate_accepted": False,
            "source_sha256": bound_pins(source_bindings(paths))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, nargs="?", const=OUTPUT)
    parser.add_argument("--check", type=Path, nargs="?", const=OUTPUT)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n")
    if args.check:
        check_receipt(args.check, result)
    if not args.write and not args.check:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))
    else:
        print("Source-executed C/E shared-event clocks and four checkpoints verified.")


if __name__ == "__main__":
    main()
