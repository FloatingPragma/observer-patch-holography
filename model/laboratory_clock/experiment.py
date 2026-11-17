"""Replay exact clock attachment controls without acquiring new timestamps."""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
import hashlib
from pathlib import Path

from tools import strict_json
from ..native_clock import clock, host_clock
from ..primitive_observer import transfer
from .attachment import Endpoint, analyse, fit_rate, validate
from . import resource_clock, event_geometry, regional_clock
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_PATHS = tuple(sorted(set((*host_clock.SOURCE_PATHS, *clock.SOURCE_PATHS,
    *event_geometry.SOURCE_PATHS, *regional_clock.SOURCE_PATHS,
    *resource_clock.target.SOURCE_PATHS, *resource_clock.compiled.source.SOURCE_PATHS,
    *resource_clock.exact_core.SOURCE_PATHS,
    "model/compiled_clock_history/receipt.py", "model/compiled_clock_history/history.py",
    "model/repair_identification/__init__.py", "model/record_laws/__init__.py",
    "model/predictions/__init__.py",
    "model/laboratory_clock/resource_clock.py",
    "model/native_clock/host_clock_records.json", "model/laboratory_clock/attachment.py",
    "model/laboratory_clock/experiment.py", "model/laboratory_clock/__main__.py",
    "model/laboratory_clock/__init__.py", "model/laboratory_clock/endpoints_example.json"))))


def deadline_survival(reduction, seconds_per_attempt, deadline_seconds, *, start="stationary"):
    """Push the native first-current law through a conditional constant cadence.

    This concerns elapsed duration from an actual attempt boundary. It adds
    no independent waiting-time distribution or exponential-clock premise.
    """
    lo, hi = seconds_per_attempt
    if (type(lo) is not F or type(hi) is not F or not 0 < lo <= hi
            or type(deadline_seconds) is not F or deadline_seconds < 0
            or start not in ("stationary", "palm")):
        raise ValueError("positive exact rate bounds and nonnegative elapsed deadline required")
    fast = deadline_seconds//lo
    slow = deadline_seconds//hi
    if fast > 4096:
        raise ValueError("the exact deadline certificate is bounded to 4096 attempts")
    no = reduction["marked"][0]
    if any(no[i][j] for i in range(3) for j in range(3) if i != j):
        raise ValueError("the source no-current reduction changed")
    prior = reduction[start]

    def survival(n):
        return sum((prior[i]*no[i][i]**n for i in range(3)), F(0))

    return {"preparation": start, "deadline_seconds": str(deadline_seconds),
            "attempts_completed_bounds": [slow, fast],
            "probability_first_current_after_deadline": [str(survival(fast)), str(survival(slow))],
            "uses_mean_wait_as_a_deterministic_conversion": False}


def _point(name, n, seconds, error=F(0), count_error=0):
    return Endpoint(name, max(0, n-count_error), n+count_error, seconds-error, seconds+error)


def _host_attachment(raw, replay):
    # Existing raw tape integrity is checked by host_clock.analyse first. It
    # supplies no accuracy bound for these independent physical timestamps.
    input_data = {"schema": "oph.affine_clock_input.v1", "data_kind": "host_measurement",
        "counter": "native_attempts_including_zero_currents",
        "provenance": {"event_mapping": "replayed original native seam/coin/current tape",
            "preparation_and_scheduler": "retained native_clock capture policy and seeded software service",
            "detection_bounds": "complete recorded attempt tape, checked against every endpoint",
            "timestamp_bounds": None, "traceability": None, "split_precommitment": None},
        "endpoints": [{"id": f"host_{i}", "role": "calibration" if i < 2 else "validation",
            "attempt_count": point["attempts"], "count_error": 0,
            "seconds": str(F(point["elapsed_ns"], 10**9)), "seconds_error": None}
            for i, point in enumerate(raw["runs"][0]["block_endpoints"])]}
    baseline = tuple(_point(f"baseline_{i}", row["attempts"], F(row["elapsed_ns"], 10**9))
                     for i, row in enumerate(raw["runs"][0]["block_endpoints"]))
    # The first timestamp defines the origin. Controls use the same elapsed
    # origin, so it need not be inserted twice with different endpoint ids.
    control = tuple(_point(f"control_{i}", row["attempts"], F(row["elapsed_ns"], 10**9))
                    for i, row in enumerate(raw["runs"][1]["block_endpoints"]) if i)
    exact_baseline = validate(baseline[:2], baseline[2:])
    exact_control = validate(baseline[:2], control)
    expected = "incompatible_with_every_calibration"
    if exact_baseline["joint_validation"] != expected or exact_control["joint_validation"] != expected:
        raise ArithmeticError("retained exact-cadence negative controls changed")
    return {"retained_raw_sha256": hashlib.sha256(host_clock.RAW_PATH.read_bytes()).hexdigest(),
        "all_native_records_replayed": replay["replay_integrity_passed"],
        "native_attempts_each_run": [r["native_attempts"] for r in replay["native_record_replays"]],
        "same_record_tape": replay["timing_control"]["same_complete_native_records"],
        "actual_host_measurement_with_unsupplied_accuracy": analyse(input_data),
        "zero_error_counterfactual": {
            "interpretation": "exact affine cadence fails even for the baseline; zero error is not a measured accuracy bound",
            "baseline": exact_baseline, "same_record_workload_control": exact_control},
        "new_timestamps_acquired": False,
        "previous_movement_clock_engineering_tolerance_changed": False}


def run():
    fixture = strict_json.load(HERE/"endpoints_example.json")
    fixture_result = analyse(fixture)
    synthetic_control = deepcopy(fixture)
    for row in synthetic_control["endpoints"]:
        if row["role"] == "validation":
            row["seconds"] = str(F(row["seconds"])+F(1, 10))
    changed = analyse(synthetic_control)
    if (fixture_result["numerical_test"]["joint_validation"] != "compatible_with_some_calibration"
            or changed["numerical_test"]["joint_validation"] != "incompatible_with_every_calibration"
            or changed["numerical_test"]["calibration"] != fixture_result["numerical_test"]["calibration"]):
        raise ArithmeticError("the unchanged-calibration synthetic control failed")
    calibration = (_point("c0", 0, F(0), F(1)), _point("c1", 10, F(10), F(1)))
    joint = validate(calibration, (_point("v0", 20, F(19)), _point("v1", 30, F(32))))
    if (not all(row["compatible_with_calibration"] for row in joint["per_endpoint"])
            or joint["joint_validation"] != "incompatible_with_every_calibration"):
        raise ArithmeticError("the common-rate/offset joint control failed")
    exact_cal = (_point("c0", 0, F(0)), _point("c1", 1000, F(1)))
    missed = {"no_missing_events": validate(exact_cal, (_point("v", 2000, F(401, 200)),)),
              "up_to_five_missing_events": validate(exact_cal, (_point("v", 2000, F(401, 200), count_error=5),))}

    kernel = transfer.native_counts()
    reduced = clock.lumping(kernel)  # checks all native marked rows again
    cal = fixture_result["numerical_test"]["calibration"]
    rates = tuple(F(cal[f"seconds_per_attempt_{side}"]) for side in ("lower", "upper"))
    word = (0, 1, 0, 1)
    word_probability = sum(clock.word_law(reduced, word), F(0))
    scales = [F(1, 1000), F(1, 100), F(1)]
    scale_rows = [{"seconds_per_attempt": str(rate), "native_word": list(word),
                   "native_word_probability": str(word_probability),
                   "stationary_moving_events_per_second": str(reduced["move_mass"]/rate)}
                  for rate in scales]
    native = {"parameter": "all native attempts, including zero currents",
        "carrier_ports": 12, "native_configurations": len(kernel.states),
        "all_marked_source_rows_rechecked": True,
        "synthetic_calibration_seconds_per_attempt": list(map(str, rates)),
        "conditional_deadline_predictions": [deadline_survival(reduced, rates, deadline, start=start)
            for start in ("stationary", "palm") for deadline in (F(1, 100), F(1, 50))],
        "scale_nonidentifiability": scale_rows,
        "nonaffine_retiming_control": {
            "same_native_attempt_endpoints": [0, 1, 2, 3, 4],
            "linear_seconds": [str(F(n, 1000)) for n in range(5)],
            "quadratic_seconds": [str(F(n*n, 1000)) for n in range(5)],
            "same_complete_word_law_when_feedback_uses_records_only": True,
            "lab_time_triggered_feedback_is_outside_this_invariance": True}}
    raw = strict_json.load(host_clock.RAW_PATH)
    replay = host_clock.analyse(raw)
    return {"schema": "oph.laboratory_clock_attachment.v1",
            "sources": bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS}),
            "synthetic_endpoint_fixture": fixture_result,
            "synthetic_held_out_rate_change": changed,
            "joint_validation_control": joint,
            "missed_event_control": missed,
            "native_waiting_law_attachment": native,
            "compiled_operation_clock_attachment": resource_clock.run(),
            "visible_event_clock_geometry": event_geometry.run(),
            "collective_regional_clock_laws": regional_clock.run(),
            "retained_host_attachment": _host_attachment(raw, replay),
            "scope": {"federation_required": False, "new_laboratory_data": False,
                "physical_proper_time_derived": False, "physical_clock_scale_selected_by_repair": False,
                "native_modular_phase_identified_with_attempts": False,
                "traceability_authenticated": False, "conditional_affine_clock_test_complete": True}}
