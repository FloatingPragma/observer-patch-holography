"""Exact and source-executed inference of shared native incidence from joint readouts.

This successor keeps the C/E two-attempt law fixed. It coarse-grains complete
source clock and Ramsey records only after their native event count is known.
The inferred quantity is a dimensionless count, not calibrated clock time.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import joint_clock, ramsey_service
from .joins_contract import source_bindings, check_receipt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUTPUT = HERE / "clock_count_inference_results.json"
CAMPAIGN = HERE / "clock_count_inference_campaign.json"
COUNTS = (1, 2)
CURRENTS = (-1, 0, 1)
SIGNS = (1, -1)
Y = tuple((current, sign) for current in CURRENTS for sign in SIGNS)


def _pred(name):
    receipt = strict_json.load(HERE / name)
    if source_bindings(receipt["source_sha256"]) != receipt["source_sha256"]:
        raise ValueError(f"{name} source custody changed")
    return receipt


def _sum_coeff(rows):
    return tuple(sum((v[i] for v in rows), F()) for i in range(3))


def _interval(coeff):
    return joint_clock.evaluate_interval(coeff)


def _risk(a, b):
    """The optimal squared-error risk for N in {1,2}, with rational bounds."""
    if min((*a, *b)) < 0:
        raise ArithmeticError("negative probability enclosure")
    lo = F() if a[0] + b[0] == 0 else a[0]*b[0]/(a[0]+b[0])
    hi = F() if a[1] + b[1] == 0 else a[1]*b[1]/(a[1]+b[1])
    return lo, hi


def _sum_interval(rows):
    rows = tuple(rows)
    return sum((x[0] for x in rows), F()), sum((x[1] for x in rows), F())


def _show(interval, places=12):
    scale = 10**places
    return [str(F(math.floor(interval[0]*scale), scale)),
            str(F(math.ceil(interval[1]*scale), scale))]


def _float_interval(value):
    """Nearest-grid center plus one unit, only for checked float display."""
    if not math.isfinite(value) or value < -1e-9:
        raise ArithmeticError("invalid nonnegative source diagnostic")
    scale = 10**8
    center = round(value*scale)
    return [str(F(max(0, center-1), scale)),
            str(F(center+1, scale))]


def _report_source_intervals(row):
    for key in ("source_read_mass_before_conditioning", "source_complete_mass"):
        row[key+"_diagnostic_interval"] = _float_interval(row.pop(key))
    for key in ("finite_preparation_flag_probabilities",
                "source_clock_export_abort_mass_by_wave_prefix",
                "source_conditional_squared_error_risk"):
        row[key+"_diagnostic_intervals"] = {name: _float_interval(v)
                                                for name, v in row.pop(key).items()}
    for output in row["source_conditional_probabilities"]:
        output["probability_diagnostic_interval"] = _float_interval(output.pop("probability"))
    return row


def exact_law(law):
    """Independent rational postprocessing of complete inherited exact rows."""
    grouped = defaultdict(list)
    rows = joint_clock.exact_clock_rows(law)
    for (_history, _coin, _label, current, count), coeff in rows.items():
        for sign, p in zip(SIGNS, joint_clock.field_probabilities(count), strict=True):
            grouped[count, current, sign].append(tuple(x*p for x in coeff))
    table = {key: _interval(_sum_coeff(values)) for key, values in grouped.items()}
    for count in COUNTS:
        for current, sign in Y:
            table.setdefault((count, current, sign), (F(), F()))
    if not (sum(table[count, current, sign][0] for count in COUNTS for current, sign in Y)
            <= 1 <= sum(table[count, current, sign][1] for count in COUNTS for current, sign in Y)):
        raise ArithmeticError("complete conditional ideal law is not normalized")
    def risk(groups):
        return _sum_interval(_risk(
            _sum_interval(table[1, current, sign] for current, sign in group),
            _sum_interval(table[2, current, sign] for current, sign in group))
            for group in groups)
    prior = risk((Y,))
    field = risk(tuple((r, s) for r in CURRENTS) for s in SIGNS)
    clock = risk(tuple((r, s) for s in SIGNS) for r in CURRENTS)
    joint = risk(((y,) for y in Y))
    gain = clock[0]-joint[1], clock[1]-joint[0]
    if gain[0] <= 0 or not (joint[1] < clock[0] < prior[0]):
        raise ArithmeticError("joint readout gave no certified uncertainty reduction")
    odds = {}
    for current, sign in Y:
        a, b = table[1, current, sign], table[2, current, sign]
        if b[1] == 0:
            odds[f"{current},{sign}"] = ["0", "0"]
        else:
            odds[f"{current},{sign}"] = _show((b[0]/(a[1]+b[0]), b[1]/(a[0]+b[1])))
    return {"law": law, "retained_exact_read_rows": len(rows),
            "hidden_source_variable": "N in {1,2}, count of zero-inclusive port-0 incidents",
            "visible_readout": "terminal detector current and source-compiled field Ramsey sign",
            "ideal_conditional_probabilities": [
                {"N": count, "detector_current": current, "field_sign": sign,
                 "probability_interval": _show(table[count, current, sign])}
                for count in COUNTS for current, sign in Y],
            "posterior_probability_N2_by_readout": odds,
            "optimal_squared_error_risk_count_units": {
                "prior": _show(prior), "field_only": _show(field),
                "clock_only": _show(clock), "joint": _show(joint),
                "joint_gain_over_clock": _show(gain)},
            "matched_marginal_independent_field_control_risk": _show(clock),
            "readout_complete_on_ideal_target": True}


def source_law(law, field_mass, wait_mass):
    """Run source arrays on complete clock frames and source field factors."""
    frames = joint_clock.clock_frames(law, raw=True)
    prep = joint_clock.clock.preparation.population_maps(law=law)
    density = np.diag(prep["success"][:, 0].astype(complex))
    grouped = defaultdict(float)
    abort = defaultdict(float)
    read_slots = abort_slots = 0
    for frame in frames:
        p = float(np.trace(frame.matrix@density@frame.matrix.conj().T).real)
        if frame.abort:
            abort[frame.ticks] += p*wait_mass[frame.ticks]
            abort_slots += 1
        else:
            current = frame.key[-1]
            for sign, q in zip(SIGNS, field_mass[frame.ticks], strict=True):
                grouped[frame.ticks, current, sign] += p*q
                read_slots += 1
    flags = {name: float(np.sum(array[:, 0])) for name, array in prep.items()
             if name != "success"}
    total = sum(grouped.values())+sum(abort.values())+sum(flags.values())
    if abs(total-1) > 2e-7:
        raise ArithmeticError("complete finite-preparation source mass lost")
    conditional = {}
    read_mass = sum(grouped.values())
    for count in COUNTS:
        for current, sign in Y:
            conditional[count, current, sign] = grouped[count, current, sign]/read_mass
    def risk(groups):
        answer = 0.
        for group in groups:
            a = sum(conditional[1, *y] for y in group)
            b = sum(conditional[2, *y] for y in group)
            if a+b: answer += a*b/(a+b)
        return answer
    clock = risk(tuple((r,s) for s in SIGNS) for r in CURRENTS)
    joint = risk(((y,) for y in Y))
    if not joint < clock:
        raise ArithmeticError("executed source comparison lost joint information")
    return {"law": law, "complete_source_clock_frames": len(frames),
            "source_read_slots_including_both_field_signs": read_slots,
            "source_clock_export_abort_slots": abort_slots,
            "finite_preparation_flag_probabilities": flags,
            "source_clock_export_abort_mass_by_wave_prefix": {str(n): abort[n] for n in (0,1,2)},
            "source_read_mass_before_conditioning": read_mass,
            "source_complete_mass": total,
            "source_conditional_squared_error_risk": {
                "clock_only": clock, "joint": joint,
                "joint_gain_over_clock": clock-joint},
            "source_conditional_probabilities": [
                {"N": n, "detector_current": r, "field_sign": s,
                 "probability": conditional[n,r,s]}
                for n in COUNTS for r,s in Y]}


def run():
    contract = strict_json.load(CAMPAIGN)
    joint_receipt = _pred("joint_clock_results.json")
    ramsey_receipt = _pred("ramsey_service_results.json")
    wave_receipt = _pred("../observer_time_fields/wave_reference_results.json")
    if (contract["predecessor_campaign"] != ramsey_receipt["campaign_id"]
            or contract["scope"] != "FINITE_COMPONENT_DIAGNOSTICS"
            or contract["boundary_status"] != "FROZEN_DIAGNOSTIC_CONTRACT"
            or contract["native_laws"] != ["C", "E"]
            or contract["frozen_before_execution"] is not True
            or contract["whole_program_frozen"] is not False
            or contract["native_family_changed"] is not False
            or contract["physical_comparison_count_added"] != 0
            or contract["physical_seconds_calibrated"] is not False
            or contract["integrated_gate_accepted"] is not False):
        raise ValueError("wrong frozen clock-count campaign")
    for key, name in (("joint_result_sha256", "joint_clock_results.json"),
                      ("ramsey_result_sha256", "ramsey_service_results.json"),
                      ("wave_reference_result_sha256", "../observer_time_fields/wave_reference_results.json")):
        if contract[key] != hashlib.sha256((HERE/name).read_bytes()).hexdigest():
            raise ValueError("predecessor result bytes changed")
    del wave_receipt  # source custody, no capped reference substituted into the readout
    exact = [exact_law(law) for law in ("C", "E")]
    fields = {n: ramsey_service.field_columns(n) for n in COUNTS}
    field_mass = {n: tuple(float(np.vdot(v, v).real) for v in fields[n])
                  for n in COUNTS}
    waits = {n: ramsey_service.field_columns(n, read=False) for n in (0, 1, 2)}
    wait_mass = {n: float(np.vdot(v, v).real) for n, v in waits.items()}
    source = [source_law(law, field_mass, wait_mass) for law in ("C", "E")]
    source_error = F(joint_receipt["source_and_service_budget"]["complete_preparation_and_joint_halftrace_upper"])
    per_law_wave = {row["law"]: F(row["wave_halftrace_bound"])
                    for row in ramsey_receipt["joint_laws"]}
    for row, compiled in zip(exact, source, strict=True):
        failure = source_error+per_law_wave[row["law"]]
        p = F(5, 33)
        conditional_tv = 2*failure/(p-failure)
        gain = F(row["optimal_squared_error_risk_count_units"]["joint_gain_over_clock"][0])
        if failure >= p or gain <= 2*conditional_tv:
            raise ArithmeticError("source/error budget cannot certify strict gain")
        row["complete_finite_preparation_and_field_halftrace_upper"] = str(failure)
        row["conditional_read_TV_upper"] = str(conditional_tv)
        row["source_robust_joint_gain_lower"] = str(gain-2*conditional_tv)
        ideal = {(r["N"],r["detector_current"],r["field_sign"]):
                 tuple(map(F,r["probability_interval"]))
                 for r in row["ideal_conditional_probabilities"]}
        for output in compiled["source_conditional_probabilities"]:
            key = output["N"],output["detector_current"],output["field_sign"]
            lo, hi = ideal[key]
            if not float(lo-conditional_tv)-1e-10 <= output["probability"] <= float(hi+conditional_tv)+1e-10:
                raise ArithmeticError("executed source readout exceeds retained bound")
    source = [_report_source_intervals(row) for row in source]
    paths = {
        "model/clock_defect_energy/clock_count_inference.py",
        "model/clock_defect_energy/test_clock_count_inference.py",
        "model/clock_defect_energy/CLOCK_COUNT_INFERENCE.md",
        "model/clock_defect_energy/clock_count_inference_campaign.json",
        "model/clock_defect_energy/joint_clock_results.json",
        "model/clock_defect_energy/ramsey_service_results.json",
        "model/observer_time_fields/wave_reference_results.json",
        "Lean/Time/ClockReadout.lean", "Lean/Time/ClockComparison.lean",
        "Lean/Dynamics/ChoiCPTP.lean", "Lean/Dynamics/PrivateInner.lean"}
    return {"schema": "oph.clock_defect_energy.clock_count_inference.v1",
            "campaign_id": contract["campaign_id"],
            "contract_sha256": hashlib.sha256(CAMPAIGN.read_bytes()).hexdigest(),
            "native_laws": ["C", "E"], "exact": exact, "compiled": source,
            "modular_parameter_per_successful_incident": "4",
            "physical_seconds_calibrated": False,
            "observer_clock_interface_inhabited_by_stochastic_detector": False,
            "affine_three_event_comparison_claimed": False,
            "native_family_changed": False,
            "integrated_gate_accepted": False,
            "source_sha256": source_bindings(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, nargs="?", const=OUTPUT)
    parser.add_argument("--check", type=Path, nargs="?", const=OUTPUT)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    if args.check:
        check_receipt(args.check, result)
    if not args.write and not args.check:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    else:
        print("Complete C/E source and exact clock-count uncertainty verified.")


if __name__ == "__main__":
    main()
