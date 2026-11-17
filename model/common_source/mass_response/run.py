"""Reproduce the unreleased clocked mass-response experiment and its controls."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy

from . import dispersion, experiment as e, gap_certificate, hubbard, model, native_attachment
from ...observer.source_scalar import compiler
from tools import strict_json


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CASE_NAMES = frozenset(("interacting_n32", "interacting_n64", "free_n64", "interacting_later_read",
                        "no_kick", "pinned_work", "late_interaction_off"))


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_bindings():
    """Bind loaded in-repository Python and the compiler's actual data inputs."""
    paths = {path for path in HERE.glob("*.py") if not path.name.startswith("test_")}
    for module in tuple(sys.modules.values()):
        location = getattr(module, "__file__", None)
        if location:
            path = Path(location).resolve()
            if (path.is_file() and path.suffix == ".py" and path.is_relative_to(ROOT)
                    and path.relative_to(ROOT).parts[0] in ("model", "code", "tools")):
                paths.add(path)
    paths.update(HERE/name for name in ("README.md", "DERIVATION.md", "HUBBARD_GAP.md"))
    bindings = {path.relative_to(ROOT).as_posix(): _hash(path) for path in sorted(paths)}
    bindings.update(compiler.source_summary()["source_sha256"])
    return dict(sorted(bindings.items()))


def _summary(result):
    probabilities = result.probabilities()
    abort_probability = math.fsum(map(e.norm, result.aborts.values()))
    prefix_error = max(abs(old-math.fsum(p for record, p in probabilities.items() if record[:2] == prefix))
                       for prefix, old in result.prefix_probabilities.items())
    forbidden = np.array([(mask&255).bit_count() != 4 for mask in range(1024)])
    charge_leakage = math.fsum(e.norm(state[forbidden]) for state in result.branches.values())
    return {
        "specification": result.specification,
        "record_probabilities": {"".join(map(str, key)): value for key, value in probabilities.items()},
        "source_abort_probabilities": {json.dumps(key): e.norm(value) for key, value in result.aborts.items()},
        "normal_terminal_records": len(probabilities), "retained_abort_records": len(result.aborts),
        "total_record_mass": math.fsum(probabilities.values())+abort_probability,
        "abort_probability": abort_probability,
        "old_record_marginal_residual": prefix_error,
        "endpoint_total_charge_four_leakage": charge_leakage,
        "receiver_occupation_probability": e.marginal(result, 2),
        "mass_return_probability": e.marginal(result, 1),
        "observer_gap_inverse": e.recover_gap_ratio(result),
        "resources_per_realized_nonabort_history": result.resources,
    }


def _case(name, **kwargs):
    started = time.monotonic()
    print(f"Executing {name}", flush=True)
    actual = e.execute(backend="source", **kwargs)
    target = e.execute(backend="continuous", **kwargs)
    duration = kwargs.get("interval", F(1, 4))
    steps = kwargs.get("steps", 64)
    bound = 2*model.split_operator_bound(duration/2, steps)
    if kwargs.get("propagation", True):
        bound += model.split_operator_bound(duration, steps)
    summary = _summary(actual)
    summary.update({
        "independent_continuous": {
            "receiver_occupation_probability": e.marginal(target, 2),
            "mass_return_probability": e.marginal(target, 1),
            "observer_gap_inverse": e.recover_gap_ratio(target),
        },
        "complete_quantum_output_half_trace_distance": e.quantum_distance(actual, target),
        "ideal_split_history_error_upper_exact": str(bound),
        "ideal_split_history_error_upper": float(bound),
        "bound_includes_floating_execution_error": False,
        "bound_includes_physical_services": False,
        "elapsed_host_seconds": time.monotonic()-started,
    })
    if summary["complete_quantum_output_half_trace_distance"] > float(bound)+2e-8:
        raise ArithmeticError("source/target diagnostic exceeds the ideal bound plus declared numerical tolerance")
    print(f"  inferred gap/clock {summary['observer_gap_inverse']['finite_single_charged_gap_over_clock']:.12g}; "
          f"complete-output distance {summary['complete_quantum_output_half_trace_distance']:.3g}", flush=True)
    return summary, actual


def run():
    initial_bindings = source_bindings()
    cases, histories = {}, {}
    settings = (
        ("interacting_n32", dict(steps=32)),
        ("interacting_n64", dict(steps=64)),
        ("free_n64", dict(steps=64, interaction=False)),
        ("interacting_later_read", dict(steps=64, interval=F(1, 2))),
        ("no_kick", dict(steps=64, kick=False)),
        ("pinned_work", dict(steps=64, propagation=False)),
        ("late_interaction_off", dict(steps=64, second_interval_interaction=False)),
    )
    for name, kwargs in settings:
        cases[name], histories[name] = _case(name, **kwargs)
    attachment = native_attachment.attach(histories["interacting_n64"], histories["no_kick"])
    # Adaptive record control is a channel on the complete output. The
    # distance between feedback and fixed-E laws changes by at most the sum
    # of their source/continuous splitting budgets.
    attachment["feedback_comparison_ideal_split_budget"] = 2*sum(
        float(weight)*cases[name]["ideal_split_history_error_upper"]
        for weight, name in ((F(5, 33), "interacting_n64"), (F(28, 33), "no_kick")))
    attachment["public_feedback_separation_after_ideal_split_budgets"] = (
        attachment["public_native_receiver_transcript_total_variation"]
        - attachment["feedback_comparison_ideal_split_budget"])
    if (attachment["actual"]["normal_paths"] != 8448
            or attachment["final_native_classical_total_variation"] <= 1e-8
            or attachment["public_feedback_separation_after_ideal_split_budgets"] <= 0):
        raise ArithmeticError("the complete native attachment did not produce measurable feedback")
    primary = cases["interacting_n64"]
    receiver = primary["receiver_occupation_probability"]
    epsilon = primary["ideal_split_history_error_upper"]
    signals = {}
    for label in ("no_kick", "pinned_work", "late_interaction_off"):
        difference = abs(receiver-cases[label]["receiver_occupation_probability"])
        signals[label] = {
            "absolute_record_probability_difference": difference,
            "ideal_split_budget_sum": epsilon+cases[label]["ideal_split_history_error_upper"],
            "separation_after_ideal_split_budgets": difference-epsilon-cases[label]["ideal_split_history_error_upper"],
            "floating_and_physical_service_errors_still_separate": True,
        }
        if signals[label]["separation_after_ideal_split_budgets"] <= 0:
            raise ArithmeticError(f"{label} response is not separated by the declared ideal split budgets")
    p = model.parameters()
    finite = dispersion.ring_gap_predictions(intra=p.intra, inter=p.inter, coupling=p.interaction)
    free = dispersion.ring_gap_predictions(intra=p.intra, inter=p.inter, coupling=0.)
    theorem = gap_certificate.certificate()
    finite_hubbard = hubbard.finite_volume_controls()
    momenta = ((0., 0., 0.), (.2, -.3, .1), (1., .5, -.7))
    refinement = []
    for level in (1, 2, 4, 8, 16, 32):
        for momentum in momenta:
            discrete = dispersion.refined_dispersion_squared(momentum, level)
            continuum = dispersion.continuum_dispersion_squared(momentum)
            bound = dispersion.refinement_remainder_bound(momentum, level)
            if not -2e-13 <= continuum-discrete <= bound+2e-13:
                raise ArithmeticError("free propagation remainder failed")
            refinement.append({"level": level, "momentum": momentum, "discrete_energy_squared": discrete,
                               "continuum_energy_squared": continuum, "error_upper": bound})
    document = {
        "schema": "oph.common_source.mass_response.v1", "publication_status": "unreleased candidate experiment",
        "base_release_commit": subprocess.check_output(["git", "rev-parse", "r2053^{commit}"], cwd=ROOT, text=True).strip(),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "native_parameters": p.record(), "cases": cases, "matched_intervention_signals": signals,
        "same_history_native_feedback": attachment,
        "independent_analytic_finite_spectra": {"interacting": finite, "free": free},
        "free_three_axis_extension": {
            "contract": "separate free eight-component continuation; not the interacting ring's bulk pole",
            "positive_band_gap_over_clock": dispersion.finite_torus_gap(1)/p.clock,
            "all_volume_gap_formula": "sqrt(3)*abs(a-b), every positive periodic side",
            "extra_regulator_components_retained": 8, "refinement_checks": refinement,
        },
        "alternative_hubbard_law": {"exact_charge_threshold_certificate": theorem,
                                    "finite_volume_diagnostics": finite_hubbard,
                                    "same_as_primary_interacting_ring": False},
        "claim_scope": {
            "complete_native_pulse_experiment": True,
            "same_history_clock_charge_interaction_and_remote_response": True,
            "no_empirical_mass_inputs": True, "architecture_map_is_declared": True,
            "structural_mass_is_spontaneous_mass_generation": False,
            "finite_cell_gap_is_thermodynamic_mass": False,
            "hubbard_charge_gap_is_full_many_body_gap": False,
            "standard_model_particle_identification": False,
            "physical_electromagnetic_charge_identification": False,
            "autonomous_physical_clock_derived": False,
            "interacting_three_spatial_dimensional_continuum_proved": False,
            "responsive_gravity_derived": False,
        },
    }
    document["source_bindings"] = source_bindings()
    if document["source_bindings"] != initial_bindings:
        raise RuntimeError("source files changed during execution; no receipt is valid")
    return document


def verify(path):
    document = strict_json.load(Path(path))
    def finite(value):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError("nonfinite receipt value")
        if isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(document)
    if document.get("schema") != "oph.common_source.mass_response.v1":
        raise ValueError("wrong receipt schema")
    required = {f"model/common_source/mass_response/{name}.py" for name in
                ("model", "history", "experiment", "run", "gap_certificate", "dispersion", "hubbard", "native_attachment")}
    required.update(compiler.source_summary()["source_sha256"])
    if not required <= set(document["source_bindings"]):
        raise ValueError("incomplete computational source custody")
    mismatches = []
    for name, expected in document["source_bindings"].items():
        current = (ROOT/name).resolve()
        if not current.is_relative_to(ROOT):
            raise ValueError("source binding escaped the repository")
        if not current.is_file() or _hash(current) != expected:
            mismatches.append(name)
    if mismatches:
        raise ValueError(f"source bindings changed: {mismatches}")
    for name in ("standard_model_particle_identification", "physical_electromagnetic_charge_identification",
                 "autonomous_physical_clock_derived", "interacting_three_spatial_dimensional_continuum_proved",
                 "responsive_gravity_derived", "finite_cell_gap_is_thermodynamic_mass"):
        if document["claim_scope"][name]:
            raise ValueError("unsupported physical promotion in retained receipt")
    if set(document["cases"]) != CASE_NAMES:
        raise ValueError("the complete seven-case experiment is required")
    for case in document["cases"].values():
        probabilities = case["record_probabilities"]
        if set(probabilities) != {format(i, "05b") for i in range(32)} or any(p < 0 for p in probabilities.values()):
            raise ValueError("incomplete or negative normal record probabilities")
        aborts = case["source_abort_probabilities"]
        if len(aborts) != 12 or any(p < 0 for p in aborts.values()):
            raise ValueError("incomplete or negative abort probabilities")
        recomputed_mass = math.fsum(probabilities.values())+math.fsum(aborts.values())
        if abs(recomputed_mass-case["total_record_mass"]) > 1e-12:
            raise ValueError("retained probability sum differs from the reported mass")
        if abs(math.fsum(aborts.values())-case["abort_probability"]) > 1e-12:
            raise ValueError("reported abort mass differs from the retained records")
        for position, field in ((1, "mass_return_probability"), (2, "receiver_occupation_probability")):
            value = math.fsum(p for record, p in probabilities.items() if record[position] == "0")
            if abs(value-case[field]) > 1e-12:
                raise ValueError("reported readout differs from the retained records")
        if abs(case["total_record_mass"]-1) > 1e-9 or case["old_record_marginal_residual"] > 1e-9:
            raise ValueError("incomplete recorded history in receipt")
        if case["endpoint_total_charge_four_leakage"] > 1e-9:
            raise ValueError("charge leakage in retained logical history")
        if not 0 <= case["complete_quantum_output_half_trace_distance"] <= case["ideal_split_history_error_upper"]+2e-8:
            raise ValueError("declared numerical comparison exceeds retained ideal split bound")
        # The point decoder needs only probabilities and apparatus calibration.
        compressed = e.Result({tuple(map(int, key)): np.array([[math.sqrt(value)]])
                               for key, value in probabilities.items()},
                              {key: np.array([[math.sqrt(value)]]) for key, value in aborts.items()},
                              {}, {}, case["specification"])
        decoded = e.recover_gap_ratio(compressed)["finite_single_charged_gap_over_clock"]
        if abs(decoded-case["observer_gap_inverse"]["finite_single_charged_gap_over_clock"]) > 1e-12:
            raise ValueError("reported gap differs from the retained observer records")
    signals = document["matched_intervention_signals"]
    if set(signals) != {"no_kick", "pinned_work", "late_interaction_off"}:
        raise ValueError("all three matched intervention controls are required")
    primary = document["cases"]["interacting_n64"]
    for name, signal in signals.items():
        control = document["cases"][name]
        difference = abs(primary["receiver_occupation_probability"]-control["receiver_occupation_probability"])
        budget = primary["ideal_split_history_error_upper"]+control["ideal_split_history_error_upper"]
        if (abs(difference-signal["absolute_record_probability_difference"]) > 1e-12
                or abs(budget-signal["ideal_split_budget_sum"]) > 1e-12
                or abs(difference-budget-signal["separation_after_ideal_split_budgets"]) > 1e-12):
            raise ValueError("intervention summary differs from its retained controls")
        if signal["separation_after_ideal_split_budgets"] <= 0:
            raise ValueError("retained intervention signal is not separated")
    attachment = document["same_history_native_feedback"]
    if (attachment["actual"]["normal_paths"] != 8448 or attachment["actual"]["abort_paths"] != 1584
            or abs(attachment["actual"]["total_probability"]-1) > 1e-9
            or attachment["kick_probability_exact"] != "5/33"
            or attachment["final_native_classical_total_variation"] <= 0
            or attachment["public_feedback_separation_after_ideal_split_budgets"] <= 0):
        raise ValueError("incomplete or unresponsive native feedback attachment")
    print(f"Verified custody ({len(document['source_bindings'])} source bindings) and retained seven-case checks; "
          "this command does not rerun the experiment")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        if args.output is not None:
            parser.error("verification does not write a receipt")
        verify(args.verify)
        return
    if args.output is None:
        parser.error("choose a fresh --output path; retained receipts are never overwritten")
    if args.output.exists():
        parser.error("output exists; use a fresh replay path")
    document = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(document, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
