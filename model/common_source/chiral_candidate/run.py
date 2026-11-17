"""Replay the declared chiral action, clocked source history and mass tests."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import action, hierarchy_controls, history, kinetic, spacetime_attachment
from ..particle_identification import run as particle_run
from ..mass_response import experiment as e
from tools import strict_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASELINE = HERE.parent/"particle_identification/results.json"
NOTE = ROOT/"docs/research/OPH_PARTICLE_MASS_RESPONSE.md"
SPECTRAL_RESPONSE = ROOT/"model/observer_particles/spectral_response.py"


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _check_bindings(bindings):
    for name, expected in bindings.items():
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file() or _hash(path) != expected:
            raise ValueError(f"chiral candidate source custody changed: {name}")


def _bindings(action_report,geometry_report):
    baseline = strict_json.load(BASELINE)
    result = dict(baseline["source_bindings"])
    for report in (action_report,geometry_report):
        for name,value in report["source_bindings"].items():
            if name in result and result[name] != value:
                raise ValueError(f"inconsistent source custody: {name}")
            result[name] = value
    paths = [p for p in HERE.glob("*.py") if not p.name.startswith("test_")]
    paths += list(HERE.glob("*.md"))+[BASELINE,BASELINE.with_suffix(".npz"),NOTE,SPECTRAL_RESPONSE]
    result.update({path.relative_to(ROOT).as_posix():_hash(path) for path in paths})
    return dict(sorted(result.items()))


def run():
    particle_run.verify(BASELINE)
    print("Checking the exact chiral/Higgs action and Weyl spectrum",flush=True)
    matter = action.report()
    geometry = spacetime_attachment.report()
    initial = _bindings(matter,geometry)
    _check_bindings(initial)
    kinematics = kinetic.report()
    experiment,arrays = history.execute()
    print("Computing fixed nonlinear hierarchy controls before comparison",flush=True)
    predictions = hierarchy_controls.predictions()
    comparison = hierarchy_controls.lepton_comparison(predictions)
    result = {
        "schema":"oph.chiral_candidate.v1","publication_status":"unreleased trial architecture",
        "action":matter,"kinematics":kinematics,"continuing_source_history":experiment,
        "frozen_metric_attachment":geometry,
        "source_hierarchy_predictions":predictions,"separate_empirical_comparison":comparison,
        "verdict":{
            "decoded_internal_module_joined_to_chiral_continuum_symbol":True,
            "verified_Higgs_lines_used_for_mass_pairing":True,
            "full90_mode_CAR_and_clocked_source_probe_executed":True,
            "native_family_nonscalar_feature_used":True,
            "common_frozen_metric_principal_symbol_constructed":True,
            "same_reconstructed_native_spacetime_field_history_proved":False,
            "chosen_architecture_has_no_empirical_mass_fit":True,
            "linear_trial_matches_observed_mass_hierarchy":False,
            "any_tested_fixed_nonlinear_texture_matches_both_lepton_ratios":any(
                row["both_ratios_pass"] for row in comparison["comparisons"]),
            "entire_OPH_excluded_by_these_trials":False,
            "quantum_gauge_and_Higgs_dynamics_executed":False,
            "local_nonperturbative_chiral_gauge_regulator_established":False,
            "Higgs_background_or_physical_units_source_selected":False,
            "observed_neutrino_masses_derived":False,
            "observed_standard_model_mass_spectrum_derived":False,
            "physical_particle_identification_complete":False,
        },
        "source_bindings":_bindings(matter,geometry),
    }
    if result["source_bindings"] != initial:
        raise RuntimeError("candidate files changed during execution")
    _check_bindings(initial)
    return result,arrays


def verify(path):
    path = Path(path)
    document = strict_json.load(path)
    if document.get("schema") != "oph.chiral_candidate.v1":
        raise ValueError("wrong chiral receipt schema")
    required = {p.relative_to(ROOT).as_posix() for p in HERE.glob("*.py") if not p.name.startswith("test_")}
    required.update(p.relative_to(ROOT).as_posix() for p in HERE.glob("*.md"))
    required.add(BASELINE.relative_to(ROOT).as_posix())
    required.add(NOTE.relative_to(ROOT).as_posix())
    required.add(SPECTRAL_RESPONSE.relative_to(ROOT).as_posix())
    if not required <= set(document["source_bindings"]):
        raise ValueError("incomplete chiral source custody")
    _check_bindings(document["source_bindings"])
    artifact = document["artifact"]
    if Path(artifact["filename"]).name != artifact["filename"]:
        raise ValueError("quantum artifact must be adjacent to receipt")
    array_path = path.with_name(artifact["filename"])
    if _hash(array_path) != artifact["sha256"]:
        raise ValueError("retained quantum artifact changed")
    report = document["continuing_source_history"]
    with np.load(array_path,allow_pickle=False) as arrays:
        if (arrays["probe_outputs"].shape != (64,256,4)
                or arrays["gaussian_covariances"].shape != (16,180,180)
                or arrays["gaussian_weights"].shape != (16,)
                or set(arrays["records"]) != {format(i,"06b") for i in range(64)}
                or set(arrays["field_records"]) != {format(i,"04b") for i in range(16)}):
            raise ValueError("incomplete retained chiral field or probe")
        for name in ("probe_outputs","gaussian_covariances","gaussian_weights"):
            if not np.isfinite(arrays[name]).all():
                raise ValueError("nonfinite retained field output")
        probabilities = report["source_record_probabilities"]
        if set(probabilities) != set(arrays["records"]):
            raise ValueError("missing observer records")
        for record,state in zip(arrays["records"],arrays["probe_outputs"],strict=True):
            value = probabilities[str(record)]
            if not 0 <= value <= 1 or abs(e.norm(state)-value) > 1e-12:
                raise ValueError("source state and observer weight disagree")
        if abs(math.fsum(probabilities.values())-1) > 1e-9:
            raise ValueError("observer probability mass lost")
        for record,weight in zip(arrays["field_records"],arrays["gaussian_weights"],strict=True):
            actual = math.fsum(value for key,value in probabilities.items() if key[:4] == record)
            if not 0 <= weight <= 1 or abs(actual-weight) > 1e-9:
                raise ValueError("full Gaussian history and source factor disagree")
        decoded = history.infer_ratios({tuple(map(int,key)):value for key,value in probabilities.items()},
                                      report["clock_calibration"],report["phase_window_contract"])
        if max(abs(a-b) for a,b in zip(decoded["native_family_order_gap_over_clock"],
                                      report["observer_inverse"]["native_family_order_gap_over_clock"],strict=True)) > 1e-12:
            raise ValueError("reported rest gaps differ from the observer records")
    for name in ("physical_particle_identification_complete","observed_standard_model_mass_spectrum_derived",
                 "local_nonperturbative_chiral_gauge_regulator_established","observed_neutrino_masses_derived"):
        if document["verdict"][name]:
            raise ValueError("unsupported physical promotion")
    print("Verified chiral candidate custody and retained field/clock weights; this does not replay dynamics")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        if args.output is not None:
            parser.error("verification does not write outputs")
        verify(args.verify); return
    if args.output is None or args.output.exists() or args.output.with_suffix(".npz").exists():
        parser.error("choose fresh JSON and companion NPZ output paths")
    document,arrays = run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    artifact = args.output.with_suffix(".npz")
    with artifact.open("xb") as handle:
        np.savez_compressed(handle,**arrays)
    document["artifact"] = {"filename":artifact.name,"sha256":_hash(artifact),
                            "contract":"all64 literal source probe outputs and16 full90-mode Gaussian conditional states with their weights"}
    with args.output.open("x") as handle:
        json.dump(document,handle,indent=2,allow_nan=False); handle.write("\n")
    print(f"Wrote {args.output} and {artifact}")


if __name__ == "__main__":
    main()
