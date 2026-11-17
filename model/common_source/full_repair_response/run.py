"""Recompute full repair/mass-observable diagnostics and preserve their scope."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import coherent_arc, coupling, memory, quantum_response, mass_observable, current_response
from ..chiral_candidate import run as candidate_run
from ...native_modular_reference import common_caps, local_preparation
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASELINE = HERE.parent / "hierarchy_diagnosis/results.json"
NOTES = (
    ROOT / "docs/research/OPH_PARTICLE_MASS_RESPONSE.md",
)

REPLAY_ABSOLUTE_TOLERANCE = 5e-11
REPLAY_RELATIVE_TOLERANCE = 2e-9


def compare_replay(retained, current, path="report"):
    """Exact discrete claims and algebra; tolerant floating diagnostics only."""
    if type(retained) is not type(current):
        raise ValueError(f"replay type changed: {path}")
    if isinstance(retained, dict):
        if retained.keys() != current.keys():
            raise ValueError(f"replay keys changed: {path}")
        for key in retained:
            compare_replay(retained[key], current[key], f"{path}.{key}")
    elif isinstance(retained, list):
        if len(retained) != len(current):
            raise ValueError(f"replay length changed: {path}")
        for i, (a, b) in enumerate(zip(retained, current, strict=True)):
            compare_replay(a, b, f"{path}[{i}]")
    elif isinstance(retained, float):
        if not (math.isfinite(retained) and math.isfinite(current)
                and math.isclose(retained, current, rel_tol=REPLAY_RELATIVE_TOLERANCE,
                                 abs_tol=REPLAY_ABSOLUTE_TOLERANCE)):
            raise ValueError(f"floating replay diagnostic changed: {path}")
    elif retained != current:
        raise ValueError(f"replay claim or exact value changed: {path}")


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _bindings():
    previous = strict_json.load(BASELINE)
    bindings = dict(bound_pins(previous["source_bindings"]))
    candidate_run._check_bindings(bindings)
    paths = [BASELINE, *NOTES, ROOT / "requirements.txt"]
    paths += sorted(HERE.glob("*.py"))+sorted(HERE.glob("*.md"))
    names = set(common_caps.SOURCE_PATHS) | set(local_preparation.SOURCE_PATHS)
    names |= set(memory.SOURCE_FILES)
    names |= {
        "model/observer/coherent_continuation.py",
        "model/observer/native.py", "model/observer/quantum_mixing.py",
        "model/repair_scheduler/coins.py",
        "model/common_source/mass_response/history.py",
        "model/observer/source_scalar/compiler.py",
    }
    paths += [ROOT/name for name in names]
    for path in paths:
        name, digest = path.relative_to(ROOT).as_posix(), _hash(path)
        if name in bindings and bindings[name] != digest:
            raise ValueError(f"frozen source changed: {name}")
        bindings[name] = digest
    return dict(sorted(bindings.items()))


def run():
    candidate_run.verify(HERE.parent / "chiral_candidate/results.json")
    before = _bindings()
    print("Computing full stopped repair response, complement and incident clock", flush=True)
    population = memory.report()
    print("Computing harmonic kernel and conditional Higgs attachment", flush=True)
    harmonic = coupling.report()
    print("Checking all native quantum coherence sectors and modular flow", flush=True)
    quantum = quantum_response.report()
    print("Executing the specified coherent continuation and source coin words", flush=True)
    coherent, arrays = coherent_arc.report(with_arrays=True)
    print("Calibrating charged poles, dispersion and Euclidean decay", flush=True)
    observable = mass_observable.report()
    print("Computing exact native current fluctuations and record-preparation response", flush=True)
    current = current_response.report()
    for packet in (population, harmonic, quantum, coherent, observable, current):
        if any(before.get(name) != digest for name, digest in bound_pins(packet.get("source_bindings", {})).items()):
            raise ValueError("a response dependency is absent from source custody")
    result = {
        "schema": "oph.full_repair_response.v1",
        "population_memory_and_clock": population,
        "harmonic_coupling_control": harmonic,
        "native_quantum_response": quantum,
        "specified_coherent_continuation": coherent,
        "mass_observable_calibration": observable,
        "native_current_response": current,
        "verdict": {
            "full_complement_and_termination_outcomes_accounted_for": True,
            "neutral_population_response_identified_as_physical_charged_propagator": False,
            "new_native_to_chiral_field_intertwiner_proved": False,
            "source_has_more_family_information_than_stationary_mean": True,
            "mass_observable_calibrated_in_declared_chiral_control": True,
            "measured_Standard_Model_mass_ratios_recovered": False,
            "failed_conditional_attachments_exclude_all_OPH": False,
            "empirical_masses_used_to_fit_a_new_source_law": False,
            "probability_assigned_to_OPH_being_correct": False,
        },
        "source_bindings": bound_pins(before),
        "replay_comparison": {
            "rational_strings_discrete_claims_and_source_hashes": "exact",
            "floating_absolute_tolerance": REPLAY_ABSOLUTE_TOLERANCE,
            "floating_relative_tolerance": REPLAY_RELATIVE_TOLERANCE,
            "comparison_tolerances_are_physical_error_certificates": False,
        },
    }
    if _bindings() != before:
        raise RuntimeError("response sources changed during computation")
    candidate_run._check_bindings(before)
    # Normalize JSON container types and reject nonfinite or unsupported data.
    return json.loads(json.dumps(result, allow_nan=False)), arrays


def verify(path):
    path = Path(path)
    retained = strict_json.load(path)
    if retained.get("schema") != "oph.full_repair_response.v1":
        raise ValueError("wrong full-response receipt schema")
    artifact = retained.get("artifact", {})
    name = artifact.get("filename", "")
    if not name or Path(name).name != name:
        raise ValueError("array artifact must be adjacent to receipt")
    array_path = path.with_name(name)
    if _hash(array_path) != artifact.get("sha256"):
        raise ValueError("full-response work artifact changed")
    current, arrays = run()
    current["artifact"] = artifact
    compare_replay(retained, current)
    with np.load(array_path, allow_pickle=False) as stored:
        if set(stored.files) != set(arrays):
            raise ValueError("missing retained work or terminal detector arrays")
        for key, value in arrays.items():
            if (stored[key].shape != value.shape or not np.isfinite(stored[key]).all()
                    or not np.allclose(stored[key], value, rtol=REPLAY_RELATIVE_TOLERANCE,
                                       atol=REPLAY_ABSOLUTE_TOLERANCE)):
                raise ValueError(f"retained work array differs from recomputation: {key}")
    print("Verified by recomputing the full-response report and retained coherent work", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--verify", "--check", dest="verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        verify(args.verify)
        return
    if args.output.exists() or args.output.with_suffix(".npz").exists():
        parser.error("choose fresh JSON and NPZ output paths")
    result, arrays = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    artifact = args.output.with_suffix(".npz")
    with artifact.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    result["artifact"] = {"filename": artifact.name, "sha256": _hash(artifact),
                          "contract": "complete coherent arc trajectory including padding, untouched reference columns and both final Ramsey work outputs"}
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"Wrote {args.output} and {artifact}", flush=True)


if __name__ == "__main__":
    main()
