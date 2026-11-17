"""Replay the candidate field readout and audit its particle identification."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import clifford, family_attachment, kinetic_response, response
from ..mass_response import experiment as e, model
from ..mass_response import run as mass_run
from tools import strict_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MASS_RECEIPT = HERE.parent/"mass_response/results.json"
PUBLIC_RESULTS = ROOT/"model/observer_particles/results"
PUBLIC_NAMES = ("representation_identification", "family_naturality", "spin_matter_join",
                "shared_spatial_matter", "unified_yukawa_action", "gauge_particle_identification")


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_public():
    """Read frozen algebra witnesses; never execute their writing producers."""
    packets, bindings = {}, {}
    for name in PUBLIC_NAMES:
        path = PUBLIC_RESULTS/(name+".json")
        packet = strict_json.load(path)
        producer = path.parent.parent/(name+".py")
        if _hash(producer) != packet["producer_sha256"]:
            raise ValueError(f"public producer custody changed: {name}")
        bindings[path.relative_to(ROOT).as_posix()] = _hash(path)
        bindings[producer.relative_to(ROOT).as_posix()] = packet["producer_sha256"]
        for location, expected in packet["source_sha256"].items():
            source = (ROOT/location).resolve()
            if not source.is_relative_to(ROOT) or _hash(source) != expected:
                raise ValueError(f"public input custody changed: {location}")
            bindings[location] = expected
        for location, expected in packet.get("artifacts", {}).items():
            artifact = (path.parent/location).resolve()
            if not artifact.is_relative_to(ROOT) or _hash(artifact) != expected:
                raise ValueError(f"public artifact custody changed: {location}")
            bindings[artifact.relative_to(ROOT).as_posix()] = expected
        packets[name] = packet
    return packets, bindings


def append_transfer(result, *, backend="source", quadrature="X"):
    """Continue every work/reference branch with an even CAR field read."""
    if backend not in ("source", "continuous"):
        raise ValueError("source or independent continuous backend required")
    branches, aborts = {}, dict(result.aborts)
    resources = dict(result.resources)
    cost = {}
    target = response.transfer_kraus(quadrature=quadrature) if backend == "continuous" else None
    for record, state in result.branches.items():
        if backend == "source":
            read = response.transfer_read(state, quadrature=quadrature)
            outputs, cost = read.outputs, read.resources
            for failure, output in zip(read.failure_records, read.failure_outputs, strict=True):
                aborts[(*record, "field_abort", *failure)] = output
        else:
            outputs = tuple(model.apply_subsystem(state, kraus, range(8), 10) for kraus in target)
        for bit, output in enumerate(outputs):
            branches[(*record, bit)] = output
    for key in ("native_pulses", "proper_code_service_events", "fresh_pointer_blanks", "committed_classical_bits"):
        resources[key] = resources.get(key, 0)+cost.get(key, 0)
    if backend == "source":
        resources["proper_code_load_events"] = resources["native_pulses"]
        resources["proper_code_export_events"] = resources["native_pulses"]
    specification = {**result.specification,
                     "record_order": [*result.specification["record_order"], "even_CAR_transfer_"+quadrature],
                     "additional_work_preparations": 0,
                     "field_read_charge_meaning": "global occupation U1; no electromagnetic dressing supplied"}
    return e.Result(branches, aborts, result.prefix_probabilities, resources, specification)


def _bindings(mass, family, public):
    result = dict(mass["source_bindings"])
    for additional in (family["source_bindings"], public):
        for key, value in additional.items():
            if key in result and result[key] != value:
                raise ValueError(f"inconsistent source custody for {key}")
            result[key] = value
    paths = [p for p in HERE.glob("*.py") if not p.name.startswith("test_")]
    paths += list(HERE.glob("*.md"))+[MASS_RECEIPT]
    result.update({path.relative_to(ROOT).as_posix(): _hash(path) for path in paths})
    return dict(sorted(result.items()))


def run():
    mass_run.verify(MASS_RECEIPT)
    mass = strict_json.load(MASS_RECEIPT)
    family = family_attachment.report()
    public, public_bindings = _read_public()
    initial_bindings = _bindings(mass, family, public_bindings)
    print("Executing continuing source history and appending the even CAR field read", flush=True)
    before = e.execute(steps=64, backend="source")
    baseline = mass["cases"]["interacting_n64"]["record_probabilities"]
    replay_residual = max(abs(value-baseline["".join(map(str, record))])
                          for record, value in before.probabilities().items())
    if replay_residual > 1e-10:
        raise ArithmeticError("the prior five-record experiment did not replay")
    actual = append_transfer(before)
    target = append_transfer(e.execute(steps=64, backend="continuous"), backend="continuous")
    old_marginal_error = max(abs(value-math.fsum(e.norm(state) for final, state in actual.branches.items()
                                                if final[:-1] == record))
                             for record, value in before.probabilities().items())
    comparison = e.quantum_distance(actual, target)
    budget = mass["cases"]["interacting_n64"]["ideal_split_history_error_upper"]
    if old_marginal_error > 1e-10 or comparison > budget+2e-8:
        raise ArithmeticError("the continuing field instrument failed its retained-history comparison")
    source_records = sorted(actual.branches)
    abort_records = sorted(actual.aborts, key=lambda value: json.dumps(value))
    arrays = {"records": np.array(["".join(map(str, key)) for key in source_records]),
              "outputs": np.stack([actual.branches[key] for key in source_records]),
              "abort_records": np.array([json.dumps(key) for key in abort_records]),
              "abort_outputs": np.stack([actual.aborts[key] for key in abort_records])}
    record_probabilities = {"".join(map(str, key)): e.norm(actual.branches[key]) for key in source_records}
    field_mean = math.fsum((1 if key[-1] == "0" else -1)*value for key, value in record_probabilities.items())
    representation = public["representation_identification"]
    result = {
        "schema": "oph.particle_identification.v1", "publication_status": "unreleased response and representation audit",
        "continuing_field_read": {
            "specification": actual.specification, "source_record_probabilities": record_probabilities,
            "record_probability_mass": math.fsum(record_probabilities.values()),
            "abort_probability_mass": math.fsum(map(e.norm, actual.aborts.values())),
            "normal_records": len(actual.branches), "abort_records": len(actual.aborts),
            "full_work_reference_shape": list(next(iter(actual.branches.values())).shape),
            "all_five_old_records_marginal_residual": old_marginal_error,
            "earlier_retained_replay_residual": replay_residual,
            "complete_quantum_output_half_trace_distance": comparison,
            "ideal_split_history_error_upper": budget,
            "floating_execution_and_physical_services_certified": False,
            "even_transfer_expectation_after_the_existing_five_reads": field_mean,
            "this_final_mean_is_a_separate_pole_fit": False,
            "observer_gap_inverse_after_extra_read": e.recover_gap_ratio(actual),
            "resources_per_realized_nonabort_history": actual.resources,
        },
        "finite_interacting_CAR_response": response.report(),
        "finite_free_control_CAR_response": response.report(interaction=False),
        "finite_inverse_response_kinetics": kinetic_response.report(),
        "separate_free_continuum_spinor": clifford.report(),
        "native_family_feature": family,
        "existing_conditional_standard_model_module": {
            "status": representation["status"],
            "decoded_internal_multiplets": representation["decoded_even_internal_multiplets"],
            "scope": representation["scope"],
            "unbroken_family_mass_conclusion": public["family_naturality"]["conditional_mass_conclusion"],
        },
        "verdict": {
            "same_history_charged_gap_has_nonzero_CAR_response": True,
            "even_transfer_is_source_measurable_without_work_reset": True,
            "spin_half_double_Dirac_continuum_candidate_classified": True,
            "native_covariant_nonscalar_family_feature_constructed": True,
            "entire_OPH_excluded_by_small_candidate": False,
            "observed_standard_model_species_identified": [],
            "physical_particle_identification_complete": False,
            "observed_mass_spectrum_source_derived": False,
            "gap_to_close": "Bind the decoded chiral internal current and a mass-generating action to the same propagating Spin/CAR response; the native nonscalar family feature alone is not that coupling.",
        },
        "source_bindings": _bindings(mass, family, public_bindings),
    }
    if result["source_bindings"] != initial_bindings:
        raise RuntimeError("source changed during particle audit")
    return result, arrays


def verify(path):
    path = Path(path)
    document = strict_json.load(path)
    if document.get("schema") != "oph.particle_identification.v1":
        raise ValueError("wrong particle audit schema")
    required = {p.relative_to(ROOT).as_posix() for p in HERE.glob("*.py") if not p.name.startswith("test_")}
    required.update(p.relative_to(ROOT).as_posix() for p in HERE.glob("*.md"))
    required.add(MASS_RECEIPT.relative_to(ROOT).as_posix())
    if not required <= set(document["source_bindings"]):
        raise ValueError("incomplete particle-audit source custody")
    for key, expected in document["source_bindings"].items():
        current = (ROOT/key).resolve()
        if not current.is_relative_to(ROOT) or _hash(current) != expected:
            raise ValueError(f"particle audit custody changed: {key}")
    artifact = document["artifact"]
    if Path(artifact["filename"]).name != artifact["filename"]:
        raise ValueError("output artifact must be adjacent to the receipt")
    array_path = path.with_name(artifact["filename"])
    if _hash(array_path) != artifact["sha256"]:
        raise ValueError("full quantum output artifact changed")
    with np.load(array_path, allow_pickle=False) as arrays:
        if (arrays["outputs"].shape != (64, 1024, 4) or arrays["abort_outputs"].shape != (12, 1024, 4)
                or set(arrays["records"]) != {format(i, "06b") for i in range(64)}):
            raise ValueError("incomplete retained quantum output")
        if not all(np.isfinite(arrays[name]).all() for name in ("outputs", "abort_outputs")):
            raise ValueError("nonfinite retained quantum state")
        probabilities = document["continuing_field_read"]["source_record_probabilities"]
        if set(probabilities) != set(arrays["records"]):
            raise ValueError("incomplete retained classical records")
        for record, state in zip(arrays["records"], arrays["outputs"], strict=True):
            probability = probabilities[str(record)]
            if not 0 <= probability <= 1 or abs(e.norm(state)-probability) > 1e-12:
                raise ValueError("retained state and record probability disagree")
        total = math.fsum(e.norm(state) for state in arrays["outputs"])
        total += math.fsum(e.norm(state) for state in arrays["abort_outputs"])
        if abs(total-1) > 1e-9:
            raise ValueError("quantum output lost probability mass")
    if (document["verdict"]["physical_particle_identification_complete"]
            or document["verdict"]["observed_standard_model_species_identified"]
            or document["verdict"]["observed_mass_spectrum_source_derived"]):
        raise ValueError("unsupported physical identification")
    kinetic = document["finite_inverse_response_kinetics"]
    if (kinetic.get("schema") != "oph.finite_car_inverse_kinetic.v1"
            or any(kinetic.get(flag) is not False for flag in (
                "response_inferred_from_observer_records", "spatial_momentum_derivatives_supplied",
                "positive_physical_kinetic_metric_established", "physical_mass_or_species_identification",
                "finite_shot_or_outward_numerical_certificate"))):
        raise ValueError("unsupported finite inverse-response identification")
    print("Verified particle-audit custody and all retained quantum output weights; no dynamics replay performed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        if args.output is not None:
            parser.error("verification does not write outputs")
        verify(args.verify)
        return
    if args.output is None or args.output.exists() or args.output.with_suffix(".npz").exists():
        parser.error("choose fresh --output JSON and companion NPZ paths")
    result, arrays = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    artifact = args.output.with_suffix(".npz")
    with artifact.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    result["artifact"] = {"filename": artifact.name, "sha256": _hash(artifact),
                          "contract": "all source work/reference outputs, normal records and terminal aborts"}
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"Wrote {args.output} and {artifact}")


if __name__ == "__main__":
    main()
