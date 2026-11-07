"""Replay the finite interfaces supporting the particle-identification route."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import capacity, dense_effects, native_capacity, passive_clock, uncertainty
from ..native_causal_dataflow import dataflow
from ..native_lorentz_transition import certificate as lorentz
from ..field_cone_compatibility import certificate as field_cone
from ..source_state_bridge import bridge as state_bridge
from ..protected_reference import reference as protected_reference
from ..gauge_defect_transport import transport as gauge_transport
from ..native_capacity_restriction import packet as capacity_restriction
from ..native_record_response import response as record_response
from ..native_record_causality import packet as record_causality
from ..native_gauge_strings import certificate as gauge_strings
from ..native_clock_observables import packet as clock_observables

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "completion_results.json"


def source_bindings():
    paths = integration_source_paths(ROOT, {
        "domain_and_proof": "model/particle_identification/CAPACITY.md",
        "source_receipts": ["model/interacting_extension/higgs_join/results.json"],
        "independent_checks": list(capacity.SOURCE_FILES + passive_clock.SOURCE_FILES
                                   + native_capacity.SOURCE_FILES
                                   + state_bridge.SOURCE_FILES
                                   + protected_reference.SOURCE_FILES
                                   + gauge_transport.SOURCE_FILES
                                   + dataflow.SOURCE_FILES + lorentz.SOURCE_FILES
                                   + field_cone.SOURCE_FILES),
    })
    paths.update(uncertainty.source_bindings())
    for producer in (capacity_restriction, record_response, record_causality,
                     gauge_strings, clock_observables):
        paths.update(producer.source_bindings())
        paths.add(str(producer.OUTPUT.relative_to(ROOT)))
    paths.update("model/particle_identification/" + name for name in (
        "DENSE_EFFECTS.md", "UNCERTAINTY.md", "SCATTERING.md", "ROADMAP.md"))
    paths.update(("Lean/EventAlgebra/FiniteBuschGleason.lean",
                  "Lean/EventAlgebra/OperationalBornBridge.lean",
                  "Lean/EventAlgebra/Robertson.lean",
                  "model/native_lorentz_transition/results.json",
                  "model/gauge_defect_transport/results.json",
                  "model/protected_reference/results.json",
                  "model/field_cone_compatibility/results.json"))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def report():
    before = source_bindings()
    frames, cones = lorentz.build(), field_cone.build()
    gauge = gauge_transport.certificate()
    protected = protected_reference.certificate()
    restrictions = capacity_restriction.produce()
    responses = record_response.report()
    causal_records = record_causality.produce()
    strings = gauge_strings.produce()
    observables = clock_observables.produce()
    for producer, value in ((capacity_restriction, restrictions),
                            (record_response, responses),
                            (record_causality, causal_records),
                            (gauge_strings, strings),
                            (clock_observables, observables)):
        if value != strict_json.load(producer.OUTPUT):
            raise ValueError(f"retained {producer.OUTPUT.parent.name} receipt differs")
    if frames != strict_json.load(lorentz.OUTPUT):
        raise ValueError("retained Lorentz transition receipt differs")
    if cones != strict_json.load(field_cone.OUTPUT):
        raise ValueError("retained field-cone receipt differs")
    if gauge != strict_json.load(ROOT / "model/gauge_defect_transport/results.json"):
        raise ValueError("retained gauge-defect transport receipt differs")
    if protected != strict_json.load(ROOT / "model/protected_reference/results.json"):
        raise ValueError("retained protected-reference receipt differs")
    result = {
        "schema": "oph.particle_identification.finite_interfaces.v1",
        "dense_rational_effects": dense_effects.build(),
        "compensated_probe_uncertainty": uncertainty.certificate(),
        "local_capacity_comparison": capacity.certificate(),
        "native_FIFO_capacity_equivalence": native_capacity.certificate(),
        "direct_native_state_and_carrier_bounds": state_bridge.certificate(),
        "protected_native_reference_and_gauge_join": protected,
        "complete_gauge_defect_transport": gauge,
        "fair_source_canonical_restriction": restrictions,
        "weighted_native_complete_record_response": responses,
        "disjoint_native_complete_record_equivalence": causal_records,
        "declared_gauge_tree_propagation": strings,
        "native_coherent_clock_observables": observables,
        "passive_source_clock_boundary": passive_clock.certificate(),
        "causal_dataflow": dataflow.certificate(),
        "native_lorentz_transitions": frames,
        "field_cone_compatibility": cones,
        "scope": {
            "mathematical_prose_proofs": ["dense rational-effect extension",
                                           "complete local-history error composition",
                                           "passive source unitary obstruction",
                                           "native state and retained-carrier bounds",
                                           "protected stationary states and modular repair covariance",
                                           "complete finite gauge transport and conditional common reference",
                                           "fair canonical restriction and exact nesting",
                                           "actual weighted-source complete record response",
                                           "complete disjoint-command exchange with references",
                                           "declared tree Gauss-sector exhaustion and propagation",
                                           "reflected native record transfer and clock distinction",
                                           "adaptive independent-product dynamic restriction",
                                           "native coherent Pauli instruments and finite uncertainty",
                                           "finite null-transition and field-cone criteria"],
            "finite_fixtures_are_not_all_domain_proofs": True,
            "native_dense_effect_menu_derived": False,
            "active_record_assisted_encodings_excluded": False,
            "physical_capacity_to_local_action_map_derived": False,
            "fair_federation_and_weighted_clock_are_one_source_family": False,
            "native_record_relaxation_factors_are_particle_masses": False,
            "record_transfer_and_working_modular_clocks_identified": False,
            "coherent_clock_services_lifted_to_growing_family": False,
            "common_native_to_gauge_Higgs_intertwiner_derived": False,
            "physical_position_momentum_identified": False,
            "relativistic_interacting_particle_spectrum_identified": False,
            "physical_collider_predictions_added": 0,
            "canonical_physical_premises_discharged": [],
        },
        "source_sha256": before,
    }
    if source_bindings() != before:
        raise RuntimeError("particle interface sources changed during replay")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    mode.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.write is not None:
        args.write.write_text(json.dumps(value, sort_keys=True, indent=2,
                                        allow_nan=False) + "\n")
    elif strict_json.load(args.check) != value:
        raise ValueError("particle interface receipt differs from complete replay")
    print("Verified finite effect, uncertainty and local-history interfaces; physical joins remain explicit.")


if __name__ == "__main__":
    main()
