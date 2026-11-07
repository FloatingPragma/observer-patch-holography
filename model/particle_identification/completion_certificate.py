"""Replay the finite interfaces supporting the particle-identification route."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools.build_completion_program import integration_source_paths
from . import capacity, dense_effects, native_capacity, passive_clock, uncertainty
from ..native_causal_dataflow import dataflow
from ..native_lorentz_transition import certificate as lorentz
from ..field_cone_compatibility import certificate as field_cone
from ..source_state_bridge import bridge as state_bridge

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
                                   + dataflow.SOURCE_FILES + lorentz.SOURCE_FILES
                                   + field_cone.SOURCE_FILES),
    })
    paths.update(uncertainty.source_bindings())
    paths.update("model/particle_identification/" + name for name in (
        "DENSE_EFFECTS.md", "UNCERTAINTY.md", "SCATTERING.md", "ROADMAP.md"))
    paths.update(("Lean/EventAlgebra/FiniteBuschGleason.lean",
                  "Lean/EventAlgebra/OperationalBornBridge.lean",
                  "Lean/EventAlgebra/Robertson.lean",
                  "model/native_lorentz_transition/results.json",
                  "model/field_cone_compatibility/results.json"))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def report():
    before = source_bindings()
    frames, cones = lorentz.build(), field_cone.build()
    if frames != json.loads(lorentz.OUTPUT.read_text()):
        raise ValueError("retained Lorentz transition receipt differs")
    if cones != json.loads(field_cone.OUTPUT.read_text()):
        raise ValueError("retained field-cone receipt differs")
    result = {
        "schema": "oph.particle_identification.finite_interfaces.v1",
        "dense_rational_effects": dense_effects.build(),
        "compensated_probe_uncertainty": uncertainty.certificate(),
        "local_capacity_comparison": capacity.certificate(),
        "native_FIFO_capacity_equivalence": native_capacity.certificate(),
        "direct_native_state_and_carrier_bounds": state_bridge.certificate(),
        "passive_source_clock_boundary": passive_clock.certificate(),
        "causal_dataflow": dataflow.certificate(),
        "native_lorentz_transitions": frames,
        "field_cone_compatibility": cones,
        "scope": {
            "mathematical_prose_proofs": ["dense rational-effect extension",
                                           "complete local-history error composition",
                                           "passive source unitary obstruction",
                                           "native state and retained-carrier bounds",
                                           "finite null-transition and field-cone criteria"],
            "finite_fixtures_are_not_all_domain_proofs": True,
            "native_dense_effect_menu_derived": False,
            "active_record_assisted_encodings_excluded": False,
            "physical_capacity_to_local_action_map_derived": False,
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
    elif json.loads(args.check.read_text()) != value:
        raise ValueError("particle interface receipt differs from complete replay")
    print("Verified finite effect, uncertainty and local-history interfaces; physical joins remain explicit.")


if __name__ == "__main__":
    main()
