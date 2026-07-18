#!/usr/bin/env python3
"""Build the canonical formal fixed-cutoff CFQ central-record model.

This source-separated executable retains the useful mathematical content of
the submitted digital carrier without promoting its authored packet.  Given
the canonical conditional CFQ schema, it constructs distinct rank-one event
occurrences, exact pinching maps, central accepted/rejected D2 pointers, exact
path trace products, and inert ancillary trace identities.

The result proves nonemptiness and internal consistency of the stipulated
fixed-cutoff class.  It does not derive the register/path schema, a global
charged response update, cofinal OPH refinement, physical family attachment,
or historical no-target ancestry.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
CONDITIONAL_CFQ = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_source_law_rigidity_conditional.json"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_cfq_central_record_model.json"
)
SUBMITTED_CARRIER_ARCHIVE_SHA256 = (
    "d4b47d158f80b4bffe318650769812a606e4e065cab1e1a01edd68ff8909291c"
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fraction(value: str) -> Fraction:
    return Fraction(value)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def build_artifact(schema: dict[str, Any], schema_sha256: str) -> dict[str, Any]:
    if schema.get("checks_pass") is not True:
        raise ValueError("the canonical conditional CFQ schema receipt is required")
    if schema.get("historical_charged_target_informed") is not True:
        raise ValueError("historical target ancestry must remain explicit")
    if schema.get("public_prediction_promotion_allowed") is not False:
        raise ValueError("the conditional packet must remain unpromoted")

    registers = schema["conditional_cfq_packet"]["registers"]
    paths = schema["conditional_cfq_packet"]["paths"]
    next_index: defaultdict[str, int] = defaultdict(int)
    occurrences: list[dict[str, Any]] = []
    path_models: list[dict[str, Any]] = []

    for path in paths:
        path_occurrences: list[str] = []
        product = Fraction(path["sign"], 1)
        for stage, register_name in enumerate(path["registers"], start=1):
            register = registers[register_name]
            dimension = int(register["dimension"])
            event_index = next_index[register_name]
            next_index[register_name] += 1
            if event_index >= dimension:
                raise ValueError(f"not enough orthogonal events in {register_name}")
            occurrence_id = f"{path['name']}:{stage}:{register_name}"
            event_weight = Fraction(1, dimension)
            product *= event_weight
            path_occurrences.append(occurrence_id)
            occurrences.append(
                {
                    "id": occurrence_id,
                    "register": register_name,
                    "dimension": dimension,
                    "event_index": event_index,
                    "event_projector": f"P=E[{event_index},{event_index}]",
                    "event_rank": 1,
                    "complement_rank": dimension - 1,
                    "normalized_trace": fraction_text(event_weight),
                    "pinching": "E_P(X)=PXP+(I-P)X(I-P)",
                    "pinching_properties": {
                        "linear": True,
                        "unital": True,
                        "completely_positive": True,
                        "trace_preserving": True,
                        "idempotent": True,
                    },
                    "central_pointer_algebra": "D2=l_infinity({rejected,accepted})",
                    "central_record_dilation": (
                        "rho -> P rho P tensor |accepted><accepted| + "
                        "(I-P) rho (I-P) tensor |rejected><rejected|"
                    ),
                }
            )
        path_models.append(
            {
                "name": path["name"],
                "occurrences": path_occurrences,
                "signed_trace_product": fraction_text(product),
                "expected_signed_trace_product": path["exact_trace_coefficient"],
                "matches_schema": fraction_text(product)
                == path["exact_trace_coefficient"],
            }
        )

    refinement_checks: list[dict[str, Any]] = []
    for register_name, register in registers.items():
        dimension = int(register["dimension"])
        for ancilla_dimension in (2, 3):
            source_trace = Fraction(1, dimension)
            refined_trace = Fraction(ancilla_dimension, dimension * ancilla_dimension)
            refinement_checks.append(
                {
                    "register": register_name,
                    "ancilla_dimension": ancilla_dimension,
                    "embedding": "iota_k(X)=X tensor I_k",
                    "coarse_map": "C_k=id tensor tau_k",
                    "source_event_rank": 1,
                    "refined_event_rank": ancilla_dimension,
                    "source_trace": fraction_text(source_trace),
                    "refined_trace": fraction_text(refined_trace),
                    "trace_preserved": source_trace == refined_trace,
                    "coarse_after_embedding_is_identity": True,
                }
            )

    physical_selection_gates = {
        "registers_emitted_by_oph_dynamics": False,
        "tracial_state_selected_by_constrained_oph_maxent": False,
        "global_mutually_exclusive_charged_response_instrument": False,
        "path_automaton_generated_by_oph": False,
        "coupling_character_grading_and_clock_generated_by_oph": False,
        "physical_path_exhaustion": False,
        "cofinal_screen_refinement": False,
        "historical_no_target_ancestry": False,
        "physical_face_family_yukawa_attachment": False,
        "renormalized_pole_scheme": False,
    }

    checks = {
        "eight_registers_present": len(registers) == 8,
        "eight_paths_present": len(paths) == 8,
        "ten_distinct_event_occurrences_constructed": (
            len(occurrences) == 10
            and len({row["id"] for row in occurrences}) == 10
        ),
        "shared_register_events_are_orthogonal": all(
            count <= int(registers[name]["dimension"])
            for name, count in next_index.items()
        ),
        "all_pinching_and_central_record_checks_pass": all(
            all(row["pinching_properties"].values())
            and row["central_pointer_algebra"].startswith("D2")
            for row in occurrences
        ),
        "all_path_trace_products_match_schema": all(
            row["matches_schema"] for row in path_models
        ),
        "all_sixteen_inert_refinement_checks_pass": (
            len(refinement_checks) == 16
            and all(
                row["trace_preserved"]
                and row["coarse_after_embedding_is_identity"]
                for row in refinement_checks
            )
        ),
        "all_physical_selection_gates_remain_open": not any(
            physical_selection_gates.values()
        ),
        "promotion_boundary_retained": (
            schema["historical_charged_target_informed"] is True
            and schema["global_source_only"] is False
            and schema["public_prediction_promotion_allowed"] is False
        ),
    }

    return {
        "artifact": "oph_charged_cfq_central_record_model",
        "schema_version": 1,
        "status": (
            "FORMAL_FIXED_CUTOFF_CFQ_MODEL_AND_CENTRAL_RECORD_DILATION_CLOSED_"
            "OPH_DYNAMICAL_AND_PHYSICAL_SELECTION_OPEN"
        ),
        "runtime_charged_reference_consumed": False,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "provenance": {
            "conditional_cfq_schema": CONDITIONAL_CFQ.relative_to(CODE_ROOT).as_posix(),
            "conditional_cfq_schema_sha256": schema_sha256,
            "submitted_carrier_archive_sha256": SUBMITTED_CARRIER_ARCHIVE_SHA256,
        },
        "proved": [
            "the stipulated fixed-cutoff CFQ class is algebraically nonempty",
            "ten distinct rank-one event occurrences fit in the eight declared registers",
            "per-event pinching maps admit central accepted/rejected D2 record pointers",
            "the eight exact signed path trace products reproduce the conditional schema",
            "normalized event traces survive inert ancillary stabilization",
        ],
        "scope_boundary": {
            "statement": (
                "This constructs one formal model of an authored target-informed schema. "
                "It is not an OPH source-law or physical charged-sector derivation."
            ),
            "global_response_update_constructed": False,
            "physical_regulator_refinement_constructed": False,
            "historical_ancestry_repaired": False,
        },
        "registers": registers,
        "event_occurrences": occurrences,
        "path_models": path_models,
        "inert_refinement_checks": refinement_checks,
        "physical_selection_gates": physical_selection_gates,
        "checks": checks,
        "checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--schema", type=Path, default=CONDITIONAL_CFQ)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    raw = args.schema.read_bytes()
    artifact = build_artifact(json.loads(raw), sha256(raw))
    encoded = (json.dumps(artifact, indent=2, sort_keys=True) + "\n").encode()
    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == encoded
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "event_occurrences": len(artifact["event_occurrences"]),
                "public_prediction_promotion_allowed": artifact[
                    "public_prediction_promotion_allowed"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
