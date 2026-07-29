#!/usr/bin/env python3
"""Issue-506 fixed-verdict builder for the same-scheme alpha/HVP audit.

The builder evaluates the four-class acceptance matrix under the frozen
class protocol.  The tabulated dispersive class replays the recorded
instrument: the pinned endpoint payload and anchor-scheme bridge are
loaded by byte hash, the reference deficit is recomputed from the
recorded on-shell decomposition, and containment in the certified
same-scheme gap interval is re-decided with exact decimal arithmetic.
The three classes without frozen repository ingests receive explicit
not-evaluable verdicts carrying their exact ingest requirements, and
the cross-class agreement row is not evaluable with a single class.

Nothing here promotes an empirical input to a source output or emits a
physical alpha prediction; the guard block records both refusals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

PACKET_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKET_DIR.parents[2]

PROTOCOL_PATH = PACKET_DIR / "data" / "frozen_class_protocol_v1.json"
ENDPOINT_PATH = (
    REPO_ROOT / "code" / "P_derivation" / "runtime"
    / "empirical_thomson_endpoint_current.json"
)
BRIDGE_PATH = (
    REPO_ROOT / "code" / "P_derivation" / "runtime"
    / "anchor_scheme_bridge_current.json"
)
OUT_PATH = PACKET_DIR / "outputs" / "alpha_hvp_class_verdict.json"

SCHEMA = "oph.alpha_hvp_class_verdict.v1"
ISSUE = 506

getcontext().prec = 60


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def load_pinned(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    return json.loads(raw.decode("utf-8")), sha256_bytes(raw)


def evaluate_tabulated_class(
    endpoint: dict[str, Any], bridge: dict[str, Any]
) -> dict[str, Any]:
    """Replay the recorded tabulated-class comparison with exact decimals."""

    gap_interval = bridge["certified_gap_from_endpoint"][
        "same_scheme_anchor_gap_interval"
    ]
    lower = Decimal(str(gap_interval[0]))
    upper = Decimal(str(gap_interval[1]))

    reference = bridge["reference_decomposition_compare_only"]
    alpha_inv_mz_phys = Decimal(str(reference["alpha_inv_mz_phys_on_shell"]))
    anchor = Decimal(str(bridge["anchor_provenance"]["a0_oph"]))
    recomputed_deficit = alpha_inv_mz_phys - anchor
    recorded_deficit = Decimal(str(reference["gap_phys_minus_oph"]))
    deficit_agrees = abs(recomputed_deficit - recorded_deficit) < Decimal(
        "1e-9"
    )

    endpoint_interval = endpoint["endpoint"]["alpha_inv_interval"]
    endpoint_central = Decimal(endpoint["endpoint"]["alpha_inv_central"])
    endpoint_contained = (
        Decimal(endpoint_interval[0])
        <= endpoint_central
        <= Decimal(endpoint_interval[1])
    )

    inside = lower <= recomputed_deficit <= upper
    recorded_flag = bool(
        bridge["verdict"]["reference_deficit_inside_certified_gap"]
    )
    return {
        "class": "tabulated_dispersive",
        "evaluation_mode": "replay_of_recorded_instrument",
        "endpoint_alpha_inv_central": str(endpoint_central),
        "endpoint_interval": [str(v) for v in endpoint_interval],
        "endpoint_central_inside_interval": bool(endpoint_contained),
        "certified_gap_interval": [str(lower), str(upper)],
        "recomputed_reference_deficit": str(recomputed_deficit),
        "recorded_reference_deficit": str(recorded_deficit),
        "deficit_recomputation_agrees": bool(deficit_agrees),
        "deficit_inside_certified_gap": bool(inside),
        "recorded_flag_agrees": bool(inside == recorded_flag),
        "class_verdict": (
            "CONFIRMED_WITHIN_FROZEN_BAND"
            if inside and deficit_agrees and endpoint_contained
            else "REFUTED_OR_INCONSISTENT"
        ),
    }


def build_verdict() -> dict[str, Any]:
    protocol, protocol_sha = load_pinned(PROTOCOL_PATH)
    endpoint, endpoint_sha = load_pinned(ENDPOINT_PATH)
    bridge, bridge_sha = load_pinned(BRIDGE_PATH)

    guards_ok = bool(
        endpoint["guards"]["promotable_as_oph_source_theorem"] is False
        and endpoint["guards"]["measured_alpha_in_solve_path"] is False
        and bridge["guards"]["public_promotion_allowed"] is False
        and bridge["guards"]["measured_values_in_any_oph_solve_path"] is False
    )

    tabulated = evaluate_tabulated_class(endpoint, bridge)

    class_matrix: dict[str, Any] = {"tabulated_dispersive": tabulated}
    for name in ("raw_dispersive", "independent_code", "lattice_hvp"):
        row = protocol["classes"][name]
        class_matrix[name] = {
            "class": name,
            "evaluation_mode": row["evaluation_mode"],
            "class_verdict": "NOT_EVALUABLE_MISSING_FROZEN_INGEST",
            "ingest_requirement": row["ingest_requirement"],
            "fabrication_excluded": True,
        }

    evaluated = [
        row
        for row in class_matrix.values()
        if row["class_verdict"] == "CONFIRMED_WITHIN_FROZEN_BAND"
    ]
    cross_class = {
        "evaluated_class_count": len(evaluated),
        "verdict": (
            "NOT_EVALUABLE_SINGLE_CLASS"
            if len(evaluated) == 1
            else "NOT_EVALUABLE_NO_CLASS"
        ),
        "note": (
            "cross-class agreement requires at least two independently "
            "evaluated classes under the shared protocol"
        ),
    }

    overall = (
        "TABULATED_CLASS_CONFIRMED__OTHER_CLASSES_NOT_EVALUABLE"
        if tabulated["class_verdict"] == "CONFIRMED_WITHIN_FROZEN_BAND"
        and guards_ok
        else "PROTOCOL_INCONSISTENT"
    )

    payload = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "row_class": "empirical_same_scheme_falsification_instrument",
        "pins": {
            "protocol_sha256": protocol_sha,
            "endpoint_sha256": endpoint_sha,
            "bridge_sha256": bridge_sha,
        },
        "guards": {
            "payload_guards_verified": guards_ok,
            "empirical_input_promoted_to_source_output": False,
            "physical_alpha_prediction_emitted": False,
        },
        "class_matrix": class_matrix,
        "cross_class_agreement": cross_class,
        "verdict": overall,
        "reactivation": {
            "raw_dispersive": protocol["classes"]["raw_dispersive"][
                "ingest_requirement"
            ],
            "independent_code": protocol["classes"]["independent_code"][
                "ingest_requirement"
            ],
            "lattice_hvp": protocol["classes"]["lattice_hvp"][
                "ingest_requirement"
            ],
            "rule": (
                "landing any frozen ingest reopens only its class row; the "
                "frozen band and shared protocol stay unchanged"
            ),
        },
        "claim_boundary": (
            "Fixed multi-class verdict of the empirical same-scheme "
            "alpha/HVP audit under the frozen protocol: the tabulated "
            "dispersive class replays confirmed inside the certified band "
            "that was declared before this evaluation, the raw-dispersive, "
            "independent-code, and lattice classes are explicitly not "
            "evaluable for want of frozen repository ingests, and "
            "cross-class agreement is not evaluable with a single class. "
            "No empirical input is promoted to a source output and no "
            "physical alpha prediction is emitted."
        ),
    }
    payload["verdict_sha256"] = sha256_bytes(
        canonical_json(
            {k: v for k, v in payload.items() if k != "verdict_sha256"}
        ).encode("utf-8")
    )
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    payload = build_verdict()
    if args.verify:
        stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        if stored != payload:
            print("VERDICT_DRIFT", file=sys.stderr)
            return 1
        print("VERDICT_VERIFIED")
        return 0
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n", encoding="utf-8"
    )
    print(json.dumps({"verdict": payload["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
