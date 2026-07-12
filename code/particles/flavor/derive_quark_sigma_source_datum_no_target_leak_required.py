#!/usr/bin/env python3
"""Project the target-free quark-spread obstruction onto the promotion gate.

The detailed obstruction proves that the source-compatible spread fiber is
``(R_{>0})^2``.  This compatibility artifact keeps the historical gate path
stable for downstream builders while removing target values and fitted
correction residuals from the source-boundary packet.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SIGMA_OBSTRUCTION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_sigma_source_nonidentifiability_obstruction.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_sigma_source_datum_no_target_leak_required.json"

MISSING_FOR_PROMOTION = [
    "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM",
    "QUARK_SOURCE_SPREAD_SECTOR_ATTACHMENT_AND_REFINEMENT",
    "NO_TARGET_LEAK_DAG_QUARK_SOURCE_SPREAD",
]

def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(obstruction: dict[str, Any]) -> dict[str, Any]:
    dependency_audit = dict(obstruction["dependency_audit"])
    classification = dict(obstruction["exact_ray_classification"])
    future_extension = dict(obstruction["minimal_future_extension"])
    return {
        "artifact": "oph_quark_sigma_source_datum_no_target_leak_required",
        "generated_utc": _timestamp(),
        "status": "closed_current_corpus_nonidentifiability_obstruction",
        "proof_status": "closed_exact_current_corpus_obstruction",
        "closure_kind": "theorem_grade_sharper_obstruction",
        "claim_tier": "selected_class_conditional_support_source_sigma_nonidentifiable",
        "promotion_allowed": False,
        "public_promotion_allowed": False,
        "source_only_sigma_emitted": False,
        "downstream_algebra_closed": True,
        "selected_class": "f_P",
        "required_identity": "P -> f_P -> Sigma_ud_phys -> (sigma_u,sigma_d,sigma_seed_ud,eta_ud)",
        "resolved_github_issues": [377, 379, 380],
        "obstruction_artifact": obstruction["artifact"],
        "obstruction_theorem_statement": obstruction["theorem_statement"],
        "compatible_spread_fiber": classification["fiber"],
        "compatible_spread_fiber_dimension": classification["fiber_dimension"],
        "independent_unselected_coordinates": classification["independent_coordinates"],
        "affine_downstream_injective": obstruction["affine_downstream_injectivity"]["injective"],
        "template_ancestry": obstruction["template_ancestry"],
        "conditional_downstream_if_closed": {
            "sigma_seed_ud": "(sigma_u + sigma_d)/2",
            "eta_ud": "(sigma_u - sigma_d)/2",
            "g_u": "g_ch*exp(-(A_ud*sigma_seed_ud - B_ud*eta_ud))",
            "g_d": "g_ch*exp(-(A_ud*sigma_seed_ud + B_ud*eta_ud))",
            "then": "ordered_three_point_mass_coordinate_readout_then_declared_RG_chart_and_dimensionless_Yukawa_normalization",
        },
        "missing_for_promotion": MISSING_FOR_PROMOTION,
        "reopen_requirements": future_extension,
        "forbidden_ancestors": dependency_audit["forbidden_ancestors"],
        "dependency_audit": dependency_audit,
        "required_theorem_package": {
            "QΣ-0": "sigma descent is representative independence, not source selection",
            "QΣ-1": "the target-free compatible spread fiber is (R_{>0})^2",
            "QΣ-2": "conditional quark masses follow once a no-target source sigma datum is supplied",
            "QΣ-A": "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM",
            "QΣ-B": "QUARK_SOURCE_SPREAD_SECTOR_ATTACHMENT_AND_REFINEMENT",
            "QΣ-C": "NO_TARGET_LEAK_DAG_QUARK_SOURCE_SPREAD",
        },
        "notes": [
            "The obstruction consumes no running-quark target, exact-witness, spread-map, or fitted-correction artifact.",
            "The upstream edge packet leaves two coefficients free even when granted as exact.",
            "Issue closure in obstruction mode does not promote numeric quark or physical Yukawa rows.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark source-spread obstruction gate artifact.")
    parser.add_argument("--obstruction", default=str(SIGMA_OBSTRUCTION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.obstruction)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
