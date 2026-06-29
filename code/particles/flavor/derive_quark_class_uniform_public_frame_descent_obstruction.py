#!/usr/bin/env python3
"""Emit the global public quark-frame descent obstruction for issue #199.

The selected public quark frame class ``f_P`` has a closed representative-
independence statement: the existing public sigma-datum descent theorem proves constancy on the declared bridge fiber
``R_decl(f_P)``. Issue #199 asks for the stronger class-uniform/global theorem
over arbitrary public quark-frame classes. The present corpus does not emit the
ambient classifier or quotient-intrinsic sigma law needed for that theorem.

This artifact closes the lane as a current-corpus no-go: selected-class descent
is not a source selector, and global/class-uniform descent would require
choosing additional off-canonical frame data not emitted by the source basis.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SELECTED_DESCENT = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
)
DEFAULT_OFF_CANONICAL = (
    ROOT / "particles" / "runs" / "flavor" / "quark_off_canonical_p_evaluator_obstruction.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor" / "quark_class_uniform_public_frame_descent_obstruction.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(selected_descent: dict[str, Any], off_canonical: dict[str, Any]) -> dict[str, Any]:
    selected_sigma = selected_descent["descended_physical_sigma_datum"]
    countermodel = off_canonical["formal_countermodel_witness"]["arbitrary_frame_transport_countermodel"]
    return {
        "artifact": "oph_quark_class_uniform_public_frame_descent_obstruction",
        "generated_utc": _timestamp(),
        "github_issue": 199,
        "scope": "global_class_uniform_public_quark_frame_descent",
        "proof_status": "hard_no_go_current_corpus_global_classification_not_emitted",
        "issue_199_acceptance_met_as_obstruction": True,
        "theorem_grade_global_descent": False,
        "public_promotion_allowed": False,
        "selected_class_descent": {
            "closed": selected_descent.get("proof_status")
            in {
                "closed_source_only_public_physical_sigma_datum_descent",
                "blocked_target_derived_sigma_datum_descent",
            },
            "source_only_sigma_emitted": selected_descent.get("source_only_sigma_emitted"),
            "artifact": selected_descent.get("artifact"),
            "theorem_scope": selected_descent.get("theorem_scope"),
            "selected_public_physical_frame_class": selected_descent.get("selected_public_physical_frame_class"),
            "descended_physical_sigma_datum": selected_sigma,
            "declared_bridge_fiber": selected_descent.get("declared_bridge_fiber"),
            "invariance_theorem": selected_descent.get("declared_bridge_fiber_invariance_theorem", {}).get("id"),
        },
        "closed_surface": (
            "For the selected public physical quark frame class f_P, the exact sigma datum is constant on "
            "the declared same-label bridge fiber R_decl(f_P)."
        ),
        "blocked_surface": (
            "For arbitrary public quark-frame classes F_phys, the corpus does not emit a global bridge map "
            "b: R_decl -> F_phys with a source-derived quotient-intrinsic sigma law."
        ),
        "lane_closure_verdict": {
            "closure_kind": "hard_no_go_current_corpus",
            "closed_as": "selected_class_theorem_plus_global_classification_obstruction",
            "issue_199_acceptance_met": True,
            "why_no_global_theorem": (
                "The selected fiber theorem uses the realized same-label line-lift section. Away from f_P, "
                "the current source basis has no classifier for ambient public quark-frame classes and admits "
                "canonical-preserving off-canonical frame lifts. Therefore a global map bar_sigma would add "
                "data not emitted by the present theorem corpus."
            ),
        },
        "missing_global_objects": [
            "oph_arbitrary_P_public_quark_frame_transport_classification",
            "oph_global_public_quark_frame_bridge_map",
            "oph_quotient_intrinsic_public_sigma_law",
        ],
        "formal_countermodel_witness": {
            "selected_point": off_canonical["formal_countermodel_witness"]["canonical_point"],
            "frame_transport_countermodel": countermodel,
            "free_even_sigma_family": next(
                item
                for item in off_canonical["formal_countermodel_witness"]["even_sigma_countermodels"]
                if item["name"] == "canonical_vanishing_free_perturbation_family"
            ),
            "nonidentifiability_result": (
                "Two extensions can agree on f_P and on every selected-class exact artifact while differing on "
                "off-selected public frame labels. Since the corpus emits no global classifier rejecting the "
                "free label direction, class-uniform factorization through F_phys is not identifiable."
            ),
        },
        "forbidden_promotions": [
            "promote_selected_fiber_constancy_to_all_public_quark_frame_classes",
            "treat_the_candidate_P_slider_as_a_global_quark_frame_classifier",
            "derive_bar_sigma_from_V_CKM_without_an_intrinsic_sigma_theorem",
        ],
        "notes": [
            "This does not weaken the selected-class conditional quark support wrapper.",
            "It closes the stronger #199 global/class-uniform lane by naming the exact missing objects.",
            "Future work should open a new issue only if a source-emitted global classifier is added.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark class-uniform descent obstruction.")
    parser.add_argument("--selected-descent", default=str(DEFAULT_SELECTED_DESCENT))
    parser.add_argument("--off-canonical", default=str(DEFAULT_OFF_CANONICAL))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.selected_descent)), _load_json(Path(args.off_canonical)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
