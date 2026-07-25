#!/usr/bin/env python3
"""Emit the standard-math-fixed neutrino transport-load selector.

Chain role: fix the one-dimensional positive load selector that sits between
the emitted cocycle loads `chi = 1 + eps` and `1 + gamma_half`.

Mathematics: on an affine one-dimensional segment, the balanced barycentric
selector and the least-distortion selector for any positive
translation-invariant quadratic form coincide at the unique midpoint. This
uses only the live endpoint loads and no oscillation targets or external mass
anchors.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the neutrino transport-load segment selector.")
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    cocycle = _load_json(Path(args.cocycle))
    gamma = float(cocycle["theorem_gap_gamma"])
    eps = float(cocycle["defect_gap_ratio"])
    gamma_half = float(cocycle["hermitian_descendant_riesz_margin"]["gamma_half"])
    chi = 1.0 + eps
    right_load = 1.0 + gamma_half
    selected_d_nu = 0.5 * (chi + right_load)

    payload = {
        "artifact": "oph_neutrino_transport_load_segment_selector",
        "generated_utc": _timestamp(),
        "source_artifact": str(Path(args.cocycle)),
        "status": "standard_math_fixed_balanced_segment_selector_closed",
        "closure_scope": "standard_math_fixed",
        "segment_kind": "positive_affine_load_segment",
        "segment_definition": "D_tau = (1-tau_nu) * chi + tau_nu * (1 + gamma_half)",
        "endpoints": {
            "left_name": "chi",
            "left_value": chi,
            "right_name": "1 + gamma_half",
            "right_value": right_load,
        },
        "selector_family": ["balanced", "least_distortion"],
        "selector_equivalence_statement": (
            "On a one-dimensional affine segment, the balanced barycentric selector and the "
            "least-distortion selector for any positive translation-invariant quadratic form "
            "coincide uniquely at the midpoint."
        ),
        "selected_selector": "balanced_equals_least_distortion_midpoint",
        "selected_tau_nu": 0.5,
        "selected_D_nu": selected_d_nu,
        "selectors": {
            "balanced": {
                "tau_nu": 0.5,
                "D_nu": selected_d_nu,
            },
            "least_distortion": {
                "tau_nu": 0.5,
                "D_nu": selected_d_nu,
                "ambient_metric_kind": "positive_translation_invariant_quadratic_form_on_segment",
            },
        },
        "derived_quantities": {
            "gamma": gamma,
            "gamma_half": gamma_half,
            "eps": eps,
            "chi": chi,
            "weight_exponent_formula": "p = 1 + gamma + eps / D_nu",
            "weight_exponent_value": 1.0 + gamma + eps / selected_d_nu,
        },
        "notes": [
            "The selector uses only the emitted positive endpoint loads chi and 1 + gamma_half.",
            "No PMNS target, oscillation target, or external mass anchor enters the selector.",
            "In one dimension the midpoint is independent of the positive overall scale of the quadratic form, so the balanced and least-distortion selectors agree exactly.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
