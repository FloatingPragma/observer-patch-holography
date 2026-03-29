#!/usr/bin/env python3
"""Emit a current-family ordered charged-lepton shape candidate."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_mean_eigenvalue_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_ordered_shape_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _softmax(values: list[float]) -> list[float]:
    anchor = max(values)
    unnormalized = [math.exp(value - anchor) for value in values]
    denom = sum(unnormalized)
    return [value / denom for value in unnormalized]


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a current-family ordered charged-lepton shape readout.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    eigenvalues = [float(value) for value in payload["current_family_eigenvalues"]]
    mean_value = sum(eigenvalues) / len(eigenvalues)
    delta_e = [value - mean_value for value in eigenvalues]
    shape_weights = _softmax(delta_e)

    artifact = {
        "artifact": "oph_lepton_ordered_shape_readout",
        "generated_utc": _timestamp(),
        "parent_artifact": "charged_shared_absolute_scale_writeback",
        "source_mean_artifact": payload.get("artifact"),
        "proof_status": "diagnostic_only_softmax_consumer",
        "predictive_promotion_allowed": False,
        "theorem_scope": "shared_budget_only",
        "writeback_scope": "current_family_only",
        "labels": ["f1", "f2", "f3"],
        "center_rule": "delta_e = current_family_eigenvalues - mean(current_family_eigenvalues)",
        "shape_rule": "shape_weights_i = exp(delta_e_i) / sum_j exp(delta_e_j)",
        "current_family_eigenvalues": eigenvalues,
        "delta_e": delta_e,
        "shape_weights_current_order": shape_weights,
        "expected_singular_values_shape_sorted": sorted(shape_weights),
        "Y_e_shape": {
            "real": [
                [shape_weights[0], 0.0, 0.0],
                [0.0, shape_weights[1], 0.0],
                [0.0, 0.0, shape_weights[2]],
            ],
            "imag": [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
        },
        "shape_closed": True,
        "consumer_builder": "build_forward_charged_leptons.py",
        "expected_closure_state": "shared_writeback_current_family",
        "metadata": {
            "note": "Current-family ordered charged-lepton softmax readout built from the populated common-refinement eigenvalue triple. This is diagnostic-only: the centered-softmax primitive is too weak to carry the observed charged hierarchy and should not be promoted as a predictive mass consumer.",
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
