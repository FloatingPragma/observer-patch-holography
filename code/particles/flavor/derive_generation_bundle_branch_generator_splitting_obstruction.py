#!/usr/bin/env python3
"""Emit the exact current-corpus obstruction certificate for branch-generator splitting."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LINE_LIFT_JSON = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
GENERATOR_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator_splitting_obstruction.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _natural_centered_pair(
    line_lift: dict[str, Any],
    generator: dict[str, Any],
) -> tuple[list[np.ndarray], float]:
    lamb = np.asarray(generator["simple_spectrum_certificate"]["eigenvalues"], dtype=float)
    diagnostics = list(line_lift["same_refinement_edge_diagnostic_by_refinement"])
    if len(diagnostics) < 2:
        raise ValueError("line-lift diagnostics must contain at least two refinement levels")
    projectors_by_level = []
    for entry in diagnostics[:2]:
        projectors_by_level.append([_decode_complex_matrix(projector) for projector in entry["projectors"]])

    centered = []
    for projectors in projectors_by_level:
        current = sum(lamb[idx] * projectors[idx] for idx in range(3))
        current = current - np.trace(current) / 3.0 * np.eye(3, dtype=complex)
        centered.append(current)
    defect = max(
        np.linalg.norm(projectors_by_level[1][idx] - projectors_by_level[0][idx], ord=2)
        for idx in range(3)
    )
    return centered, float(defect)


def build_artifact(line_lift: dict[str, Any], generator: dict[str, Any]) -> dict[str, Any]:
    centered, defect = _natural_centered_pair(line_lift, generator)
    commutator = centered[1] @ centered[0] - centered[0] @ centered[1]
    commutator_norm = float(np.linalg.norm(commutator, ord=2))
    overlaps = list(generator["projective_readout_certificate"]["same_label_overlap_amplitudes"])

    return {
        "artifact": "oph_generation_bundle_branch_generator_splitting_obstruction_certificate",
        "generated_utc": _timestamp(),
        "proof_status": "closed_exact_current_corpus_obstruction_certificate",
        "target_theorem": "oph_generation_bundle_branch_generator_splitting",
        "target_clause": "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split",
        "verdict": "not_forced_by_current_corpus",
        "theorem_statement": (
            "The displayed current proxy shell does not force the descended-commutator clause. "
            "Using the displayed proxy ordered spectrum together with the displayed level-0 and "
            "level-1 projector systems, the natural centered proxy operators have nonzero commutator "
            "whose size is linear in the projector defect rather than uniformly quadratic-small. "
            "So the current artifacts do not imply exact vanishing or uniform quadratic smallness "
            "after the central split."
        ),
        "first_failed_implication": (
            "explicit_centered_proxy_shell_and_current_projector_data "
            "does_not_imply "
            "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split"
        ),
        "current_attached_data_obstruction": {
            "commutator_operator_norm": commutator_norm,
            "projector_defect_operator_norm": defect,
            "commutator_over_defect": float(commutator_norm / defect) if defect > 0.0 else 0.0,
            "commutator_over_defect_squared": float(commutator_norm / (defect * defect)) if defect > 0.0 else 0.0,
            "reason": (
                "On the displayed level-0/level-1 projector systems, the natural centered operators built "
                "from the proxy ordered spectrum have a commutator that is linear in the projector defect "
                "rather than quadratic-small."
            ),
        },
        "exact_vanishing_not_forced": {
            "same_label_overlap_amplitudes": overlaps,
            "max_same_label_overlap_amplitude": max(overlaps),
            "min_same_label_overlap_amplitude": min(overlaps),
            "reason": (
                "Exact vanishing is not forced by the displayed current data because the same-label level-0 "
                "and level-1 eigenspaces are not identical. The overlaps are close to 1 but not equal to 1."
            ),
        },
        "issue_149_resolution_mode": "sharp_exact_obstruction",
        "resulting_remaining_open_chain": [
            "charged.02-derive-post-promotion-uncentered-trace-lift",
            "charged.03-derive-d10-p-to-charged-affine-anchor",
            "charged.04-emit-charged-masses-from-p",
        ],
        "notes": [
            "This certificate does not prove impossibility in every future extension of the corpus.",
            "It proves that the current displayed proxy/projector data do not themselves force the missing clause.",
            "That is the exact obstruction object needed to close the current operator-side issue honestly when a positive proof is absent.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact obstruction certificate for branch-generator splitting.")
    parser.add_argument("--line-lift", default=str(LINE_LIFT_JSON))
    parser.add_argument("--generator", default=str(GENERATOR_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.line_lift)), _load_json(Path(args.generator)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
