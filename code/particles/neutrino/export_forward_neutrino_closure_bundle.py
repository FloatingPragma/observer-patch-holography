#!/usr/bin/env python3
"""Bundle the active neutrino forward lane into one closure-surface artifact.

Chain role: collect the emitted neutrino matrix, splitting, and selector-law
artifacts into the forward bundle used by audits and public-surface gating.

Mathematics: packaging only; this file does not derive new masses, but it keeps
the theorem tier, selector certificates, and PMNS blocker explicit.

OPH-derived inputs: the local Majorana scalar, matrix, splitting, phase, and
holonomy artifacts.

Output: the forward neutrino closure bundle for downstream reporting.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCALAR = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"
DEFAULT_MATRIX = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
DEFAULT_SPLITTINGS = ROOT / "particles" / "runs" / "neutrino" / "forward_splittings.json"
DEFAULT_PHASE = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"
DEFAULT_HOLONOMY = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the forward neutrino closure-bundle artifact.")
    parser.add_argument("--scalar", default=str(DEFAULT_SCALAR))
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--splittings", default=str(DEFAULT_SPLITTINGS))
    parser.add_argument("--phase", default=str(DEFAULT_PHASE))
    parser.add_argument("--holonomy", default=str(DEFAULT_HOLONOMY))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    scalar = json.loads(Path(args.scalar).read_text(encoding="utf-8"))
    matrix = json.loads(Path(args.matrix).read_text(encoding="utf-8"))
    splittings = json.loads(Path(args.splittings).read_text(encoding="utf-8"))
    phase = json.loads(Path(args.phase).read_text(encoding="utf-8"))
    holonomy = json.loads(Path(args.holonomy).read_text(encoding="utf-8"))

    payload = {
        "artifact": "oph_forward_neutrino_closure_bundle",
        "generated_utc": _timestamp(),
        "closure_tier": "selector_law_standard_math_forward",
        "theorem_candidate_id": scalar.get("theorem_candidate_id"),
        "sublemma_candidate_id": scalar.get("sublemma_candidate_id"),
        "bundle_descent_candidate_id": scalar.get("bundle_descent_candidate_id"),
        "required_overlap_certificate": scalar.get("required_overlap_certificate"),
        "selector_center": scalar.get("selector_center"),
        "selector_point_absolute": scalar.get("selector_point_absolute"),
        "formula_affine_candidate": scalar.get("formula_affine_candidate"),
        "norm_identification_formula": scalar.get("norm_identification_formula"),
        "phase_mode": matrix.get("phase_mode") or "canonical_selector",
        "selector_point_certified": bool(matrix.get("selector_point_certified")),
        "selector_law_certified": bool(splittings.get("selector_law_certified", matrix.get("selector_law_certified"))),
        "certification_status": splittings.get("certification_status", "selector_law_standard_math_forward"),
        "phase_certificate_source": str(Path(args.phase)),
        "projector_certificate_source": str(Path(args.holonomy)),
        "majorana_matrix_real": matrix.get("majorana_matrix_real"),
        "majorana_matrix_imag": matrix.get("majorana_matrix_imag"),
        "u_nu_left_real": matrix.get("u_nu_left_real"),
        "u_nu_left_imag": matrix.get("u_nu_left_imag"),
        "collective_mode_dominance": splittings.get("collective_mode_dominance"),
        "masses_gev_sorted": splittings.get("masses_gev_sorted"),
        "delta_m21_sq_gev2": splittings.get("delta_m21_sq_gev2"),
        "delta_m31_sq_gev2": splittings.get("delta_m31_sq_gev2"),
        "splitting_ratio_r": splittings.get("splitting_ratio_r"),
        "ordering_phase_certified": splittings.get("ordering_phase_certified"),
        "oph_only_promotion_gate": {
            "proof_status": scalar.get("proof_status"),
            "quadraticity_clause_status": scalar.get("quadraticity_clause_status"),
            "overlap_nonvanishing_status": scalar.get("overlap_nonvanishing_status"),
            "overlap_nonvanishing_witness_hint": scalar.get("overlap_nonvanishing_witness_hint"),
        },
        "pmns_status": "blocked_pending_shared_charged_lepton_left_basis",
        "notes": [
            "This wrapper promotes the current centered selector-law branch to the forward-usable neutrino sandbox surface under standard-math certification while keeping OPH-only promotion blocked behind the overlap certificate.",
            "The live forward numerical defect remains the exact 1-2 near-degeneracy from isotropic mu_nu, so the next mass-moving object is still the defect-weighted mu_e family.",
            "Flavor-labeled neutrino rows in the public table remain qualitative until the shared charged-lepton left basis and PMNS orientation close.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
