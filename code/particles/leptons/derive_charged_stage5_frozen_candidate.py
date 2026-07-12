#!/usr/bin/env python3
"""Emit the reference-free numerical Stage-5 charged-lepton candidate.

This builder intentionally has no charged-mass reference input.  It freezes the
historical D12/Stage-5 continuation in one auditable artifact while preserving
the distinction between a numerical candidate and a theorem-grade prediction.
"""

from __future__ import annotations

import argparse
import json
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "P_derivation"))

from paper_math import PaperMathContext  # noqa: E402


D10 = (
    ROOT / "particles" / "runs" / "calibration"
    / "d10_ew_forward_transmutation_certificate.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_stage5_frozen_candidate.json"
)
FROZEN_PIXEL_P = Decimal("1.63094")
FROZEN_V_GEV = Decimal("246.76711732749683")


def build_artifact(d10: dict[str, Any]) -> dict[str, Any]:
    v = Decimal(str(d10["source_only_reconstruction"]["v_from_source_transmutation_gev"]))
    p = Decimal(str(d10["oph_inputs"]["p"]))
    if p != FROZEN_PIXEL_P or v != FROZEN_V_GEV:
        raise ValueError(
            "The legacy Stage-5 receipt is frozen at "
            f"P={FROZEN_PIXEL_P}, v={FROZEN_V_GEV} GeV; received P={p}, v={v}. "
            "Create a new candidate rather than mutating this one."
        )
    context = PaperMathContext(precision=80, su2_cutoff=1, su3_cutoff=1, n_c=3, n_g=3)
    masses = context.charged_lepton_masses(v)
    vectors = context.stage5_vectors
    roots = [str(value) for value in vectors["roots"]]
    mass_rows = {
        "electron": {"value_gev": str(masses["e"]), "status": "frozen_retrodictive_candidate"},
        "muon": {"value_gev": str(masses["mu"]), "status": "frozen_retrodictive_candidate"},
        "tau": {"value_gev": str(masses["tau"]), "status": "frozen_retrodictive_candidate"},
    }
    return {
        "artifact": "oph_charged_stage5_frozen_candidate",
        "status": "FROZEN_ACCURATE_NUMERICAL_CANDIDATE_NOT_THEOREM_GRADE",
        "public_prediction_promotion_allowed": False,
        "reference_data_consumed_by_builder": False,
        "inputs": {
            "frozen_pixel_P": str(p),
            "v_from_source_transmutation_gev": str(v),
            "v_source_artifact": "code/particles/runs/calibration/d10_ew_forward_transmutation_certificate.json",
            "v_source_family": d10.get("family_source_id", "unknown"),
            "N_c": 3,
            "N_g": 3,
            "beta_EW": 4,
        },
        "frozen_construction": {
            "epsilon_Z6": str(vectors["epsilon"]),
            "delta_e": str(vectors["delta"]),
            "delta_formula": "(N_c+1)/(2*N_c*N_g)",
            "balanced_carrier": "C_delta = I + (exp(i delta) R + exp(-i delta) R^2)/sqrt(2)",
            "ordered_roots": roots,
            "mass_direction": "m_i proportional to r_i^2",
            "integer_exponents_n_e": list(vectors["n_e"]),
            "determinant_law": "det(M_e) = v(P)^3 / (2 * 6^14)",
            "geometric_mean_mass_gev": str(v / (Decimal(2) * Decimal(6) ** 14) ** (Decimal(1) / Decimal(3))),
        },
        "candidate_mass_rows": mass_rows,
        "theorem_grade_audit": {
            "closed_inputs": [
                "N_c=3",
                "N_g=3",
                "beta_EW=N_c+1",
                "D10 source-side v(P)",
            ],
            "unproved_or_historically_post_hoc_inputs": [
                "balanced circulant charged-family carrier",
                "phase transport law delta=(N_c+1)/(2*N_c*N_g)",
                "integer exponent prescription n_e=(7,4,3)",
                "historical normalization factor 2^(1/6)",
                "normalized determinant attachment det(M_e)=v(P)^3/(2*6^14)",
                "ordered-root attachment to electron, muon, tau",
            ],
            "reason_not_promoted": (
                "The builder is reference-free, but its frozen formula family was developed "
                "with the charged spectrum already known. Reference-free evaluation is not "
                "the same as a prospective blind prediction."
            ),
        },
        "branch_and_scheme_boundary": {
            "frozen_branch": "legacy D10 running-tree branch at P=1.63094",
            "current_public_pixel_used": False,
            "mass_scheme": "unspecified historical Stage-5 mass coordinate",
            "pole_mass_comparison_allowed_here": False,
            "reason": (
                "The formula predates the current public pixel branch and supplies no charged-lepton "
                "RG or threshold theorem mapping its mass coordinate to pole masses. Branch and pole-mass "
                "comparisons belong only in the downstream audit."
            ),
        },
        "falsification_contract": (
            "The formula is now frozen. Future changes to improve agreement require a new "
            "version and must not overwrite this candidate."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--d10", type=Path, default=D10)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact(json.loads(args.d10.read_text(encoding="utf-8")))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(artifact["candidate_mass_rows"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
