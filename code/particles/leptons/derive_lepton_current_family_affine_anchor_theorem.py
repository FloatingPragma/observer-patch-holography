#!/usr/bin/env python3
"""Emit the exact same-family charged affine-anchor theorem artifact."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"
QUADRATIC_THEOREM_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_quadratic_readout_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(exact_readout: dict, quadratic_theorem: dict) -> dict:
    masses = [float(value) for value in exact_readout["predicted_singular_values_abs"]]
    centered_logs = [float(value) for value in exact_readout["centered_log_shape_exact"]]
    affine_anchor = math.log(math.prod(masses)) / 3.0
    geometric_mean = math.exp(affine_anchor)
    determinant = math.prod(masses)
    reconstructed = [math.exp(affine_anchor + value) for value in centered_logs]

    return {
        "artifact": "oph_lepton_current_family_affine_anchor_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_affine_anchor",
        "theorem_scope": "current_family_only",
        "public_promotion_allowed": False,
        "supporting_exact_readout_artifact": exact_readout["artifact"],
        "supporting_quadratic_readout_theorem": quadratic_theorem["artifact"],
        "theorem_statement": (
            "On the fixed exact same-family charged witness, the unique affine coordinate is "
            "A_ch^(cf) = (1/3) log det(Y_e^(cf)) = mean(log m_i). "
            "Equivalently g_e^(cf) = exp(A_ch^(cf)), and the exact charged triple is recovered by "
            "m_i = exp(A_ch^(cf) + ell_i^ctr) from the centered log shape."
        ),
        "current_family_affine_anchor": {
            "name": "A_ch_current_family",
            "formula": "(1/3) * log(det(Y_e_current_family)) = mean(log(m_e), log(m_mu), log(m_tau))",
            "value": affine_anchor,
        },
        "current_family_geometric_mean": {
            "name": "g_e_current_family",
            "formula": "exp(A_ch_current_family)",
            "value": geometric_mean,
        },
        "current_family_determinant": {
            "formula": "det(Y_e_current_family) = m_e * m_mu * m_tau",
            "value": determinant,
        },
        "centered_log_shape_exact": centered_logs,
        "predicted_singular_values_abs": masses,
        "reconstructed_from_affine_anchor": reconstructed,
        "exact_fit_residuals_abs": [
            reconstructed[idx] - masses[idx]
            for idx in range(3)
        ],
        "do_not_claim_now": [
            "theorem-grade A_ch on the live charged theorem lane",
            "theorem-grade mu_phys(Y_e) on the live charged theorem lane",
        ],
        "notes": [
            "This closes the exact same-family affine coordinate on the target-anchored charged witness.",
            "It does not promote the theorem-grade charged absolute lane, which remains behind C_hat_e^{cand} promotion and the descended scalar mu_phys(Y_e).",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact same-family charged affine-anchor theorem artifact.")
    parser.add_argument("--exact-readout", default=str(EXACT_READOUT_JSON))
    parser.add_argument("--quadratic-theorem", default=str(QUADRATIC_THEOREM_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    exact_readout = json.loads(Path(args.exact_readout).read_text(encoding="utf-8"))
    quadratic_theorem = json.loads(Path(args.quadratic_theorem).read_text(encoding="utf-8"))
    artifact = build_artifact(exact_readout, quadratic_theorem)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
