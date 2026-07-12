#!/usr/bin/env python3
"""Emit the current-family target-audit affine-coordinate identity."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_affine_anchor_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(exact_readout: dict) -> dict:
    g_u = float(exact_readout["g_u"])
    g_d = float(exact_readout["g_d"])
    masses_u = [float(value) for value in exact_readout["predicted_singular_values_u"]]
    masses_d = [float(value) for value in exact_readout["predicted_singular_values_d"]]
    exact_u = [float(value) for value in exact_readout["E_u_log_exact"]]
    exact_d = [float(value) for value in exact_readout["E_d_log_exact"]]

    affine_anchor = 0.5 * (math.log(g_u) + math.log(g_d))
    sector_split = 0.5 * (math.log(g_u) - math.log(g_d))
    shared_sector_mean = math.exp(affine_anchor)
    reconstructed_g_u = math.exp(affine_anchor + sector_split)
    reconstructed_g_d = math.exp(affine_anchor - sector_split)
    reconstructed_u = [reconstructed_g_u * math.exp(value) for value in exact_u]
    reconstructed_d = [reconstructed_g_d * math.exp(value) for value in exact_d]

    return {
        "artifact": "oph_quark_current_family_affine_anchor_theorem",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_affine_anchor",
        "theorem_scope": "current_family_only",
        "public_promotion_allowed": False,
        "target_anchored": True,
        "source_only_prediction": False,
        "single_running_quark_sextet_claim_allowed": False,
        "supporting_exact_readout_artifact": exact_readout["artifact"],
        "theorem_statement": (
            "On the fixed target-anchored current-family quark witness, the common affine coordinate is "
            "A_q^(cf) = (1/2) * (log g_u + log g_d) = log(sqrt(g_u g_d)). Equivalently "
            "delta_q^(cf) = (1/2) * (log g_u - log g_d), so g_u = exp(A_q^(cf) + delta_q^(cf)) and "
            "g_d = exp(A_q^(cf) - delta_q^(cf)). Together with the exact centered ordered shapes "
            "E_u^log and E_d^log, this reconstructs the chosen mixed-convention target coordinates."
        ),
        "current_family_affine_anchor": {
            "name": "A_q_current_family",
            "formula": "(1/2) * (log(g_u) + log(g_d)) = log(sqrt(g_u * g_d))",
            "value": affine_anchor,
        },
        "current_family_sector_split": {
            "name": "delta_q_current_family",
            "formula": "(1/2) * (log(g_u) - log(g_d))",
            "value": sector_split,
        },
        "current_family_shared_sector_mean": {
            "name": "g_q_current_family",
            "formula": "exp(A_q_current_family) = sqrt(g_u * g_d)",
            "value": shared_sector_mean,
        },
        "current_family_sector_means": {
            "g_u_exact": g_u,
            "g_d_exact": g_d,
            "reconstructed_g_u": reconstructed_g_u,
            "reconstructed_g_d": reconstructed_g_d,
        },
        "centered_log_shapes_exact": {
            "E_u_log_exact": exact_u,
            "E_d_log_exact": exact_d,
        },
        "predicted_singular_values_exact": {
            "u_sector": masses_u,
            "d_sector": masses_d,
        },
        "reconstructed_from_affine_anchor": {
            "u_sector": reconstructed_u,
            "d_sector": reconstructed_d,
        },
        "exact_fit_residuals": {
            "u_sector": [
                reconstructed_u[idx] - masses_u[idx]
                for idx in range(3)
            ],
            "d_sector": [
                reconstructed_d[idx] - masses_d[idx]
                for idx in range(3)
            ],
        },
        "do_not_claim_now": [
            "theorem-grade target-free A_q^phys on the live quark theorem lane",
            "theorem-grade physical-sheet promotion of the current-family affine anchor",
        ],
        "notes": [
            "This closes the exact affine coordinate on the target-anchored current-family quark witness.",
            "It does not promote Theta_ud^abs, which still requires a target-free physical-sheet carrier A_q^phys.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact current-family quark affine-anchor theorem artifact.")
    parser.add_argument("--exact-readout", default=str(EXACT_READOUT_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    exact_readout = json.loads(Path(args.exact_readout).read_text(encoding="utf-8"))
    artifact = build_artifact(exact_readout)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
