#!/usr/bin/env python3
"""Emit a compare-only exact Higgs/top reference adapter on the D11 Jacobian.

Chain role: expose the exact inverse-slice adapter that hits the canonical
Higgs/top reference targets exactly on the linear D11 readout, while keeping
the live predictive forward seed reference-free.

Mathematics: solve the linear D11 Jacobian readout against the canonical
reference pair `(m_t, m_H)` to get one exact compare-only adapter.

OPH-derived inputs: the emitted D11 critical-surface core and Jacobian plus the
repo-local canonical reference-value store.

Output: a machine-readable exact-reference adapter that is explicitly
non-promotable and compare-only.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
D11_CRITICAL_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_critical_surface_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d11_surface: dict, references: dict) -> dict:
    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])
    higgs_ref = dict(references["higgs"])
    top_ref = dict(references["top_quark"])

    m_h_core = float(core["mH_core_gev"])
    m_t_core = float(core["mt_pole_core_gev"])
    d_mh_d_lambda = float(jacobian["d_mH_d_lambda"])
    d_mt_d_y = float(jacobian["d_mt_pole_d_y_t"])

    m_h_target = float(higgs_ref["value_gev"])
    m_t_target = float(top_ref["value_gev"])

    delta_lambda = (m_h_target - m_h_core) / d_mh_d_lambda
    delta_y = (m_t_target - m_t_core) / d_mt_d_y

    predicted_m_h = m_h_core + d_mh_d_lambda * delta_lambda
    predicted_m_t = m_t_core + d_mt_d_y * delta_y

    pi_y = delta_y / float(core["y_t_core_mt"])
    pi_lambda = -(9.0 / 16.0) * delta_lambda / float(core["lambda_core_mt"])

    return {
        "artifact": "oph_d11_reference_exact_adapter",
        "generated_utc": _timestamp(),
        "proof_status": "compare_only_exact_reference_adapter",
        "scope": "compare_only_inverse_slice",
        "promotable": False,
        "source_artifact": str(D11_CRITICAL_JSON),
        "reference_source": str(REFERENCE_JSON),
        "live_predictive_branch_artifact": str(
            ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"
        ),
        "exact_reference_targets": {
            "mH_gev": m_h_target,
            "mt_pole_gev": m_t_target,
        },
        "inverse_slice_coordinates": {
            "delta_lambda_mt": delta_lambda,
            "delta_y_t_mt": delta_y,
        },
        "predicted_outputs": {
            "mH_gev": predicted_m_h,
            "mt_pole_gev": predicted_m_t,
        },
        "exact_fit_residuals_gev": {
            "mH_gev": predicted_m_h - m_h_target,
            "mt_pole_gev": predicted_m_t - m_t_target,
        },
        "normalized_readback": {
            "pi_y": pi_y,
            "pi_lambda": pi_lambda,
            "eta_HT": 0.5 * (pi_y - pi_lambda),
            "w_HT": pi_y - pi_lambda,
        },
        "formulas": {
            "delta_y_t_mt": "(mt_reference - mt_pole_core_gev) / d_mt_pole_d_y_t",
            "delta_lambda_mt": "(mH_reference - mH_core_gev) / d_mH_d_lambda",
            "mt_pole_gev": "mt_pole_core_gev + d_mt_pole_d_y_t * delta_y_t_mt",
            "mH_gev": "mH_core_gev + d_mH_d_lambda * delta_lambda_mt",
        },
        "strictly_not_claimed": [
            "live_forward_seed_equality",
            "predictive_reference_free_closure",
            "theorem_grade_higgs_top_exact_hit",
        ],
        "notes": [
            "This is a compare-only inverse slice on the linear D11 Jacobian, not the live predictive forward seed.",
            "It exists only to surface an exact Higgs/top reference hit on the current emitted D11 core.",
            "The live public Higgs/top rows remain the reference-free forward-seed outputs, not this inverse adapter.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the compare-only exact D11 Higgs/top reference adapter.")
    parser.add_argument("--d11-surface", default=str(D11_CRITICAL_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    artifact = build_artifact(d11_surface, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
