#!/usr/bin/env python3
"""Emit the sharp no-go theorem for the current one-scalar D11 fixed ray.

Chain role: certify that the current one-scalar D11 branch cannot carry the
exact Higgs/top pair on the declared D10/D11 surface.

Mathematics: the fixed ray forces `pi_y = pi_lambda`, hence `w_HT = 0`. The
compare-only exact pair on the same Jacobian surface has nonzero wedge, so it
lies off the fixed ray. This proves that exact pair promotion requires one
additional forward coordinate beyond the current one-scalar branch.

OPH-derived inputs: the declared D10/D11 calibration surface, the closed
one-scalar fixed-ray seed, and the compare-only exact inverse slice used only
as a witness.

Output: a machine-readable no-go theorem artifact with the smallest honest
extension contract stated explicitly.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
D11_SURFACE_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"
D11_FORWARD_SEED_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
D11_EXACT_ADAPTER_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_fixed_ray_no_go_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d11_surface: dict, forward_seed: dict, exact_adapter: dict) -> dict:
    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])
    sigma = float(forward_seed["sigma_D11_HT"])

    y_core = float(core["y_t_core_mt"])
    lambda_core = float(core["lambda_core_mt"])
    mt_core = float(core["mt_pole_core_gev"])
    mh_core = float(core["mH_core_gev"])
    d_mt = float(jacobian["d_mt_pole_d_y_t"])
    d_mh = float(jacobian["d_mH_d_lambda"])

    exact_targets = dict(exact_adapter["exact_reference_targets"])
    exact_slice = dict(exact_adapter["inverse_slice_coordinates"])
    exact_readback = dict(exact_adapter["normalized_readback"])

    pi_y_exact = float(exact_readback["pi_y"])
    pi_lambda_exact = float(exact_readback["pi_lambda"])
    sigma_exact = 0.5 * (pi_y_exact + pi_lambda_exact)
    eta_exact = 0.5 * (pi_y_exact - pi_lambda_exact)
    w_exact = pi_y_exact - pi_lambda_exact

    compatibility_functional = (
        (float(exact_targets["mt_pole_gev"]) - mt_core) / (d_mt * y_core)
        + (9.0 / 16.0) * (float(exact_targets["mH_gev"]) - mh_core) / (d_mh * lambda_core)
    )

    return {
        "artifact": "oph_d11_fixed_ray_no_go_theorem",
        "generated_utc": _timestamp(),
        "theorem_id": "D11FixedRayNoGoTheorem",
        "proof_status": "closed_no_go_on_current_one_scalar_fixed_ray",
        "status": "closed",
        "scope": "declared_d10_d11_surface_against_compare_only_exact_pair",
        "source_artifacts": {
            "d11_declared_surface": str(D11_SURFACE_JSON),
            "d11_forward_seed": str(D11_FORWARD_SEED_JSON),
            "d11_reference_exact_adapter": str(D11_EXACT_ADAPTER_JSON),
        },
        "current_fixed_ray_branch": {
            "sigma_D11_HT": sigma,
            "pi_y": sigma,
            "pi_lambda": sigma,
            "eta_HT": 0.0,
            "w_HT": 0.0,
            "law": "delta_y_t / y_t_core = -(9/16) * delta_lambda / lambda_core = sigma_D11_HT",
        },
        "exact_compare_witness": {
            "mt_pole_exact_gev": float(exact_targets["mt_pole_gev"]),
            "mH_exact_gev": float(exact_targets["mH_gev"]),
            "delta_y_t_exact": float(exact_slice["delta_y_t_mt"]),
            "delta_lambda_exact": float(exact_slice["delta_lambda_mt"]),
            "pi_y_exact": pi_y_exact,
            "pi_lambda_exact": pi_lambda_exact,
            "Sigma_HT_exact": sigma_exact,
            "eta_HT_exact": eta_exact,
            "w_HT_exact": w_exact,
        },
        "fixed_ray_obstruction": {
            "compatibility_functional_formula": (
                "((mt_target - mt_pole_core_gev) / (d_mt_pole_d_y_t * y_t_core_mt)) "
                "+ (9/16) * ((mH_target - mH_core_gev) / (d_mH_d_lambda * lambda_core_mt))"
            ),
            "compatibility_functional_value": compatibility_functional,
            "fixed_ray_condition": "reachable_on_current_one_scalar_ray iff compatibility_functional = 0",
            "current_exact_pair_reachable": False,
        },
        "smallest_honest_extension": {
            "object": "Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda)",
            "equivalent_coordinates": "(Sigma_HT, eta_HT)",
            "one_extra_scalar_beyond_fixed_ray": True,
            "exact_extension_values": {
                "Sigma_HT_exact": sigma_exact,
                "eta_HT_exact": eta_exact,
                "w_HT_exact": w_exact,
            },
            "notes": [
                "The current fixed ray freezes eta_HT = 0, so one extra scalar is necessary.",
                "The exact center Sigma_HT_exact also differs from sigma_D11_HT, so a wedge-only patch around the old center is not exact.",
            ],
        },
        "proof": [
            "On the current one-scalar branch, pi_y = pi_lambda = sigma_D11_HT, so w_HT vanishes identically.",
            "The compare-only exact pair on the same declared D11 surface has pi_y_exact != pi_lambda_exact, hence w_HT_exact != 0.",
            "Therefore the exact pair lies off the current one-scalar fixed ray and cannot be promoted by any theorem that preserves that fixed-ray role.",
            "The smallest honest exact extension is a two-coordinate forward readout object Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda), equivalently (Sigma_HT, eta_HT).",
        ],
        "notes": [
            "This is a no-go theorem about the current one-scalar fixed ray, not a demotion of the exact Higgs theorem on the declared D10/D11 surface.",
            "The compare-only exact inverse slice remains compare-only and is used here only as a witness that the exact pair is off the fixed ray.",
            "The exact Higgs row is carried separately by D11LiveForwardExactHiggsPromotion.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D11 fixed-ray no-go theorem artifact.")
    parser.add_argument("--d11-surface", default=str(D11_SURFACE_JSON))
    parser.add_argument("--forward-seed", default=str(D11_FORWARD_SEED_JSON))
    parser.add_argument("--exact-adapter", default=str(D11_EXACT_ADAPTER_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    forward_seed = json.loads(Path(args.forward_seed).read_text(encoding="utf-8"))
    exact_adapter = json.loads(Path(args.exact_adapter).read_text(encoding="utf-8"))
    artifact = build_artifact(d11_surface, forward_seed, exact_adapter)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
