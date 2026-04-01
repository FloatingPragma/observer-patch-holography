#!/usr/bin/env python3
"""Record the scalarized D12 quark continuation bundle.

Chain role: gather the strongest current D12 mass-side and mixing-side
continuation diagnostics in one place without promoting them as OPH-forward
closure.

Mathematics: two-scalar mass branch on the ordered family plus gauge-fixed
physical invariant decomposition of the same-label left-transport generator.

OPH-derived inputs: the D12 overlap transport law, the quadratic-even scalar
shell, and the forward-emitted D12 same-label transport closure.

Output: one diagnostic bundle summarizing the strongest current D12
continuation math and the exact remaining mass-side value laws.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_QUADRATIC = ROOT / "particles" / "runs" / "flavor" / "quark_quadratic_even_transport_scalar.json"
DEFAULT_PHYSICAL = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_same_label_physical_invariant_bundle.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_scalarized_continuation_bundle.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the scalarized D12 quark continuation bundle.")
    parser.add_argument("--quadratic", default=str(DEFAULT_QUADRATIC))
    parser.add_argument("--physical", default=str(DEFAULT_PHYSICAL))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    quadratic = json.loads(Path(args.quadratic).read_text(encoding="utf-8"))
    physical = json.loads(Path(args.physical).read_text(encoding="utf-8"))
    payload = {
        "artifact": "oph_quark_scalarized_continuation_bundle",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "mixing_closed_mass_value_laws_open_on_d12_continuation",
        "mass_side": {
            "selector_scalar_name": "Delta_ud_overlap",
            "odd_transport_scalar_name": "Lambda_ud_B_transport",
            "quadratic_scalar_name": "eta_Q_centered",
            "B_ord": [-1.0, 0.0, 1.0],
            "Q_ord": quadratic["Q_ord"],
            "x2": quadratic["x2"],
            "Q_ord_formula": quadratic["Q_ord_formula"],
            "eta_Q_projector_formula_from_centered_residual": quadratic["eta_Q_projector_formula"],
            "kappa_Q_formula_from_eta_Q": quadratic["kappa_Q_formula"],
            "quadratic_even_log_formula_direct": quadratic["quadratic_even_log_formula_direct"],
            "odd_transport_formula": "(sigma_u_total_log_per_side * sigma_d_total_log_per_side / (2 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))) * Delta_ud_overlap",
            "tau_u_formula": "(Delta_ud_overlap / 2) * sigma_d_total_log_per_side / (sigma_u_total_log_per_side + sigma_d_total_log_per_side)",
            "tau_d_formula": "(Delta_ud_overlap / 2) * sigma_u_total_log_per_side / (sigma_u_total_log_per_side + sigma_d_total_log_per_side)",
            "best_theorem_mean_surface": {
                "Delta_ud_overlap": 0.14049998075678738,
                "Lambda_ud_B_transport": 0.14591158065891177,
                "eta_Q_centered": -0.0181047533864264,
                "kappa_Q": -0.012364871560017752,
                "tau_u_log_per_side": 0.026099916046670564,
                "tau_d_log_per_side": 0.04415007433172313,
                "m_u_GeV": [0.002158247734638124, 1.2719518494594315, 172.23296806851607],
                "m_d_GeV": [0.004660608415841266, 0.09271974673707856, 4.147651052785027],
                "rms_log_mass_error": 0.005984702111797554,
            },
            "best_exact_mean_specialization": {
                "alpha_u": 1.0007763698011345,
                "alpha_d": 1.008463281557513,
                "Delta_ud_overlap": 0.14049991320632976,
                "Lambda_ud_B_transport": 0.1459115105066314,
                "eta_Q_centered": -0.018104730181494357,
                "kappa_Q": -0.01236485571191432,
                "tau_u_log_per_side": 0.026099903498190035,
                "tau_d_log_per_side": 0.04415005310497486,
                "m_u_GeV": [0.0021599234039183404, 1.272939334772138, 172.36668155105923],
                "m_d_GeV": [0.004700052692857295, 0.09350445861546268, 4.182753646232867],
                "rms_log_mass_error": 5.2130513964914564e-05,
            },
            "sample_same_family_ray_point": {
                "ray_modulus": 0.6695617711471163,
                "t1_sample": 0.6695617711471163,
                "Delta_ud_overlap_formula": "ray_modulus / 5",
                "Delta_ud_overlap": 0.13391235422942327,
                "kappa_Q_formula": "-ray_modulus / 54",
                "eta_Q_centered_formula": "2 * (1 - x2^2) * kappa_Q = -((1 - x2^2) / 27) * ray_modulus",
                "eta_Q_centered": -0.018155152181872827,
                "theorem_mean_rms": 0.006597329560808134,
                "exact_mean_rms": 0.0027768355284593363,
                "status": "sample_only_not_theorem",
            },
        },
        "mixing_side": {
            "full_matrix_artifact": physical["artifact"],
            "same_label_generator_formula": "V_CKM^fwd = U_u^dagger @ U_d, K_CKM = Log_pr(V_CKM^fwd)",
            "gauge_fixed_physical_decomposition": dict(
                physical["physical_invariants"],
                **physical["generator_invariants"],
            ),
            "physical_invariant_formula": physical["physical_invariant_formula"],
            "transport_closure_residual_fro_norm": physical["closure_residual_fro_norm"],
        },
        "honest_remaining_value_laws": [
            "quark_d12_t1_value_law",
            "quark_exact_mean_split_value_law_or_carrier_repair",
        ],
        "notes": [
            "The D12 quark mass branch reduces to two scalar laws, Delta_ud_overlap and eta_Q_centered, once the current ordered-family carrier is fixed.",
            "On the stricter same-family branch, those two scalars collapse further to the emitted D12_ud_mass_ray with unresolved ray_modulus; the retained numerical point is sample-only.",
            "The CKM/CP side closes on the D12 continuation branch because the same-label transport unitary is emitted directly by the forward Yukawa step and its principal logarithm exists uniquely on the standard gauge representative.",
            "This bundle is diagnostic only and does not alter the live public quark rows.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
