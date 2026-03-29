#!/usr/bin/env python3
"""Build the D10 electroweak observable family from the OPH calibration core.

Chain role: lift the paper D10 core into the baseline electroweak observables
that every downstream calibration branch consumes.

Mathematics: paper D10 fixed-point solve plus tree-level electroweak algebra for
`MW`, `MZ`, `alpha_em`, `sin^2(theta_W)`, and the shared Higgs vev.

OPH-derived inputs: `P`, `alpha_u`, `mu_u`, and `n_c` from the D10 core.

Output: the baseline calibration artifact that the reduced two-scalar D10
transport chain factorizes into source, selector, and readout objects.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _find_ancillary_particle_code_dir() -> Path:
    for base in Path(__file__).resolve().parents:
        candidate = base / "arXiv" / "RC1" / "ancillary" / "code" / "particles"
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate arXiv/RC1/ancillary/code/particles")


PARTICLE_CODE_DIR = _find_ancillary_particle_code_dir()
for candidate in [PARTICLE_CODE_DIR / "core", PARTICLE_CODE_DIR]:
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from particle_masses_paper_d10_d11 import (  # type: ignore
    P_DEFAULT,
    alpha_em_from_alpha1_alpha2,
    build_paper_d10,
    sin2_theta_w,
    solve_mz_fixed_point_tree,
)


DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(p_value: float) -> dict[str, object]:
    d10 = build_paper_d10(pix_area=p_value)
    alpha_em = alpha_em_from_alpha1_alpha2(d10.alpha1_mz, d10.alpha2_mz)
    sin2 = sin2_theta_w(d10.alpha1_mz, d10.alpha2_mz)
    m_w_run = 0.5 * d10.v * math.sqrt(4.0 * math.pi * d10.alpha2_mz)
    mz_run_check, v_check, a1_check, a2_check, a3_check = solve_mz_fixed_point_tree(d10.alpha_u, d10.p, d10.n_c, d10.mu_u)
    delta_rho_stage3 = 3.0 / (32.0 * math.pi ** 2)
    m_z_pole_stage3 = d10.mz_run / math.sqrt(1.0 + delta_rho_stage3)
    running_mass_ratio_residual = (1.0 - d10.m_w_run**2 / d10.mz_run**2) - d10.sin2w_mz
    stage3_mass_ratio_residual = (1.0 - d10.m_w_run**2 / d10.m_z_pole_stage3**2) - d10.sin2w_mz

    running_outputs = {
        "m_w_run": d10.m_w_run,
        "m_z_run": d10.mz_run,
        "alpha_em_inv_mz": 1.0 / d10.alpha_em_mz,
        "sin2w_mz": d10.sin2w_mz,
        "v": d10.v,
    }

    return {
        "artifact": "oph_d10_ew_observable_family",
        "generated_utc": _timestamp(),
        "observable_family_id": "d10_running_tree",
        "p": d10.p,
        "core_source": {
            "p": d10.p,
            "mu_u": d10.mu_u,
            "alpha_u": d10.alpha_u,
            "mz_run": d10.mz_run,
            "v": d10.v,
            "alpha1_mz": d10.alpha1_mz,
            "alpha2_mz": d10.alpha2_mz,
            "alpha3_mz": d10.alpha3_mz,
        },
        "generator_keys_unrounded": [
            "p",
            "mu_u",
            "alpha_u",
            "mz_run",
            "v",
            "alpha1_mz",
            "alpha2_mz",
            "alpha3_mz",
        ],
        "generator_values_unrounded": {
            "p": d10.p,
            "mu_u": d10.mu_u,
            "alpha_u": d10.alpha_u,
            "mz_run": d10.mz_run,
            "v": d10.v,
            "alpha1_mz": d10.alpha1_mz,
            "alpha2_mz": d10.alpha2_mz,
            "alpha3_mz": d10.alpha3_mz,
        },
        "identity_residuals": {
            "alpha_em_from_alpha1_alpha2": alpha_em - d10.alpha_em_mz,
            "sin2_from_alpha1_alpha2": sin2 - d10.sin2w_mz,
            "mw_from_v_alpha2": m_w_run - d10.m_w_run,
            "mz_fixed_point_residual": mz_run_check - d10.mz_run,
            "v_fixed_point_residual": v_check - d10.v,
            "alpha1_fixed_point_residual": a1_check - d10.alpha1_mz,
            "alpha2_fixed_point_residual": a2_check - d10.alpha2_mz,
            "alpha3_fixed_point_residual": a3_check - d10.alpha3_mz,
        },
        "mixed_family_residuals": {
            "mz_stage3_surrogate_residual": m_z_pole_stage3 - d10.m_z_pole_stage3,
            "sin2_from_mass_ratio_running_minus_reported": running_mass_ratio_residual,
            "sin2_from_mass_ratio_stage3_minus_reported": stage3_mass_ratio_residual,
        },
        "neutral_transport": {
            "kind": "none_for_running_family",
            "status": "not_applied",
        },
        "charged_transport": {
            "kind": "none_for_running_family",
            "status": "not_applied",
        },
        "stage3_surrogate_transport": {
            "kind": "delta_rho_stage3_z_only",
            "delta_rho_stage3": delta_rho_stage3,
            "status": "outside_common_family",
        },
        "coherence_witness": {
            "running_mass_ratio_residual": running_mass_ratio_residual,
            "stage3_mass_ratio_residual": stage3_mass_ratio_residual,
            "mixed_sources_detected": True,
        },
        "reported_outputs": running_outputs,
        "mixed_reporting_surface": {
            "w_row_current_source": "d10_running_tree",
            "z_row_current_source": "stage3_pole_surrogate",
            "alpha_em_row_current_source": "d10_running_tree",
            "sin2w_row_current_source": "d10_running_tree",
            "mixed_scheme": True,
            "family_purity_violation": True,
            "stage3_surrogate_kind": "delta_rho_stage3_z_only",
            "exact_missing_transport_object_for_pole_family": "EWTransportKernel_D10",
        },
        "promotion_gate": {
            "mixed_scheme": False,
            "common_origin_certified": True,
            "descendant_targets_only": list(running_outputs.keys()),
            "exact_missing_transport_object_for_pole_family": "EWTransportKernel_D10",
            "smaller_exact_missing_clause": "EWTransportReadoutCoherence_D10",
            "family_purity_readout_criterion": "all_running_or_all_common_pole_effective",
        },
        "notes": [
            "This artifact freezes the exact single-P running electroweak family that the current D10 code already realizes algebraically.",
            "It should be treated separately from any pole/effective reporting surface.",
            "The current W/Z mismatch is a mixed-family reporting defect rather than a single-P precision defect: W, alpha_em, and sin^2(theta_W) already live on the running family while Z is still reported from a Stage-3 Z-only surrogate.",
            "If a pole/effective family is required, the remaining exact missing object is one common EWTransportKernel_D10 rather than more digits of P.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current D10 electroweak observable-family artifact.")
    parser.add_argument("--p", type=float, default=P_DEFAULT)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(args.p)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
