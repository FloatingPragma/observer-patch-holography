#!/usr/bin/env python3
"""Declare the honest D10 repair frontier beyond the present current carrier.

Chain role: separate the builder-local current-carrier selector shell from the
broader exact-PDG electroweak repair program.

Mathematics: frontier bookkeeping on top of the already-closed current-carrier
chart, keeping track of which residual is merely builder-local and which one is
the honest global repair burden.

OPH-derived inputs: the selected D10 carrier readout, the exact current-carrier
mass chart, and the current-carrier obstruction artifacts already emitted on the
live calibration path.

Output: a machine-readable repair-branch artifact stating that exact PDG W/Z
requires a branch beyond the present current carrier.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
DEFAULT_EXACT_CHART = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_mass_pair_chart_current_carrier.json"
DEFAULT_TAU2_OBSTRUCTION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_branch_beyond_current_carrier.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(readout: dict, exact_chart: dict, tau2_obstruction: dict | None = None) -> dict:
    public_quintet = dict(readout.get("public_emitted_quintet") or {})
    return {
        "artifact": "oph_d10_ew_repair_branch_beyond_current_carrier",
        "generated_utc": _timestamp(),
        "object_id": "D10RepairBranchBeyondCurrentCarrier",
        "status": "open",
        "proof_status": "current_carrier_chart_closed_exact_pdg_wz_still_requires_family_escape",
        "required_closure_kind": "single_family_single_P_no_mixed_readout",
        "builder_local_frontier": exact_chart.get("next_single_residual_object") or "EWExactMassPairSelector_D10",
        "replaces_builder_local_frontier": "EWExactMassPairSelector_D10",
        "stronger_residual_object": "EWSinglePostTransportTreeIdentity_D10",
        "current_carrier_exact_mass_pair": {
            "MW_pole": float(public_quintet.get("MW_pole", 0.0)),
            "MZ_pole": float(public_quintet.get("MZ_pole", 0.0)),
        },
        "current_carrier_chart_artifact": exact_chart.get("artifact"),
        "current_carrier_chart_status": exact_chart.get("status"),
        "current_carrier_tau2_obstruction_artifact": None if tau2_obstruction is None else tau2_obstruction.get("artifact"),
        "current_carrier_tau2_obstruction_status": None if tau2_obstruction is None else tau2_obstruction.get("status"),
        "operative_primitive": {
            "object_id": "EWSinglePostTransportTreeIdentity_D10",
            "status": "closed_smaller_primitive",
            "proof_status": "coherent_tree_law_closed_value_open",
            "family_source_id": "d10_running_tree",
            "scheme_id": "freeze_once",
            "origin_kernel_id": "EWTransportKernel_D10",
            "base_running_quintet": {
                "alphaY_prime": "alphaY_*",
                "alpha2_prime": "alpha2_*",
                "alpha_em_eff_inv": "(alphaY_* + alpha2_*) / (alphaY_* * alpha2_*)",
                "sin2w_eff": "alphaY_* / (alphaY_* + alpha2_*)",
                "MW_pole": "v_report * sqrt(pi * alpha2_*)",
                "MZ_pole": "v_report * sqrt(pi * (alphaY_* + alpha2_*))",
            },
            "coupling_update": {
                "alphaY_prime": "alphaY_* + delta_alphaY_tree",
                "alpha2_prime": "alpha2_* + delta_alpha2_tree",
            },
            "coherent_emitter_formula": {
                "alpha_em_eff_inv": "(alphaY_prime + alpha2_prime) / (alphaY_prime * alpha2_prime)",
                "sin2w_eff": "alphaY_prime / (alphaY_prime + alpha2_prime)",
                "MW_pole": "v_report * sqrt(pi * alpha2_prime)",
                "MZ_pole": "v_report * sqrt(pi * (alphaY_prime + alpha2_prime))",
            },
            "tree_identity_residuals": {
                "R_alpha": "(alphaY_prime + alpha2_prime) / (alphaY_prime * alpha2_prime) - alpha_em_eff_inv",
                "R_sin": "1 - (MW_pole / MZ_pole)^2 - sin2w_eff",
            },
            "proof_gate": "single_family_single_P_no_mixed_readout",
            "forbidden_inverse_witnesses": [
                "inverse_reference_WZ_slice",
                "mixed_running_alpha_plus_reference_mass_readout",
                "mixed_running_sin2_plus_reference_mass_readout",
            ],
        },
        "target_scope": [
            "MW_pole",
            "MZ_pole",
            "alpha_em_eff_inv",
            "sin2w_eff",
            "v_report",
        ],
        "notes": [
            "The present selected/current carrier closes its own exact W/Z chart.",
            "The builder-local next object on that chart is `EWExactMassPairSelector_D10`.",
            "The strongest strictly smaller constructive primitive beneath the broader repair branch is `EWSinglePostTransportTreeIdentity_D10`.",
            "Exact PDG W/Z should still be treated as a broader D10 repair problem beyond the present current carrier rather than as a hidden selector rename on the same carrier.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 repair-branch artifact beyond the current carrier.")
    parser.add_argument("--readout", default=str(DEFAULT_READOUT))
    parser.add_argument("--exact-chart", default=str(DEFAULT_EXACT_CHART))
    parser.add_argument("--tau2-obstruction", default=str(DEFAULT_TAU2_OBSTRUCTION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    exact_chart = json.loads(Path(args.exact_chart).read_text(encoding="utf-8"))
    tau2_path = Path(args.tau2_obstruction)
    tau2_obstruction = json.loads(tau2_path.read_text(encoding="utf-8")) if tau2_path.exists() else None
    artifact = build_artifact(readout, exact_chart, tau2_obstruction)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
