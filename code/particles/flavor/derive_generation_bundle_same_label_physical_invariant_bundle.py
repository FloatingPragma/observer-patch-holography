#!/usr/bin/env python3
"""Record the gauge-fixed physical-invariant shell for the D12 quark route.

Chain role: separate the honest same-sheet CKM/CP invariant shell from the
remaining physical-branch and mass-side burdens on the D12 continuation branch.

Mathematics: read the honest forward same-label transport unitary and its
principal logarithm from the D12 mass/transport closure artifact, then expose
the standard CKM angles, CP phase, Jarlskog, and generator-side invariant
package.

OPH-derived inputs: the D12 mass branch together with the forward-emitted
same-label transport closure.

Output: a D12 continuation bundle of physical invariants with zero remaining
transport-generator residual on the current D12 sheet.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TRANSPORT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_same_label_physical_invariant_bundle.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(transport: dict[str, object]) -> dict[str, object]:
    standard = dict(transport["standard_ckm_parameters"])
    generator = dict(transport["same_label_transport_generator"])
    invariants = dict(generator["generator_invariants"])
    selector_value = transport.get("quark_relative_sheet_selector")
    return {
        "artifact": "oph_generation_bundle_same_label_physical_invariant_bundle",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "forward_same_label_transport_closed_on_current_D12_sheet",
        "full_matrix_artifact": transport["artifact"],
        "gauge_convention": (
            "standard rephasing makes V_ud, V_us, V_cs, V_cb, V_tb real-positive; "
            "a diagonal conjugation then makes K_12 and K_23 real-positive"
        ),
        "physical_invariants": {
            "theta_12": standard["theta_12"],
            "theta_23": standard["theta_23"],
            "theta_13": standard["theta_13"],
            "phi_cp": standard["delta_ckm"],
            "jarlskog": standard["jarlskog"],
        },
        "diagonal_phase_bookkeeping": {
            "chi_diagonal_imag": invariants["chi_diagonal_imag"],
        },
        "physical_invariant_formula": (
            "V_CKM^fwd = U_u^dagger U_d, "
            "K_CKM = Log_pr(V_CKM^fwd), "
            "K_gf = i diag(chi_1,chi_2,chi_3) + "
            "theta_12(E12-E21) + theta_23(E23-E32) + "
            "theta_13(e^{i phi_cp} E13 - e^{-i phi_cp} E31)"
        ),
        "generator_invariants": invariants,
        "closure_residual_fro_norm": transport["closure_residual"]["fro_norm"],
        "current_sheet_status": transport.get("current_sheet_status"),
        "physical_branch_status": transport.get("physical_branch_status"),
        "quark_relative_sheet_selector": selector_value,
        "next_single_residual_object": "quark_d12_t1_value_law" if selector_value else "quark_relative_sheet_selector",
        "notes": [
            "This bundle records the gauge-fixed CKM/CP shell emitted by the forward same-label transport unitary on the current D12 sheet.",
            "The transport-generator residual is closed to machine precision on that sheet, but the current D12 sheet is not the physical quark branch.",
            (
                "The solver-side selector is now fixed to sigma_ref, so the next exact object is the one-scalar D12 value law on the emitted mass ray."
                if selector_value
                else "The exact next object is one discrete relative-sheet selector; mass-side scale fixing remains separate after branch selection."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the gauge-fixed D12 quark physical-invariant bundle.")
    parser.add_argument("--transport", default=str(DEFAULT_TRANSPORT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    transport = json.loads(Path(args.transport).read_text(encoding="utf-8"))
    payload = build_artifact(transport)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
