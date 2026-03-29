#!/usr/bin/env python3
"""Record the gauge-fixed physical-invariant shell for the D12 quark route.

Chain role: separate the physical CKM/CP invariant shell from the still-open
same-label left-transport value law.

Mathematics: gauge-fix the anti-Hermitian left-transport generator so the
off-diagonal physical invariants become `(theta_12, theta_23, theta_13, phi_cp)`.

OPH-derived inputs: the current same-label left-transport continuation analysis.

Output: a diagnostic-only bundle of physical invariants and the next open value
law for the matrix-valued transport generator.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_same_label_physical_invariant_bundle.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, object]:
    return {
        "artifact": "oph_generation_bundle_same_label_physical_invariant_bundle",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "diagnostic_invariants_written_from_target_ckm_not_OPH_derived",
        "full_matrix_artifact": "oph_generation_bundle_same_label_left_transport_candidate",
        "gauge_convention": "diagonal conjugation chosen so that K_12 and K_23 are real-positive",
        "physical_invariants": {
            "theta_12": 0.2256228751541542,
            "theta_23": 0.04380491981435288,
            "theta_13": 0.003473054182841341,
            "phi_cp": -2.6946364809227097,
        },
        "diagonal_phase_bookkeeping": {
            "chi_diagonal_imag": [
                -0.006095179275809967,
                -0.7300387679202932,
                -0.3260719358139316,
            ],
        },
        "physical_invariant_formula": "K_gf = i diag(chi_1,chi_2,chi_3) + theta_12(E12-E21) + theta_23(E23-E32) + theta_13(e^{i phi_cp} E13 - e^{-i phi_cp} E31)",
        "relative_branch_rigidity": 0.00023461092057035416,
        "next_single_residual_object": "same_label_left_transport_physical_invariant_value_laws",
        "notes": [
            "This bundle records the gauge-fixed physical CKM/CP shell on the current D12 continuation branch.",
            "The recorded invariant values come from the current compare-derived continuation analysis and are not promoted as OPH-forward laws.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the gauge-fixed D12 quark physical-invariant bundle.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    payload = build_artifact()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
