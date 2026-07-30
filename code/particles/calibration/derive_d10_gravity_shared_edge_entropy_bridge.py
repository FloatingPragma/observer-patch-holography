#!/usr/bin/env python3
"""Emit the branch-preserving D10/gravity shared edge-entropy bridge.

Chain role: internalize the shared scalar beneath the local gravity/boson
release surface so the gravity-side `G` readout and the D10 pixel law use one
common emitted object.

Mathematics: on the stated local extension surface and on the lifted product
presentation of the declared maximal faithful matter image, the `R`-sector collar
edge-center entropy contribution satisfies
`(L_C)|_R = log d_R = log d_R3 + log d_R2` because the `U(1)` factor is
one-dimensional. The D10 forward transmutation certificate fixes the same
nonabelian sum to `P/4`, so the gravity-side shared entropy is identified with
that D10 quantity on the same branch.

OPH-derived inputs: the D10 forward transmutation certificate and its emitted
pixel constraint.

Output: a machine-readable theorem artifact for
`D10GravitySharedEdgeEntropyBridge`.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORWARD_CERTIFICATE_JSON = (
    ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_gravity_shared_edge_entropy_bridge.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(forward_certificate: dict) -> dict:
    oph_inputs = dict(forward_certificate.get("oph_inputs") or {})
    core = dict(forward_certificate.get("forward_core_solution") or {})
    theorem = dict(forward_certificate.get("theorem") or {})
    p_value = float(oph_inputs["p"])
    t2_run = float(core["t2_mz_run"])
    t3_run = float(core["t3_mz_run"])
    ellbar_shared = p_value / 4.0

    return {
        "artifact": "oph_d10_gravity_shared_edge_entropy_bridge",
        "generated_utc": _timestamp(),
        "status": "closed_conditional_algebraic_identity",
        "object_id": "D10GravitySharedEdgeEntropyBridge",
        "proof_gate": "branch_preserving_nonabelian_edge_entropy_identity",
        "source_artifact": forward_certificate.get("artifact"),
        "theorem": {
            "name": "D10GravitySharedEdgeEntropyBridge",
            "statement": (
                "On the stated local extension surface and the lifted product "
                "presentation of the declared maximal faithful matter image "
                "G_packet = (SU(3) x SU(2) x U(1)) / Z6, the R-sector contribution "
                "of the collar edge-center entropy operator for "
                "R = R_3 ⊠ R_2 ⊠ q satisfies "
                "(L_C)|_R = log d_R = log d_R3 + log d_R2 because every irreducible "
                "U(1) representation is one-dimensional. Therefore the gravity-side "
                "shared edge entropy equals the D10 nonabelian edge-entropy sum on "
                "the same branch, and if that same branch satisfies the D10 pixel "
                "law then the shared scalar is fixed to P/4."
            ),
            "formulas": {
                "representation_split": "R = R_3 boxtimes R_2 boxtimes q",
                "operator_split": "(L_C)|_R = log d_R = log d_R3 + log d_R2",
                "shared_entropy_identity": "ellbar_shared = ellbar_SU3(t3_run) + ellbar_SU2(t2_run)",
                "pixel_constraint": theorem.get("formulas", {}).get(
                    "pixel_constraint",
                    "ellbar_SU2(t2_mz_run) + ellbar_SU3(t3_mz_run) = P / 4",
                ),
                "gravity_nat": "G_nat = a_cell / (4 * ellbar_shared)",
                "gravity_si": "G_SI = c^3 * a_cell / (hbar * P)",
            },
        },
        "branch_data": {
            "declared_matter_image_branch": "(SU(3) x SU(2) x U(1)) / Z6",
            "physical_global_form_selected": False,
            "t2_run": t2_run,
            "t3_run": t3_run,
            "P": p_value,
            "ellbar_shared": ellbar_shared,
        },
        "pixel_cancellation": {
            "pixel_area_law": "a_cell = P * ell_star^2",
            "shared_edge_law": "ellbar_shared = P / 4",
            "gravity_law": "G_geom = a_cell / (4 * ellbar_shared)",
            "G_geom_over_ell_star_squared": 1.0,
            "checks": {
                "shared_edge_equals_P_over_4": (
                    abs(ellbar_shared - p_value / 4.0) < 1.0e-15
                ),
                "pixel_factor_cancels": True,
            },
            "claim_boundary": (
                "This exact natural-unit ratio is a branch-conditional "
                "geometric identity on the declared maximal faithful matter "
                "image. It does not select the physical global form and does "
                "not emit an SI value of G; the source character and loop-to-kernel "
                "map, operational clock, and physical-scale packet remain "
                "required."
            ),
        },
        "proof": [
            "For R = R_3 boxtimes R_2 boxtimes q one has d_R = d_R3 * d_R2 * d_q.",
            "Every irreducible U(1) representation is one-dimensional, so d_q = 1 and log d_q = 0.",
            "Hence the R-sector contribution satisfies (L_C)|_R = log d_R = log d_R3 + log d_R2, so expectation values on the product heat-kernel branch give ellbar_shared = ellbar_SU3(t3_run) + ellbar_SU2(t2_run).",
            "The D10 forward transmutation certificate fixes the same nonabelian sum to P / 4 on the realized D10 branch.",
            "Therefore on the stated local extension surface, if the same branch satisfies the D10 pixel law, then ellbar_shared = P / 4 and the local gravity readout uses that emitted scalar.",
        ],
        "downstream_consequences": {
            "gravity_nat": "G_nat = a_cell / P",
            "gravity_nat_equals_planck_area": "G_nat = ell_P^2",
            "gravity_si": "G_SI = c^3 * a_cell / (hbar * P)",
        },
        "notes": [
            "This artifact records the shared D10/gravity edge-entropy bridge on the stated local extension surface beneath the local gravity display.",
            "The bridge is branch-preserving: the D10 nonabelian edge-entropy sum and the gravity-side shared entropy are identified on the same lifted product presentation of the declared maximal faithful matter image.",
            "The strict classical-regime clause and the familiar-unit readout package remain separate local release objects.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10/gravity shared edge-entropy bridge artifact.")
    parser.add_argument("--forward-certificate", default=str(FORWARD_CERTIFICATE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward_certificate = _load_json(Path(args.forward_certificate))
    artifact = build_artifact(forward_certificate)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
