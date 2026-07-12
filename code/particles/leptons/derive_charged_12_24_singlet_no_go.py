#!/usr/bin/env python3
"""Prove that the invariant 12/24 hierarchy lock cannot split three families."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
HIERARCHY = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_local_global_hierarchy_resonance_closeout_335.json"
)
SCALARS = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_sector_local_support_extension_source_scalar_pair_readback.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_12_24_singlet_no_go.json"
)


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def centered(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    return [value - mean for value in values]


def build_artifact(hierarchy: dict[str, Any], scalars: dict[str, Any]) -> dict[str, Any]:
    sieve = hierarchy["screen_sieve_certificate"]
    register = hierarchy["round_count_certificate"]
    ports = int(sieve["orbit_size"])
    slots = int(register["m_rep"])
    if ports != 12 or slots != 24:
        raise ValueError("the declared hierarchy branch must carry the 12/24 lock")

    # A transitive invariant load has one value on every port. Any family-blind
    # equivariant readout therefore gives the same value on all three families.
    test_family_image = [1.0, 1.0, 1.0]
    centered_image = centered(test_family_image)
    return {
        "artifact": "oph_charged_12_24_singlet_no_go",
        "status": "CLOSED_NO_GO_INVARIANT_12_24_LOAD_INSUFFICIENT",
        "public_charged_mass_promotion_allowed": False,
        "hierarchy_inputs": {
            "screen_ports": ports,
            "repair_register_slots": slots,
            "screen_action": "transitive_A5_over_C5_orbit",
            "load_rule": "invariant_scalar_X_maps_to_X_over_12_per_port",
            "register_role": "common_hierarchy_exponent_normalization",
        },
        "theorem": {
            "id": "invariant_12_24_singlet_cannot_emit_charged_family_shape",
            "statement": (
                "Let the twelve screen ports carry the transitive A5/C5 action used by "
                "the hierarchy sieve, and let the 24-slot write/check register act only "
                "through its invariant common repair count. Any family-blind equivariant "
                "readout of that invariant load into a transitive three-family carrier has "
                "image proportional to (1,1,1). Its centered projection is zero. Hence the "
                "invariant 12/24 hierarchy lock can fix a common electroweak scale but cannot "
                "emit either independent charged-family shape scalar."
            ),
            "proof": [
                "Transitivity makes the invariant subspace of the twelve-port permutation representation one-dimensional, spanned by the all-ones vector.",
                "The declared screen sieve retains exactly this invariant component through X -> X/12.",
                "The oriented 24-slot count changes the common exponent normalization but supplies no family-labelled non-singlet record.",
                "Equivariance and absence of a marked family force the three-family image to be proportional to (1,1,1).",
                "Centering removes that image, so both independent coordinates of the three-family sum-zero plane vanish.",
            ],
            "centered_equal_family_image": centered_image,
            "forced_shape_on_invariant_branch": {
                "eta_source_support_extension_log_per_side": 0.0,
                "sigma_source_support_extension_total_log_per_side": 0.0,
                "mass_ratio_class": "1:1:1",
            },
        },
        "comparison_with_live_gap": {
            "eta_live": scalars.get("eta_source_support_extension_log_per_side"),
            "sigma_live": scalars.get("sigma_source_support_extension_total_log_per_side"),
            "live_values_still_open": True,
            "conclusion": "the invariant 12/24 certificate cannot populate the two open scalar slots",
        },
        "minimal_new_theorem": {
            "id": "charged_family_non_singlet_port_attachment",
            "required_objects": [
                "a refinement-stable non-singlet moment of the twelve-port or oriented-repair record",
                "an ID-independent equivariant attachment from that moment to three physical charged-family lines",
                "two independent centered readouts spanning the charged-family sum-zero plane",
                "a no-target-leak receipt proving the attachment was fixed before charged-mass comparison",
            ],
            "representation_options_not_yet_selected": ["A5 irrep 3", "A5 irrep 3-prime", "other derived three-family carrier"],
            "warning": "choosing an irrep or attachment because it reproduces lepton masses would be target fitting",
        },
        "claim_boundary": (
            "The 12/24 resonance can participate in a future charged theorem only after a "
            "derived non-singlet family attachment is supplied. The existing invariant "
            "hierarchy theorem alone predicts degeneracy, not the charged-lepton hierarchy."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hierarchy", type=Path, default=HIERARCHY)
    parser.add_argument("--scalars", type=Path, default=SCALARS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact(_load(args.hierarchy), _load(args.scalars))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(artifact["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
