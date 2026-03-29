#!/usr/bin/env python3
"""Emit the honest D12 light-quark isospin selector law.

Chain role: make the paper-facing D12 continuation selector explicit without
pretending that the selector value itself is already emitted from recovered
core.

Mathematics: odd-budget neutrality on the light-quark pure-`B` sector, written
as one D12 selector scalar whose value would fix the light-quark payload pair.

OPH-derived inputs: the already-closed pure-`B` source-readback law on the
ordered three-point family.

Output: a machine-readable D12 selector-law artifact beneath the larger quark
mass/CKM/CP continuation program.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_READBACK = (
    ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor" / "light_quark_isospin_overlap_defect_selector_law.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_readback: dict) -> dict:
    b_ord = list(source_readback.get("B_ord") or [-1.0, 0.0, 1.0])
    return {
        "artifact": "oph_light_quark_isospin_overlap_defect_selector_law",
        "generated_utc": _timestamp(),
        "status": "closed_smaller_primitive",
        "proof_status": "selector_law_closed_value_open",
        "scope": "D12_continuation_only",
        "recovered_core_promotion_allowed": False,
        "recovered_core_no_go_for_nonzero_light_quark_pure_b_selector": True,
        "selector_scalar_name": "Delta_ud_overlap",
        "next_single_residual_object": "Delta_ud_overlap",
        "odd_budget_neutrality_formula": "beta_u_diag_B_source + beta_d_diag_B_source = 0",
        "selector_equivalence_formula": "Delta_ud_overlap = beta_u_diag_B_source - beta_d_diag_B_source",
        "beta_u_diag_B_source_formula": "Delta_ud_overlap / 2",
        "beta_d_diag_B_source_formula": "-Delta_ud_overlap / 2",
        "B_ord": b_ord,
        "source_readback_u_log_per_side_formula": "(Delta_ud_overlap / 2) * B_ord",
        "source_readback_d_log_per_side_formula": "(-Delta_ud_overlap / 2) * B_ord",
        "tau_u_log_per_side_formula": "beta_u_diag_B_source",
        "tau_d_log_per_side_formula": "beta_d_diag_B_source",
        "Delta_ud_overlap": None,
        "beta_u_diag_B_source": None,
        "beta_d_diag_B_source": None,
        "source_readback_u_log_per_side": None,
        "source_readback_d_log_per_side": None,
        "tau_u_log_per_side": None,
        "tau_d_log_per_side": None,
        "notes": [
            "This is the honest D12 continuation-level selector law for light-quark isospin splitting.",
            "It does not override the recovered-core no-go: the selector value itself is still open and therefore not promotable as a recovered-core nonzero pure-B source selector.",
            "Once Delta_ud_overlap is emitted on the honest D12 route, the light-sector pure-B payload pair and tau pair follow algebraically on the ordered three-point family.",
        ],
        "source_artifact": source_readback.get("artifact"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 light-quark isospin overlap-defect selector-law artifact.")
    parser.add_argument("--source-readback", default=str(DEFAULT_SOURCE_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_readback = json.loads(Path(args.source_readback).read_text(encoding="utf-8"))
    artifact = build_artifact(source_readback)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
