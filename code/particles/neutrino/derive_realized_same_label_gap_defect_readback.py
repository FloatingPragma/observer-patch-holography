#!/usr/bin/env python3
"""Export the realized same-label gap/defect readback artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "neutrino" / "defect_weighted_mu_e_family.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the realized same-label gap/defect readback artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    arrows = payload.get("realized_same_label_arrows", ["psi12", "psi23", "psi31"])
    null_map = {arrow: None for arrow in arrows}
    strict_missing = "realized_arrow_pullback_from_flavor_gap_and_defect_certificates"
    artifact = {
        "artifact": "oph_realized_same_label_gap_defect_readback",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only",
        "payload_status": "reduced_to_strict_smallest_missing_object",
        "upstream_exact_clause": payload.get("upstream_exact_clause"),
        "same_label_readback_origin": payload.get("same_label_readback_origin"),
        "selector_center": payload.get("selector_center"),
        "kernel_choice": payload.get("kernel_choice"),
        "realized_same_label_arrows": arrows,
        "same_label_overlap_sq": dict(null_map),
        "gap_e": dict(null_map),
        "defect_e": dict(null_map),
        "same_label_defect_rule": payload.get("same_label_defect_rule"),
        "q_e_rule": payload.get("raw_edge_score_rule"),
        "eta_e_rule": payload.get("centered_log_rule"),
        "base_mu_nu": payload.get("base_mu_nu"),
        "mu_e_rule": payload.get("weight_rule"),
        "derived_q_e": dict(null_map),
        "q_mean": None,
        "q_spread": None,
        "q_e": dict(null_map),
        "eta_e": dict(null_map),
        "mu_e": dict(null_map),
        "neutrino_only_isotropy_obstruction": payload.get("neutrino_only_isotropy_obstruction"),
        "smallest_constructive_missing_object": strict_missing,
        "strict_smallest_exact_missing_object": strict_missing,
        "missing_fields_by_arrow": {arrow: ["gap_e", "defect_e"] for arrow in arrows},
        "complete_by_arrow": {arrow: False for arrow in arrows},
        "metadata": {
            "note": "This artifact is now the exact builder-facing shell for the realized-arrow pullback bundle. The downstream neutrino chain is fixed once builder-grade gap_e and defect_e values are populated on the realized same-label arrows.",
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
