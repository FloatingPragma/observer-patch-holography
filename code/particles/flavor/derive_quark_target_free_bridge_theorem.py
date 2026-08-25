#!/usr/bin/env python3
"""Emit the conditional D12 light-ratio bridge identity package.

Chain role: make the exact missing public bridge theorem explicit as a theorem
package with statement, proof skeleton, and corollaries, rather than leaving it
only as a named frontier object.

Mathematics: on the minimal light branch
    y_u = c_u * epsilon^6,  y_d = c_d * epsilon^6,  epsilon = 1/6,
the common sixth-order factor cancels in the light ratio:
    log(y_d / y_u) = log(c_d / c_u).
The D12 overlap-defect scalar is the normalized one-step defect on that
six-step branch, so
    Delta_ud_overlap = (1/6) * log(y_d / y_u) = (1/6) * log(c_d / c_u).
On the emitted D12 mass ray one also has
    Delta_ud_overlap = t1 / 5,
therefore
    t1 = (5/6) * log(c_d / c_u).

Output: exact conditional identities and induced corollaries.  The current
source corpus does not emit ``c_d/c_u`` or ``t1``; an optional target-anchored
sidecar may be attached only by an explicit command-line argument and never
changes promotion status.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
MASS_RAY_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"
OVERLAP_LAW_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
CURRENT_FAMILY_TARGET_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_current_family_target_anchored_value_package.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_target_free_bridge_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(
    mass_ray: dict[str, Any],
    overlap_law: dict[str, Any],
    spread_map: dict[str, Any],
    current_family_target: dict[str, Any] | None = None,
) -> dict[str, Any]:
    x2 = float(mass_ray["sample_same_family_point"]["x2"])
    sigma_u = float(spread_map["sigma_u_total_log_per_side"])
    sigma_d = float(spread_map["sigma_d_total_log_per_side"])
    target_values = None
    if current_family_target is not None:
        target_values = dict(current_family_target["computed_candidate_theorem_values"]["quark_d12_t1_value_law"])
    return {
        "artifact": "oph_quark_target_free_bridge_theorem",
        "generated_utc": _timestamp(),
        "scope": "public_bridge_theorem_package",
        "proof_status": "conditional_bridge_identity_value_source_open",
        "public_promotion_allowed": False,
        "bridge_closure_status": "conditional_algebra_closed_value_emission_open",
        "exact_missing_object": "source_derived_c_d_over_c_u_or_t1_value",
        "value_emission_status": "open_viable_scalar_emission_route",
        "theorem_ids": [
            "light_quark_overlap_defect_value_law",
            "quark_d12_t1_value_law",
        ],
        "principal_theorem_statement": (
            "Conditional on a source-derived positive coefficient pair c_u,c_d, assume the minimal light branch "
            "y_u = c_u * epsilon^6, y_d = c_d * epsilon^6, epsilon = 1/6, and the emitted same-family "
            "D12 mass ray D12_ud^mass = R_D12^ud. Then the target-free light overlap-defect scalar is "
            "Delta_ud_overlap = (1/6) * log(c_d / c_u). On the emitted D12 mass ray one has "
            "Delta_ud_overlap = t1 / 5 and ray_modulus = t1. Therefore "
            "t1 = (5/6) * log(c_d / c_u), equivalently "
            "Delta_ud_overlap = t1 / 5 = (1/6) * log(c_d / c_u)."
        ),
        "equivalent_wrappers": {
            "light_quark_overlap_defect_value_law": (
                "Delta_ud_overlap = (1/6) * log(c_d / c_u)"
            ),
            "quark_d12_t1_value_law": (
                "t1 = ray_modulus = (5/6) * log(c_d / c_u)"
            ),
        },
        "proof_skeleton": [
            {
                "step": 1,
                "statement": (
                    "On the minimal light branch, y_u = c_u * epsilon^6 and y_d = c_d * epsilon^6 "
                    "with epsilon = 1/6. The common sixth-order factor cancels, so "
                    "log(y_d / y_u) = log(c_d / c_u)."
                ),
            },
            {
                "step": 2,
                "statement": (
                    "The light-quark overlap-defect scalar is the normalized one-step defect on the "
                    "six-step light branch. Therefore Delta_ud_overlap = (1/6) * log(y_d / y_u) = "
                    "(1/6) * log(c_d / c_u)."
                ),
            },
            {
                "step": 3,
                "statement": (
                    "Use the emitted D12 mass ray identity Delta_ud_overlap = t1 / 5. Substituting step 2 gives "
                    "t1 = 5 * Delta_ud_overlap = (5/6) * log(c_d / c_u)."
                ),
            },
        ],
        "induced_d12_mass_side_corollaries": {
            "Delta_ud_overlap": "(1/6) * log(c_d / c_u)",
            "t1": "(5/6) * log(c_d / c_u)",
            "ray_modulus": "t1",
            "eta_Q_centered": "-((1 - x2^2) / 27) * t1",
            "kappa_Q": "-t1 / 54",
            "x2": x2,
        },
        "induced_source_corollaries": {
            "beta_u_diag_B_source": "t1 / 10",
            "beta_d_diag_B_source": "-t1 / 10",
            "source_readback_u_log_per_side": "(-t1/10, 0, +t1/10)",
            "source_readback_d_log_per_side": "(+t1/10, 0, -t1/10)",
        },
        "induced_transport_corollaries": {
            "sigma_u_total_log_per_side": sigma_u,
            "sigma_d_total_log_per_side": sigma_d,
            "tau_u_log_per_side": "sigma_d_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            "tau_d_log_per_side": "sigma_u_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
            "Lambda_ud_B_transport": "sigma_u_total_log_per_side * sigma_d_total_log_per_side * t1 / (10 * (sigma_u_total_log_per_side + sigma_d_total_log_per_side))",
        },
        "single_bridge_gap_to_internalize": "source_derived_c_d_over_c_u_or_t1_value",
        "computed_current_family_target_check": target_values,
        "target_check_attachment": {
            "attached": target_values is not None,
            "classification": "explicit_comparison_only_target_anchored_sidecar",
            "used_as_proof_premise": False,
            "changes_public_promotion": False,
        },
        "notes": [
            "This artifact internalizes only the conditional bridge identities; it does not emit their scalar input.",
            "The D12 mass ray remains a viable open source-emission route, not a no-go.",
            "A current-family target check is absent by default and, when explicitly attached, remains comparison-only.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the synthesized quark target-free bridge theorem package.")
    parser.add_argument("--mass-ray", default=str(MASS_RAY_JSON))
    parser.add_argument("--overlap-law", default=str(OVERLAP_LAW_JSON))
    parser.add_argument("--spread-map", default=str(SPREAD_MAP_JSON))
    parser.add_argument(
        "--current-family-target",
        default=None,
        help=(
            "Optional explicit comparison-only target sidecar. It is never "
            "loaded by default and cannot change promotion status."
        ),
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    current_family_target_path = (
        Path(args.current_family_target) if args.current_family_target else None
    )
    payload = build_payload(
        _load_json(Path(args.mass_ray)),
        _load_json(Path(args.overlap_law)),
        _load_json(Path(args.spread_map)),
        _load_json(current_family_target_path)
        if current_family_target_path is not None and current_family_target_path.exists()
        else None,
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
