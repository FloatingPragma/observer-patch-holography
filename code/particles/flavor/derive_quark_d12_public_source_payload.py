#!/usr/bin/env python3
"""Emit the selected-class D12 pure-B conditional payload audit.

Chain role: materialize the pure-B source readback pair from the public exact
quark theorem on the selected frame class.

Mathematics: the public exact quark theorem emits the selected-class light
ratio. The closed D12 value law gives

    Delta_ud_overlap = (1/6) log(m_d / m_u)
    t1 = 5 Delta_ud_overlap.

The pure-B source payload is then forced by the D12 source corollary:

    beta_u = t1 / 10, beta_d = -t1 / 10,
    source_u = beta_u * B_ord, source_d = beta_d * B_ord.

The current upstream artifact contains target-attached mixed-scheme GeV mass
textures, not source-emitted dimensionless Yukawas, and the D12 value package
does not emit ``t1``.  The numerical projection is retained as a comparison
audit but promotion is inherited fail-closed from both upstreams.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PUBLIC_YUKAWA = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
DEFAULT_T1_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_d12_t1_value_law.json"
DEFAULT_SOURCE_READBACK = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_public_source_payload.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(public_yukawa: dict[str, Any], t1_law: dict[str, Any], source_readback: dict[str, Any]) -> dict[str, Any]:
    exact_values = dict(public_yukawa["public_exact_outputs"]["exact_running_values_gev"])
    m_u = float(exact_values["u"])
    m_d = float(exact_values["d"])
    ell_ud = math.log(m_d / m_u)
    delta_ud_overlap = ell_ud / 6.0
    t1 = 5.0 * delta_ud_overlap
    beta_u = t1 / 10.0
    beta_d = -t1 / 10.0
    b_ord = [float(value) for value in source_readback["B_ord"]]
    source_u = [beta_u * value for value in b_ord]
    source_d = [beta_d * value for value in b_ord]
    b_norm_sq = float(sum(value * value for value in b_ord))

    public_yukawa_ready = (
        public_yukawa.get("public_promotion_allowed") is True
        and public_yukawa.get("proof_status")
        == "closed_source_only_public_exact_yukawa_end_to_end_theorem"
        and public_yukawa.get("physical_yukawa_construction_closed") is True
        and (public_yukawa.get("non_circularity_status") or {}).get(
            "promotion_allowed"
        )
        is True
    )
    t1_value_ready = (
        t1_law.get("public_promotion_allowed") is True
        and t1_law.get("exact_missing_object") in {None, ""}
    )
    promotion_blockers: list[str] = []
    if not public_yukawa_ready:
        promotion_blockers.append(
            "PUBLIC_COMMON_SCALE_DIMENSIONLESS_SOURCE_YUKAWA_THEOREM"
        )
    if not t1_value_ready:
        promotion_blockers.append("SOURCE_DERIVED_C_D_OVER_C_U_OR_T1_VALUE")
    promotion_allowed = not promotion_blockers

    return {
        "artifact": "oph_quark_d12_public_source_payload",
        "generated_utc": _timestamp(),
        "scope": "selected_class_conditional_payload_and_target_attached_audit_only",
        "proof_status": (
            "closed_public_selected_class_pure_B_source_payload"
            if promotion_allowed
            else "blocked_unpromotable_yukawa_input_and_open_d12_value_source"
        ),
        "public_promotion_allowed": promotion_allowed,
        "off_canonical_promotion_allowed": False,
        "value_classification": (
            "source_emitted_selected_class_payload"
            if promotion_allowed
            else "comparison_only_projection_from_unpromotable_inputs"
        ),
        "promotion_blockers": promotion_blockers,
        "input_status_audit": {
            "public_yukawa_proof_status": public_yukawa.get("proof_status"),
            "public_yukawa_promotion_allowed": public_yukawa.get(
                "public_promotion_allowed"
            ),
            "physical_yukawa_construction_closed": public_yukawa.get(
                "physical_yukawa_construction_closed"
            ),
            "d12_t1_proof_status": t1_law.get("proof_status"),
            "d12_t1_promotion_allowed": t1_law.get("public_promotion_allowed"),
            "d12_t1_exact_missing_object": t1_law.get("exact_missing_object"),
            "promotion_inherited_fail_closed": True,
        },
        "input_artifacts": {
            "public_exact_yukawa_theorem": public_yukawa.get("artifact"),
            "d12_t1_value_law": t1_law.get("artifact"),
            "source_readback_law": source_readback.get("artifact"),
        },
        "selected_public_physical_frame_class": public_yukawa["selected_public_physical_frame_class"],
        "theorem_statement": (
            "Conditional on a source-emitted common-scale dimensionless light ratio and a closed D12 value law, "
            "Delta_ud_overlap = (1/6) log(y_d/y_u), t1 = 5 Delta_ud_overlap, and the pure-B payload is "
            "beta_u = t1/10, beta_d = -t1/10. The current stored numbers come from an unpromotable "
            "target-attached mixed-scheme mass-texture audit and do not satisfy those premises."
        ),
        "comparison_only_light_ratio_input": {
            "m_u_gev": m_u,
            "m_d_gev": m_d,
            "ell_ud": ell_ud,
            "ell_ud_formula": "log(m_d / m_u)",
        },
        "d12_scalars": {
            "Delta_ud_overlap": delta_ud_overlap,
            "Delta_ud_overlap_formula": "(1/6) * log(m_d / m_u)",
            "t1": t1,
            "t1_formula": "(5/6) * log(m_d / m_u)",
        },
        "B_ord": b_ord,
        "B_ord_norm_sq": b_norm_sq,
        "beta_u_diag_B_source": beta_u,
        "beta_d_diag_B_source": beta_d,
        "J_B_source_u": beta_u,
        "J_B_source_d": beta_d,
        "source_readback_u_log_per_side": source_u,
        "source_readback_d_log_per_side": source_d,
        "J_B_source_u_formula": "t1 / 10",
        "J_B_source_d_formula": "-t1 / 10",
        "beta_u_diag_B_source_formula": "t1 / 10",
        "beta_d_diag_B_source_formula": "-t1 / 10",
        "source_readback_u_log_per_side_formula": "beta_u_diag_B_source * B_ord",
        "source_readback_d_log_per_side_formula": "beta_d_diag_B_source * B_ord",
        "pure_B_certificates": {
            "center_entry_u": source_u[1],
            "center_entry_d": source_d[1],
            "endpoint_sum_u": source_u[0] + source_u[2],
            "endpoint_sum_d": source_d[0] + source_d[2],
            "J_B_from_endpoint_u": (source_u[2] - source_u[0]) / 2.0,
            "J_B_from_endpoint_d": (source_d[2] - source_d[0]) / 2.0,
            "dot_u_over_norm": sum(u * b for u, b in zip(source_u, b_ord, strict=True)) / b_norm_sq,
            "dot_d_over_norm": sum(d * b for d, b in zip(source_d, b_ord, strict=True)) / b_norm_sq,
        },
        "off_canonical_boundary": {
            "arbitrary_P_transport_closed": False,
            "reason": (
                "This payload is evaluated on the selected public class emitted by OPH axioms + P. "
                "It does not classify all public quark frame classes or prove an arbitrary-P source-payload family."
            ),
        },
        "notes": [
            "The arithmetic payload is retained only as a comparison audit while either upstream promotion gate is open.",
            "A future source payload requires common-scale dimensionless Yukawas with scheme, threshold, and running-v(mu) provenance, plus a source-derived light ratio or t1.",
            "The D12 scalar-emission route remains open and viable; this blocked payload is not a no-go for that route.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the selected-public-class quark D12 source payload.")
    parser.add_argument("--public-yukawa", default=str(DEFAULT_PUBLIC_YUKAWA))
    parser.add_argument("--t1-law", default=str(DEFAULT_T1_LAW))
    parser.add_argument("--source-readback", default=str(DEFAULT_SOURCE_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        _load_json(Path(args.public_yukawa)),
        _load_json(Path(args.t1_law)),
        _load_json(Path(args.source_readback)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
