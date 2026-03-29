#!/usr/bin/env python3
"""Compress same-label pullback data to a proof-facing scalar certificate.

Chain role: reduce the proof-facing neutrino input from raw same-label matrix
data to the scalar family the intrinsic eta-chain actually consumes.

Mathematics: same-label overlap/gap/defect compression followed by the exact
`q_e -> eta_e -> mu_e` reduction.

OPH-derived inputs: a full overlap edge-sector bundle or an existing same-label
pullback payload from the live neutrino/flavor chain.

Output: either a complete scalar certificate or, if the pullback values are
still missing, an honest shell stating the exact remaining readback object.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
EDGE_ORDER = ("psi12", "psi23", "psi31")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _save_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_matrix(obj: Any) -> np.ndarray:
    if isinstance(obj, dict) and "real" in obj:
        real = np.asarray(obj["real"], dtype=float)
        imag = np.asarray(obj.get("imag", np.zeros_like(real)), dtype=float)
        return real + 1j * imag
    return np.asarray(obj, dtype=float)


def _frob_norm_sq(mat: np.ndarray) -> float:
    return float(np.real(np.trace(mat.conjugate().T @ mat)))


def _is_complete_scalar_map(value: Any) -> bool:
    return isinstance(value, dict) and all(value.get(edge) is not None for edge in EDGE_ORDER)


def _compute_edge_from_full_bundle(edge_name: str, edge: dict[str, Any]) -> dict[str, Any]:
    label = str(edge["same_label"])
    tail = edge["tail"]
    head = edge["head"]
    p_tail = _parse_matrix(tail["same_label_projector"])
    p_head = _parse_matrix(head["same_label_projector"])
    lift = _parse_matrix(edge["transport"])

    rank_tail = float(np.real(np.trace(p_tail)))
    rank_head = float(np.real(np.trace(p_head)))
    if abs(rank_tail - rank_head) > 1.0e-8 or rank_tail <= 0.0:
        raise ValueError(f"{edge_name}: same-label projectors must have equal positive rank")

    lp_t_l = lift @ p_tail @ lift.conjugate().T
    overlap_raw = float(np.real(np.trace(p_head @ lp_t_l)))
    overlap_sq = overlap_raw / rank_tail
    defect_overlap = float(1.0 - overlap_sq)
    defect_projector = float(0.5 * _frob_norm_sq(p_head - lp_t_l) / rank_tail)
    projector_residual = float(abs(defect_overlap - defect_projector))

    eig_tail = {str(k): float(v) for k, v in dict(tail["sector_eigenvalues"]).items()}
    eig_head = {str(k): float(v) for k, v in dict(head["sector_eigenvalues"]).items()}
    gap_tail = min(abs(v - eig_tail[label]) for k, v in eig_tail.items() if k != label)
    gap_head = min(abs(v - eig_head[label]) for k, v in eig_head.items() if k != label)
    gap_e = math.sqrt(gap_tail * gap_head)

    return {
        "same_label": label,
        "rank": rank_tail,
        "same_label_overlap_sq": float(overlap_sq),
        "defect_e": defect_overlap,
        "gap_e": gap_e,
        "same_label_gap_witness": {
            "tail_gap": float(gap_tail),
            "head_gap": float(gap_head),
            "gap_rule": "g_e = sqrt(gap_tail * gap_head)",
        },
        "same_label_defect_witness": {
            "projector_trace_rank": rank_tail,
            "defect_from_projector_identity": defect_projector,
            "projector_identity_residual": projector_residual,
            "defect_rule": "d_e = (1/(2 r_e)) ||P_head - L_e P_tail L_e^*||_HS^2 = 1 - overlap_sq_e",
        },
    }


def _build_shell(source: dict[str, Any], source_path: Path) -> dict[str, Any]:
    smallest = source.get(
        "smallest_constructive_missing_object",
        "realized_arrow_pullback_from_flavor_gap_and_defect_certificates",
    )
    return {
        "artifact": "oph_neutrino_same_label_scalar_certificate",
        "generated_utc": _timestamp(),
        "source_path": str(source_path),
        "source_artifacts": [source.get("artifact")],
        "source_kind": "same_label_pullback_shell",
        "proof_status": "shell_waiting_live_same_label_scalars",
        "compression_status": "proof_facing_compression_defined_values_open",
        "sufficient_for_intrinsic_mass_eigenstates": False,
        "insufficient_for_pmns_or_public_flavor_rows": True,
        "exact_downstream_factorization_object": "same_label_scalar_certificate_(gap_e,overlap_sq_e)_mod_common_scale",
        "builder_facing_exact_object": "centered_log_pullback_class_[log_q_e]",
        "realized_same_label_arrows": list(source.get("realized_same_label_arrows", EDGE_ORDER)),
        "same_label": dict(source.get("same_label", {edge: None for edge in EDGE_ORDER})),
        "rank_e": {edge: None for edge in EDGE_ORDER},
        "same_label_overlap_sq": dict(source.get("same_label_overlap_sq", {edge: None for edge in EDGE_ORDER})),
        "defect_e": dict(source.get("defect_e", {edge: None for edge in EDGE_ORDER})),
        "gap_e": dict(source.get("gap_e", {edge: None for edge in EDGE_ORDER})),
        "q_e": {edge: None for edge in EDGE_ORDER},
        "eta_e": {edge: None for edge in EDGE_ORDER},
        "base_mu_nu": source.get("base_mu_nu"),
        "mu_e": {edge: None for edge in EDGE_ORDER},
        "certificate_residuals": {
            "max_projector_identity_residual": None,
        },
        "smallest_constructive_missing_object": smallest,
        "notes": [
            "The scalar certificate contract is fixed, but the live same-label overlap/gap/defect scalars are not yet populated on all realized arrows.",
            "Once those scalars are emitted, the downstream intrinsic neutrino mass-eigenstate branch factors through this certificate without needing the raw matrix payloads.",
        ],
    }


def _build_complete_certificate(source: dict[str, Any], source_path: Path, *, base_mu: float | None) -> dict[str, Any]:
    artifact_name = str(source.get("artifact", ""))
    edge_info: dict[str, dict[str, Any]] = {}
    if artifact_name == "oph_overlap_edge_sector_data_bundle":
        raw_edges = dict(source["edges"])
        for edge in EDGE_ORDER:
            edge_info[edge] = _compute_edge_from_full_bundle(edge, raw_edges[edge])
        source_kind = "full_overlap_edge_sector_bundle"
    else:
        for edge in EDGE_ORDER:
            edge_info[edge] = {
                "same_label": (source.get("same_label") or {}).get(edge),
                "rank": float(((source.get("rank_e") or {}).get(edge)) or 1.0),
                "same_label_overlap_sq": float(source["same_label_overlap_sq"][edge]),
                "defect_e": float(source["defect_e"][edge]),
                "gap_e": float(source["gap_e"][edge]),
                "same_label_gap_witness": ((source.get("same_label_gap_witness") or {}).get(edge)),
                "same_label_defect_witness": ((source.get("same_label_defect_witness") or {}).get(edge)),
            }
        source_kind = "existing_same_label_pullback_payload"

    overlap_sq = {edge: float(edge_info[edge]["same_label_overlap_sq"]) for edge in EDGE_ORDER}
    defect_e = {edge: float(edge_info[edge]["defect_e"]) for edge in EDGE_ORDER}
    gap_e = {edge: float(edge_info[edge]["gap_e"]) for edge in EDGE_ORDER}
    rank_e = {edge: float(edge_info[edge]["rank"]) for edge in EDGE_ORDER}
    same_label = {edge: edge_info[edge]["same_label"] for edge in EDGE_ORDER}
    q_e = {edge: float(math.sqrt(gap_e[edge] * defect_e[edge])) for edge in EDGE_ORDER}
    log_q = np.asarray([math.log(q_e[edge]) for edge in EDGE_ORDER], dtype=float)
    mean_log_q = float(np.mean(log_q))
    eta_e = {edge: float(log_q[idx] - mean_log_q) for idx, edge in enumerate(EDGE_ORDER)}
    base_mu_value = float(base_mu if base_mu is not None else source.get("base_mu_nu", 1.0))
    exp_eta = {edge: math.exp(eta_e[edge]) for edge in EDGE_ORDER}
    mean_exp_eta = float(sum(exp_eta.values()) / len(exp_eta))
    mu_e = {edge: float(base_mu_value * exp_eta[edge] / mean_exp_eta) for edge in EDGE_ORDER}

    max_projector_identity_residual = 0.0
    for edge in EDGE_ORDER:
        witness = edge_info[edge].get("same_label_defect_witness")
        if isinstance(witness, dict) and witness.get("projector_identity_residual") is not None:
            max_projector_identity_residual = max(max_projector_identity_residual, float(witness["projector_identity_residual"]))

    return {
        "artifact": "oph_neutrino_same_label_scalar_certificate",
        "generated_utc": _timestamp(),
        "source_path": str(source_path),
        "source_artifacts": [source.get("artifact")],
        "source_kind": source_kind,
        "proof_status": "fixed_cutoff_scalar_sufficient_downstream_certificate",
        "compression_status": "matrix_data_compressed_to_scalar_certificate",
        "sufficient_for_intrinsic_mass_eigenstates": True,
        "insufficient_for_pmns_or_public_flavor_rows": True,
        "exact_downstream_factorization_object": "same_label_scalar_certificate_(gap_e,overlap_sq_e)_mod_common_scale",
        "builder_facing_exact_object": "centered_log_pullback_class_[log_q_e]",
        "realized_same_label_arrows": list(EDGE_ORDER),
        "same_label": same_label,
        "rank_e": rank_e,
        "same_label_overlap_sq": overlap_sq,
        "defect_e": defect_e,
        "gap_e": gap_e,
        "q_e": q_e,
        "eta_e": eta_e,
        "base_mu_nu": base_mu_value,
        "mu_e": mu_e,
        "same_label_gap_witness": {edge: edge_info[edge].get("same_label_gap_witness") for edge in EDGE_ORDER},
        "same_label_defect_witness": {edge: edge_info[edge].get("same_label_defect_witness") for edge in EDGE_ORDER},
        "certificate_residuals": {
            "max_projector_identity_residual": float(max_projector_identity_residual),
        },
        "rules": {
            "overlap_scalar": "omega_e = same_label_overlap_sq_e = (1/r_e) Tr(P_head L_e P_tail L_e^*)",
            "defect_rule": "d_e = 1 - omega_e",
            "gap_rule": "g_e = sqrt(gap_tail * gap_head)",
            "q_rule": "q_e = sqrt(g_e * d_e)",
            "eta_rule": "eta_e = log(q_e) - mean_f(log q_f)",
            "mu_rule": "mu_e = base_mu_nu * exp(eta_e) / mean_f(exp(eta_f))",
        },
        "notes": [
            "This is the smallest proof-facing scalar artifact the downstream intrinsic neutrino lane needs.",
            "The intrinsic neutrino mass-eigenstate branch factors through this scalar certificate once the same-label scalars are populated.",
            "PMNS and flavor-labeled rows still require the shared charged-lepton left basis.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the neutrino same-label scalar certificate.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--base-mu", type=float, default=None)
    args = parser.parse_args()

    source_path = Path(args.input)
    source = _load_json(source_path)

    if _is_complete_scalar_map(source.get("same_label_overlap_sq")) and _is_complete_scalar_map(source.get("gap_e")) and _is_complete_scalar_map(source.get("defect_e")):
        payload = _build_complete_certificate(source, source_path, base_mu=args.base_mu)
    else:
        payload = _build_shell(source, source_path)

    out_path = Path(args.output)
    _save_json(out_path, payload)
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
