#!/usr/bin/env python3
"""Emit exact compare-only neutrino bridge coordinates on the exact segment adapter.

Chain role: make the exact two-parameter neutrino adapter more explicit by
recording the induced compare-only bridge coordinates on that same branch.

Mathematics: the exact segment adapter already fixes one selector coordinate
``tau_nu`` and one positive rescaling ``lambda_nu`` so that the representative
central splittings are hit exactly on the existing positive selector segment.
On that same compare-only branch one can evaluate

    B_nu = lambda_nu * q_mean^p_nu / m_star_eV
    P_nu = sqrt(I_nu) * sqrt(ratio_hat) / sum_defect
    C_nu = B_nu / P_nu

with the exact-adapter ratio ``ratio_hat`` and the same emitted phase/defect
stack. This remains a compare-only sidecar and must not be promoted into the
theorem lane.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ADAPTER_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
NORMALIZER_JSON = ROOT / "particles" / "runs" / "neutrino" / "same_label_overlap_defect_weight_normalizer.json"
READBACK_JSON = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"
HESSIAN_JSON = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_hessian.json"
REPAIR_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
SCALE_ANCHOR_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
CORRECTION_SCAFFOLD_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_bridge_correction_invariant_scaffold.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_exact_adapter_bridge_coordinate.json"

EDGE_ORDER = ("psi12", "psi23", "psi31")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _repo_ref(path: Path) -> str:
    return str(path.relative_to(ROOT.parent))


def build_payload(
    *,
    adapter: dict[str, Any],
    normalizer: dict[str, Any],
    readback: dict[str, Any],
    hessian: dict[str, Any],
    repair: dict[str, Any],
    scale_anchor: dict[str, Any],
    correction_scaffold: dict[str, Any],
) -> dict[str, Any]:
    lambda_nu = float(adapter["exact_solution"]["lambda_nu"])
    p_nu = float(adapter["exact_solution"]["p_nu"])
    tau_nu = float(adapter["exact_solution"]["tau_nu"])
    d_nu = float(adapter["exact_solution"]["D_nu"])
    q_mean = float(normalizer["q_mean"])
    ratio_hat = float(adapter["exact_outputs"]["ratio_21_over_32"])
    m_star_eV = float(scale_anchor["anchors"]["m_star_gev"]) * 1.0e9
    defect_sum = float(sum(float(value) for value in readback["defect_e"].values()))
    qbar = {edge: float(normalizer["qbar_e"][edge]) for edge in EDGE_ORDER}
    selector = {edge: float(repair["selector_phases_absolute"][edge]) for edge in EDGE_ORDER}
    center = {edge: float(hessian["selector_point"][edge]) for edge in EDGE_ORDER}
    i_nu = sum(qbar[edge] * (1.0 - math.cos(selector[edge] - center[edge])) for edge in EDGE_ORDER)

    b_nu = lambda_nu * (q_mean**p_nu) / m_star_eV
    p_nu_proxy = math.sqrt(i_nu) * math.sqrt(ratio_hat) / defect_sum
    c_nu = b_nu / p_nu_proxy

    correction_window = correction_scaffold["strongest_compare_only_correction_window"]["interval"]
    bridge_window = correction_scaffold["induced_target_containing_bridge_scalar_window"]["interval"]
    c_lower, c_upper = float(correction_window[0]), float(correction_window[1])
    b_lower, b_upper = float(bridge_window[0]), float(bridge_window[1])

    return {
        "artifact": "oph_neutrino_exact_adapter_bridge_coordinate",
        "generated_utc": _timestamp(),
        "proof_status": "exact_compare_only_bridge_coordinate_on_exact_segment_adapter",
        "scope": adapter["scope"],
        "promotable": False,
        "source_artifacts": {
            "exact_adapter": _repo_ref(ADAPTER_JSON),
            "same_label_overlap_defect_weight_normalizer": _repo_ref(NORMALIZER_JSON),
            "realized_same_label_gap_defect_readback": _repo_ref(READBACK_JSON),
            "majorana_overlap_defect_hessian": _repo_ref(HESSIAN_JSON),
            "neutrino_weighted_cycle_repair": _repo_ref(REPAIR_JSON),
            "neutrino_scale_anchor": _repo_ref(SCALE_ANCHOR_JSON),
            "neutrino_bridge_correction_invariant_scaffold": _repo_ref(CORRECTION_SCAFFOLD_JSON),
        },
        "exact_adapter_solution": {
            "tau_nu": tau_nu,
            "D_nu": d_nu,
            "p_nu": p_nu,
            "lambda_nu": lambda_nu,
            "masses_eV": list(adapter["exact_outputs"]["masses_eV"]),
            "delta_m_sq_eV2": dict(adapter["exact_outputs"]["delta_m_sq_eV2"]),
            "ratio_21_over_32": ratio_hat,
        },
        "bridge_coordinates": {
            "definition": "B_nu = P_nu * C_nu",
            "paper_facing_amplitude": {
                "symbol": "B_nu",
                "formula": "lambda_nu * q_mean^p_nu / m_star_eV",
                "value": b_nu,
            },
            "internal_proxy": {
                "symbol": "P_nu",
                "formula": "I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1",
                "value": p_nu_proxy,
                "ingredients": {
                    "I_nu": i_nu,
                    "ratio_hat": ratio_hat,
                    "sum_defect": defect_sum,
                },
            },
            "reduced_correction_invariant": {
                "symbol": "C_nu",
                "formula": "B_nu / P_nu",
                "value": c_nu,
            },
        },
        "consistency_checks": {
            "reconstruction_error_abs": abs(b_nu - p_nu_proxy * c_nu),
            "within_current_compare_only_c_window": c_lower <= c_nu <= c_upper,
            "within_current_compare_only_b_window": b_lower <= b_nu <= b_upper,
            "current_compare_only_c_window": [c_lower, c_upper],
            "current_compare_only_b_window": [b_lower, b_upper],
        },
        "theorem_boundary": {
            "status": "non_promotable_compare_only_bridge_coordinate_sidecar",
            "statement": (
                "This artifact evaluates exact compare-only bridge coordinates on the already-emitted positive selector "
                "segment after the exact two-parameter fit. It does not emit theorem-grade C_nu or B_nu and must not "
                "feed back into the theorem lane."
            ),
            "forbidden_feedback": "exact_compare_only_bridge_coordinate_must_not_feed_back_into_theorem_state_or_C_nu_emission",
        },
        "notes": [
            "The exact adapter moves only one selector coordinate and one positive rescaling on the already-emitted positive selector segment.",
            "The same emitted phase and defect stack is reused here; the only exact-adapter-dependent bridge ingredient is ratio_hat, taken from the exact two-parameter fit.",
            "This sidecar strengthens the explicit exact branch description without changing the open theorem object C_nu.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact compare-only neutrino bridge-coordinate sidecar.")
    parser.add_argument("--adapter", default=str(ADAPTER_JSON))
    parser.add_argument("--normalizer", default=str(NORMALIZER_JSON))
    parser.add_argument("--readback", default=str(READBACK_JSON))
    parser.add_argument("--hessian", default=str(HESSIAN_JSON))
    parser.add_argument("--repair", default=str(REPAIR_JSON))
    parser.add_argument("--scale-anchor", default=str(SCALE_ANCHOR_JSON))
    parser.add_argument("--correction-scaffold", default=str(CORRECTION_SCAFFOLD_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        adapter=_load_json(Path(args.adapter)),
        normalizer=_load_json(Path(args.normalizer)),
        readback=_load_json(Path(args.readback)),
        hessian=_load_json(Path(args.hessian)),
        repair=_load_json(Path(args.repair)),
        scale_anchor=_load_json(Path(args.scale_anchor)),
        correction_scaffold=_load_json(Path(args.correction_scaffold)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
