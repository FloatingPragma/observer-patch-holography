#!/usr/bin/env python3
"""Build the freeze-once D10 electroweak repair law.

Chain role: emit the exact coherent D10 repair package once one authoritative
target `W/Z` pair is frozen.

Mathematics: exact algebra on the selected D10 carrier. First solve the unique
charged anchor `tau2` that matches the frozen `W`; then factor the remaining
neutral motion into the fiber-parallel hypercharge leg plus one orthogonal
neutral-shear scalar. The same repaired coupling pair then emits one coherent
electroweak quintet.

OPH-derived inputs: the current D10 source pair together with the
machine-readable local frozen target values for `W` and `Z`.

Output: the closed freeze-once coherent repair law, the target-box dominance
diagnostic, and the explicit target-freeze / subobject split.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_FACTORIZATION_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
DEFAULT_BOX_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_box_dominance.json"
DEFAULT_SPLIT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_freeze_and_subobject_split.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _float(value: Any, default: float | None = None) -> float:
    if value is None:
        if default is None:
            raise ValueError("missing required float")
        return default
    return float(value)


@dataclass(frozen=True)
class TargetPoint:
    MW_target_gev: float
    MZ_target_gev: float
    tau2_w_anchor: float
    tauY_fiber_dagger: float
    n_fiber_dagger: float
    delta_alpha2_tree: float
    delta_alpha2_dagger: float
    MZ_fiber_after_exact_W_anchor_gev: float
    delta_MZ_after_exact_W_anchor_gev: float
    delta_MZ_after_exact_W_anchor_mev: float
    coeff_delta_n_per_gev_delta_MZ: float
    coeff_delta_alphaY_perp_per_gev_delta_MZ: float
    delta_n_tree_exact: float
    delta_n_dagger: float
    delta_alphaY_parallel: float
    delta_alphaY_perp: float
    delta_alphaY_total: float
    alpha2_star: float
    alphaY_star: float
    alpha2_dagger: float
    alphaY_dagger: float
    alpha_em_eff_inv_dagger: float
    sin2w_eff_dagger: float
    neutral_shear_share_of_total: float
    fiber_parallel_share_of_total: float


class D10Basis:
    def __init__(self, *, mw_current: float, mz_current: float, v_report: float, alphaY_mz: float, alpha2_mz: float, eta_source: float):
        self.mw_current = mw_current
        self.mz_current = mz_current
        self.v_report = v_report
        self.alphaY_mz = alphaY_mz
        self.alpha2_mz = alpha2_mz
        self.alpha_sum = alphaY_mz + alpha2_mz
        self.eta_source = eta_source
        self.alpha2_star = alpha2_mz
        self.alphaY_star = alphaY_mz * (1.0 - 2.0 * eta_source)

    def tau2_exact_from_mw(self, mw_target: float) -> float:
        return mw_target * mw_target / (self.mw_current * self.mw_current) - 1.0

    def delta_alpha2_from_mw(self, mw_target: float) -> float:
        return mw_target * mw_target / (math.pi * self.v_report * self.v_report) - self.alpha2_mz

    def tauY_fiber(self, tau2: float) -> float:
        return -(tau2 + 2.0 * self.eta_source) / (1.0 + 4.0 * tau2 * tau2)

    def n_fiber(self, tau2: float) -> float:
        return 1.0 + (self.alphaY_mz * self.tauY_fiber(tau2) + self.alpha2_mz * tau2) / self.alpha_sum

    def delta_alphaY_parallel_from_mw(self, mw_target: float) -> float:
        tau2 = self.tau2_exact_from_mw(mw_target)
        return self.alphaY_mz * (8.0 * self.eta_source * tau2 * tau2 - tau2) / (1.0 + 4.0 * tau2 * tau2)

    def mz_fiber_after_exact_w_anchor(self, mw_target: float) -> float:
        return self.v_report * math.sqrt(math.pi * self.alpha_sum * self.n_fiber(self.tau2_exact_from_mw(mw_target)))

    def delta_mz_after_exact_w_anchor(self, mw_target: float, mz_target: float) -> float:
        return mz_target - self.mz_fiber_after_exact_w_anchor(mw_target)

    def coeff_delta_n_per_gev_delta_mz(self, mw_target: float, mz_target: float) -> float:
        return (mz_target + self.mz_fiber_after_exact_w_anchor(mw_target)) / (math.pi * self.v_report * self.v_report * self.alpha_sum)

    def coeff_delta_alphaY_perp_per_gev_delta_mz(self, mw_target: float, mz_target: float) -> float:
        return (mz_target + self.mz_fiber_after_exact_w_anchor(mw_target)) / (math.pi * self.v_report * self.v_report)

    def point(self, mw_target: float, mz_target: float) -> TargetPoint:
        tau2 = self.tau2_exact_from_mw(mw_target)
        tauY = self.tauY_fiber(tau2)
        n_fiber = self.n_fiber(tau2)
        da2 = self.delta_alpha2_from_mw(mw_target)
        mz_fiber = self.mz_fiber_after_exact_w_anchor(mw_target)
        dmz_gev = mz_target - mz_fiber
        coeff_dn = self.coeff_delta_n_per_gev_delta_mz(mw_target, mz_target)
        coeff_dy_perp = self.coeff_delta_alphaY_perp_per_gev_delta_mz(mw_target, mz_target)
        delta_n = coeff_dn * dmz_gev
        delta_alphaY_parallel = self.delta_alphaY_parallel_from_mw(mw_target)
        delta_alphaY_perp = coeff_dy_perp * dmz_gev
        delta_alphaY_total = delta_alphaY_parallel + delta_alphaY_perp
        alpha2_dagger = self.alpha2_star + da2
        alphaY_dagger = self.alphaY_star + delta_alphaY_parallel + delta_alphaY_perp
        alpha_sum_dagger = alpha2_dagger + alphaY_dagger
        return TargetPoint(
            MW_target_gev=mw_target,
            MZ_target_gev=mz_target,
            tau2_w_anchor=tau2,
            tauY_fiber_dagger=tauY,
            n_fiber_dagger=n_fiber,
            delta_alpha2_tree=da2,
            delta_alpha2_dagger=da2,
            MZ_fiber_after_exact_W_anchor_gev=mz_fiber,
            delta_MZ_after_exact_W_anchor_gev=dmz_gev,
            delta_MZ_after_exact_W_anchor_mev=1000.0 * dmz_gev,
            coeff_delta_n_per_gev_delta_MZ=coeff_dn,
            coeff_delta_alphaY_perp_per_gev_delta_MZ=coeff_dy_perp,
            delta_n_tree_exact=delta_n,
            delta_n_dagger=delta_n,
            delta_alphaY_parallel=delta_alphaY_parallel,
            delta_alphaY_perp=delta_alphaY_perp,
            delta_alphaY_total=delta_alphaY_total,
            alpha2_star=self.alpha2_star,
            alphaY_star=self.alphaY_star,
            alpha2_dagger=alpha2_dagger,
            alphaY_dagger=alphaY_dagger,
            alpha_em_eff_inv_dagger=alpha_sum_dagger / (alphaY_dagger * alpha2_dagger),
            sin2w_eff_dagger=alphaY_dagger / alpha_sum_dagger,
            neutral_shear_share_of_total=(delta_alphaY_perp / delta_alphaY_total) if delta_alphaY_total != 0.0 else None,
            fiber_parallel_share_of_total=(delta_alphaY_parallel / delta_alphaY_total) if delta_alphaY_total != 0.0 else None,
        )


def _load_basis(source_pair: dict[str, Any]) -> D10Basis:
    pair = dict(source_pair.get("source_pair") or {})
    slice_payload = dict(source_pair.get("compact_hypercharge_only_mass_slice") or {})
    compact_quintet = dict(slice_payload.get("coherent_output_quintet") or {})
    return D10Basis(
        mw_current=_float(compact_quintet.get("MW_pole")),
        mz_current=_float(compact_quintet.get("MZ_pole")),
        v_report=_float(compact_quintet.get("v_report")),
        alphaY_mz=_float(pair.get("alphaY_mz")),
        alpha2_mz=_float(pair.get("alpha2_mz")),
        eta_source=_float(slice_payload.get("eta_EW")),
    )


def build_factorization_report(reference: dict[str, Any], basis: D10Basis) -> dict[str, Any]:
    mw_target = _float(reference["w_boson"]["value_gev"])
    mz_target = _float(reference["z_boson"]["value_gev"])
    point = basis.point(mw_target, mz_target)
    return {
        "artifact": "oph_d10_ew_w_anchor_neutral_shear_factorization",
        "generated_utc": _timestamp(),
        "status": "closed_freeze_once_coherent_repair_law",
        "forward_claim_allowed": False,
        "public_surface_candidate_allowed": True,
        "target_free_predictive_emission_closed": False,
        "exact_missing_law": {
            "object_id": "FreezeOnceCoherentD10ElectroweakRepairLaw_D10",
            "status": "closed",
            "kind": "freeze_once_authoritative_target_repair",
            "proof_gate": "single_family_single_P_no_mixed_readout",
            "stricter_still_open_object": "EWTargetFreeRepairValueLaw_D10",
        },
        "target_spec": {
            "kind": "local_machine_readable_reference_surface",
            "source_label": reference["w_boson"]["source"]["label"],
            "epoch": "2025",
            "MW_pole_central_gev": mw_target,
            "MW_pole_err_gev": _float(reference["w_boson"].get("error_plus_gev"), 0.0),
            "MZ_pole_central_gev": mz_target,
            "MZ_pole_err_gev": _float(reference["z_boson"].get("error_plus_gev"), 0.0),
            "note": "This law is exact on the frozen authoritative D10 target surface. The stronger target-free emission from P alone remains open.",
        },
        "new_smaller_primitive": {
            "proposed_object_id": "EWAnchoredNeutralShearPrimitive_D10",
            "meaning": "Once the target pair is frozen, exact D10 repair factors into an exact W anchor plus one remaining neutral-shear scalar.",
            "strictly_smaller_than": "D10RepairBranchBeyondCurrentCarrier",
            "depends_on": [
                "EWSinglePostTransportTreeIdentity_D10",
                "oph_d10_target_spec_freeze",
            ],
        },
        "current_basis": {
            "MW_current_gev": basis.mw_current,
            "MZ_current_gev": basis.mz_current,
            "v_report_gev": basis.v_report,
            "alpha2_mz": basis.alpha2_mz,
            "alphaY_mz": basis.alphaY_mz,
            "alpha_sum_mz": basis.alpha_sum,
            "eta_source": basis.eta_source,
            "alpha2_star": basis.alpha2_star,
            "alphaY_star": basis.alphaY_star,
        },
        "theorem": {
            "name": "Freeze-once coherent D10 electroweak repair law",
            "statement": (
                "Fix the current D10 carrier basis and one authoritative frozen target pair (MW_target, MZ_target) with 0 < MW_target < MZ_target. "
                "Then there exists a unique coherent repair package "
                "(delta_alpha2_dagger, delta_alphaY_parallel, delta_alphaY_perp), equivalently "
                "(tau2_w_anchor, delta_n_dagger). The corresponding repaired coupling pair "
                "(alpha2_dagger, alphaY_dagger) emits one coherent electroweak quintet on a single family, "
                "with no mixed readout and no hidden selector."
            ),
            "formulas": {
                "tau2_w_anchor": "MW_target^2 / MW_current^2 - 1",
                "delta_alpha2_dagger": "MW_target^2 / (pi * v_report^2) - alpha2_mz",
                "tauY_fiber": "-(tau2_w_anchor + 2*eta_source) / (1 + 4*tau2_w_anchor^2)",
                "n_fiber": "1 + (alphaY_mz * tauY_fiber + alpha2_mz * tau2_w_anchor) / (alphaY_mz + alpha2_mz)",
                "MZ_fiber_after_exact_W_anchor": "v_report * sqrt(pi * (alphaY_mz + alpha2_mz) * n_fiber)",
                "delta_MZ_after_exact_W_anchor": "MZ_target - MZ_fiber_after_exact_W_anchor",
                "delta_n_dagger": "((MZ_target + MZ_fiber_after_exact_W_anchor) * delta_MZ_after_exact_W_anchor) / (pi * v_report^2 * (alphaY_mz + alpha2_mz))",
                "delta_alphaY_parallel": "alphaY_mz * (8*eta_source*tau2_w_anchor^2 - tau2_w_anchor) / (1 + 4*tau2_w_anchor^2)",
                "delta_alphaY_perp": "((MZ_target + MZ_fiber_after_exact_W_anchor) * delta_MZ_after_exact_W_anchor) / (pi * v_report^2)",
                "alpha2_dagger": "alpha2_mz + delta_alpha2_dagger = MW_target^2 / (pi * v_report^2)",
                "alphaY_dagger": "alphaY_mz * (1 - 2*eta_source) + delta_alphaY_parallel + delta_alphaY_perp",
                "alpha_em_eff_inv_dagger": "(alphaY_dagger + alpha2_dagger) / (alphaY_dagger * alpha2_dagger)",
                "sin2w_eff_dagger": "alphaY_dagger / (alphaY_dagger + alpha2_dagger)",
            },
            "equivalent_single_remaining_scalars": [
                "delta_MZ_after_exact_W_anchor",
                "delta_n_dagger",
                "delta_alphaY_perp",
            ],
            "future_forward_split_suggested_by_theorem": [
                "EWChargedAnchorLaw_D10",
                "EWNeutralShearLaw_D10",
            ],
        },
        "repair_law": {
            "chart_coordinates": {
                "tau2_w_anchor": point.tau2_w_anchor,
                "delta_n_dagger": point.delta_n_dagger,
            },
            "coherent_repair_package": {
                "delta_alpha2_dagger": point.delta_alpha2_dagger,
                "delta_alphaY_parallel": point.delta_alphaY_parallel,
                "delta_alphaY_perp": point.delta_alphaY_perp,
            },
            "equivalent_two_coordinate_chart": "(tau2_w_anchor, delta_n_dagger)",
        },
        "central_target_point": asdict(point),
        "coherent_repaired_couplings": {
            "alpha2_star": point.alpha2_star,
            "alphaY_star": point.alphaY_star,
            "alpha2_dagger": point.alpha2_dagger,
            "alphaY_dagger": point.alphaY_dagger,
        },
        "coherent_repaired_quintet": {
            "MW_pole": point.MW_target_gev,
            "MZ_pole": point.MZ_target_gev,
            "alpha_em_eff_inv": point.alpha_em_eff_inv_dagger,
            "sin2w_eff": point.sin2w_eff_dagger,
            "v_report": basis.v_report,
        },
        "conclusion": {
            "meaning": "On the frozen authoritative D10 target surface, the exact repair package and coherent repaired quintet are fully emitted.",
            "still_compare_only": False,
            "stricter_still_open_object": "EWTargetFreeRepairValueLaw_D10",
            "stricter_still_open_statement": "Emit the same nonzero repair directly from P alone with no frozen target input.",
        },
    }


def build_box_report(reference: dict[str, Any], basis: D10Basis) -> dict[str, Any]:
    mw_c = _float(reference["w_boson"]["value_gev"])
    mz_c = _float(reference["z_boson"]["value_gev"])
    mw_err = _float(reference["w_boson"].get("error_plus_gev"), 0.0)
    mz_err = _float(reference["z_boson"].get("error_plus_gev"), 0.0)
    corners = [
        basis.point(mw_c - mw_err, mz_c - mz_err),
        basis.point(mw_c - mw_err, mz_c + mz_err),
        basis.point(mw_c + mw_err, mz_c - mz_err),
        basis.point(mw_c + mw_err, mz_c + mz_err),
    ]

    def bounds(getter: str) -> list[float]:
        values = [float(getattr(point, getter)) for point in corners]
        return [min(values), max(values)]

    return {
        "artifact": "oph_d10_ew_w_anchor_neutral_shear_box_dominance",
        "generated_utc": _timestamp(),
        "status": "diagnostic_compare_only",
        "forward_claim_allowed": False,
        "target_box": {
            "kind": "local_machine_readable_reference_1sigma_box",
            "MW_box_gev": [mw_c - mw_err, mw_c + mw_err],
            "MZ_box_gev": [mz_c - mz_err, mz_c + mz_err],
        },
        "theorem": {
            "name": "Reference-box neutral-shear dominance theorem",
            "statement": (
                "Over the current local reference 1-sigma W/Z box, the W-anchored remaining D10 repair scalar stays positive and the neutral-shear component dominates the total hypercharge repair."
            ),
        },
        "corner_samples": [asdict(point) for point in corners],
        "bounds": {
            "delta_MZ_after_exact_W_anchor_mev": bounds("delta_MZ_after_exact_W_anchor_mev"),
            "delta_n_tree_exact": bounds("delta_n_tree_exact"),
            "delta_alphaY_parallel": bounds("delta_alphaY_parallel"),
            "delta_alphaY_perp": bounds("delta_alphaY_perp"),
            "neutral_shear_share_of_total": bounds("neutral_shear_share_of_total"),
            "fiber_parallel_share_of_total": bounds("fiber_parallel_share_of_total"),
        },
        "verdict": {
            "delta_MZ_after_exact_W_anchor_positive_everywhere": all(point.delta_MZ_after_exact_W_anchor_gev > 0.0 for point in corners),
            "neutral_shear_dominates_total_hypercharge_repair_everywhere": all(
                point.neutral_shear_share_of_total is not None and point.neutral_shear_share_of_total > 0.5
                for point in corners
            ),
        },
    }


def build_target_freeze_split(reference: dict[str, Any], basis: D10Basis) -> dict[str, Any]:
    factorization = build_factorization_report(reference, basis)
    central = factorization["central_target_point"]
    return {
        "artifact": "oph_d10_ew_target_freeze_and_subobject_split",
        "generated_utc": _timestamp(),
        "status": "closed_freeze_once_subobject_split",
        "forward_claim_allowed": False,
        "target_freeze_required": True,
        "reason": "Exact D10 repair is a value-law problem only after one authoritative W/Z target pair is frozen.",
        "frozen_reference_target": {
            "MW_pole_target_gev": central["MW_target_gev"],
            "MZ_pole_target_gev": central["MZ_target_gev"],
        },
        "subobject_split": {
            "charged_anchor_object": "EWChargedAnchorLaw_D10",
            "neutral_shear_object": "EWNeutralShearLaw_D10",
            "current_frozen_target_values": {
                "tau2_w_anchor": central["tau2_w_anchor"],
                "delta_n_dagger": central["delta_n_dagger"],
                "delta_n_tree_exact": central["delta_n_tree_exact"],
                "delta_alpha2_dagger": central["delta_alpha2_dagger"],
                "delta_alpha2_tree": central["delta_alpha2_tree"],
                "delta_alphaY_parallel": central["delta_alphaY_parallel"],
                "delta_alphaY_perp": central["delta_alphaY_perp"],
            },
        },
        "note": "This file records the exact freeze-once subobject split. The only stricter remaining D10 step is the target-free repair value law from P alone.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 W-anchor / neutral-shear factorization diagnostics.")
    parser.add_argument("--references", default=str(REFERENCE_JSON))
    parser.add_argument("--source-pair", default=str(SOURCE_PAIR_JSON))
    parser.add_argument("--readout", default="", help="Unused legacy compatibility flag.")
    parser.add_argument("--factorization-output", default=str(DEFAULT_FACTORIZATION_OUT))
    parser.add_argument("--box-output", default=str(DEFAULT_BOX_OUT))
    parser.add_argument("--split-output", default=str(DEFAULT_SPLIT_OUT))
    args = parser.parse_args()

    references = json.loads(Path(args.references).read_text(encoding="utf-8"))["entries"]
    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    basis = _load_basis(source_pair)

    outputs = [
        (Path(args.factorization_output), build_factorization_report(references, basis)),
        (Path(args.box_output), build_box_report(references, basis)),
        (Path(args.split_output), build_target_freeze_split(references, basis)),
    ]
    for path, payload in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"saved: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
