#!/usr/bin/env python3
"""Export the current D10 pole/effective transport-kernel boundary artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_transport_kernel.json"
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(family: dict[str, object], source_pair_payload: dict[str, object] | None = None) -> dict[str, object]:
    core_source = dict(family.get("core_source", {}))
    reported_outputs = dict(family.get("reported_outputs", {}))
    generator_keys = list(family.get("generator_keys_unrounded", []))
    coherence_witness = dict(family.get("coherence_witness", {}))
    base_w = float(reported_outputs.get("m_w_run", 0.0))
    base_z = float(reported_outputs.get("m_z_run", 0.0))
    base_alpha = float(reported_outputs.get("alpha_em_inv_mz", 0.0))
    base_sin2w = float(reported_outputs.get("sin2w_mz", 0.0))
    base_v = float(reported_outputs.get("v", 0.0))
    shared_source = {
        "family_source_id": "d10_running_tree",
        "scheme_id": "freeze_once",
        "origin_kernel_id": "EWTransportKernel_D10",
    }
    frozen_candidate_outputs = {
        "MW_pole": 80.3692,
        "MZ_pole": 91.1880,
        "alpha_em_eff_inv": 127.930,
        "sin2w_eff": 0.23161,
        "v_report": base_v,
    }
    frozen_candidate_scalars = {
        "delta_alpha": float(frozen_candidate_outputs["alpha_em_eff_inv"] / base_alpha - 1.0),
        "delta_kappa": float(frozen_candidate_outputs["sin2w_eff"] / base_sin2w - 1.0),
        "delta_rho": float(frozen_candidate_outputs["MZ_pole"] / base_z - 1.0),
        "delta_rW": float(frozen_candidate_outputs["MW_pole"] / base_w - frozen_candidate_outputs["MZ_pole"] / base_z),
    }
    frozen_candidate_mass_ratio_residual = float(
        1.0 - (frozen_candidate_outputs["MW_pole"] / frozen_candidate_outputs["MZ_pole"]) ** 2 - frozen_candidate_outputs["sin2w_eff"]
    )

    source_pair_payload = source_pair_payload or {}
    return {
        "artifact": "oph_d10_ew_transport_kernel",
        "generated_utc": _timestamp(),
        "kernel_id": "EWTransportKernel_D10",
        "active_readout_selector_candidate_id": "EWActiveReadoutSelector_D10",
        "active_readout_selector_status": "candidate_only",
        "post_selector_missing_object": "EWSharedScalarReadoutPackage_D10",
        "smallest_populated_missing_object": "EWSharedScalarReadoutPackage_D10",
        "strict_missing_object_beneath_transport_entry_preimage": "EWGaugeSourceTransportPair_D10",
        "source_pair_population_status": source_pair_payload.get("two_scalar_population_status"),
        "shared_scalar_readout_map_status": "candidate_only",
        "current_emission_status": "zero_benchmark_only__frozen_nonzero_calibration_candidate_available",
        "next_live_mover": "populated_shared_scalar_package",
        "shared_scalar_package_id": "Sigma_EW_D10",
        "shared_scalar_package_theorem_candidate": "EWSharedScalarReadoutPackage_D10",
        "shared_scalar_emitter_candidate": "EWSharedScalarEmitter_D10",
        "readout_coherence_clause_id": "EWTransportReadoutCoherence_D10",
        "family_source_artifact": str(family.get("artifact")),
        "family_source_id": str(family.get("observable_family_id")),
        "observable_family_id": "d10_pole_effective",
        "scheme_id": "freeze_once",
        "observable_quartet": ["W", "Z", "alpha_em", "sin2w"],
        "core_source_ref": {
            "generator_keys_unrounded": generator_keys,
            "generator_values_unrounded": core_source,
        },
        "neutral_transport": {
            "basis": ["A", "Z"],
            "entries": ["Pi_AA", "Pi_AZ", "Pi_ZZ"],
            "proof_status": "candidate_only",
        },
        "charged_transport": {
            "basis": ["W"],
            "entries": ["Pi_WW"],
            "proof_status": "candidate_only",
        },
        "derived_readouts": {
            "delta_alpha_from": ["Pi_AA", "Pi_AZ"],
            "delta_kappa_from": ["Pi_AA", "Pi_AZ", "Pi_ZZ"],
            "delta_rho_from": ["Pi_WW", "Pi_ZZ"],
            "delta_rW_from": ["Pi_WW", "Pi_AA", "Pi_AZ", "Pi_ZZ"],
        },
        "transport_source_entries": {
            "neutral": ["Pi_AA", "Pi_AZ", "Pi_ZZ"],
            "charged": ["Pi_WW"],
        },
        "transport_entry_values": {
            "Pi_AA": None,
            "Pi_AZ": None,
            "Pi_ZZ": None,
            "Pi_WW": None,
        },
        "source_transport_pair_artifact": source_pair_payload.get("artifact"),
        "source_transport_pair_symbol": source_pair_payload.get("source_pair_symbol"),
        "source_transport_pair": source_pair_payload.get("source_pair"),
        "shared_scalar_package": {
            "symbol": "Sigma_EW_D10 = (delta_alpha, delta_kappa, delta_rho, delta_rW)",
            "zero_normalization_rule": "F_i(0; x) = x for every electroweak readout row",
            "common_provenance_required": True,
        },
        "shared_scalar_values": {
            "delta_alpha": 0.0,
            "delta_kappa": 0.0,
            "delta_rho": 0.0,
            "delta_rW": 0.0,
        },
        "frozen_nonzero_calibration_candidate": {
            "status": "calibration_readout_candidate_not_oph_only",
            "observable_quartet": {
                "MW_pole_gev": frozen_candidate_outputs["MW_pole"],
                "MZ_pole_gev": frozen_candidate_outputs["MZ_pole"],
                "alpha_em_eff_inv": frozen_candidate_outputs["alpha_em_eff_inv"],
                "sin2w_eff": frozen_candidate_outputs["sin2w_eff"],
            },
            "shared_scalar_values": dict(frozen_candidate_scalars),
            "coherent_output_quintet": dict(frozen_candidate_outputs),
            "mass_ratio_residual": frozen_candidate_mass_ratio_residual,
            "source_kind": "frozen_external_calibration_readout",
            "common_provenance_fields": dict(shared_source),
            "transport_support_required": ["Pi_AA", "Pi_AZ", "Pi_ZZ", "Pi_WW"],
            "not_yet_oph_only_reason": "transport_entry_preimage_unpopulated",
        },
        "external_nonzero_quartet_image_test": source_pair_payload.get("external_nonzero_quartet_image_test"),
        "first_nonzero_oph_seed_trial": source_pair_payload.get("first_nonzero_oph_seed_trial"),
        "source_pair_two_scalar_family": source_pair_payload.get("family_coordinates"),
        "source_pair_special_slices": source_pair_payload.get("special_slices"),
        "source_pair_first_nonzero_trial": source_pair_payload.get("first_nonzero_oph_seed_trial"),
        "emitter_family": {
            "kind": "zero_normalized_affine_relative_emitter",
            "delta_alpha_formula": "alpha_em_eff_inv / a0 - 1",
            "delta_kappa_formula": "sin2w_eff / s0 - 1",
            "delta_rho_formula": "MZ_pole / Z0 - 1",
            "delta_rW_formula": "MW_pole / W0 - MZ_pole / Z0",
            "inverse_readout": {
                "alpha_em_eff_inv": "a0 * (1 + delta_alpha)",
                "sin2w_eff": "s0 * (1 + delta_kappa)",
                "MZ_pole": "Z0 * (1 + delta_rho)",
                "MW_pole": "W0 * (1 + delta_rho + delta_rW)",
                "v_report": "v0",
            },
        },
        "reported_outputs": {
            "MW_pole": "derived",
            "MZ_pole": "derived",
            "alpha_em_eff_inv": "derived",
            "sin2w_eff": "derived",
            "v_report": {
                "kind": "inherit_running_core",
            },
        },
        "readout_assignments": {
            "MW_pole": {
                "origin_kernel_id": "EWTransportKernel_D10",
                "via": ["delta_rW", "delta_rho"],
                "zero_normalized": True,
            },
            "MZ_pole": {
                "origin_kernel_id": "EWTransportKernel_D10",
                "via": ["delta_rho"],
                "zero_normalized": True,
            },
            "alpha_em_eff_inv": {
                "origin_kernel_id": "EWTransportKernel_D10",
                "via": ["delta_alpha"],
                "zero_normalized": True,
            },
            "sin2w_eff": {
                "origin_kernel_id": "EWTransportKernel_D10",
                "via": ["delta_kappa"],
                "zero_normalized": True,
            },
            "v_report": {
                "kind": "inherit_running_core",
            },
        },
        "active_readout_family_id": "d10_running_tree",
        "active_readout_reason": "fallback_until_common_readout_certified",
        "selector_rule": {
            "if_common_readout_certified": "d10_pole_effective",
            "if_not_common_readout_certified": "d10_running_tree",
        },
        "coherent_quintet_when_running_family": {
            "W": base_w,
            "Z": base_z,
            "alpha_em_inv": base_alpha,
            "sin2w": base_sin2w,
            "v": base_v,
        },
        "base_quintet_running": {
            "MW_pole": base_w,
            "MZ_pole": base_z,
            "alpha_em_eff_inv": base_alpha,
            "sin2w_eff": base_sin2w,
            "v_report": base_v,
        },
        "base_quintet_source_path": {
            "MW_pole": "reported_outputs.m_w_run",
            "MZ_pole": "reported_outputs.m_z_run",
            "alpha_em_eff_inv": "reported_outputs.alpha_em_inv_mz",
            "sin2w_eff": "reported_outputs.sin2w_mz",
            "v_report": "reported_outputs.v",
        },
        "coherent_output_quintet": {
            "MW_pole": base_w,
            "MZ_pole": base_z,
            "alpha_em_eff_inv": base_alpha,
            "sin2w_eff": base_sin2w,
            "v_report": base_v,
        },
        "zero_scalar_benchmark": {
            "rule": "R_EW_D10(Q_run_D10, 0) = Q_run_D10",
            "W": base_w,
            "Z": base_z,
            "alpha_em_inv": base_alpha,
            "sin2w": base_sin2w,
            "v": base_v,
        },
        "quartet_source_lock": {
            "W": "d10_running_tree",
            "Z": "d10_running_tree",
            "alpha_em": "d10_running_tree",
            "sin2w": "d10_running_tree",
            "v": "d10_running_tree",
        },
        "provenance_equality_fields": [
            "family_source_id",
            "scheme_id",
            "origin_kernel_id",
        ],
        "scalar_provenance": {
            "delta_alpha": dict(shared_source),
            "delta_kappa": dict(shared_source),
            "delta_rho": dict(shared_source),
            "delta_rW": dict(shared_source),
        },
        "family_purity_gate": {
            "no_run_pole_mix": True,
            "z_only_surrogate_forbidden": True,
            "orphan_scalar_corrections_forbidden": True,
            "independently_rounded_targets_forbidden": True,
            "common_readout_certified": False,
        },
        "common_provenance_witness": {
            "all_equal_family_source_id": True,
            "all_equal_scheme_id": True,
            "all_equal_origin_kernel_id": True,
            "shared_source": dict(shared_source),
        },
        "coherence_witness": {
            "running_mass_ratio_residual": coherence_witness.get("running_mass_ratio_residual"),
            "stage3_mass_ratio_residual": coherence_witness.get("stage3_mass_ratio_residual"),
            "mixed_stage3_mass_ratio_residual": coherence_witness.get("stage3_mass_ratio_residual"),
            "mixed_sources_detected": coherence_witness.get("mixed_sources_detected", False),
        },
        "promotion_gate": {
            "mixed_scheme": False,
            "z_only_surrogate_forbidden": True,
            "provenance_equality_required": True,
            "independently_rounded_targets_forbidden": True,
            "common_readout_certified": False,
            "populated_shared_scalar_package_required": True,
            "no_orphan_scalar_corrections": True,
            "smaller_exact_missing_clause": "EWTransportReadoutCoherence_D10",
            "strictly_smaller_next_subclause": "EWScalarProvenanceEquality_D10",
        },
        "notes": [
            "This artifact is the exact boundary for any pole/effective electroweak family built on top of the current exact D10 running-family core.",
            "The selector phase is already behind this lane: until one coherent pole/effective family is certified, the active public quintet should remain the running-family quintet instead of mixing W_run with a Stage-3 Z-only surrogate.",
            "The next live mover is the populated shared scalar package Sigma_EW_D10 with one zero-normalized readout family beneath EWTransportKernel_D10 rather than another per-observable family split.",
            "The running-family base quintet is now exported explicitly so a future populated Sigma_EW_D10 package can read directly from one fixed running source rather than a partially implicit benchmark.",
            "The smallest constructive object beneath the public package boundary is now the zero-normalized affine relative emitter EWSharedScalarEmitter_D10.",
            "At the current boundary, the emitted scalar tuple is the zero benchmark Sigma_EW_D10 = (0,0,0,0), which keeps the coherent emitted quintet equal to the running-family quintet until the first non-zero populated package is justified.",
            "A frozen non-zero calibration-readout candidate is now recorded on the same emitter formulas, but it is not promoted to OPH-only status because the local transport-entry preimage on (Pi_AA, Pi_AZ, Pi_ZZ, Pi_WW) is still unpopulated.",
            "Any non-zero future scalar package must still be emitted on the same fixed running-family base quintet and the same transport-entry support (Pi_AA, Pi_AZ, Pi_ZZ, Pi_WW).",
            "The four-entry A/Z/W shell is no longer the smallest exact source-side object; beneath it the reduced missing source is the D10 gauge source transport pair Tau_EW_D10 = (tau_Y, tau_2).",
            "The strict smallest live missing object is therefore the source-side pair EWGaugeSourceTransportPair_D10 rather than the full Pi-shell itself.",
            "The current source-pair artifact now records the simplest augmentation beyond the tested anti-diagonal seed: the two-scalar family (sigma_EW, eta_EW) on that same pair together with its first nonzero OPH seed trial.",
            "The smaller exact missing clause is EWTransportReadoutCoherence_D10: either W, Z, alpha_em, and sin^2(theta_W) all stay on the running family or all move together to one common pole/effective family with one shared kernel and scheme.",
            "Inside that criterion, the strictly smaller exact subclause is EWScalarProvenanceEquality_D10: delta_alpha, delta_kappa, delta_rho, and delta_rW must expose one common family source, one frozen scheme, and one origin kernel.",
            "A singleton Z-only surrogate is formally non-promotable under this artifact.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 electroweak transport-kernel boundary artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    family = json.loads(Path(args.input).read_text(encoding="utf-8"))
    source_pair_path = Path(args.source_pair)
    source_pair_payload = json.loads(source_pair_path.read_text(encoding="utf-8")) if source_pair_path.exists() else None
    artifact = build_artifact(family, source_pair_payload)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
