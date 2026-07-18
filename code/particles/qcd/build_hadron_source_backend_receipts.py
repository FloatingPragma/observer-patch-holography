#!/usr/bin/env python3
"""Build the hadronic precision backend receipt scaffold.

This is a simulator support surface, not a lattice-QCD solver. It emits the
minimum receipt tree required by the hadronic precision endpoint audit so
future numerical backends have a frozen schema, claim tier, and no-target-leak
contract before any comparison with alpha(0), g-2, e+e- data, rare decays, or
hadron masses.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "qcd" / "hadron_source_backend"
DEFAULT_LAMBDA = ROOT / "particles" / "runs" / "qcd" / "lambda_msbar_descendant.json"

CLAIMS = {
    "CONVENTIONAL_QCD_REFERENCE",
    "SOURCE_PROTOTYPE_NOT_PROMOTED",
    "SOURCE_INTERVAL_PROMOTED",
    "EMPIRICAL_CLOSURE_ONLY",
    "COMPARISON_ONLY",
}

CLAIM_TIERS = ("H0", "H1", "H2", "H3", "H4", "H5", "H6", "H7")

FORBIDDEN_TARGETS = (
    "CODATA_ALPHA",
    "MUON_G_MINUS_2",
    "EE_TO_HADRONS",
    "RARE_DECAY_DATA",
    "HADRON_MASS_TARGETS",
    "PDG_QCD_FITS",
)

REQUIRED_FILES = (
    "manifest.json",
    "source_dag.json",
    "qcd_ensemble/quotient_schema.json",
    "qcd_ensemble/gamma_groupoid.json",
    "qcd_ensemble/base_measure.json",
    "qcd_ensemble/source_action.json",
    "qcd_ensemble/source_parameter_map.json",
    "qcd_ensemble/coarse_maps.json",
    "vacuum/euclidean_slab.json",
    "vacuum/transfer_operator.json",
    "vacuum/reflection_positivity.json",
    "vacuum/vacuum_promotion.json",
    "currents/ward_current_definition.json",
    "currents/current_normalization_ZV.json",
    "currents/contact_terms.json",
    "currents/ward_residuals.csv",
    "correlators/vector_current_2pt_raw.json",
    "correlators/vector_current_2pt_covariance.json",
    "correlators/disconnected_diagrams.json",
    "correlators/autocorrelation_report.json",
    "spectral/moments.json",
    "spectral/hankel_positivity.json",
    "spectral/stieltjes_bounds.json",
    "spectral/J24Q.json",
    "spectral/omegaQ.json",
    "spectral/spectral_interval.json",
    "endpoint/kernel_definition.json",
    "endpoint/Xi_same_scheme.json",
    "endpoint/Delta_had_interval.json",
    "endpoint/ATh_interval.json",
    "endpoint/pixel_contraction_interval.json",
    "higher_point/Q4_HLbL_receipt.json",
    "higher_point/transition_B_to_K_receipt.json",
    "higher_point/transition_Sigma_to_p_receipt.json",
    "controls/no_target_leak_dag.json",
    "controls/empirical_data_exclusion_manifest.json",
    "controls/frozen_code_hashes.json",
    "controls/replay_receipts.json",
    "controls/comparison_data_manifest.json",
    "claim.md",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_text(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_optional_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"exists": False, "path": path.relative_to(ROOT).as_posix()}
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.setdefault("exists", True)
    payload.setdefault("path", path.relative_to(ROOT).as_posix())
    return payload


def base_receipt(name: str, *, tier: str, promoted: bool = False) -> dict[str, Any]:
    return {
        "artifact": name,
        "claim_tier": tier,
        "promotion_allowed": promoted,
        "status": "MISSING_SOURCE_CERTIFICATE" if not promoted else "SOURCE_INTERVAL_PROMOTED",
        "external_targets_used": [],
        "forbidden_targets": list(FORBIDDEN_TARGETS),
    }


def source_dag(*, claim: str, tier: str, lambda_payload: dict[str, Any]) -> dict[str, Any]:
    lambda_inputs = lambda_payload.get("source_inputs", {}) if lambda_payload.get("exists") else {}
    return {
        "artifact": "hadron_source_backend_source_dag",
        "claim": claim,
        "claim_tier": tier,
        "source_inputs": {
            "OPH_PIXEL_SOURCE": "P",
            "D10_QCD_DESCENDANT": lambda_payload.get("path"),
            "lambda_msbar_source_inputs": lambda_inputs,
        },
        "forbidden_targets": list(FORBIDDEN_TARGETS),
        "comparison_data_excluded": True,
        "no_target_leak_status": "PASS_EMPTY_COMPARISON_DAG",
        "required_missing_source_objects": [
            "QCD_SOURCE_PARAMETER_MAP",
            "OPH_QCD_QUOTIENT_ENSEMBLE_LAW",
            "EUCLIDEAN_QCD_SLAB_TRANSFER",
            "WARD_CURRENT_NORMALIZATION_LEDGER",
            "POSITIVE_STIELTJES_JACOBI_EXPORT",
            "SAME_SCHEME_XI_REMAINDER",
        ],
    }


def build_payloads(*, claim: str, tier: str, lambda_payload: dict[str, Any]) -> dict[str, str | dict[str, Any]]:
    promoted = claim == "SOURCE_INTERVAL_PROMOTED"
    source_action_status = "DECLARED_PLACEHOLDER_NOT_SOURCE_CLOSED"
    if claim == "CONVENTIONAL_QCD_REFERENCE":
        source_action_status = "CONVENTIONAL_QCD_REFERENCE"

    payloads: dict[str, str | dict[str, Any]] = {
        "source_dag.json": source_dag(claim=claim, tier=tier, lambda_payload=lambda_payload),
        "qcd_ensemble/quotient_schema.json": {
            **base_receipt("qcd_quotient_schema", tier=tier, promoted=False),
            "Sigma_QCD_r": "finite gauge links, quark fields, boundary conditions, currents, regulator metadata",
            "Gamma_QCD_r": "gauge representatives, mesh labels, port labels, worker ids, seed presentation artifacts",
            "Q_QCD_r": "Sigma_QCD_r/Gamma_QCD_r",
        },
        "qcd_ensemble/gamma_groupoid.json": {
            **base_receipt("qcd_gamma_groupoid", tier=tier),
            "quotiented_labels": [
                "gauge_representatives",
                "mesh_labels",
                "port_labels",
                "worker_ids",
                "repair_schedule_ids",
                "random_seed_presentations",
                "inert_ancillas",
            ],
        },
        "qcd_ensemble/base_measure.json": {
            **base_receipt("qcd_base_measure", tier=tier),
            "base_measure": "m_r^QCD",
            "status": "REQUIRED_NOT_POPULATED",
        },
        "qcd_ensemble/source_action.json": {
            **base_receipt("qcd_source_action", tier=tier),
            "law": "mu_r(q;P)=Z_r^-1 m_r(q) exp[-S_r(q;P)]",
            "status": source_action_status,
        },
        "qcd_ensemble/source_parameter_map.json": {
            **base_receipt("qcd_source_parameter_map", tier=tier),
            "required_map": "P -> (g3, theta_QCD, m_u, m_d, m_s, m_c, m_b, m_t, Z_scheme)",
            "lambda_msbar_descendant": lambda_payload,
            "status": "PARTIAL_DESCENDANT_PRESENT_SOURCE_MAP_OPEN",
        },
        "qcd_ensemble/coarse_maps.json": {
            **base_receipt("qcd_coarse_maps", tier=tier),
            "required": "c_sr refinement maps",
        },
        "vacuum/euclidean_slab.json": {
            **base_receipt("qcd_euclidean_slab", tier=tier),
            "required_tuple": ["Q_r^QCD", "m_r^0", "J_r", "V_r", "a_t_r", "Theta_RP"],
        },
        "vacuum/transfer_operator.json": {
            **base_receipt("qcd_transfer_operator", tier=tier),
            "required_operator": "T_r=exp[-a_t,r(H_r^E-E_0,r)]",
        },
        "vacuum/reflection_positivity.json": {
            **base_receipt("qcd_reflection_positivity", tier=tier),
            "status": "REFLECTION_POSITIVITY_CERTIFICATE_REQUIRED",
        },
        "vacuum/vacuum_promotion.json": {
            **base_receipt("qcd_vacuum_promotion", tier=tier),
            "status": "HMC_STATIONARY_DISTRIBUTION_NOT_ENOUGH",
        },
        "currents/ward_current_definition.json": {
            **base_receipt("ward_current_definition", tier=tier),
            "current": "J_Q^{W,R,mu}",
            "required": ["Ward_projection", "singlet_projection", "confinement_projection", "Z_V", "contact_terms"],
        },
        "currents/current_normalization_ZV.json": {
            **base_receipt("current_normalization_ZV", tier=tier),
            "required_ledger": ["Z_V", "J_conserved", "J_local", "J_contact", "Pi_Ward", "Omega_Q", "scheme"],
        },
        "currents/contact_terms.json": {
            **base_receipt("ward_contact_terms", tier=tier),
            "status": "CONTACT_TERM_LEDGER_REQUIRED",
        },
        "currents/ward_residuals.csv": "momentum,residual,bound,status\n",
        "correlators/vector_current_2pt_raw.json": {
            **base_receipt("vector_current_2pt_raw", tier=tier),
            "observable": "C_QQ(t)",
            "status": "NOT_COMPUTED",
        },
        "correlators/vector_current_2pt_covariance.json": {
            **base_receipt("vector_current_2pt_covariance", tier=tier),
            "status": "NOT_COMPUTED",
        },
        "correlators/disconnected_diagrams.json": {
            **base_receipt("disconnected_diagram_ledger", tier=tier),
            "status": "NOT_COMPUTED",
        },
        "correlators/autocorrelation_report.json": {
            **base_receipt("autocorrelation_report", tier=tier),
            "status": "NOT_COMPUTED",
        },
        "spectral/moments.json": {
            **base_receipt("spectral_moments", tier=tier),
            "derivation_required": "moments from Euclidean current correlators, not endpoint matching",
        },
        "spectral/hankel_positivity.json": {
            **base_receipt("hankel_positivity", tier=tier),
            "required": "Hankel matrices H^(k) positive semidefinite",
        },
        "spectral/stieltjes_bounds.json": {
            **base_receipt("stieltjes_bounds", tier=tier),
            "required": "positive Stieltjes interval bounds",
        },
        "spectral/J24Q.json": {
            **base_receipt("J24Q", tier=tier),
            "status": "MUST_DERIVE_FROM_MOMENTS_OR_LANCZOS",
        },
        "spectral/omegaQ.json": {
            **base_receipt("omegaQ", tier=tier),
            "status": "CURRENT_NORMALIZATION_REQUIRED",
        },
        "spectral/spectral_interval.json": {
            **base_receipt("spectral_interval", tier=tier),
            "status": "NO_SOURCE_INTERVAL",
        },
        "endpoint/kernel_definition.json": {
            **base_receipt("hadronic_endpoint_kernel", tier=tier),
            "Delta_had": "mZ^2/(3*pi) integral d rho_Q^(2)(s)/(s*(s+mZ^2)) + Xi_Q",
        },
        "endpoint/Xi_same_scheme.json": {
            **base_receipt("Xi_same_scheme", tier=tier),
            "required_terms": ["sub", "EW", "QED", "isospin", "heavy", "FV", "cont", "contact"],
        },
        "endpoint/Delta_had_interval.json": {
            **base_receipt("Delta_had_interval", tier=tier, promoted=promoted),
            "interval": None,
        },
        "endpoint/ATh_interval.json": {
            **base_receipt("ATh_interval", tier=tier, promoted=promoted),
            "interval": None,
        },
        "endpoint/pixel_contraction_interval.json": {
            **base_receipt("pixel_contraction_interval", tier=tier, promoted=promoted),
            "interval": None,
        },
        "higher_point/Q4_HLbL_receipt.json": {
            **base_receipt("Q4_HLbL_receipt", tier=tier),
            "status": "TWO_POINT_MEASURE_INSUFFICIENT",
        },
        "higher_point/transition_B_to_K_receipt.json": {
            **base_receipt("B_to_K_transition_spectral_receipt", tier=tier),
            "status": "NOT_STARTED",
        },
        "higher_point/transition_Sigma_to_p_receipt.json": {
            **base_receipt("Sigma_to_p_transition_spectral_receipt", tier=tier),
            "status": "NOT_STARTED",
        },
        "controls/no_target_leak_dag.json": {
            "artifact": "no_target_leak_dag",
            "claim": claim,
            "forbidden_targets": list(FORBIDDEN_TARGETS),
            "edges_from_forbidden_targets": [],
            "status": "PASS_EMPTY_COMPARISON_DAG",
        },
        "controls/empirical_data_exclusion_manifest.json": {
            "artifact": "empirical_data_exclusion_manifest",
            "excluded_from_source_mode": list(FORBIDDEN_TARGETS),
            "comparison_allowed_after_freeze": True,
        },
        "controls/frozen_code_hashes.json": {
            "artifact": "frozen_code_hashes",
            "builder": Path(__file__).relative_to(ROOT).as_posix(),
            "builder_sha256": sha256_text(Path(__file__).read_text(encoding="utf-8")),
        },
        "controls/replay_receipts.json": {
            **base_receipt("replay_receipts", tier=tier),
            "status": "SCHEMA_ONLY_REPLAY_NOT_RUN",
        },
        "controls/comparison_data_manifest.json": {
            "artifact": "comparison_data_manifest",
            "comparison_data": [],
            "status": "NO_COMPARISON_DATA_ATTACHED",
        },
        "claim.md": claim + "\n",
    }
    return payloads


def write_payload(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
        return
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_bundle(out_dir: Path, *, claim: str, tier: str, lambda_path: Path) -> dict[str, Any]:
    if claim not in CLAIMS:
        raise ValueError(f"unknown claim: {claim}")
    if tier not in CLAIM_TIERS:
        raise ValueError(f"unknown claim tier: {tier}")

    lambda_payload = read_optional_json(lambda_path)
    payloads = build_payloads(claim=claim, tier=tier, lambda_payload=lambda_payload)
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        target = out_dir / rel_path
        write_payload(target, payload)
        file_hashes[rel_path] = "sha256:" + hashlib.sha256(target.read_bytes()).hexdigest()

    manifest = {
        "artifact": "hadron_source_backend_manifest",
        "generated_utc": now_utc(),
        "milestone": "HVP_ALPHA_SOURCE_PROTOTYPE",
        "claim": claim,
        "claim_tier": tier,
        "promotion_allowed": claim == "SOURCE_INTERVAL_PROMOTED",
        "source_open": claim != "SOURCE_INTERVAL_PROMOTED",
        "required_files": list(REQUIRED_FILES),
        "missing_files": [
            rel for rel in REQUIRED_FILES if rel != "manifest.json" and not (out_dir / rel).is_file()
        ],
        "forbidden_targets": list(FORBIDDEN_TARGETS),
        "file_hashes": file_hashes,
    }
    write_payload(out_dir / "manifest.json", manifest)
    manifest["file_hashes"]["manifest.json"] = "sha256:" + hashlib.sha256(
        (out_dir / "manifest.json").read_bytes()
    ).hexdigest()
    write_payload(out_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build hadronic precision backend receipt scaffold.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--claim", choices=sorted(CLAIMS), default="SOURCE_PROTOTYPE_NOT_PROMOTED")
    parser.add_argument("--tier", choices=CLAIM_TIERS, default="H2")
    parser.add_argument("--lambda-msbar", default=str(DEFAULT_LAMBDA))
    args = parser.parse_args()

    manifest = build_bundle(
        Path(args.output),
        claim=args.claim,
        tier=args.tier,
        lambda_path=Path(args.lambda_msbar),
    )
    print(Path(args.output) / "manifest.json")
    print(f"claim={manifest['claim']} tier={manifest['claim_tier']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
