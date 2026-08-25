#!/usr/bin/env python3
"""Build the JWST compact-object source-release receipt scaffold.

This is a paper-stack mirror for the simulator workbench in oph-physics-sim.
It is deliberately fail-closed: it can emit receipt files, source-leak guards,
and the strongest allowed claim tier, but it does not turn JWST compact-object
catalog rows into an OPH confirmation or a physical mass/age prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "jwst" / "compact_object_source_release"

CLAIM_TIERS = (
    "J0_DIAGNOSTIC_PROXY",
    "J1_CATALOG_RECORD",
    "J2_SPECTROSCOPIC_OR_PHOTOMETRIC_OBJECT",
    "J3_DEGENERACY_AUDITED_OBJECT",
    "J4_CONDITIONAL_PHYSICAL_OBJECT",
    "J5_SOURCE_RELEASE_CANDIDATE",
    "J6_SOURCE_ONLY_OBJECT_ABUNDANCE",
    "J7_FORWARD_MOCK_PHYSICAL_SPECTRUM",
    "J8_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION",
)

CLAIM_REQUIREMENTS = {
    "J1_CATALOG_RECORD": ("CATALOG_INGESTION_RECEIPT", "CATALOG_PROVENANCE_RECEIPT"),
    "J2_SPECTROSCOPIC_OR_PHOTOMETRIC_OBJECT": (
        "OBJECT_ID_RECEIPT",
        "APERTURE_RECEIPT",
        "REDSHIFT_POSTERIOR_RECEIPT",
        "PHOTOMETRY_OR_SPECTRUM_RECEIPT",
        "MORPHOLOGY_RECORD_RECEIPT",
    ),
    "J3_DEGENERACY_AUDITED_OBJECT": (
        "REDSHIFT_RECEIPT",
        "DUST_RECEIPT",
        "AGN_RECEIPT",
        "NEBULAR_RECEIPT",
        "STELLAR_POPULATION_RECEIPT",
        "LENSING_RECEIPT",
        "MORPHOLOGY_PSF_RECEIPT",
        "SELECTION_RECEIPT",
        "DEGENERACY_AUDIT_RECEIPT",
    ),
    "J4_CONDITIONAL_PHYSICAL_OBJECT": (
        "OBJECT_PARENT_RECEIPT",
        "PACKET_MASS_SHELL_RECEIPT",
        "FINITE_PACKET_STRESS_READOUT_RECEIPT",
        "TOTAL_STRESS_CLOSURE_RECEIPT",
        "RADIATIVE_TRANSFER_RECEIPT",
        "LENSING_SOURCE_PLANE_RECEIPT",
        "CHEMICAL_SFH_RECEIPT",
    ),
    "J5_SOURCE_RELEASE_CANDIDATE": (
        "OBJECT_RELEASE_STATE_RECEIPT",
        "SOURCE_RELEASE_RESIDUAL_RECEIPT",
    ),
    "J6_SOURCE_ONLY_OBJECT_ABUNDANCE": (
        "OBJECT_QUOTIENT_ENSEMBLE_RECEIPT",
        "OBJECT_SOURCE_LAW_RECEIPT",
        "OBJECT_ABUNDANCE_SOURCE_RECEIPT",
        "NO_TARGET_LEAKAGE_RECEIPT",
        "LOAD_REFINEMENT_COMPATIBILITY_RECEIPT",
    ),
    "J7_FORWARD_MOCK_PHYSICAL_SPECTRUM": (
        "JWST_FORWARD_OPERATOR_RECEIPT",
        "JWST_SELECTION_RECEIPT",
        "FROZEN_FORWARD_OPERATOR_RECEIPT",
    ),
    "J8_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION": (
        "FROZEN_SOURCE_RECEIPT",
        "FROZEN_CATALOG_DATA_RECEIPT",
        "FROZEN_CATALOG_LIKELIHOOD_RECEIPT",
        "JWST_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION_RECEIPT",
    ),
}

FORBIDDEN_SOURCE_TOKENS = (
    "jwst_catalog",
    "catalog_count",
    "catalog_counts",
    "observed_count",
    "observed_counts",
    "anomaly_label",
    "interesting_object",
    "posterior_summary",
    "likelihood_residual",
    "model_comparison_score",
    "mass_age_tension_label",
)

NONCLAIMS = (
    "compact record surface is not physical assembly",
    "red color is not old age",
    "luminosity is not stellar mass",
    "broad lines are not automatically large black-hole mass",
    "source release is not galaxy formation",
    "record density is not stellar-mass density",
    "JWST catalog rows are not OPH confirmations",
)

REQUIRED_FILES = (
    "manifest.json",
    "source_artifact.json",
    "quotient_ensemble.json",
    "source_law.json",
    "release_state.json",
    "compact_record_surface.json",
    "object_parent.json",
    "jwst_forward_operator.json",
    "degeneracy_audit.json",
    "abundance_selector.json",
    "frozen_likelihood.json",
    "no_target_leak_dag.json",
    "claim_ladder.json",
    "claim.md",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def target_leak_hits(config: Path | None) -> list[str]:
    if config is None or not config.is_file():
        return []
    text = config.read_text(encoding="utf-8")
    haystack = text.lower()
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        parsed = None
    if parsed is not None:
        haystack = json.dumps(parsed, sort_keys=True).lower()
    return sorted(token for token in FORBIDDEN_SOURCE_TOKENS if token in haystack)


def path_meta(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"path": None, "exists": False}
    meta: dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if path.is_file():
        raw = path.read_bytes()
        meta.update({"sha256": sha256_bytes(raw), "byte_count": len(raw)})
    return meta


def base_receipt(artifact: str, receipts: dict[str, bool], *, status: str) -> dict[str, Any]:
    missing = [name for name, passed in receipts.items() if not passed]
    return {
        "artifact": artifact,
        "generated_utc": now_utc(),
        "status": status,
        "readiness_gates": receipts,
        "blocking_receipts": missing,
        "nonclaims": list(NONCLAIMS),
        "claim_boundary": (
            "Fail-closed JWST compact-object source-release paper-stack mirror. "
            "This artifact is a receipt surface, not a physical promotion."
        ),
    }


def strongest_allowed_claim(receipts: dict[str, bool]) -> tuple[str, str | None, list[str]]:
    strongest = "J0_DIAGNOSTIC_PROXY"
    for tier in CLAIM_TIERS[1:]:
        missing = [name for name in CLAIM_REQUIREMENTS[tier] if not receipts.get(name, False)]
        if missing:
            return strongest, missing[0], missing
        strongest = tier
    return strongest, None, []


def build_payloads(*, config: Path | None) -> dict[str, str | dict[str, Any]]:
    leak_hits = target_leak_hits(config)
    receipts: dict[str, bool] = {
        "NO_TARGET_LEAKAGE_RECEIPT": not leak_hits,
        "OBJECT_QUOTIENT_ENSEMBLE_RECEIPT": False,
        "OBJECT_SOURCE_LAW_RECEIPT": False,
        "OBJECT_RELEASE_STATE_RECEIPT": False,
        "COMPACT_RECORD_SURFACE_RECEIPT": False,
        "SYNC_RESIDUAL_RECEIPT": False,
        "RECORD_DENSITY_RECEIPT": False,
        "OBJECT_PARENT_RECEIPT": False,
        "PACKET_MASS_SHELL_RECEIPT": False,
        "FINITE_PACKET_STRESS_READOUT_RECEIPT": False,
        "TOTAL_STRESS_CLOSURE_RECEIPT": False,
        "RADIATIVE_TRANSFER_RECEIPT": False,
        "LENSING_SOURCE_PLANE_RECEIPT": False,
        "CHEMICAL_SFH_RECEIPT": False,
        "DEGENERACY_AUDIT_RECEIPT": False,
        "MASS_AGE_TENSION_PROMOTION_RECEIPT": False,
        "OBJECT_ABUNDANCE_SOURCE_RECEIPT": False,
        "LOAD_REFINEMENT_COMPATIBILITY_RECEIPT": False,
        "JWST_FORWARD_OPERATOR_RECEIPT": False,
        "JWST_SELECTION_RECEIPT": False,
        "FROZEN_FORWARD_OPERATOR_RECEIPT": False,
        "FROZEN_SOURCE_RECEIPT": False,
        "FROZEN_CATALOG_DATA_RECEIPT": False,
        "FROZEN_CATALOG_LIKELIHOOD_RECEIPT": False,
        "JWST_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION_RECEIPT": False,
    }
    claim, first_blocked, missing = strongest_allowed_claim(receipts)
    source_status = "TARGET_LEAK_INVALID" if leak_hits else "SOURCE_OPEN_NOT_PROMOTED"

    payloads: dict[str, str | dict[str, Any]] = {
        "source_artifact.json": {
            **base_receipt("jwst_compact_object_source_artifact", receipts, status=source_status),
            "config": path_meta(config),
            "target_leak_hits": leak_hits,
            "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
        },
        "quotient_ensemble.json": {
            **base_receipt(
                "jwst_object_quotient_ensemble",
                {"OBJECT_QUOTIENT_ENSEMBLE_RECEIPT": False},
                status="QUOTIENT_LAW_REQUIRED",
            ),
            "presentation_space": "Sigma_obj_r",
            "quotient": "Q_obj_r = Sigma_obj_r/Gamma_obj_r",
            "quotiented_labels": [
                "gauge_representatives",
                "mesh_labels",
                "port_labels",
                "packet_labels",
                "basis_choices",
                "worker_ids",
                "queue_order",
                "retry_counters",
                "inert_ancillas",
            ],
        },
        "source_law.json": {
            **base_receipt(
                "jwst_object_source_law",
                {"OBJECT_SOURCE_LAW_RECEIPT": False, "NO_TARGET_LEAKAGE_RECEIPT": not leak_hits},
                status=source_status,
            ),
            "required_law": "mu_obj_r(q)=Z_obj_r^-1 m_obj_r(q) exp[-S_obj_r(q)]",
            "hard_rule": "NORMAL_FORM_MAP != SOURCE_LAW",
        },
        "release_state.json": {
            **base_receipt(
                "jwst_object_release_state",
                {"OBJECT_RELEASE_STATE_RECEIPT": False, "SOURCE_RELEASE_RESIDUAL_RECEIPT": False},
                status="RELEASE_STATE_REQUIRED",
            ),
            "nonclaim": "source release is record availability, not galaxy formation",
        },
        "compact_record_surface.json": {
            **base_receipt(
                "jwst_compact_record_surface",
                {
                    "COMPACT_RECORD_SURFACE_RECEIPT": False,
                    "SYNC_RESIDUAL_RECEIPT": False,
                    "RECORD_DENSITY_RECEIPT": False,
                },
                status="RECORD_SURFACE_REQUIRED",
            ),
            "nonclaim": "compactness/readout indices do not imply assembly, age, or mass",
        },
        "object_parent.json": {
            **base_receipt(
                "jwst_finite_object_parent",
                {
                    "OBJECT_PARENT_RECEIPT": False,
                    "PACKET_MASS_SHELL_RECEIPT": False,
                    "FINITE_PACKET_STRESS_READOUT_RECEIPT": False,
                    "TOTAL_STRESS_CLOSURE_RECEIPT": False,
                    "RADIATIVE_TRANSFER_RECEIPT": False,
                    "LENSING_SOURCE_PLANE_RECEIPT": False,
                    "CHEMICAL_SFH_RECEIPT": False,
                },
                status="FINITE_PARENT_REQUIRED",
            ),
            "nonclaim": "scalar rows cannot be promoted to stress, gas mass, dynamical mass, or maturity",
        },
        "jwst_forward_operator.json": {
            **base_receipt(
                "jwst_forward_operator",
                {
                    "JWST_FORWARD_OPERATOR_RECEIPT": False,
                    "JWST_SELECTION_RECEIPT": False,
                    "FROZEN_FORWARD_OPERATOR_RECEIPT": False,
                },
                status="FORWARD_OPERATOR_REQUIRED",
            ),
            "freeze_rule": "operator and thresholds must be frozen before catalog comparison",
        },
        "degeneracy_audit.json": {
            **base_receipt(
                "jwst_degeneracy_audit",
                {
                    "REDSHIFT_RECEIPT": False,
                    "DUST_RECEIPT": False,
                    "AGN_RECEIPT": False,
                    "NEBULAR_RECEIPT": False,
                    "STELLAR_POPULATION_RECEIPT": False,
                    "LENSING_RECEIPT": False,
                    "MORPHOLOGY_PSF_RECEIPT": False,
                    "SELECTION_RECEIPT": False,
                    "DEGENERACY_AUDIT_RECEIPT": False,
                    "MASS_AGE_TENSION_PROMOTION_RECEIPT": False,
                },
                status="DEGENERACY_OPEN",
            ),
            "recommended_open_label": "DEGENERACY_OPEN",
        },
        "abundance_selector.json": {
            **base_receipt(
                "jwst_object_abundance_selector",
                {
                    "OBJECT_ABUNDANCE_SOURCE_RECEIPT": False,
                    "NO_TARGET_LEAKAGE_RECEIPT": not leak_hits,
                    "LOAD_REFINEMENT_COMPATIBILITY_RECEIPT": False,
                },
                status=source_status,
            ),
            "nonclaim": "selection is not physical abundance",
        },
        "frozen_likelihood.json": {
            **base_receipt(
                "jwst_frozen_catalog_likelihood",
                {
                    "FROZEN_SOURCE_RECEIPT": False,
                    "FROZEN_CATALOG_DATA_RECEIPT": False,
                    "FROZEN_CATALOG_LIKELIHOOD_RECEIPT": False,
                    "JWST_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION_RECEIPT": False,
                },
                status="LIKELIHOOD_NOT_EVALUATED",
            ),
            "nonclaim": "catalog comparison without frozen hashes is diagnostic only",
        },
        "no_target_leak_dag.json": {
            "artifact": "jwst_no_target_leak_dag",
            "generated_utc": now_utc(),
            "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
            "target_leak_hits": leak_hits,
            "edges_from_forbidden_targets": [] if not leak_hits else leak_hits,
            "NO_TARGET_LEAKAGE_RECEIPT": not leak_hits,
            "status": "PASS_EMPTY_COMPARISON_DAG" if not leak_hits else "FAIL_FORBIDDEN_SOURCE_INPUT",
        },
        "claim_ladder.json": {
            "artifact": "jwst_compact_object_claim_ladder",
            "generated_utc": now_utc(),
            "claim_tiers": list(CLAIM_TIERS),
            "claim_requirements": {key: list(value) for key, value in CLAIM_REQUIREMENTS.items()},
            "strongest_allowed_claim": claim,
            "first_blocked_gate": first_blocked,
            "missing_receipts": missing,
            "nonclaims": list(NONCLAIMS),
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


def build_bundle(out_dir: Path, *, config: Path | None = None) -> dict[str, Any]:
    payloads = build_payloads(config=config)
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        target = out_dir / rel_path
        write_payload(target, payload)
        file_hashes[rel_path] = sha256_bytes(target.read_bytes())

    claim_ladder = payloads["claim_ladder.json"]
    assert isinstance(claim_ladder, dict)
    manifest = {
        "artifact": "jwst_compact_object_source_release_manifest",
        "generated_utc": now_utc(),
        "milestone": "Q3_JWST_COMPACT_OBJECT_SOURCE_RELEASE_AUDIT",
        "strongest_allowed_claim": claim_ladder["strongest_allowed_claim"],
        "first_blocked_gate": claim_ladder["first_blocked_gate"],
        "source_open": True,
        "promotion_allowed": False,
        "required_files": list(REQUIRED_FILES),
        "missing_files": [
            rel for rel in REQUIRED_FILES if rel != "manifest.json" and not (out_dir / rel).is_file()
        ],
        "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
        "nonclaims": list(NONCLAIMS),
        "paper_problem_note": "physics-problems/jwst_compact_object_source_release.md",
        "simulator_surface": "oph-physics-sim/oph_fpe/jwst",
        "file_hashes": file_hashes,
    }
    write_payload(out_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build JWST compact-object source-release receipts.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--config", default=None, type=Path)
    args = parser.parse_args()

    manifest = build_bundle(Path(args.output), config=args.config)
    print(Path(args.output) / "manifest.json")
    print(f"strongest_allowed_claim={manifest['strongest_allowed_claim']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
