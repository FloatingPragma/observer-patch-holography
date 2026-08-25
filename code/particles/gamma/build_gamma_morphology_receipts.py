#!/usr/bin/env python3
"""Build the gamma morphology receipt scaffold.

This is a paper-stack mirror for the simulator gamma workbench. It freezes the
objects required for OPH gamma-ray morphology claims, but it does not analyze
Fermi data and does not promote any gamma excess into an OPH detection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "gamma" / "morphology_claims"

CLAIM_TIERS = (
    "DIAGNOSTIC_GAMMA_MAP",
    "SOURCE_DERIVED_GAMMA_TEMPLATE",
    "INSTRUMENT_CONVOLVED_GAMMA_TEMPLATE",
    "IDENTIFIABLE_GAMMA_TEMPLATE",
    "LIKELIHOOD_EVALUATED_GAMMA_MORPHOLOGY",
    "CROSS_TRACER_VALIDATED_GAMMA_MORPHOLOGY",
    "OPH_GAMMA_MORPHOLOGY_CANDIDATE",
)

RECEIPTS = (
    "GAMMA_SOURCE_ARTIFACT_RECEIPT",
    "GAMMA_ROUTE_DECLARATION_RECEIPT",
    "GAMMA_SOURCE_LAW_RECEIPT",
    "GAMMA_STRESS_PARENT_RECEIPT",
    "ANOMALY_SM_CURRENT_NULL_RECEIPT",
    "STRESS_TO_GAMMA_CONTRACTION_RECEIPT",
    "PHOTON_RESPONSE_KERNEL_RECEIPT",
    "BOUNDARY_RECORD_SOURCE_RECEIPT",
    "BOUNDARY_DIPOLE_AXIS_FREEZE_RECEIPT",
    "ASTRO_BRIDGE_FREEZE_RECEIPT",
    "LINE_OF_SIGHT_PROJECTION_RECEIPT",
    "INSTRUMENT_RESPONSE_RECEIPT",
    "SIGNED_TEMPLATE_POSITIVITY_RECEIPT",
    "GAMMA_TEMPLATE_FREEZE_RECEIPT",
    "GAMMA_NO_DATA_USE_RECEIPT",
    "FOREGROUND_ALTERNATIVE_SET_RECEIPT",
    "GAMMA_IDENTIFIABILITY_RECEIPT",
    "FROZEN_GAMMA_LIKELIHOOD_RECEIPT",
    "HELDOUT_GAMMA_VALIDATION_RECEIPT",
    "GAMMA_CROSS_TRACER_RECEIPT",
    "OPH_GAMMA_MORPHOLOGY_PREDICTION_RECEIPT",
)

CLAIM_REQUIREMENTS = {
    "SOURCE_DERIVED_GAMMA_TEMPLATE": (
        "GAMMA_SOURCE_ARTIFACT_RECEIPT",
        "GAMMA_ROUTE_DECLARATION_RECEIPT",
        "GAMMA_SOURCE_LAW_RECEIPT",
        "GAMMA_NO_DATA_USE_RECEIPT",
        "GAMMA_TEMPLATE_FREEZE_RECEIPT",
    ),
    "INSTRUMENT_CONVOLVED_GAMMA_TEMPLATE": (
        "PHOTON_RESPONSE_KERNEL_RECEIPT",
        "LINE_OF_SIGHT_PROJECTION_RECEIPT",
        "INSTRUMENT_RESPONSE_RECEIPT",
        "SIGNED_TEMPLATE_POSITIVITY_RECEIPT",
    ),
    "IDENTIFIABLE_GAMMA_TEMPLATE": (
        "FOREGROUND_ALTERNATIVE_SET_RECEIPT",
        "GAMMA_IDENTIFIABILITY_RECEIPT",
    ),
    "LIKELIHOOD_EVALUATED_GAMMA_MORPHOLOGY": (
        "FROZEN_GAMMA_LIKELIHOOD_RECEIPT",
        "HELDOUT_GAMMA_VALIDATION_RECEIPT",
    ),
    "CROSS_TRACER_VALIDATED_GAMMA_MORPHOLOGY": ("GAMMA_CROSS_TRACER_RECEIPT",),
    "OPH_GAMMA_MORPHOLOGY_CANDIDATE": ("OPH_GAMMA_MORPHOLOGY_PREDICTION_RECEIPT",),
}

FORBIDDEN_SOURCE_TOKENS = (
    "gamma_residual_map",
    "gamma_residual_maps",
    "lat_count_residual",
    "lat_count_residuals",
    "likelihood_value",
    "likelihood_values",
    "posterior_summary",
    "foreground_nuisance_fit",
    "foreground_nuisance_fits",
    "template_covariance_from_target_data",
    "human_picked_axis_after_gamma_inspection",
    "gamma_derived_axis",
    "residual_amplitude",
)

REQUIRED_FILES = (
    "manifest.json",
    "gamma_source_artifact.json",
    "source_dag.json",
    "stress_parent_adapter.json",
    "boundary_record_adapter.json",
    "astro_bridge.json",
    "stress_contraction.json",
    "photon_response.json",
    "los_projection.json",
    "instrument_response.json",
    "template_compiler.json",
    "foreground_registry.json",
    "identifiability.json",
    "likelihood.json",
    "cross_tracer.json",
    "nulls.json",
    "receipts.json",
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


def strongest_allowed_claim(receipts: dict[str, bool]) -> tuple[str, str | None, list[str]]:
    strongest = "DIAGNOSTIC_GAMMA_MAP"
    for tier in CLAIM_TIERS[1:]:
        required: list[str] = []
        for claim in CLAIM_TIERS[1 : CLAIM_TIERS.index(tier) + 1]:
            required.extend(CLAIM_REQUIREMENTS[claim])
        missing = [name for name in required if not receipts.get(name, False)]
        if missing:
            return strongest, missing[0], missing
        strongest = tier
    return strongest, None, []


def build_receipts(*, route: str, config: Path | None, direct_anomaly_gamma: bool) -> dict[str, bool]:
    leak_hits = target_leak_hits(config)
    receipts = {name: False for name in RECEIPTS}
    receipts["GAMMA_ROUTE_DECLARATION_RECEIPT"] = route in {"TRANSPORTED_STRESS", "BOUNDARY_RECORD_DIPOLE", "BOTH"}
    receipts["GAMMA_NO_DATA_USE_RECEIPT"] = not leak_hits
    receipts["ANOMALY_SM_CURRENT_NULL_RECEIPT"] = not direct_anomaly_gamma
    receipts["FOREGROUND_ALTERNATIVE_SET_RECEIPT"] = True
    return receipts


def base_payload(name: str, *, route: str, claim: str, receipts: dict[str, bool]) -> dict[str, Any]:
    return {
        "artifact": name,
        "route": route,
        "strongest_allowed_claim": claim,
        "promotion_allowed": claim == "OPH_GAMMA_MORPHOLOGY_CANDIDATE",
        "readiness_gates": receipts,
    }


def build_payloads(*, route: str, config: Path | None, direct_anomaly_gamma: bool) -> dict[str, str | dict[str, Any]]:
    receipts = build_receipts(route=route, config=config, direct_anomaly_gamma=direct_anomaly_gamma)
    claim, first_blocked, missing = strongest_allowed_claim(receipts)
    leak_hits = target_leak_hits(config)
    base = {"route": route, "claim": claim}
    return {
        "gamma_source_artifact.json": {
            **base_payload("gamma_source_artifact", route=route, claim=claim, receipts=receipts),
            "required_object": "mathfrak_G_gamma_r",
            "status": "SOURCE_ARTIFACT_REQUIRED",
        },
        "source_dag.json": {
            **base,
            "artifact": "gamma_source_dag",
            "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
            "target_leak_hits": leak_hits,
            "status": "PASS_EMPTY_COMPARISON_DAG" if not leak_hits else "FAIL_FORBIDDEN_SOURCE_INPUT",
        },
        "stress_parent_adapter.json": {
            **base,
            "artifact": "stress_parent_adapter",
            "scalar_only_status": "SCALAR_ONLY_GAMMA_DIAGNOSTIC",
            "required_receipts": [
                "FINITE_COVARIANT_COLLAR_PACKET_PARENT_RECEIPT",
                "PACKET_MASS_SHELL_RECEIPT",
                "TOTAL_STRESS_CLOSURE_RECEIPT",
                "SM_CURRENT_NULL_RECEIPT",
                "RETARDED_RESPONSE_RECEIPT",
                "RESPONSE_STABILITY_RECEIPT",
                "REFINEMENT_CONVERGENCE_RECEIPT",
                "CDM_LIMIT_RECOVERY_RECEIPT",
            ],
        },
        "boundary_record_adapter.json": {**base, "artifact": "boundary_record_adapter", "status": "BOUNDARY_RECORD_SOURCE_REQUIRED"},
        "astro_bridge.json": {**base, "artifact": "astro_bridge", "status": "ASTRO_BRIDGE_FREEZE_REQUIRED"},
        "stress_contraction.json": {
            **base,
            "artifact": "stress_contraction",
            "default_contraction": "sigma_T = T_A^{mu nu} u_mu u_nu",
            "extended_branch_requires": "EXTENDED_STRESS_CONTRACTION_RECEIPT",
        },
        "photon_response.json": {
            **base,
            "artifact": "photon_response",
            "direct_anomaly_gamma_default": 0.0,
            "nonzero_direct_requires": "ANOMALY_EM_CURRENT_RECEIPT",
        },
        "los_projection.json": {**base, "artifact": "line_of_sight_projection", "status": "PROJECTION_REQUIRED"},
        "instrument_response.json": {**base, "artifact": "instrument_response", "operator": "intensity_to_binned_counts"},
        "template_compiler.json": {**base, "artifact": "template_compiler", "status": "TEMPLATE_FREEZE_REQUIRED"},
        "foreground_registry.json": {
            **base,
            "artifact": "foreground_registry",
            "minimum": [
                "isotropic",
                "galactic_diffuse",
                "point_sources",
                "gas_rings",
                "dust",
                "inverse_compton",
                "fermi_bubbles",
                "bulge",
                "disk_msp",
                "annihilation_rho_squared",
                "decay_rho",
            ],
        },
        "identifiability.json": {**base, "artifact": "identifiability", "metric": "eta_id"},
        "likelihood.json": {**base, "artifact": "frozen_gamma_likelihood", "likelihood": "binned_poisson_counts"},
        "cross_tracer.json": {**base, "artifact": "cross_tracer", "status": "ROUTE_SPECIFIC_VALIDATION_REQUIRED"},
        "nulls.json": {**base, "artifact": "gamma_nulls", "status": "SYSTEMATICS_NULLS_REQUIRED"},
        "receipts.json": {**base, "artifact": "gamma_receipts", "readiness_gates": receipts},
        "claim_ladder.json": {
            **base,
            "artifact": "gamma_claim_ladder",
            "claim_ladder": list(CLAIM_TIERS),
            "first_blocked_gate": first_blocked,
            "missing_for_next_claim": missing,
        },
        "claim.md": claim + "\n",
    }


def write_payload(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_bundle(out_dir: Path, *, route: str, config: Path | None, direct_anomaly_gamma: bool) -> dict[str, Any]:
    if route not in {"TRANSPORTED_STRESS", "BOUNDARY_RECORD_DIPOLE", "BOTH"}:
        raise ValueError(f"unknown route: {route}")
    payloads = build_payloads(route=route, config=config, direct_anomaly_gamma=direct_anomaly_gamma)
    out_dir.mkdir(parents=True, exist_ok=True)
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        target = out_dir / rel_path
        write_payload(target, payload)
        file_hashes[rel_path] = sha256_bytes(target.read_bytes())

    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    ladder = json.loads((out_dir / "claim_ladder.json").read_text(encoding="utf-8"))
    manifest = {
        "artifact": "gamma_morphology_manifest",
        "generated_utc": now_utc(),
        "milestone": "GAMMA_MORPHOLOGY_AUDIT",
        "route": route,
        "strongest_allowed_claim": claim,
        "first_blocked_gate": ladder["first_blocked_gate"],
        "promotion_allowed": claim == "OPH_GAMMA_MORPHOLOGY_CANDIDATE",
        "required_files": list(REQUIRED_FILES),
        "missing_files": [rel for rel in REQUIRED_FILES if rel != "manifest.json" and not (out_dir / rel).is_file()],
        "file_hashes": file_hashes,
    }
    write_payload(out_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build gamma morphology receipt scaffold.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--route", choices=["TRANSPORTED_STRESS", "BOUNDARY_RECORD_DIPOLE", "BOTH"], default="TRANSPORTED_STRESS")
    parser.add_argument("--config")
    parser.add_argument("--direct-anomaly-gamma", action="store_true")
    args = parser.parse_args()

    manifest = build_bundle(
        Path(args.output),
        route=args.route,
        config=Path(args.config) if args.config else None,
        direct_anomaly_gamma=bool(args.direct_anomaly_gamma),
    )
    print(Path(args.output) / "manifest.json")
    print(f"claim={manifest['strongest_allowed_claim']} route={manifest['route']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
