#!/usr/bin/env python3
"""Build the high-energy messenger coefficient-emission algebra scaffold.

This is a paper-stack mirror for the simulator coefficient-emitter.  The
default values are planted examples: no source artifact is loaded, no moment
interiority or feature-minimality calculation is performed, and no coefficient
solve is run.  Therefore this script cannot emit a source receipt or authorize
promotion.  It also checks that an optional configuration does not smuggle UHE
target data into the source side.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "uhe" / "coefficient_emission"

REQUIRED_RECEIPTS = (
    "BASELINE_FULL_SUPPORT",
    "FEATURE_MINIMALITY",
    "MOMENT_INTERIOR",
    "SOURCE_LOAD_QUOTIENT_VISIBLE",
    "NO_UHE_DATA_USE",
    "REFINEMENT_COMPATIBILITY",
    "COEFFICIENT_SOLVE_CONVERGED",
    "COMMON_SOURCE_LOCK",
)

FORBIDDEN_SOURCE_TOKENS = (
    "uhe_event_coordinates",
    "event_coordinates",
    "arrival_direction",
    "arrival_directions",
    "event_energy",
    "event_energies",
    "association_failure",
    "association_failures",
    "catalog_match_after_events",
    "post_event_catalog_match",
    "likelihood_value",
    "likelihood_values",
    "posterior_summary",
    "diagnostic_overlay",
    "diagnostic_overlays",
    "residual_map",
    "residual_maps",
    "human_selected_event_pattern",
    "human_picked_event_pattern",
)

REQUIRED_FILES = (
    "manifest.json",
    "source_release_quotient.json",
    "source_law.json",
    "source_loads.json",
    "baseline_measure.json",
    "feature_map.json",
    "moment_targets.json",
    "coefficient_solver.json",
    "emitted_coefficients.json",
    "source_dag.json",
    "claim_ladder.json",
    "claim.md",
)

NONCLAIMS = (
    "source-only coefficient emission is not a detected UHE source",
    "a finite MaxEnt solve is not an event-map fit",
    "neutrino, cosmic-ray, and gamma maps must share the source coefficient",
    "propagation and detector kernels are downstream likelihood plumbing",
    "target-data leakage invalidates a source-only label",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def stable_hash(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return sha256_bytes(raw)


def logit(p: float) -> float:
    if not 0.0 < p < 1.0:
        raise ValueError("probability must be strictly between 0 and 1")
    return math.log(p / (1.0 - p))


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


def example_coefficients() -> dict[str, Any]:
    binary_eta = logit(0.70) - logit(0.40)
    poisson_eta = math.log(6.0 / 2.0)
    finite_eta = [0.25, 0.35, -0.20, 0.10, 0.05]
    return {
        "finite_maxent_eta": finite_eta,
        "binary_logit_example": {
            "p_OPH": 0.70,
            "p0": 0.40,
            "eta": binary_eta,
            "formula": "logit(p_OPH)-logit(p0)",
        },
        "poisson_opportunity_example": {
            "lambda_OPH": 6.0,
            "lambda0": 2.0,
            "eta": poisson_eta,
            "formula": "log(lambda_OPH/lambda0)",
        },
    }


def build_payloads(*, config: Path | None) -> dict[str, str | dict[str, Any]]:
    leak_hits = target_leak_hits(config)
    receipts = {
        "BASELINE_FULL_SUPPORT": False,
        "FEATURE_MINIMALITY": False,
        "MOMENT_INTERIOR": False,
        "SOURCE_LOAD_QUOTIENT_VISIBLE": False,
        "NO_UHE_DATA_USE": not leak_hits,
        "REFINEMENT_COMPATIBILITY": False,
        "COEFFICIENT_SOLVE_CONVERGED": False,
        "COMMON_SOURCE_LOCK": False,
    }
    claim = (
        "SYNTHETIC_COEFFICIENT_SCAFFOLD_ONLY"
        if receipts["NO_UHE_DATA_USE"]
        else "INVALIDATED_COEFFICIENT_DAG"
    )
    coeffs = example_coefficients()
    base = {
        "generated_utc": now_utc(),
        "artifact_type": "UHE_COEFFICIENT_EMISSION_SCAFFOLD",
        "claim": claim,
        "physical_claim": False,
        "synthetic_fixture": True,
        "source_artifact_bound": False,
        "promotion_allowed": False,
        "readiness_gates": receipts,
        "readiness_gate_status": "UNVERIFIED_EXCEPT_NO_UHE_DATA_USE",
        "nonclaims": list(NONCLAIMS),
    }
    source_dag = {
        **base,
        "artifact": "uhe_source_dag",
        "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
        "target_leak_hits": leak_hits,
        "status": (
            "PASS_NO_FORBIDDEN_INPUTS_SCAFFOLD_ONLY"
            if not leak_hits
            else "FAIL_FORBIDDEN_SOURCE_INPUT"
        ),
    }
    feature_map = {
        "feature_names": [
            "Ahat_r",
            "C_compact",
            "H_hidden",
            "Ahat_r*C_compact",
            "Ahat_r*H_hidden",
        ],
        "minimality_receipt": receipts["FEATURE_MINIMALITY"],
        "minimality_check_run": False,
        "hash": stable_hash(["Ahat_r", "C_compact", "H_hidden", "Ahat_r*C_compact", "Ahat_r*H_hidden"]),
    }
    return {
        "source_release_quotient.json": {
            **base,
            "artifact": "uhe_source_release_quotient",
            "quotient": "Q_r^rel = Sigma_r^rel/Gamma_r^rel",
            "quotiented_labels": [
                "hidden_representatives",
                "port_labels",
                "repair_schedules",
                "worker_metadata",
                "mesh_labels",
                "inert_coordinates",
            ],
        },
        "source_law.json": {
            **base,
            "artifact": "uhe_source_law",
            "law": "mu_r^rel(q)=Z_r^-1 m_r(q) exp[-S_r(q)]",
            "source_only": False,
            "source_law_instantiated": False,
        },
        "source_loads.json": {
            **base,
            "artifact": "compact_engine_source_loads",
            "observable": "L_alpha,r^CE(q;g)=Pi_g,r sum_C V_C,r^phys omega_alpha,r S_nu,r,C",
            "quotient_visible": False,
            "source_loads_supplied": False,
            "forbidden": "event coordinates, event energies, residual maps, likelihood values",
        },
        "baseline_measure.json": {
            **base,
            "artifact": "uhe_baseline_measure",
            "full_support": False,
            "baseline_supplied": False,
            "fixed_before_uhe_comparison": False,
        },
        "feature_map.json": {**base, "artifact": "uhe_feature_map", **feature_map},
        "moment_targets.json": {
            **base,
            "artifact": "uhe_source_moment_targets",
            "definition": "c_alpha,r(g)=N_alpha,r E_mu[L_alpha,r^CE(Q;g)]",
            "target_data_inputs_allowed": False,
            "source_moments_supplied": False,
        },
        "coefficient_solver.json": {
            **base,
            "artifact": "uhe_coefficient_solver",
            "solver": "not run; deterministic damped Newton is the intended production solver",
            "converged": False,
            "solver_status": "PLANTED_EXAMPLE_VALUES_ONLY",
            "receipts_checked": list(REQUIRED_RECEIPTS),
        },
        "emitted_coefficients.json": {
            **base,
            "artifact": "uhe_emitted_coefficients",
            "coefficients": coeffs,
            "coefficient_status": "PLANTED_EXAMPLES_NOT_SOURCE_DERIVED",
            "common_source_lock": False,
            "species_downstream_only": ["neutrino", "cosmic_ray", "gamma"],
        },
        "source_dag.json": source_dag,
        "claim_ladder.json": {
            **base,
            "artifact": "uhe_claim_ladder",
            "claim_tiers": [
                "SYNTHETIC_COEFFICIENT_SCAFFOLD_ONLY",
                "SOURCE_ONLY_COEFFICIENT_EMITTED",
                "CONDITIONAL_SOURCE_MODEL",
                "FITTED_OPH_COEFFICIENT",
                "INVALIDATED_COEFFICIENT_DAG",
            ],
            "required_receipts": list(REQUIRED_RECEIPTS),
            "first_blocked_gate": (
                "SOURCE_ARTIFACT_AND_COMPUTED_READINESS_RECEIPTS"
                if receipts["NO_UHE_DATA_USE"]
                else "NO_UHE_DATA_USE"
            ),
        },
        "claim.md": claim + "\n",
    }


def write_payload(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def build_bundle(out_dir: Path, *, config: Path | None = None) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    payloads = build_payloads(config=config)
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        path = out_dir / rel_path
        write_payload(path, payload)
        file_hashes[rel_path] = sha256_bytes(path.read_bytes())
    missing = [
        rel_path
        for rel_path in REQUIRED_FILES
        if rel_path != "manifest.json" and not (out_dir / rel_path).is_file()
    ]
    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    ladder = json.loads((out_dir / "claim_ladder.json").read_text(encoding="utf-8"))
    source_dag = json.loads((out_dir / "source_dag.json").read_text(encoding="utf-8"))
    manifest = {
        "artifact": "uhe_coefficient_emission_manifest",
        "generated_utc": now_utc(),
        "milestone": "UHE_COEFFICIENT_EMISSION_AUDIT",
        "artifact_type": "UHE_COEFFICIENT_EMISSION_SCAFFOLD",
        "strongest_allowed_claim": claim,
        "first_blocked_gate": ladder.get("first_blocked_gate"),
        "promotion_allowed": False,
        "physical_claim": False,
        "synthetic_fixture": True,
        "source_artifact_bound": False,
        "required_files": list(REQUIRED_FILES),
        "missing_files": missing,
        "target_leak_hits": source_dag.get("target_leak_hits", []),
        "file_hashes": file_hashes,
    }
    write_payload(out_dir / "manifest.json", manifest)
    # A file cannot contain its own byte hash without a separately specified
    # canonicalization/fixed-point convention.  Hash only the other bundle
    # members; the old write-then-rewrite sequence stored a necessarily stale
    # manifest self-hash.
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=DEFAULT_OUT, type=Path)
    parser.add_argument("--config", default=None, type=Path)
    args = parser.parse_args(argv)
    manifest = build_bundle(args.output, config=args.config)
    print(json.dumps(manifest, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
