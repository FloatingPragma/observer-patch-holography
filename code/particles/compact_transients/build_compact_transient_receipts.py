#!/usr/bin/env python3
"""Build the compact-transient receipt scaffold.

This is a paper-stack mirror for the compact-transient simulator workbench. It
freezes the CR0-CR4 claim ladder and receipt files for FRBs, old-host compact
sources, and black-hole recycling, but it does not analyze transient catalogs
or promote any event class into an OPH confirmation.

NON_CLAIMS:
The checker verifies presence, byte-level hash integrity, schema conformance,
declared provenance, freshness against recorded inputs, and non-self-authorship.
It does not verify that any value is physically meaningful, correctly derived,
or scientifically true. Requiring an external signature or a reference to an
independently retrievable dataset is an obvious hardening direction that would
make forgery cost more than doing the science, but this checker does not do
either.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "compact_transients" / "receipt_scaffold"
BUILDER_ID = "compact-transient-receipt-builder/v1"
EVIDENCE_SCHEMA = "compact-transient-evidence-v1"
EVIDENCE_MANIFEST_SCHEMA = "compact-transient-evidence-manifest-v1"
VERDICT_CHECKING_TIER = "SCHEMA_AND_PROVENANCE_ONLY"
SHA256_RE = re.compile(r"sha256:[0-9a-f]{64}\Z")

CLAIM_TIERS = (
    "CR0_VOCABULARY_ONLY",
    "CR1_QUOTIENT_DIAGNOSTIC",
    "CR2_CONDITIONAL_PHENOMENOLOGY",
    "CR3_FROZEN_PHYSICAL_PREDICTION",
    "CR4_SOURCE_ONLY_OPH_PREDICTION",
)

FORBIDDEN_SOURCE_TOKENS = (
    "ringdown_residual",
    "ringdown_residuals",
    "postfit_repair_tail_amplitude",
    "post_fit_repair_tail_amplitude",
    "echo_score",
    "echo_scores",
    "waveform_template_tuned_after_residual_inspection",
)

REQUIRED_FILES = (
    "manifest.json",
    "compact_history.json",
    "compact_quotient.json",
    "source_law.json",
    "repair_emission_kernel.json",
    "packet_parent.json",
    "detector_thinning.json",
    "censoring.json",
    "point_process_likelihood.json",
    "frb_controls.json",
    "bh_recycling.json",
    "refinement_accuracy.json",
    "promotion_audit.json",
    "claim.md",
)

RECEIPT_NAMES = (
    "COMPACT_HISTORY_RECEIPT",
    "COMPACT_QUOTIENT_RECEIPT",
    "COMPACT_SOURCE_LAW_RECEIPT",
    "PACKETIZED_KERNEL_RECEIPT",
    "PHYSICAL_CLOCK_RECEIPT",
    "FINITE_PACKET_PARENT_RECEIPT",
    "PACKET_CONSERVATION_RECEIPT",
    "PROPAGATION_RECEIPT",
    "DETECTION_THINNING_RECEIPT",
    "CENSORING_AND_UPPER_LIMIT_RECEIPT",
    "POINT_PROCESS_LIKELIHOOD_RECEIPT",
    "REPEATER_HISTORY_LIKELIHOOD_RECEIPT",
    "FRB_SOURCE_IDENTITY_RECEIPT",
    "FRB_CADENCE_EXPOSURE_RECEIPT",
    "BH_GENEALOGY_DAG_RECEIPT",
    "NO_GENERATION_LEAKAGE_RECEIPT",
    "CONTROL_MODEL_RECEIPT",
    "REFINEMENT_STABILITY_RECEIPT",
    "SIMULATOR_ACCURACY_RECEIPT",
    "FROZEN_HASHES_RECEIPT",
    "HELDOUT_LIKELIHOOD_RECEIPT",
    "PROMOTION_AUDIT_RECEIPT",
    "COMPACT_SOURCE_ACTION_DERIVED_RECEIPT",
    "EMISSION_MICROPHYSICS_DERIVED_RECEIPT",
    "PHYSICAL_CLOCK_DERIVED_RECEIPT",
    "OLD_HOST_FRB_SOURCE_THEOREM_RECEIPT",
    "BH_GENEALOGY_PRIOR_THEOREM_RECEIPT",
)

EVIDENCE_KINDS = {
    "COMPACT_QUOTIENT_RECEIPT": "compact_quotient_validation",
    "COMPACT_SOURCE_LAW_RECEIPT": "normalized_channel_source_law",
    "PACKETIZED_KERNEL_RECEIPT": "packetized_repair_kernel",
    "PHYSICAL_CLOCK_RECEIPT": "physical_clock_calibration",
    "FINITE_PACKET_PARENT_RECEIPT": "finite_packet_parent",
    "PROPAGATION_RECEIPT": "normalized_propagation_kernel",
    "DETECTION_THINNING_RECEIPT": "fitted_detector_thinning",
    "CENSORING_AND_UPPER_LIMIT_RECEIPT": "censoring_upper_limit_model",
    "POINT_PROCESS_LIKELIHOOD_RECEIPT": "point_process_likelihood_evaluation",
    "HELDOUT_LIKELIHOOD_RECEIPT": "heldout_likelihood_evaluation",
    "CONTROL_MODEL_RECEIPT": "frozen_control_model_comparison",
    "REFINEMENT_STABILITY_RECEIPT": "refinement_stability_evaluation",
    "FROZEN_HASHES_RECEIPT": "frozen_artifact_set",
}

CR2_RECEIPTS = (
    "COMPACT_QUOTIENT_RECEIPT",
    "COMPACT_SOURCE_LAW_RECEIPT",
    "PACKETIZED_KERNEL_RECEIPT",
    "PHYSICAL_CLOCK_RECEIPT",
    "FINITE_PACKET_PARENT_RECEIPT",
    "PROPAGATION_RECEIPT",
    "DETECTION_THINNING_RECEIPT",
    "POINT_PROCESS_LIKELIHOOD_RECEIPT",
    "CENSORING_AND_UPPER_LIMIT_RECEIPT",
    "HELDOUT_LIKELIHOOD_RECEIPT",
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


def _nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _finite_number(value: object) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def _true_fields(content: dict[str, Any], *names: str) -> bool:
    return all(content.get(name) is True for name in names)


def _valid_common_content(content: object) -> str | None:
    if not isinstance(content, dict) or not content:
        return "content must be a non-empty object"
    assumptions = content.get("assumptions")
    if not isinstance(assumptions, list) or not assumptions or not all(_nonempty_string(x) for x in assumptions):
        return "assumptions must be a non-empty string list"
    for name in ("units", "parameters"):
        value = content.get(name)
        if not isinstance(value, dict) or not value or not all(
            _nonempty_string(key) and _nonempty_string(item) for key, item in value.items()
        ):
            return f"{name} must be a non-empty string mapping"
    outputs = content.get("outputs")
    if not isinstance(outputs, dict) or not outputs or not all(
        _nonempty_string(key) and isinstance(value, str) and SHA256_RE.fullmatch(value)
        for key, value in outputs.items()
    ):
        return "outputs must be a non-empty mapping of names to SHA-256 hashes"
    return None


def _validate_receipt_content(
    receipt: str, content: object, input_records: list[dict[str, str]]
) -> str | None:
    common_error = _valid_common_content(content)
    if common_error is not None:
        return common_error
    assert isinstance(content, dict)

    if receipt == "COMPACT_QUOTIENT_RECEIPT":
        if not all(_nonempty_string(content.get(name)) for name in ("quotient_space", "equivalence_relation")):
            return "quotient space and equivalence relation are required"
        if not _true_fields(content, "canonicalizer_deterministic", "readouts_factor_through_quotient"):
            return "quotient invariance checks did not pass"
    elif receipt == "COMPACT_SOURCE_LAW_RECEIPT":
        if not all(_nonempty_string(content.get(name)) for name in ("channel", "source_law")):
            return "channel and source law are required"
        error = content.get("normalization_error")
        tolerance = content.get("normalization_tolerance")
        if not content.get("normalization_checked") is True or not all(
            _finite_number(value) and value >= 0 for value in (error, tolerance)
        ) or error > tolerance:
            return "source-law normalization check did not pass"
    elif receipt == "PACKETIZED_KERNEL_RECEIPT":
        variables = content.get("joint_variables")
        if not _nonempty_string(content.get("kernel")) or not _true_fields(
            content, "normalization_checked", "coupled_sampler_checked"
        ) or not isinstance(variables, list) or not {"state", "packet", "receipt", "clock"}.issubset(variables):
            return "packetized coupled-kernel checks did not pass"
    elif receipt == "PHYSICAL_CLOCK_RECEIPT":
        if not _nonempty_string(content.get("clock_map")) or not _true_fields(content, "calibration_checked"):
            return "physical-clock calibration is missing"
        if not all(_nonempty_string(content.get(name)) for name in ("input_unit", "output_unit")):
            return "physical-clock units are missing"
    elif receipt == "FINITE_PACKET_PARENT_RECEIPT":
        sectors = content.get("sectors")
        required = {"radio", "gamma", "GW", "optical", "neutrino", "environmental", "recipient"}
        if content.get("finite_schema_checked") is not True or not isinstance(sectors, list) or not required.issubset(sectors):
            return "finite packet-parent sectors are incomplete"
    elif receipt == "PROPAGATION_RECEIPT":
        if not all(_nonempty_string(content.get(name)) for name in ("channel", "propagation_kernel")) or not _true_fields(
            content, "normalization_checked", "domain_codomain_checked"
        ):
            return "propagation-kernel checks did not pass"
    elif receipt == "DETECTION_THINNING_RECEIPT":
        if not all(_nonempty_string(content.get(name)) for name in ("detector_model", "exposure")) or not _true_fields(
            content, "fit_checked", "probability_bounds_checked"
        ):
            return "detector-thinning checks did not pass"
    elif receipt == "CENSORING_AND_UPPER_LIMIT_RECEIPT":
        if not _nonempty_string(content.get("censoring_model")) or not _true_fields(
            content, "zero_count_exposure_checked", "known_upper_limits_checked", "no_double_counting_checked"
        ):
            return "censoring and upper-limit checks did not pass"
    elif receipt == "POINT_PROCESS_LIKELIHOOD_RECEIPT":
        count = content.get("event_count")
        if not _nonempty_string(content.get("likelihood")) or content.get("compensator_included") is not True:
            return "point-process likelihood or compensator is missing"
        if isinstance(count, bool) or not isinstance(count, int) or count < 0 or not _finite_number(content.get("log_likelihood")):
            return "point-process evaluation is not finite"
    elif receipt == "HELDOUT_LIKELIHOOD_RECEIPT":
        count = content.get("heldout_event_count")
        if len(input_records) < 2 or not all(
            _nonempty_string(content.get(name)) for name in ("proper_score", "preregistration_id")
        ) or content.get("split_checked") is not True:
            return "held-out split and preregistration are required"
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0 or not _finite_number(content.get("heldout_score")):
            return "held-out likelihood result is not finite"
    elif receipt == "CONTROL_MODEL_RECEIPT":
        controls = content.get("controls")
        if not isinstance(controls, list) or len(set(controls)) < 2 or not all(_nonempty_string(x) for x in controls):
            return "at least two distinct controls are required"
        if not _finite_number(content.get("acceptance_threshold")) or not _true_fields(
            content, "preregistered", "absolute_calibration_checked", "multiple_testing_checked"
        ):
            return "control-model checks did not pass"
    elif receipt == "REFINEMENT_STABILITY_RECEIPT":
        regulators = content.get("regulators")
        defect = content.get("maximum_defect")
        tolerance = content.get("defect_tolerance")
        if not isinstance(regulators, list) or len(set(regulators)) < 2 or not all(_nonempty_string(x) for x in regulators):
            return "at least two distinct regulators are required"
        if not all(_finite_number(value) and value >= 0 for value in (defect, tolerance)) or defect > tolerance or content.get("stability_checked") is not True:
            return "refinement stability check did not pass"
    elif receipt == "FROZEN_HASHES_RECEIPT":
        expected = {record["path"]: record["sha256"] for record in input_records}
        if not _nonempty_string(content.get("freeze_id")) or content.get("immutable_checked") is not True:
            return "freeze identity or immutability check is missing"
        if content.get("frozen_artifacts") != expected:
            return "frozen artifact hashes do not match current inputs"
    else:
        return "receipt has no content checker"
    return None


def _invalid_evidence_checks(reason: str) -> dict[str, dict[str, str]]:
    return {name: {"status": "INVALID", "reason": reason} for name in EVIDENCE_KINDS}


def _check_evidence_artifact(
    receipt: str, offer: object, *, manifest_dir: Path, output_dir: Path
) -> dict[str, str]:
    if not isinstance(offer, dict) or set(offer) != {"path", "sha256"}:
        return {"status": "INVALID", "reason": "offer must contain only path and sha256"}
    offered_path = offer.get("path")
    offered_hash = offer.get("sha256")
    if not _nonempty_string(offered_path) or not isinstance(offered_hash, str) or not SHA256_RE.fullmatch(offered_hash):
        return {"status": "INVALID", "reason": "offer path or SHA-256 is malformed"}
    path = Path(offered_path)
    if not path.is_absolute():
        path = manifest_dir / path
    path = path.resolve()
    resolved_output = output_dir.resolve()
    if path == resolved_output or resolved_output in path.parents:
        return {"status": "INVALID", "reason": "evidence is inside this builder's output directory"}
    if not path.is_file():
        return {"status": "INVALID", "reason": "evidence path is not a regular file", "path": str(path)}
    raw = path.read_bytes()
    if not raw:
        return {"status": "INVALID", "reason": "evidence artifact is empty", "path": str(path)}
    actual_hash = sha256_bytes(raw)
    if actual_hash != offered_hash:
        return {"status": "INVALID", "reason": "evidence content hash mismatch", "path": str(path)}
    try:
        artifact = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError):
        return {"status": "INVALID", "reason": "evidence is not valid JSON", "path": str(path)}
    if not isinstance(artifact, dict) or not artifact:
        return {"status": "INVALID", "reason": "evidence must be a non-empty object", "path": str(path)}
    if artifact.get("schema_version") != EVIDENCE_SCHEMA:
        return {"status": "INVALID", "reason": "evidence schema version mismatch", "path": str(path)}
    if artifact.get("receipt") != receipt or artifact.get("kind") != EVIDENCE_KINDS[receipt]:
        return {"status": "INVALID", "reason": "evidence receipt or kind mismatch", "path": str(path)}
    producer = artifact.get("producer")
    if not _nonempty_string(producer) or producer == BUILDER_ID:
        return {"status": "INVALID", "reason": "evidence producer is missing or self-authored", "path": str(path)}
    inputs = artifact.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        return {"status": "INVALID", "reason": "evidence must record at least one input", "path": str(path)}
    input_records: list[dict[str, str]] = []
    seen_inputs: set[Path] = set()
    for item in inputs:
        if not isinstance(item, dict) or set(item) != {"path", "sha256"}:
            return {"status": "INVALID", "reason": "input record must contain only path and sha256", "path": str(path)}
        input_name = item.get("path")
        input_hash = item.get("sha256")
        if not _nonempty_string(input_name) or not isinstance(input_hash, str) or not SHA256_RE.fullmatch(input_hash):
            return {"status": "INVALID", "reason": "input path or SHA-256 is malformed", "path": str(path)}
        input_path = Path(input_name)
        if not input_path.is_absolute():
            input_path = path.parent / input_path
        input_path = input_path.resolve()
        if input_path == path or input_path in seen_inputs or not input_path.is_file():
            return {"status": "INVALID", "reason": "input is missing, duplicate, or circular", "path": str(path)}
        seen_inputs.add(input_path)
        if sha256_bytes(input_path.read_bytes()) != input_hash:
            return {"status": "INVALID", "reason": "recorded input hash is stale", "path": str(path)}
        input_records.append({"path": str(input_name), "sha256": input_hash})
    try:
        content_error = _validate_receipt_content(receipt, artifact.get("content"), input_records)
    except (KeyError, OverflowError, TypeError, ValueError):
        content_error = "evidence content has invalid value types"
    if content_error is not None:
        return {"status": "INVALID", "reason": content_error, "path": str(path)}
    return {"status": "VALIDATED", "reason": "all checks passed", "path": str(path), "sha256": actual_hash}


def check_evidence_manifest(
    evidence_manifest: Path | None, *, output_dir: Path
) -> dict[str, dict[str, str]]:
    if evidence_manifest is None:
        return {name: {"status": "MISSING", "reason": "no evidence artifact supplied"} for name in EVIDENCE_KINDS}
    path = evidence_manifest.resolve()
    if not path.is_file():
        return _invalid_evidence_checks("evidence manifest is not a regular file")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return _invalid_evidence_checks("evidence manifest is not valid JSON")
    if not isinstance(manifest, dict) or set(manifest) != {"schema_version", "receipts"}:
        return _invalid_evidence_checks("manifest must contain only schema_version and receipts")
    if manifest.get("schema_version") != EVIDENCE_MANIFEST_SCHEMA or not isinstance(manifest.get("receipts"), dict):
        return _invalid_evidence_checks("evidence manifest schema mismatch")
    offers = manifest["receipts"]
    unknown = sorted(set(offers) - set(EVIDENCE_KINDS))
    if unknown:
        return _invalid_evidence_checks("unsupported receipt offers: " + ", ".join(unknown))
    return {
        name: _check_evidence_artifact(name, offers[name], manifest_dir=path.parent, output_dir=output_dir)
        if name in offers
        else {"status": "MISSING", "reason": "no evidence artifact supplied"}
        for name in EVIDENCE_KINDS
    }


def default_receipts(*, leak_hits: list[str]) -> dict[str, bool]:
    receipts = {name: False for name in RECEIPT_NAMES}
    receipts["NO_GENERATION_LEAKAGE_RECEIPT"] = not leak_hits
    return receipts


def strongest_claim(receipts: dict[str, bool], *, leak_hits: list[str]) -> tuple[str, str | None]:
    if leak_hits:
        return "CR1_QUOTIENT_DIAGNOSTIC", "NO_GENERATION_LEAKAGE_RECEIPT"
    for name in CR2_RECEIPTS:
        if not receipts.get(name, False):
            return "CR1_QUOTIENT_DIAGNOSTIC", name
    cr3_gate_labels = {
        "CONTROL_MODEL_RECEIPT": "CONTROLS",
        "REFINEMENT_STABILITY_RECEIPT": "REFINEMENT",
        "FROZEN_HASHES_RECEIPT": "FREEZE",
    }
    for name, gate_label in cr3_gate_labels.items():
        if not receipts.get(name, False):
            return "CR2_CONDITIONAL_PHENOMENOLOGY", gate_label
    for name in (
        "COMPACT_SOURCE_ACTION_DERIVED_RECEIPT",
        "EMISSION_MICROPHYSICS_DERIVED_RECEIPT",
        "PHYSICAL_CLOCK_DERIVED_RECEIPT",
        "OLD_HOST_FRB_SOURCE_THEOREM_RECEIPT",
        "BH_GENEALOGY_PRIOR_THEOREM_RECEIPT",
    ):
        if not receipts.get(name, False):
            return "CR3_FROZEN_PHYSICAL_PREDICTION", name
    return "CR4_SOURCE_ONLY_OPH_PREDICTION", None


def base_payload(name: str, receipts: dict[str, bool], claim: str) -> dict[str, Any]:
    return {
        "artifact": name,
        "generated_by": BUILDER_ID,
        "generated_utc": now_utc(),
        "claim": claim,
        "physical_claim": claim in {"CR3_FROZEN_PHYSICAL_PREDICTION", "CR4_SOURCE_ONLY_OPH_PREDICTION"},
        "readiness_gates": receipts,
    }


def build_payloads(
    *,
    config: Path | None,
    evidence_manifest: Path | None = None,
    output_dir: Path = DEFAULT_OUT,
) -> dict[str, str | dict[str, Any]]:
    leak_hits = target_leak_hits(config)
    checks = check_evidence_manifest(evidence_manifest, output_dir=output_dir)
    receipts = default_receipts(leak_hits=leak_hits)
    for name in EVIDENCE_KINDS:
        receipts[name] = checks[name]["status"] == "VALIDATED"
    claim, first_blocked = strongest_claim(receipts, leak_hits=leak_hits)
    return {
        "compact_history.json": {
            **base_payload("compact_transient_history", receipts, claim),
            "history_object": "Hist_r^CR",
            "requires_source_identity": True,
            "requires_genealogy_for_bh": True,
        },
        "compact_quotient.json": {
            **base_payload("compact_transient_quotient", receipts, claim),
            "quotient": "Q_r^CR = Sigma_r^CR/Gamma_r^CR",
            "likelihood_may_read_representative_labels": False,
        },
        "source_law.json": {
            **base_payload("compact_transient_source_law", receipts, claim),
            "law": "mu_r^CR plus packetized K_Gamma,r^hist",
            "normal_form_is_not_source_law": True,
        },
        "repair_emission_kernel.json": {
            **base_payload("repair_emission_kernel", receipts, claim),
            "kernel": "K_Gamma,r(dq',dPi,dell,dtau|q)",
            "independence_shortcuts_require_factorization_theorem": True,
        },
        "packet_parent.json": {
            **base_payload("finite_packet_parent", receipts, claim),
            "required_sectors": ["radio", "gamma", "GW", "optical", "neutrino", "environmental", "recipient"],
            "scalar_event_row_sufficient": False,
        },
        "detector_thinning.json": {
            **base_payload("detector_thinning", receipts, claim),
            "kernel": "Thin_c(dO|y,ObsWin_c)=p_det(y;ObsWin_c) R_det,c(dO|y)",
        },
        "censoring.json": {
            **base_payload("censoring_and_upper_limits", receipts, claim),
            "kernel": "Cens_c(dU|y,ObsWin_c)=(1-p_det)U_c(dU|y)",
            "score_nondetections": True,
        },
        "point_process_likelihood.json": {
            **base_payload("marked_catalog_likelihood", receipts, claim),
            "likelihood": "sum_i log Lambda(O_i)-integral_ObsWin Lambda(dO)",
            "compensator_required": True,
        },
        "frb_controls.json": {
            **base_payload("frb_repair_reload_controls", receipts, claim),
            "controls": {
                "M0": "young_only",
                "M1": "young_plus_old_gc_poisson_or_weibull_timing",
                "M2": "young_plus_old_gc_repair_reload_timing",
            },
            "first_prediction": "old/GC repeaters show fluence-conditioned recovery after cadence and exposure correction",
            "host_mixture_rank_required": True,
        },
        "bh_recycling.json": {
            **base_payload("black_hole_recycling", receipts, claim),
            "target_leak_hits": leak_hits,
            "genealogy_dag_required": True,
            "forbidden_path": "ringdown_residual -> generation_label -> claim_success",
            "repair_tail_template": "frozen damped sinusoid with independent generation prior",
        },
        "refinement_accuracy.json": {
            **base_payload("refinement_and_accuracy", receipts, claim),
            "accuracy_bound": (
                "epsilon_mu + E[N] epsilon_K + epsilon_E + epsilon_prop + "
                "epsilon_detector + epsilon_canon + epsilon_clock + epsilon_mc"
            ),
            "simulator_accuracy_receipt": receipts["SIMULATOR_ACCURACY_RECEIPT"],
        },
        "promotion_audit.json": {
            **base_payload("compact_transient_promotion_audit", receipts, claim),
            "claim_tiers": list(CLAIM_TIERS),
            "verdict_checking_tier": VERDICT_CHECKING_TIER,
            "scientific_validation_performed": False,
            "evidence_checks": checks,
            "first_blocked_gate": first_blocked,
            "promotion_allowed": claim in {"CR3_FROZEN_PHYSICAL_PREDICTION", "CR4_SOURCE_ONLY_OPH_PREDICTION"},
            "target_leak_hits": leak_hits,
        },
        "claim.md": claim + "\n",
    }


def write_payload(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def build_bundle(
    out_dir: Path, *, config: Path | None = None, evidence_manifest: Path | None = None
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    payloads = build_payloads(
        config=config,
        evidence_manifest=evidence_manifest,
        output_dir=out_dir,
    )
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        path = out_dir / rel_path
        write_payload(path, payload)
        file_hashes[rel_path] = sha256_bytes(path.read_bytes())
    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    audit = json.loads((out_dir / "promotion_audit.json").read_text(encoding="utf-8"))
    evidence_checks = audit["evidence_checks"]
    missing = [
        rel_path
        for rel_path in REQUIRED_FILES
        if rel_path != "manifest.json" and not (out_dir / rel_path).is_file()
    ]
    manifest = {
        "artifact": "compact_transient_receipt_manifest",
        "generated_by": BUILDER_ID,
        "generated_utc": now_utc(),
        "milestone": "COMPACT_TRANSIENT_RECEIPT_SCAFFOLD",
        "strongest_allowed_claim": claim,
        "first_blocked_gate": audit["first_blocked_gate"],
        "promotion_allowed": audit["promotion_allowed"],
        "physical_claim": audit["physical_claim"],
        "target_leak_hits": audit["target_leak_hits"],
        "evidence_checks": evidence_checks,
        "required_files": list(REQUIRED_FILES),
        "missing_files": missing,
        "file_hashes": file_hashes,
    }
    write_payload(out_dir / "manifest.json", manifest)
    manifest["file_hashes"]["manifest.json"] = sha256_bytes((out_dir / "manifest.json").read_bytes())
    write_payload(out_dir / "manifest.json", manifest)
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=DEFAULT_OUT, type=Path)
    parser.add_argument("--config", default=None, type=Path)
    parser.add_argument("--evidence-manifest", default=None, type=Path)
    args = parser.parse_args(argv)
    manifest = build_bundle(
        args.output,
        config=args.config,
        evidence_manifest=args.evidence_manifest,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
