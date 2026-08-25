#!/usr/bin/env python3
"""Report the exact backend-side readiness state of the hadron lane.

Chain role: sharpen the live execution-bound hadron frontier into explicit
backend-side readiness clauses: filled schedule receipt, publication-complete
manifest provenance, real correlator arrays, stable-channel evaluation, and the
final closure/publication gate.

Mathematics: readiness classification only; this script does not derive hadron
masses.

OPH-derived inputs: the runtime receipt, cfg/source payload shell, optional
backend manifest, optional normalized production dump, optional evaluation, and
optional closure report.

Output: a machine-readable hadron production-readiness report.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.validate_production_hadron_closure import _get_schedule_scalars, _is_finite_number
from particles.hadron.production_execution_support import (
    PRODUCTION_EXECUTION_CLASSES,
    PUBLIC_PROMOTION_BLOCKER,
)

DEFAULT_RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
DEFAULT_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_MANIFEST = ROOT / "particles" / "runs" / "hadron" / "oph_hadron_production_backend_manifest.json"
DEFAULT_DUMP = ROOT / "particles" / "runs" / "hadron" / "backend_correlator_dump.production.json"
DEFAULT_EVALUATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_evaluation.json"
DEFAULT_CLOSURE = ROOT / "particles" / "runs" / "hadron" / "hadron_production_closure_validation_report.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "hadron_production_readiness_report.json"
REQUIRED_PRODUCTION_CHANNELS = ("pi_iso", "N_iso_direct", "N_iso_exchange")
REQUIRED_LOCAL_PRODUCTS = (
    "backend_correlator_dump.production.json",
    "stable_channel_sequence_evaluation.json",
    "hadron_production_closure_validation_report.json",
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_optional_json(path: str | Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    path_obj = Path(path)
    if not path_obj.exists():
        return None
    return json.loads(path_obj.read_text(encoding="utf-8"))


def _value_at_path(payload: dict[str, Any], path: str) -> Any:
    current: Any = payload
    for token in path.split("."):
        if not isinstance(current, dict) or token not in current:
            return None
        current = current[token]
    return current


def _value_is_populated(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        text = value.strip()
        return bool(text) and not text.startswith("<fill-me>")
    if isinstance(value, dict):
        return bool(value)
    if isinstance(value, list):
        return bool(value)
    return True


def _manifest_provenance_status(manifest: dict[str, Any] | None) -> dict[str, Any]:
    required_paths = [
        "profile_id",
        "backend.family",
        "backend.name",
        "backend.version",
        "backend.git_commit",
        "backend.run_id",
        "backend.build_id",
        "backend.machine",
        "solvers.rhmc_strange.rational_coefficients",
        "integrator",
        "boundary_conditions",
        "sources",
        "contractions",
    ]
    if manifest is None:
        return {
            "manifest_present": False,
            "execution_class": None,
            "public_promotion_allowed": False,
            "production_execution_class": False,
            "required_paths": required_paths,
            "populated_paths": [],
            "missing_or_placeholder_paths": required_paths,
            "publication_complete": False,
        }
    populated = []
    missing = []
    for path in required_paths:
        value = _value_at_path(manifest, path)
        if _value_is_populated(value):
            populated.append(path)
        else:
            missing.append(path)
    execution_class = str(manifest.get("execution_class") or "unclassified")
    public_promotion_allowed = bool(manifest.get("public_promotion_allowed", False))
    return {
        "manifest_present": True,
        "execution_class": execution_class,
        "public_promotion_allowed": public_promotion_allowed,
        "production_execution_class": execution_class in PRODUCTION_EXECUTION_CLASSES,
        "required_paths": required_paths,
        "populated_paths": populated,
        "missing_or_placeholder_paths": missing,
        "publication_complete": not missing,
    }


def _dump_array_status(dump: dict[str, Any] | None) -> dict[str, Any]:
    if dump is None:
        return {
            "dump_present": False,
            "production_execution": False,
            "execution_class": None,
            "public_promotion_allowed": False,
            "all_required_arrays_finite": False,
            "required_channels": list(REQUIRED_PRODUCTION_CHANNELS),
            "missing_or_nonfinite_arrays": [],
        }
    missing: list[str] = []
    for ensemble_id, ensemble in (dump.get("ensembles") or {}).items():
        for cfg_id, cfg in ((ensemble.get("cfgs") or {}).items()):
            for src_id, source in ((cfg.get("sources") or {}).items()):
                expected_len = source.get("expected_t_extent")
                for channel in REQUIRED_PRODUCTION_CHANNELS:
                    values = source.get(channel)
                    field = f"{ensemble_id}.{cfg_id}.{src_id}.{channel}"
                    if not isinstance(values, list) or not values:
                        missing.append(field)
                        continue
                    if expected_len is not None and len(values) != int(expected_len):
                        missing.append(field)
                        continue
                    try:
                        if not all(math.isfinite(float(value)) for value in values):
                            missing.append(field)
                    except Exception:
                        missing.append(field)
    return {
        "dump_present": True,
        "production_execution": bool(dump.get("production_execution") is True),
        "execution_class": dump.get("execution_class"),
        "public_promotion_allowed": bool(dump.get("public_promotion_allowed", False)),
        "all_required_arrays_finite": not missing,
        "required_channels": list(REQUIRED_PRODUCTION_CHANNELS),
        "missing_or_nonfinite_arrays": missing,
    }


def build_readiness_report(
    receipt: dict[str, Any],
    payload: dict[str, Any],
    *,
    manifest: dict[str, Any] | None = None,
    dump: dict[str, Any] | None = None,
    evaluation: dict[str, Any] | None = None,
    closure_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    n_therm, n_sep = _get_schedule_scalars(receipt)
    receipt_filled = _is_finite_number(n_therm) and _is_finite_number(n_sep)
    manifest_status = _manifest_provenance_status(manifest)
    dump_status = _dump_array_status(dump)
    evaluation_complete = evaluation is not None and evaluation.get("status") == "production_measure_evaluation_complete"
    closure_public_ready = bool((closure_report or {}).get("public_unsuppression_ready"))
    publication_bundle_ready = (
        receipt_filled
        and manifest_status["publication_complete"]
        and manifest_status["production_execution_class"]
        and dump_status["production_execution"]
        and dump_status["public_promotion_allowed"]
        and dump_status["all_required_arrays_finite"]
        and evaluation_complete
        and closure_public_ready
    )

    if not receipt_filled:
        smallest_residual = "runtime_schedule_receipt_N_therm_and_N_sep with explicit external N_therm/N_sep inputs"
    elif not manifest_status["manifest_present"]:
        smallest_residual = (
            "production backend export bundle on the seeded family with publication-complete manifest provenance "
            "and real correlator arrays"
        )
    elif not manifest_status["production_execution_class"]:
        smallest_residual = "real production backend execution class, not a diagnostic or surrogate backend bundle"
    elif not manifest_status["public_promotion_allowed"] or not dump_status["public_promotion_allowed"]:
        smallest_residual = PUBLIC_PROMOTION_BLOCKER
    elif not manifest_status["publication_complete"]:
        smallest_residual = "publication-complete backend manifest provenance on the seeded family"
    elif not dump_status["dump_present"] or not dump_status["production_execution"] or not dump_status["all_required_arrays_finite"]:
        smallest_residual = "backend correlator arrays from real production RHMC/HMC execution on the theorem-emitted seeded family"
    elif not evaluation_complete:
        smallest_residual = "stable_channel_sequence_evaluation with populated forward-window and published statistical/systematic fields for pi_iso and N_iso"
    elif not closure_public_ready:
        smallest_residual = (
            (closure_report or {}).get("smallest_live_residual_object")
            or "hadron production closure report with public stable-channel unsuppression ready"
        )
    else:
        smallest_residual = None

    exact_remaining_runtime_object = None
    if smallest_residual is not None:
        exact_remaining_runtime_object = {
            "name": "production_backend_export_bundle",
            "status": "open" if not publication_bundle_ready else "closed",
            "definition": (
                "production backend export bundle on the seeded family with publication-complete manifest provenance "
                "and real correlator arrays"
            ),
            "required_channels": list(REQUIRED_PRODUCTION_CHANNELS),
            "required_local_products_after_normalization": list(REQUIRED_LOCAL_PRODUCTS),
            "runtime_receipt_contract": {
                "N_therm": n_therm,
                "N_sep": n_sep,
            },
            "runner_path": "particles/hadron/run_production_backend_writeback.py",
        }

    return {
        "artifact": "oph_hadron_production_readiness_report",
        "generated_utc": _timestamp(),
        "runner_path": "particles/hadron/run_production_backend_writeback.py",
        "runner_available": True,
        "runtime_receipt": {
            "artifact": receipt.get("artifact"),
            "status": receipt.get("status"),
            "receipt_filled": receipt_filled,
            "required_schedule_scalars": {"N_therm": n_therm, "N_sep": n_sep},
        },
        "payload_surface": {
            "artifact": payload.get("artifact"),
            "status": payload.get("status"),
            "cfg_support_realization_status": payload.get("cfg_support_realization_status"),
            "support_realization_schedule_status": (payload.get("support_realization_schedule") or {}).get("status"),
        },
        "backend_manifest_publication_status": manifest_status,
        "production_dump_status": dump_status,
        "sequence_evaluation_status": {
            "artifact": None if evaluation is None else evaluation.get("artifact"),
            "status": None if evaluation is None else evaluation.get("status"),
            "production_measure_evaluation_complete": evaluation_complete,
        },
        "closure_status": {
            "artifact": None if closure_report is None else closure_report.get("artifact"),
            "closure_grade": None if closure_report is None else closure_report.get("closure_grade"),
            "public_unsuppression_ready": closure_public_ready,
        },
        "publication_bundle_ready": publication_bundle_ready,
        "smallest_backend_residual_object": smallest_residual,
        "exact_remaining_runtime_object": exact_remaining_runtime_object,
        "notes": [
            "This report sharpens the backend-side hadron frontier beyond the older generic dump wording.",
            "Numeric stable-channel closure and publication/provenance readiness are tracked separately on purpose.",
            "A closure report can be numerically ready while publication_bundle_ready remains false if backend manifest provenance is still incomplete.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the hadron production readiness report.")
    parser.add_argument("--receipt", default=str(DEFAULT_RECEIPT))
    parser.add_argument("--payload", default=str(DEFAULT_PAYLOAD))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--dump", default=str(DEFAULT_DUMP))
    parser.add_argument("--evaluation", default=str(DEFAULT_EVALUATION))
    parser.add_argument("--closure-report", default=str(DEFAULT_CLOSURE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    receipt = _load_optional_json(args.receipt)
    payload = _load_optional_json(args.payload)
    if receipt is None or payload is None:
        raise FileNotFoundError("receipt and payload are required to build the hadron production readiness report")
    manifest = _load_optional_json(args.manifest)
    dump = _load_optional_json(args.dump)
    evaluation = _load_optional_json(args.evaluation)
    closure_report = _load_optional_json(args.closure_report)
    report = build_readiness_report(
        receipt,
        payload,
        manifest=manifest,
        dump=dump,
        evaluation=evaluation,
        closure_report=closure_report,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
