from __future__ import annotations

import copy
import importlib.util
import json
import shutil
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any


PGK_BUNDLE_RELATIVE = Path("amplituhedron/engineering-the-simulation/oph_kernel_sdk_v2")
DEFAULT_MANIFEST_RELATIVE = PGK_BUNDLE_RELATIVE / "examples/A614_manifest_v2.json"


def default_workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def default_pgk_root() -> Path:
    return default_workspace_root() / PGK_BUNDLE_RELATIVE


def default_pgk_manifest() -> Path:
    return default_workspace_root() / DEFAULT_MANIFEST_RELATIVE


def default_kernel_runtime() -> Path:
    return default_pgk_root() / "kernel_runtime.py"


def positive_geometry_kernel_report(
    manifest_path: Path | None = None,
    *,
    pgk_root: Path | None = None,
    registry: Any | None = None,
    source: str = "engineering_the_simulation_kernel_constitution_for_oph",
) -> dict[str, Any]:
    """Run the certified positive-geometry SDK v2 checker and emit OPH-FPE gates.

    This is an integration wrapper, not a new checker. The SDK v2 runtime is
    the authority: manifest status strings are descriptive only, and gates pass
    only when host-registered checkers execute.
    """

    root = Path(pgk_root) if pgk_root is not None else default_pgk_root()
    manifest = Path(manifest_path) if manifest_path is not None else default_pgk_manifest()
    runtime_path = root / "kernel_runtime.py"
    if not runtime_path.is_file() or not manifest.is_file():
        return _unavailable_plugin_report(
            runtime_path=runtime_path,
            manifest_path=manifest,
            source=source,
        )
    runtime = _load_kernel_runtime(runtime_path)
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    receipt = runtime.check_manifest(copy.deepcopy(manifest_payload), registry=registry)
    receipt_payload = asdict(receipt)

    gate_rows = receipt_payload.get("gates") or {}
    gate_status = {
        str(name): str(row.get("status"))
        for name, row in gate_rows.items()
        if isinstance(row, dict)
    }
    verdict = str(receipt_payload.get("verdict"))
    fallback = str(receipt_payload.get("fallback_action"))
    geometry_certified = verdict in {
        "CERTIFIED_NATIVE_BACKEND_NOT_ENABLED",
        "CERTIFIED_ACCELERATED",
    }
    acceleration_enabled = verdict == "CERTIFIED_ACCELERATED"
    gates = {
        "PGK_manifest_loaded": True,
        "PGK_native_checker_available": True,
        "PGK_structural_checks_passed": not bool(receipt_payload.get("errors")),
        "PGK_geometry_certified": geometry_certified,
        "PGK_sector_recognition_certified": gate_status.get("C1_recognition") == "PASS",
        "PGK_quotient_compiler_certified": gate_status.get("C2_compiler") == "PASS",
        "PGK_readout_equivalence_certified": gate_status.get("C4_semantics") == "PASS",
        "PGK_resource_accounting_certified": gate_status.get("C5_resources") == "PASS",
        "PGK_fail_closed_reproducibility_certified": gate_status.get("C6_fallback") == "PASS",
        "PGK_event_commit_certified": gate_status.get("C7_commit") == "PASS",
        "PGK_failed_gates_present": any(status == "FAIL" for status in gate_status.values()),
        "PGK_open_gates_present": any(status == "OPEN" for status in gate_status.values()),
        "PGK_acceleration_enabled": acceleration_enabled,
        "PGK_fallback_exact_generic": fallback == "EXACT_GENERIC",
        "trusted_oph_scattering_readout": acceleration_enabled,
        "generic_oph_repair_fallback_required": not acceleration_enabled,
    }
    sector = manifest_payload.get("sector") or {}
    native = (manifest_payload.get("native") or {}).get("payload") or {}
    return {
        "mode": "oph_positive_geometry_kernel_integration_v2",
        "source": source,
        "kernel_constitution": {
            "doctrine": "proof-carrying plugin with host-registered gates and fail-closed fallback",
            "runtime": str(runtime_path),
            "manifest_path": str(manifest),
            "fallback_path": "generic OPH repair",
            "manifest_status_strings_authoritative": False,
        },
        "sector": {
            "sector_id": manifest_payload.get("kernel_id"),
            "theory_tag": sector.get("theory_scope"),
            "native_object": sector.get("native_object"),
            "physical_OPH_sector_recognized": bool(sector.get("physical_OPH_sector_recognized", False)),
            "n": native.get("n"),
            "k": native.get("k"),
            "m": native.get("m"),
            "loop_order": native.get("loop_order"),
            "geometry_family": manifest_payload.get("plugin_type"),
        },
        "receipt": receipt_payload,
        "readiness_gates": gates,
        "execution_mode": _execution_mode(verdict),
        "trusted_acceleration_enabled": acceleration_enabled,
        "physical_prediction": False,
        "claim_boundary": (
            "Positive geometry is integrated as a certified accelerator backend only. "
            "The current SDK v2 A614 pilot can certify exact native geometry and fail-closed routing, "
            "but it is not promoted to a trusted OPH scattering readout unless host-registered "
            "sector recognition, quotient compiler, readout equivalence, and event-commit gates pass. "
            "Manifest status strings cannot self-promote a result."
        ),
    }


def _unavailable_plugin_report(
    *,
    runtime_path: Path,
    manifest_path: Path,
    source: str,
) -> dict[str, Any]:
    missing = [
        str(path)
        for path in (runtime_path, manifest_path)
        if not path.is_file()
    ]
    receipt = {
        "verdict": "DECLARED_PLUGIN_UNAVAILABLE",
        "fallback_action": "EXACT_GENERIC",
        "errors": [f"missing declared plugin input: {path}" for path in missing],
        "gates": {},
    }
    gates = {
        "PGK_manifest_loaded": manifest_path.is_file(),
        "PGK_native_checker_available": runtime_path.is_file(),
        "PGK_structural_checks_passed": False,
        "PGK_geometry_certified": False,
        "PGK_sector_recognition_certified": False,
        "PGK_quotient_compiler_certified": False,
        "PGK_readout_equivalence_certified": False,
        "PGK_resource_accounting_certified": False,
        "PGK_fail_closed_reproducibility_certified": False,
        "PGK_event_commit_certified": False,
        "PGK_failed_gates_present": False,
        "PGK_open_gates_present": True,
        "PGK_acceleration_enabled": False,
        "PGK_fallback_exact_generic": True,
        "trusted_oph_scattering_readout": False,
        "generic_oph_repair_fallback_required": True,
    }
    return {
        "mode": "oph_positive_geometry_kernel_integration_v2",
        "source": source,
        "kernel_constitution": {
            "doctrine": "declared external plugin with fail-closed fallback",
            "runtime": str(runtime_path),
            "manifest_path": str(manifest_path),
            "fallback_path": "generic OPH repair",
            "manifest_status_strings_authoritative": False,
        },
        "sector": {
            "sector_id": None,
            "theory_tag": None,
            "native_object": None,
            "physical_OPH_sector_recognized": False,
            "n": None,
            "k": None,
            "m": None,
            "loop_order": None,
            "geometry_family": None,
        },
        "receipt": receipt,
        "readiness_gates": gates,
        "execution_mode": "EXACT_GENERIC_FALLBACK",
        "trusted_acceleration_enabled": False,
        "physical_prediction": False,
        "claim_boundary": (
            "The positive-geometry SDK is a declared external plugin input. "
            "An absent runtime or manifest cannot promote acceleration, "
            "scattering, neutral-bulk, or physical-prediction receipts."
        ),
    }


def write_positive_geometry_kernel_report(
    out_dir: Path,
    manifest_path: Path | None = None,
    *,
    pgk_root: Path | None = None,
    source: str = "engineering_the_simulation_kernel_constitution_for_oph",
) -> dict[str, Any]:
    report = positive_geometry_kernel_report(
        manifest_path,
        pgk_root=pgk_root,
        source=source,
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "positive_geometry_kernel_report.json").write_text(
        json.dumps(report, indent=2, default=str),
        encoding="utf-8",
    )
    (out / "positive_geometry_kernel_report.md").write_text(_markdown_report(report), encoding="utf-8")
    manifest = Path(manifest_path) if manifest_path is not None else default_pgk_manifest()
    if manifest.exists():
        shutil.copyfile(manifest, out / "positive_geometry_kernel_manifest.json")
    (out / "positive_geometry_kernel_receipt.json").write_text(
        json.dumps(report["receipt"], indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return report


def _load_kernel_runtime(runtime_path: Path) -> Any:
    if not runtime_path.exists():
        raise FileNotFoundError(f"PGK runtime not found: {runtime_path}")
    module_name = "_oph_fpe_external_kernel_runtime_v2"
    spec = importlib.util.spec_from_file_location(module_name, runtime_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"could not load PGK runtime: {runtime_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _execution_mode(verdict: str) -> str:
    if verdict == "CERTIFIED_ACCELERATED":
        return "CERTIFIED_ACCELERATED"
    if verdict == "CERTIFIED_NATIVE_BACKEND_NOT_ENABLED":
        return "CERTIFIED_NATIVE"
    if str(verdict).startswith("REJECTED"):
        return "EXACT_GENERIC_FALLBACK"
    return "EXPERIMENTAL_UNTRUSTED"


def _fmt(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def _markdown_report(report: dict[str, Any]) -> str:
    receipt = report["receipt"]
    gates = report["readiness_gates"]
    return "\n".join(
        [
            "# OPH Positive Geometry Kernel",
            "",
            str(report["claim_boundary"]),
            "",
            "## Status",
            "",
            f"- verdict: `{receipt.get('verdict')}`",
            f"- execution mode: `{report['execution_mode']}`",
            f"- trusted acceleration enabled: `{str(report['trusted_acceleration_enabled']).lower()}`",
            f"- fallback action: `{receipt.get('fallback_action')}`",
            "",
            "## Gates",
            "",
            *[f"- {key}: `{_fmt(value)}`" for key, value in gates.items()],
            "",
        ]
    )
