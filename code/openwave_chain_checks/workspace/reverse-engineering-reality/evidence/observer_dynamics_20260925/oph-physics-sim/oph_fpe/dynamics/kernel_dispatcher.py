from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from oph_fpe.dynamics.positive_geometry import write_positive_geometry_kernel_report


SUPPORTED_KERNELS = ("positive_geometry",)


def dispatch_configured_kernels(
    config: dict[str, Any],
    run_dir: Path,
    *,
    engine: str,
    source: str = "configured_sim_run",
) -> dict[str, Any]:
    """Run opt-in certified-kernel hooks and write a fail-closed receipt.

    The dispatcher is intentionally a receipt/routing layer. It does not alter
    OPH repair dynamics unless a future backend supplies a certified readout
    substitution for the concrete run sector.
    """

    kernel_config = _kernel_config(config)
    requested = _requested_kernel_names(kernel_config)
    if not requested:
        return {}

    report: dict[str, Any] = {
        "mode": "oph_kernel_dispatch_v0",
        "engine": engine,
        "source": source,
        "configured": True,
        "requested_kernels": requested,
        "supported_kernels": list(SUPPORTED_KERNELS),
        "generic_repair_executed": True,
        "readout_substitution_implemented": False,
        "effective_acceleration_enabled": False,
        "physical_observables_changed": False,
        "claim_boundary": (
            "Kernel dispatch is an opt-in, fail-closed accelerator interface. "
            "Receipts may certify native mathematical structure, but this local "
            "simulator keeps generic OPH repair as the trusted execution path "
            "until an OPH sector compiler and readout-substitution backend are "
            "certified for the concrete run sector."
        ),
        "kernels": {},
        "unsupported_kernels": [name for name in requested if name not in SUPPORTED_KERNELS],
    }

    if "positive_geometry" in requested:
        report["kernels"]["positive_geometry"] = _dispatch_positive_geometry(
            kernel_config.get("positive_geometry"),
            run_dir,
            engine=engine,
            source=source,
        )

    report["effective_acceleration_enabled"] = any(
        bool(row.get("effective_acceleration_enabled", False))
        for row in (report.get("kernels") or {}).values()
    )
    report["routing_decision"] = _routing_decision(report)
    write_kernel_dispatch_report(run_dir, report)
    return report


def write_kernel_dispatch_report(out_dir: Path, report: dict[str, Any]) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "kernel_dispatch_report.json").write_text(
        json.dumps(report, indent=2, default=str),
        encoding="utf-8",
    )
    (out / "kernel_dispatch_report.md").write_text(_markdown_report(report), encoding="utf-8")


def kernel_dispatch_manifest_summary(report: dict[str, Any]) -> dict[str, Any]:
    if not report:
        return {}
    positive_geometry = (report.get("kernels") or {}).get("positive_geometry") or {}
    return {
        "mode": report.get("mode"),
        "routing_decision": report.get("routing_decision"),
        "generic_repair_executed": bool(report.get("generic_repair_executed", True)),
        "readout_substitution_implemented": bool(report.get("readout_substitution_implemented", False)),
        "effective_acceleration_enabled": bool(report.get("effective_acceleration_enabled", False)),
        "physical_observables_changed": bool(report.get("physical_observables_changed", False)),
        "positive_geometry": {
            "requested": bool(positive_geometry),
            "dispatch_status": positive_geometry.get("dispatch_status"),
            "native_verdict": positive_geometry.get("native_verdict"),
            "native_execution_mode": positive_geometry.get("native_execution_mode"),
            "native_trusted_acceleration": bool(positive_geometry.get("native_trusted_acceleration", False)),
            "effective_acceleration_enabled": bool(
                positive_geometry.get("effective_acceleration_enabled", False)
            ),
            "fallback_required": bool(positive_geometry.get("fallback_required", True)),
        },
    }


def _kernel_config(config: dict[str, Any]) -> dict[str, Any]:
    for key in ("kernels", "kernel_acceleration"):
        value = config.get(key)
        if isinstance(value, dict):
            return value
    return {}


def _requested_kernel_names(kernel_config: dict[str, Any]) -> list[str]:
    if not kernel_config or kernel_config.get("enabled") is False:
        return []
    names: list[str] = []
    for name, value in kernel_config.items():
        if name == "enabled":
            continue
        if _entry_enabled(value):
            names.append(str(name))
    return sorted(names)


def _entry_enabled(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, dict):
        return bool(value.get("enabled", True))
    return False


def _dispatch_positive_geometry(
    raw_config: Any,
    run_dir: Path,
    *,
    engine: str,
    source: str,
) -> dict[str, Any]:
    pg_config = raw_config if isinstance(raw_config, dict) else {}
    strict = bool(pg_config.get("strict", False))
    try:
        report = write_positive_geometry_kernel_report(
            run_dir,
            manifest_path=_optional_path(pg_config.get("manifest_path", pg_config.get("manifest"))),
            pgk_root=_optional_path(pg_config.get("pgk_root")),
            source=f"{source}:{engine}:positive_geometry",
        )
    except Exception as exc:
        error = {
            "mode": "positive_geometry_dispatch_error",
            "dispatch_status": "fail_closed_error",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "fallback_required": True,
            "native_trusted_acceleration": False,
            "effective_acceleration_enabled": False,
        }
        Path(run_dir, "positive_geometry_kernel_error.json").write_text(
            json.dumps(error, indent=2, default=str),
            encoding="utf-8",
        )
        if strict:
            raise
        return error

    native_acceleration = bool(report.get("trusted_acceleration_enabled", False))
    semantic_substitution_requested = bool(pg_config.get("allow_semantic_substitution", True))
    substitution_implemented = False
    effective_acceleration = native_acceleration and semantic_substitution_requested and substitution_implemented
    return {
        "dispatch_status": _positive_geometry_dispatch_status(
            report,
            native_acceleration=native_acceleration,
            semantic_substitution_requested=semantic_substitution_requested,
            substitution_implemented=substitution_implemented,
        ),
        "report_file": "positive_geometry_kernel_report.json",
        "receipt_file": "positive_geometry_kernel_receipt.json",
        "native_verdict": (report.get("receipt") or {}).get("verdict"),
        "native_execution_mode": report.get("execution_mode"),
        "native_trusted_acceleration": native_acceleration,
        "semantic_substitution_requested": semantic_substitution_requested,
        "substitution_implemented": substitution_implemented,
        "effective_acceleration_enabled": effective_acceleration,
        "fallback_required": not effective_acceleration,
        "physical_observables_changed": False,
    }


def _positive_geometry_dispatch_status(
    report: dict[str, Any],
    *,
    native_acceleration: bool,
    semantic_substitution_requested: bool,
    substitution_implemented: bool,
) -> str:
    if native_acceleration and semantic_substitution_requested and substitution_implemented:
        return "certified_accelerated_substitution"
    if native_acceleration and semantic_substitution_requested:
        return "certified_acceleration_available_no_sim_substitution_backend"
    if report.get("execution_mode") in {"CERTIFIED_GEOMETRY", "CERTIFIED_NATIVE"}:
        return "geometry_certified_diagnostic_only"
    return "exact_generic_fallback"


def _routing_decision(report: dict[str, Any]) -> str:
    kernels = report.get("kernels") or {}
    if any(bool(row.get("effective_acceleration_enabled", False)) for row in kernels.values()):
        return "certified_kernel_readout_substitution"
    if kernels:
        return "generic_oph_repair_with_kernel_receipts"
    if report.get("unsupported_kernels"):
        return "generic_oph_repair_unsupported_kernel_request"
    return "generic_oph_repair"


def _optional_path(value: Any) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(str(value)).expanduser()
    if path.is_absolute():
        return path
    return (Path.cwd() / path).resolve()


def _markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# OPH Kernel Dispatch",
        "",
        str(report.get("claim_boundary", "")),
        "",
        "## Routing",
        "",
        f"- engine: `{report.get('engine')}`",
        f"- routing decision: `{report.get('routing_decision')}`",
        f"- generic repair executed: `{str(report.get('generic_repair_executed')).lower()}`",
        f"- readout substitution implemented: `{str(report.get('readout_substitution_implemented')).lower()}`",
        f"- effective acceleration enabled: `{str(report.get('effective_acceleration_enabled')).lower()}`",
        f"- physical observables changed: `{str(report.get('physical_observables_changed')).lower()}`",
        "",
        "## Kernels",
        "",
    ]
    kernels = report.get("kernels") or {}
    if not kernels:
        lines.append("- none")
    for name, kernel in kernels.items():
        lines.append(
            "- "
            f"{name}: `{kernel.get('dispatch_status')}` "
            f"(native verdict `{kernel.get('native_verdict')}`, "
            f"effective acceleration `{str(kernel.get('effective_acceleration_enabled')).lower()}`)"
        )
    unsupported = report.get("unsupported_kernels") or []
    if unsupported:
        lines.extend(["", "## Unsupported", ""])
        lines.extend(f"- {name}" for name in unsupported)
    lines.append("")
    return "\n".join(lines)
