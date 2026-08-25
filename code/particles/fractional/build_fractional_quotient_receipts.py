#!/usr/bin/env python3
"""Build the fractional quotient-sector receipt scaffold.

This mirrors the simulator-side `oph_fractional` sandbox in the paper-stack
code tree. It records the quotient, topological, optical, identifiability, and
no-target-leak receipts needed for the fractional exciton/FQAH problem note.
It does not promote any material sample into an OPH detection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import generate_runtime_kernel_harness as runtime_harness


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "fractional" / "quotient_sector_sandbox"

# Lean-to-Python bridge for the four quotient-correctness gates.  The Lean
# module renders the certified checker verdicts on the exact sandbox instance
# into a JSON payload and pins the rendered text with `#guard_msgs`, so any
# drift in the instance data or verdicts fails the Lean build.  The checked-in
# certificate file carries the same payload; before any verdict is read from
# it, it is verified byte-for-byte (after parsing) against the pin in the Lean
# source.  Dependency direction: Lean source -> pinned payload -> checked-in
# certificate -> this builder.  Nothing here calls Lean at build time, and
# nothing in Lean depends on this scaffold.
LEAN_CERTIFICATE_SOURCE = (
    ROOT.parent / "Lean" / "ObserverPatchHolography" / "QuotientLumpability.lean"
)
DEFAULT_CERTIFICATE = (
    Path(__file__).resolve().parent / "fractional_quotient_certificate.json"
)
RUNTIME_CAPTURE = Path(__file__).resolve().parent / "runtime_kernel_capture.json"
RUNTIME_HARNESS = (
    ROOT.parent
    / "Lean"
    / "ObserverPatchHolography"
    / "QuotientLumpabilityRuntimeHarness.lean"
)
LEAN_CERTIFIED_RECEIPTS = (
    "CANONICALIZER_IDEMPOTENCE",
    "REPRESENTATIVE_INVARIANCE",
    "QUOTIENT_LUMPABILITY",
    "NO_ORBIT_SIZE_BIAS",
)
_LEAN_PIN_PATTERN = re.compile(
    r"/--\s*info:\s*(\{.*?\})\s*-/\s*#guard_msgs\s+in\s*"
    r"#eval\s+IO\.println\s+\w*CertificateJson",
    re.DOTALL,
)

CLAIM_TIERS = (
    "DIAGNOSTIC_ONLY",
    "FRACTIONAL_QUOTIENT_SANDBOX_DIAGNOSTIC",
    "MATERIAL_HAMILTONIAN_PROMOTED",
    "FROZEN_OPTICAL_PREDICTION",
    "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE",
)

RECEIPTS = (
    "MATERIAL_QUOTIENT_NORMAL_FORM_RECEIPT",
    "CANONICALIZER_IDEMPOTENCE",
    "REPRESENTATIVE_INVARIANCE",
    "QUOTIENT_LUMPABILITY",
    "NO_ORBIT_SIZE_BIAS",
    "SOURCE_HAMILTONIAN_FROZEN",
    "SOURCE_LAW_FROZEN",
    "NO_TARGET_LEAK_DAG",
    "ACTIVE_PROJECTOR_CERTIFICATE",
    "CHERN_BAND_CERTIFICATE",
    "BAND_GEOMETRY_CERTIFICATE",
    "GAP_CERTIFICATE",
    "GROUND_STATE_DEGENERACY_CERTIFICATE",
    "FLUX_INSERTION_CHARGE_PUMP_CERTIFICATE",
    "HALL_CONDUCTANCE_CERTIFICATE",
    "K_T_LEDGER_CERTIFICATE",
    "EDGE_COLLAR_CERTIFICATE",
    "BULK_EDGE_CONSISTENCY_CERTIFICATE",
    "OPTICAL_OPERATOR_CERTIFICATE",
    "OPTICAL_MODULE_CERTIFICATE",
    "LINE_FAN_FROZEN_PREDICTION",
    "BINDING_DRIFT_BOUND",
    "OPTICAL_SECTOR_IDENTIFIED",
    "FROZEN_COMPARISON_RECEIPT",
    "SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT",
    "MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT",
    "EXPERIMENT_SPECIFIC_SOURCE_LAW_RECEIPT",
    "FROZEN_SAMPLE_COMPARISON_RECEIPT",
)

SANDBOX_REQUIREMENTS = (
    "MATERIAL_QUOTIENT_NORMAL_FORM_RECEIPT",
    "CANONICALIZER_IDEMPOTENCE",
    "REPRESENTATIVE_INVARIANCE",
    "QUOTIENT_LUMPABILITY",
    "NO_ORBIT_SIZE_BIAS",
    "SOURCE_HAMILTONIAN_FROZEN",
    "SOURCE_LAW_FROZEN",
    "NO_TARGET_LEAK_DAG",
    "ACTIVE_PROJECTOR_CERTIFICATE",
    "CHERN_BAND_CERTIFICATE",
    "BAND_GEOMETRY_CERTIFICATE",
    "GAP_CERTIFICATE",
    "GROUND_STATE_DEGENERACY_CERTIFICATE",
    "FLUX_INSERTION_CHARGE_PUMP_CERTIFICATE",
    "HALL_CONDUCTANCE_CERTIFICATE",
    "K_T_LEDGER_CERTIFICATE",
    "EDGE_COLLAR_CERTIFICATE",
    "BULK_EDGE_CONSISTENCY_CERTIFICATE",
    "OPTICAL_OPERATOR_CERTIFICATE",
    "OPTICAL_MODULE_CERTIFICATE",
    "LINE_FAN_FROZEN_PREDICTION",
    "BINDING_DRIFT_BOUND",
    "OPTICAL_SECTOR_IDENTIFIED",
    "FROZEN_COMPARISON_RECEIPT",
    "SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT",
)

MATERIAL_REQUIREMENTS = (
    "MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT",
    "EXPERIMENT_SPECIFIC_SOURCE_LAW_RECEIPT",
    "FROZEN_SAMPLE_COMPARISON_RECEIPT",
)

FAIL_CLOSED_STATES = (
    "SOURCE_NOT_FROZEN",
    "NOT_QUOTIENT_INVARIANT",
    "CANONICALIZER_NOT_IDEMPOTENT",
    "KERNEL_NOT_LUMPABLE",
    "ORBIT_SIZE_BIAS_DETECTED",
    "NO_GAP_CERTIFICATE",
    "CHERN_NUMBER_UNSTABLE",
    "PHASE_CERTIFICATE_NONINJECTIVE",
    "SECTOR_AMBIGUOUS",
    "OPTICAL_OPERATOR_UNCERTIFIED",
    "BINDING_DRIFT_UNBOUNDED",
    "OPTICAL_SECTOR_AMBIGUOUS",
    "TARGET_LEAK_DETECTED",
    "REFINEMENT_DEFECT_TOO_LARGE",
    "DIAGNOSTIC_ONLY",
)

FORBIDDEN_SOURCE_TOKENS = (
    "optical_peak_measurement",
    "optical_peak_measurements",
    "fqah_optical_spectrum",
    "fqah_transport_data",
    "anyon_trion_line_positions",
    "experimental_line_fan",
    "post_hoc_gate_slope",
    "target_fractional_charge",
    "target_chern_number",
)

REQUIRED_FILES = (
    "manifest.json",
    "lean_certificate.json",
    "material_presentation.json",
    "quotient_schema.json",
    "source_law.json",
    "hamiltonian_certificate.json",
    "active_band.json",
    "manybody_certificate.json",
    "topological_ledger.json",
    "edge_collar.json",
    "optical_module.json",
    "line_fan.json",
    "identifiability.json",
    "simulator_correctness.json",
    "refinement.json",
    "no_target_leak_dag.json",
    "failure_states.json",
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
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        haystack = text.lower()
    else:
        haystack = json.dumps(parsed, sort_keys=True).lower()
    return sorted(token for token in FORBIDDEN_SOURCE_TOKENS if token in haystack)


def lean_pinned_certificates(source: Path = LEAN_CERTIFICATE_SOURCE) -> list[dict[str, Any]]:
    """Extract every #guard_msgs-pinned certificate payload from the Lean source."""
    if not source.is_file():
        return []
    text = source.read_text(encoding="utf-8")
    payloads: list[dict[str, Any]] = []
    for match in _LEAN_PIN_PATTERN.finditer(text):
        try:
            payloads.append(json.loads(match.group(1)))
        except json.JSONDecodeError:
            continue
    return payloads


def load_lean_certificate(certificate: Path | None) -> tuple[dict[str, Any] | None, str]:
    """Load the certificate and bind it to the Lean pin.

    Returns the parsed certificate (when parseable) and a status.  Only
    `CERTIFICATE_PINNED` allows verdicts to be read; every other status reads
    the four certified receipts as False.
    """
    path = certificate or DEFAULT_CERTIFICATE
    if not path.is_file():
        return None, "CERTIFICATE_MISSING"
    try:
        cert = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None, "CERTIFICATE_UNPARSEABLE"
    pins = lean_pinned_certificates()
    if not pins:
        return cert, "LEAN_PIN_NOT_FOUND"
    if cert not in pins:
        return cert, "CERTIFICATE_NOT_PINNED"
    return cert, "CERTIFICATE_PINNED"


def load_runtime_binding() -> tuple[dict[str, Any] | None, str]:
    """Bind the captured simulator revision to its generated Lean harness."""
    if not RUNTIME_CAPTURE.is_file():
        return None, "RUNTIME_CAPTURE_MISSING"
    try:
        cap = json.loads(RUNTIME_CAPTURE.read_text(encoding="utf-8"))
        rendered = runtime_harness.render(cap)
    except (
        json.JSONDecodeError,
        OSError,
        KeyError,
        TypeError,
        runtime_harness.CaptureError,
    ):
        return None, "RUNTIME_CAPTURE_INVALID"
    if not RUNTIME_HARNESS.is_file():
        return cap, "RUNTIME_HARNESS_MISSING"
    if RUNTIME_HARNESS.read_text(encoding="utf-8") != rendered:
        return cap, "RUNTIME_HARNESS_STALE"
    return cap, "RUNTIME_CAPTURE_PINNED"


DECLARED_SANDBOX_RECEIPTS = tuple(
    name
    for name in SANDBOX_REQUIREMENTS
    if name not in LEAN_CERTIFIED_RECEIPTS
    and name not in ("NO_TARGET_LEAK_DAG", "SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT")
)


def receipt_provenance() -> dict[str, str]:
    provenance: dict[str, str] = {}
    for name in RECEIPTS:
        if name in LEAN_CERTIFIED_RECEIPTS:
            provenance[name] = "COMPUTED_FROM_LEAN_CERTIFICATE"
        elif name == "NO_TARGET_LEAK_DAG":
            provenance[name] = "COMPUTED_TOKEN_SCAN"
        elif name == "SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT":
            provenance[name] = "DERIVED_CONJUNCTION_WITH_RUNTIME_BINDING"
        elif name in SANDBOX_REQUIREMENTS:
            provenance[name] = "DECLARED_SANDBOX_SCAFFOLD"
        else:
            provenance[name] = "DECLARED_UNMET"
    return provenance


def build_receipts(config: Path | None, certificate: Path | None = None) -> dict[str, bool]:
    leak_hits = target_leak_hits(config)
    cert, cert_status = load_lean_certificate(certificate)
    _, runtime_status = load_runtime_binding()
    verdicts = cert.get("verdicts", {}) if isinstance(cert, dict) else {}
    receipts = {name: False for name in RECEIPTS}
    for name in DECLARED_SANDBOX_RECEIPTS:
        receipts[name] = True
    for name in LEAN_CERTIFIED_RECEIPTS:
        receipts[name] = cert_status == "CERTIFICATE_PINNED" and verdicts.get(name) is True
    receipts["NO_TARGET_LEAK_DAG"] = not leak_hits
    receipts["SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT"] = all(
        receipts[name]
        for name in (
            "MATERIAL_QUOTIENT_NORMAL_FORM_RECEIPT",
            "CANONICALIZER_IDEMPOTENCE",
            "REPRESENTATIVE_INVARIANCE",
            "QUOTIENT_LUMPABILITY",
            "NO_ORBIT_SIZE_BIAS",
            "NO_TARGET_LEAK_DAG",
        )
    ) and runtime_status == "RUNTIME_CAPTURE_PINNED"
    return receipts


def strongest_allowed_claim(receipts: dict[str, bool]) -> tuple[str, str | None, list[str]]:
    missing_sandbox = [name for name in SANDBOX_REQUIREMENTS if not receipts.get(name, False)]
    if missing_sandbox:
        return "DIAGNOSTIC_ONLY", missing_sandbox[0], missing_sandbox
    missing_material = [name for name in MATERIAL_REQUIREMENTS if not receipts.get(name, False)]
    if missing_material:
        return "FRACTIONAL_QUOTIENT_SANDBOX_DIAGNOSTIC", missing_material[0], missing_material
    return "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE", None, []


def base_payload(artifact: str, *, claim: str, receipts: dict[str, bool]) -> dict[str, Any]:
    return {
        "artifact": artifact,
        "claim": claim,
        "promotion_allowed": claim == "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE",
        "material_claim": claim in {"MATERIAL_HAMILTONIAN_PROMOTED", "FROZEN_OPTICAL_PREDICTION", "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE"},
        "readiness_gates": receipts,
    }


def build_payloads(
    config: Path | None, certificate: Path | None = None
) -> dict[str, str | dict[str, Any]]:
    receipts = build_receipts(config, certificate)
    claim, first_blocked, missing = strongest_allowed_claim(receipts)
    leak_hits = target_leak_hits(config)
    cert, cert_status = load_lean_certificate(certificate)
    runtime_capture, runtime_status = load_runtime_binding()
    cert_path = certificate or DEFAULT_CERTIFICATE
    base = base_payload("fractional_quotient_base", claim=claim, receipts=receipts)
    return {
        "lean_certificate.json": {
            **base,
            "artifact": "fractional_lean_certificate_binding",
            "mechanism": (
                "The Lean module renders its certified checker verdicts on the exact "
                "sandbox instance into a JSON payload pinned by #guard_msgs; the "
                "checked-in certificate file must match that pin before any verdict "
                "is read. Lean depends on nothing here; this builder only reads."
            ),
            "lean_module": "ObserverPatchHolography.QuotientLumpability",
            "lean_source": "Lean/ObserverPatchHolography/QuotientLumpability.lean",
            "certificate_path": (
                str(cert_path.resolve().relative_to(ROOT.parent))
                if cert_path.resolve().is_relative_to(ROOT.parent)
                else str(cert_path)
            ),
            "certificate_sha256": (
                sha256_bytes(cert_path.read_bytes()) if cert_path.is_file() else None
            ),
            "certificate_status": cert_status,
            "certified_receipts": {
                name: receipts[name] for name in LEAN_CERTIFIED_RECEIPTS
            },
            "certificate_verdicts": (cert or {}).get("verdicts"),
            "certificate_theorems": (cert or {}).get("theorems"),
            "certificate_witnesses": (cert or {}).get("witnesses"),
            "runtime_capture_path": str(RUNTIME_CAPTURE.relative_to(ROOT.parent)),
            "runtime_capture_sha256": (
                sha256_bytes(RUNTIME_CAPTURE.read_bytes())
                if RUNTIME_CAPTURE.is_file()
                else None
            ),
            "runtime_capture_status": runtime_status,
            "runtime_capture_commit": (runtime_capture or {}).get(
                "captured_at_commit"
            ),
            "runtime_capture_tree": (runtime_capture or {}).get("captured_tree"),
            "runtime_harness_path": str(RUNTIME_HARNESS.relative_to(ROOT.parent)),
            "runtime_harness_sha256": (
                sha256_bytes(RUNTIME_HARNESS.read_bytes())
                if RUNTIME_HARNESS.is_file()
                else None
            ),
            "runtime_freshness": (
                "PINNED_STATIC_SNAPSHOT; run "
                "generate_runtime_kernel_harness.py --check --sim-repo <path> "
                "to verify a current simulator checkout"
            ),
            "fail_closed": (
                "any certificate_status other than CERTIFICATE_PINNED reads all four "
                "certified receipts as False; any runtime_capture_status other than "
                "RUNTIME_CAPTURE_PINNED makes SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT "
                "False"
            ),
        },
        "material_presentation.json": {
            **base,
            "artifact": "fractional_material_presentation",
            "material_class": "twisted_TMD_FCI_FQAH_sandbox",
            "normal_form_scope": "finite quotient-sector presentation",
            "status": "SANDBOX_PRESENTATION_ONLY",
        },
        "quotient_schema.json": {
            **base,
            "artifact": "fractional_quotient_schema",
            "canonicalizer_idempotence": receipts["CANONICALIZER_IDEMPOTENCE"],
            "representative_invariance": receipts["REPRESENTATIVE_INVARIANCE"],
            "quotient_lumpability": receipts["QUOTIENT_LUMPABILITY"],
            "no_orbit_size_bias": receipts["NO_ORBIT_SIZE_BIAS"],
            "verdict_source": "lean_certificate.json",
            "certificate_status": cert_status,
        },
        "source_law.json": {
            **base,
            "artifact": "fractional_source_law",
            "source_status": "TOY_SOURCE_FROZEN_FOR_SANDBOX",
            "material_specific_source_required": "EXPERIMENT_SPECIFIC_SOURCE_LAW_RECEIPT",
        },
        "hamiltonian_certificate.json": {
            **base,
            "artifact": "fractional_hamiltonian_promotion_certificate",
            "status": "SANDBOX_PHASE_PROMOTION_ONLY",
            "material_specific_required": "MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT",
        },
        "active_band.json": {
            **base,
            "artifact": "fractional_active_band",
            "chern_band_certificate": receipts["CHERN_BAND_CERTIFICATE"],
            "active_projector_certificate": receipts["ACTIVE_PROJECTOR_CERTIFICATE"],
        },
        "manybody_certificate.json": {
            **base,
            "artifact": "fractional_manybody_certificate",
            "gap_certificate": receipts["GAP_CERTIFICATE"],
            "ground_state_degeneracy": receipts["GROUND_STATE_DEGENERACY_CERTIFICATE"],
            "flux_insertion_charge_pump": receipts["FLUX_INSERTION_CHARGE_PUMP_CERTIFICATE"],
            "hall_conductance": receipts["HALL_CONDUCTANCE_CERTIFICATE"],
        },
        "topological_ledger.json": {
            **base,
            "artifact": "fractional_topological_ledger",
            "ledger": "Abelian K,t sandbox with neutral fractional exciton shadow sector",
            "receipt": "K_T_LEDGER_CERTIFICATE",
        },
        "edge_collar.json": {
            **base,
            "artifact": "fractional_edge_collar",
            "edge_collar_certificate": receipts["EDGE_COLLAR_CERTIFICATE"],
            "bulk_edge_consistency": receipts["BULK_EDGE_CONSISTENCY_CERTIFICATE"],
        },
        "optical_module.json": {
            **base,
            "artifact": "fractional_optical_module",
            "optical_operator_certificate": receipts["OPTICAL_OPERATOR_CERTIFICATE"],
            "module_certificate": receipts["OPTICAL_MODULE_CERTIFICATE"],
            "neutral_fractional_exciton_boundary": "neutral total charge may still carry nontrivial topological shadow tau",
        },
        "line_fan.json": {
            **base,
            "artifact": "fractional_line_fan",
            "line_fan_frozen_prediction": receipts["LINE_FAN_FROZEN_PREDICTION"],
            "binding_drift_bound": receipts["BINDING_DRIFT_BOUND"],
        },
        "identifiability.json": {
            **base,
            "artifact": "fractional_identifiability",
            "optical_sector_identified": receipts["OPTICAL_SECTOR_IDENTIFIED"],
            "ambiguous_default_fail_state": "OPTICAL_SECTOR_AMBIGUOUS",
        },
        "simulator_correctness.json": {
            **base,
            "artifact": "fractional_simulator_correctness",
            "simulator_quotient_correctness_receipt": receipts["SIMULATOR_QUOTIENT_CORRECTNESS_RECEIPT"],
            "runtime_capture_status": runtime_status,
            "runtime_capture_commit": (runtime_capture or {}).get(
                "captured_at_commit"
            ),
            "runtime_capture_tree": (runtime_capture or {}).get("captured_tree"),
            "freshness_scope": "PINNED_STATIC_SNAPSHOT",
        },
        "refinement.json": {
            **base,
            "artifact": "fractional_refinement",
            "status": "SANDBOX_REFINEMENT_COMPATIBLE",
            "fail_state": "REFINEMENT_DEFECT_TOO_LARGE",
        },
        "no_target_leak_dag.json": {
            **base,
            "artifact": "fractional_no_target_leak_dag",
            "forbidden_source_tokens": list(FORBIDDEN_SOURCE_TOKENS),
            "target_leak_hits": leak_hits,
            "status": "PASS_EMPTY_COMPARISON_DAG" if not leak_hits else "FAIL_FORBIDDEN_SOURCE_INPUT",
        },
        "failure_states.json": {
            **base,
            "artifact": "fractional_failure_states",
            "fail_closed_states": list(FAIL_CLOSED_STATES),
        },
        "receipts.json": {
            **base,
            "artifact": "fractional_receipts",
            "readiness_gates": receipts,
            "receipt_provenance": receipt_provenance(),
        },
        "claim_ladder.json": {
            **base,
            "artifact": "fractional_claim_ladder",
            "claim_ladder": list(CLAIM_TIERS),
            "first_blocked_gate": first_blocked,
            "missing_for_next_claim": missing,
            "claim_boundary": (
                "The sandbox may demonstrate OPH quotient-sector mechanics. A real material claim "
                "requires a material-specific Hamiltonian/source-law proof and frozen sample comparison."
            ),
        },
        "claim.md": claim + "\n",
    }


def write_payload(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_bundle(
    out_dir: Path, *, config: Path | None, certificate: Path | None = None
) -> dict[str, Any]:
    payloads = build_payloads(config, certificate)
    out_dir.mkdir(parents=True, exist_ok=True)
    file_hashes: dict[str, str] = {}
    for rel_path, payload in payloads.items():
        target = out_dir / rel_path
        write_payload(target, payload)
        file_hashes[rel_path] = sha256_bytes(target.read_bytes())

    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    ladder = json.loads((out_dir / "claim_ladder.json").read_text(encoding="utf-8"))
    manifest = {
        "artifact": "fractional_quotient_manifest",
        "generated_utc": now_utc(),
        "milestone": "FRACTIONAL_QUOTIENT_SECTOR_SANDBOX",
        "strongest_allowed_claim": claim,
        "first_blocked_gate": ladder["first_blocked_gate"],
        "promotion_allowed": claim == "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE",
        "material_claim": claim in {"MATERIAL_HAMILTONIAN_PROMOTED", "FROZEN_OPTICAL_PREDICTION", "EXPERIMENTAL_QUOTIENT_READOUT_CANDIDATE"},
        "required_files": list(REQUIRED_FILES),
        "missing_files": [rel for rel in REQUIRED_FILES if rel != "manifest.json" and not (out_dir / rel).is_file()],
        "file_hashes": file_hashes,
        "paper_problem_note": "physics-problems/fractional_excitons_as_oph_quotient_sector_readouts.md",
        "simulator_surface": "oph-physics-sim/oph_fractional",
    }
    write_payload(out_dir / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Build fractional quotient-sector receipt scaffold.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--config")
    parser.add_argument(
        "--certificate",
        help=(
            "Path to a Lean-pinned verdict certificate "
            f"(default: {DEFAULT_CERTIFICATE})"
        ),
    )
    args = parser.parse_args()

    manifest = build_bundle(
        Path(args.output),
        config=Path(args.config) if args.config else None,
        certificate=Path(args.certificate) if args.certificate else None,
    )
    print(Path(args.output) / "manifest.json")
    print(
        f"claim={manifest['strongest_allowed_claim']} "
        f"first_blocked_gate={manifest['first_blocked_gate']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
