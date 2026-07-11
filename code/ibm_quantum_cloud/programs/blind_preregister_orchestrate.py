#!/usr/bin/env python3
"""Offline production orchestration for the blinded IBM preregistration.

The command builds one production-random catalog, derives and binds the
hardened analysis lock, verifies all canonical logical OpenQASM 3 digests, and
writes the three artifacts as a commit-last transaction.  It does not import
credentials, connect to IBM, or submit jobs.
"""

from __future__ import annotations

import argparse
import json
import os
import stat
import tempfile
from pathlib import Path
from typing import Any, Callable, Mapping

import blind_preregister as preregister


class SealOrchestrationError(RuntimeError):
    """Raised when production sealing cannot complete without ambiguity."""


CatalogBuilder = Callable[[], tuple[dict[str, Any], dict[str, Any]]]
AnalysisBuilder = Callable[[Mapping[str, Any], Mapping[str, Any]], dict[str, Any]]
AnalysisValidator = Callable[[Mapping[str, Any]], Any]


def _default_analysis_builder(
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
) -> dict[str, Any]:
    try:
        from build_blinded_analysis_lock import build_locked_analysis
    except ImportError as exc:
        raise SealOrchestrationError(
            "the hardened blinded-analysis lock builder is unavailable"
        ) from exc
    return build_locked_analysis(public_manifest, reveal)


def _default_analysis_validator(lock: Mapping[str, Any]) -> None:
    try:
        from cayley_blind_likelihood_analysis import validate_analysis_lock
    except ImportError as exc:
        raise SealOrchestrationError("the hardened analysis validator is unavailable") from exc
    validate_analysis_lock(lock, verify_code_hash=True)


def build_and_bind_production_bundle(
    *,
    catalog_builder: CatalogBuilder = preregister.build_blinded_preregistration,
    analysis_builder: AnalysisBuilder = _default_analysis_builder,
    analysis_validator: AnalysisValidator = _default_analysis_validator,
    _allow_test_mode: bool = False,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Execute the two-phase seal exactly once and return verified artifacts."""

    preliminary_public, preliminary_reveal = catalog_builder()
    preregister.verify_bundle(
        preliminary_public,
        preliminary_reveal,
        rebuild_circuits=False,
    )
    mode = preliminary_reveal["sealed_payload"].get("mode")
    if mode != "production_random" and not (_allow_test_mode and mode == "deterministic_test_only"):
        raise SealOrchestrationError("production sealing requires a production-random catalog")

    catalog_digest = preliminary_public["catalog_precommitment_sha256"]
    circuit_rows = preliminary_public["circuits"]
    secret_hex = preliminary_reveal["secret_hex"]
    private_circuits = preliminary_reveal["sealed_payload"]["circuits"]

    analysis_lock = analysis_builder(preliminary_public, preliminary_reveal)
    if not isinstance(analysis_lock, dict):
        raise SealOrchestrationError("analysis builder did not return a JSON object")
    if analysis_lock.get("catalog_precommitment_sha256") != catalog_digest:
        raise SealOrchestrationError("analysis lock binds a different catalog precommitment")
    analysis_validator(analysis_lock)

    final_public, final_reveal = preregister.bind_analysis_document(
        preliminary_public,
        preliminary_reveal,
        analysis_lock,
    )
    if final_public["catalog_precommitment_sha256"] != catalog_digest:
        raise SealOrchestrationError("catalog changed while binding the analysis lock")
    if final_public["circuits"] != circuit_rows:
        raise SealOrchestrationError("public circuit catalog changed while binding the lock")
    if final_reveal["secret_hex"] != secret_hex:
        raise SealOrchestrationError("production secret changed while binding the lock")
    if final_reveal["sealed_payload"]["circuits"] != private_circuits:
        raise SealOrchestrationError("private circuit catalog changed while binding the lock")
    if final_reveal["sealed_payload"]["analysis_document"] != analysis_lock:
        raise SealOrchestrationError("bound reveal does not contain the exact hardened lock")

    analysis_validator(analysis_lock)
    preregister.verify_bundle(
        final_public,
        final_reveal,
        analysis_lock,
        rebuild_circuits=True,
    )
    return final_public, final_reveal, analysis_lock


def _resolved_new_target(path: Path) -> Path:
    parent = path.parent.resolve()
    return parent / path.name


def _require_external_reveal(reveal_path: Path, checkout_root: Path) -> None:
    checkout = checkout_root.resolve()
    try:
        reveal_path.relative_to(checkout)
    except ValueError:
        return
    raise SealOrchestrationError("private reveal path must be outside the git checkout")


def _stage_json(target: Path, value: Mapping[str, Any], mode: int) -> Path:
    descriptor, raw_path = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=target.parent,
    )
    staged = Path(raw_path)
    try:
        os.fchmod(descriptor, mode)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(preregister.canonical_json_bytes(value) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        staged.unlink(missing_ok=True)
        raise
    return staged


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_sealed_bundle_atomically(
    *,
    public_path: Path,
    reveal_path: Path,
    analysis_path: Path,
    checkout_root: Path,
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
    analysis_lock: Mapping[str, Any],
) -> tuple[Path, Path, Path]:
    """Write private and analysis artifacts first, with public manifest last.

    Each target is atomically created without overwriting.  If a Python-level
    commit step fails, every target created by this transaction is removed.
    The public manifest is the commit marker, so it is never visible before
    the private reveal and analysis lock exist.
    """

    public_target = _resolved_new_target(Path(public_path))
    reveal_target = _resolved_new_target(Path(reveal_path))
    analysis_target = _resolved_new_target(Path(analysis_path))
    targets = (public_target, reveal_target, analysis_target)
    if len(set(targets)) != 3:
        raise SealOrchestrationError("public, reveal, and analysis paths must be distinct")
    _require_external_reveal(reveal_target, Path(checkout_root))
    existing = [str(path) for path in targets if path.exists() or path.is_symlink()]
    if existing:
        raise FileExistsError("refusing to overwrite sealed artifact: " + ", ".join(existing))
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)

    staged: dict[Path, Path] = {}
    committed: list[Path] = []
    try:
        staged[reveal_target] = _stage_json(reveal_target, reveal, 0o600)
        staged[analysis_target] = _stage_json(analysis_target, analysis_lock, 0o644)
        staged[public_target] = _stage_json(public_target, public_manifest, 0o644)
        for target in (reveal_target, analysis_target, public_target):
            # Hard-linking a same-directory staging inode creates the final
            # name atomically and fails rather than replacing an existing file.
            os.link(staged[target], target, follow_symlinks=False)
            committed.append(target)
            staged[target].unlink()
            _fsync_directory(target.parent)
    except BaseException:
        for target in reversed(committed):
            target.unlink(missing_ok=True)
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)
        raise

    if stat.S_IMODE(reveal_target.stat().st_mode) != 0o600:
        for target in reversed(committed):
            target.unlink(missing_ok=True)
        raise PermissionError("private reveal was not committed with mode 0600")
    return public_target, reveal_target, analysis_target


def public_seal_summary(public_manifest: Mapping[str, Any]) -> dict[str, Any]:
    resources = public_manifest["resources"]
    return {
        "analysis_lock_sha256": public_manifest["analysis_lock_sha256"],
        "catalog_precommitment_sha256": public_manifest["catalog_precommitment_sha256"],
        "circuit_count": len(public_manifest["circuits"]),
        "manifest_sha256": public_manifest["manifest_sha256"],
        "total_planned_shots": resources["total_planned_shots"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and atomically write the offline blinded production seal."
    )
    parser.add_argument("--public-out", type=Path, required=True)
    parser.add_argument("--reveal-out", type=Path, required=True)
    parser.add_argument("--analysis-lock-out", type=Path, required=True)
    parser.add_argument(
        "--checkout-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    public, reveal, analysis_lock = build_and_bind_production_bundle()
    write_sealed_bundle_atomically(
        public_path=args.public_out,
        reveal_path=args.reveal_out,
        analysis_path=args.analysis_lock_out,
        checkout_root=args.checkout_root,
        public_manifest=public,
        reveal=reveal,
        analysis_lock=analysis_lock,
    )
    print(json.dumps(public_seal_summary(public), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
