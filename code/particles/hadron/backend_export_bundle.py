#!/usr/bin/env python3
"""Helpers for raw hadron backend export bundles.

This module defines a binary export contract that is practical for real RHMC/HMC
production jobs: a JSON manifest plus an HDF5 correlator store.

The normalized JSON dump consumed by the existing OPH hadron pipeline remains the
same. The only new logic here is loading a raw backend bundle and turning it into
that already-frozen normalized shape.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

_MANIFEST_BASENAMES = (
    "backend_run_manifest.json",
    "backend_export_manifest.json",
    "raw_backend_export_manifest.json",
)
_CHANNELS = ("pi_iso", "N_iso_direct", "N_iso_exchange")


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _canonical_dataset_path(ensemble_id: str, cfg_id: str, src_id: str, channel: str) -> str:
    return f"/ensembles/{ensemble_id}/cfgs/{cfg_id}/sources/{src_id}/{channel}"


def _coord4(values: Any, *, field_name: str) -> list[int]:
    seq = list(values or [])
    if len(seq) != 4:
        raise ValueError(f"{field_name} must contain exactly 4 coordinates")
    out: list[int] = []
    for idx, value in enumerate(seq):
        try:
            out.append(int(value))
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{field_name}[{idx}] is not an integer") from exc
    return out


def _manifest_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_dir():
        for basename in _MANIFEST_BASENAMES:
            candidate = path / basename
            if candidate.exists():
                return candidate
        raise FileNotFoundError(
            f"no backend export manifest found in {path}; expected one of {', '.join(_MANIFEST_BASENAMES)}"
        )
    return path


def _as_float_list(values: Any, *, expected_length: int, field_name: str) -> list[float]:
    try:
        seq = list(values)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{field_name} is not sequence-like") from exc
    if len(seq) != expected_length:
        raise ValueError(f"{field_name} has length {len(seq)}, expected {expected_length}")
    out = []
    for idx, value in enumerate(seq):
        try:
            f = float(value)
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{field_name}[{idx}] is not numeric") from exc
        if not math.isfinite(f):
            raise ValueError(f"{field_name}[{idx}] is not finite")
        out.append(f)
    return out


def _normalize_sources_for_cfg(cfg_entry: dict[str, Any]) -> list[dict[str, Any]]:
    sources = cfg_entry.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError(f"cfg {cfg_entry.get('cfg_id')!r} must provide a nonempty sources list")
    out = []
    for source in sources:
        src_id = str(source["src_id"])
        coord = _coord4(
            source.get("coord") or source.get("coords"),
            field_name=f"cfg {cfg_entry.get('cfg_id')!r} source {src_id!r} coord",
        )
        datasets = dict(source.get("datasets") or {})
        out.append(
            {
                "src_id": src_id,
                "coord": coord,
                "datasets": datasets,
            }
        )
    return out


def _normalize_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    if str(manifest.get("artifact")) != "oph_hadron_backend_raw_export":
        raise ValueError("manifest artifact must be 'oph_hadron_backend_raw_export'")
    ensembles = manifest.get("ensembles")
    if not isinstance(ensembles, list) or not ensembles:
        raise ValueError("raw backend export manifest must contain a nonempty ensembles list")
    normalized = {
        "artifact": "oph_hadron_backend_raw_export",
        "format_version": int(manifest.get("format_version") or 1),
        "profile_id": manifest.get("profile_id"),
        "backend": dict(manifest.get("backend") or {}),
        "physics": dict(manifest.get("physics") or {}),
        "solvers": dict(manifest.get("solvers") or {}),
        "integrator": dict(manifest.get("integrator") or {}),
        "boundary_conditions": dict(manifest.get("boundary_conditions") or {}),
        "sources": dict(manifest.get("sources") or {}),
        "contractions": dict(manifest.get("contractions") or {}),
        "files": dict(manifest.get("files") or {}),
        "ensembles": [],
    }
    for ensemble in ensembles:
        cfgs = ensemble.get("cfgs")
        if not isinstance(cfgs, list) or not cfgs:
            raise ValueError(f"ensemble {ensemble.get('ensemble_id')!r} must provide a nonempty cfgs list")
        normalized_cfgs = []
        for cfg in cfgs:
            normalized_cfgs.append(
                {
                    "cfg_id": str(cfg["cfg_id"]),
                    "trajectory_stop": cfg.get("trajectory_stop"),
                    "sources": _normalize_sources_for_cfg(cfg),
                }
            )
        normalized["ensembles"].append(
            {
                "ensemble_id": str(ensemble["ensemble_id"]),
                "cfgs": normalized_cfgs,
            }
        )
    return normalized


def load_backend_input_artifact(path: str | Path) -> dict[str, Any]:
    """Load either an already-normalized JSON export or a raw bundle manifest.

    Accepted inputs:
    - existing normalized inline JSON export consumed by current code
    - a raw manifest JSON with artifact='oph_hadron_backend_raw_export'
    - a directory containing a raw manifest plus correlators.h5
    """
    manifest_path = _manifest_path(path)
    payload = load_json(manifest_path)
    if str(payload.get("artifact")) != "oph_hadron_backend_raw_export":
        return payload
    return backend_input_from_raw_manifest(payload, base_dir=manifest_path.parent, manifest_path=manifest_path)


def backend_input_from_raw_manifest(
    manifest: dict[str, Any],
    *,
    base_dir: str | Path,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Convert a raw export manifest + HDF5 correlator store into the normalized tree."""
    normalized_manifest = _normalize_manifest(manifest)
    correlator_file = (normalized_manifest.get("files") or {}).get("correlators_hdf5")
    if not correlator_file:
        raise ValueError("raw backend export manifest must provide files.correlators_hdf5")
    correlator_path = Path(base_dir) / str(correlator_file)
    if not correlator_path.exists():
        raise FileNotFoundError(f"correlator store not found: {correlator_path}")

    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError("h5py is required to load raw backend export bundles") from exc

    out: dict[str, Any] = {
        "artifact": "oph_hadron_backend_raw_export_inlined",
        "format_version": normalized_manifest["format_version"],
        "profile_id": normalized_manifest.get("profile_id"),
        "raw_export_provenance": {
            "manifest_artifact": normalized_manifest.get("artifact"),
            "manifest_path": str(Path(manifest_path).resolve()) if manifest_path is not None else None,
            "correlators_hdf5": str(correlator_path),
            "backend": normalized_manifest.get("backend"),
            "physics": normalized_manifest.get("physics"),
            "solvers": normalized_manifest.get("solvers"),
            "integrator": normalized_manifest.get("integrator"),
            "boundary_conditions": normalized_manifest.get("boundary_conditions"),
            "sources": normalized_manifest.get("sources"),
            "contractions": normalized_manifest.get("contractions"),
        },
        "ensembles": {},
    }

    with h5py.File(correlator_path, "r") as h5:
        for ensemble in normalized_manifest["ensembles"]:
            ensemble_id = str(ensemble["ensemble_id"])
            cfg_tree: dict[str, Any] = {}
            for cfg in ensemble["cfgs"]:
                cfg_id = str(cfg["cfg_id"])
                source_tree: dict[str, Any] = {}
                for source in cfg["sources"]:
                    src_id = str(source["src_id"])
                    coords = list(source.get("coord") or [])
                    datasets = source.get("datasets") or {}
                    source_payload: dict[str, Any] = {}
                    for channel in _CHANNELS:
                        dataset_path = str(datasets.get(channel) or _canonical_dataset_path(ensemble_id, cfg_id, src_id, channel))
                        if dataset_path not in h5:
                            raise ValueError(
                                f"missing dataset {dataset_path!r} for {ensemble_id}.{cfg_id}.{src_id}.{channel}"
                            )
                        dset = h5[dataset_path]
                        if len(dset.shape) != 1:
                            raise ValueError(f"dataset {dataset_path!r} must be 1D")
                        source_payload[channel] = _as_float_list(
                            dset[...],
                            expected_length=int(dset.shape[0]),
                            field_name=f"{ensemble_id}.{cfg_id}.{src_id}.{channel}",
                        )
                    source_tree[src_id] = {
                        "coord": coords,
                        **source_payload,
                    }
                cfg_tree[cfg_id] = {
                    "trajectory_stop": cfg.get("trajectory_stop"),
                    "sources": source_tree,
                }
            out["ensembles"][ensemble_id] = {
                "ensemble_id": ensemble_id,
                "cfgs": cfg_tree,
            }
    return out


def build_backend_export_skeleton(
    receipt: dict[str, Any],
    payload: dict[str, Any],
    *,
    profile_id: str = "oph_reference_backend_v1",
) -> tuple[dict[str, Any], list[tuple[str, int]]]:
    """Return a raw bundle manifest and the canonical dataset list.

    The caller is responsible for creating the HDF5 file and datasets.
    """
    ensemble_payloads = {
        str(entry["ensemble_id"]): entry
        for entry in payload.get("ensemble_payloads", [])
    }
    manifest: dict[str, Any] = {
        "artifact": "oph_hadron_backend_raw_export",
        "format_version": 1,
        "profile_id": profile_id,
        "backend": {
            "family": "rhmc_hmc",
            "name": "<fill-me>",
            "version": "<fill-me>",
            "git_commit": "<fill-me>",
            "run_id": "<fill-me>",
            "build_id": "<fill-me>",
            "machine": "<fill-me>",
        },
        "physics": {
            "branch": "seeded_2p1_qed_off",
            "target_channels": ["pi_iso", "N_iso"],
            "notes": [
                "Fill this manifest with the exact backend profile and runtime provenance before publication.",
                "Do not replace NaN datasets with synthetic values. All datasets must come from real backend execution.",
            ],
        },
        "solvers": {},
        "integrator": {},
        "boundary_conditions": {},
        "sources": {},
        "contractions": {},
        "files": {
            "correlators_hdf5": "correlators.h5",
        },
        "ensembles": [],
    }
    datasets: list[tuple[str, int]] = []
    for sched in (receipt.get("execution_contract") or {}).get("ensemble_schedule", []):
        ensemble_id = str(sched["ensemble_id"])
        payload_entry = ensemble_payloads.get(ensemble_id)
        if payload_entry is None:
            raise ValueError(f"payload missing ensemble {ensemble_id!r}")
        payload_t = int(payload_entry["T"])
        schedule_t = sched.get("t_extent")
        if schedule_t is not None and int(schedule_t) != payload_t:
            raise ValueError(
                f"t_extent mismatch for {ensemble_id}: receipt={schedule_t}, payload={payload_t}"
            )
        cfgs = []
        for cfg_id in sched.get("cfg_ids", []):
            cfg_id = str(cfg_id)
            src_descriptors = (payload_entry.get("source_descriptors_by_cfg") or {}).get(cfg_id, [])
            if not src_descriptors:
                raise ValueError(f"payload missing source descriptors for {ensemble_id}.{cfg_id}")
            source_items = []
            for src_desc in src_descriptors:
                src_id = str(src_desc["src_id"])
                coord = _coord4(
                    src_desc.get("coords") or src_desc.get("coord"),
                    field_name=f"{ensemble_id}.{cfg_id}.{src_id}.coord",
                )
                datasets_by_channel = {}
                for channel in _CHANNELS:
                    dset_path = _canonical_dataset_path(ensemble_id, cfg_id, src_id, channel)
                    datasets_by_channel[channel] = dset_path
                    datasets.append((dset_path, payload_t))
                source_items.append(
                    {
                        "src_id": src_id,
                        "coord": coord,
                        "datasets": datasets_by_channel,
                    }
                )
            cfgs.append(
                {
                    "cfg_id": cfg_id,
                    "trajectory_stop": (sched.get("trajectory_stop_by_cfg") or {}).get(cfg_id),
                    "sources": source_items,
                }
            )
        manifest["ensembles"].append(
            {
                "ensemble_id": ensemble_id,
                "cfgs": cfgs,
            }
        )
    return manifest, datasets
