#!/usr/bin/env python3
"""Assemble a backend export from production array files.

This is the no-HDF5 ingestion bridge for backend teams that can emit one plain
numeric array file per required correlator path. It preserves the production
manifest metadata and writes the inline backend export accepted by
run_production_backend_writeback.py.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.production_execution_support import normalize_source_id  # noqa: E402


CHANNELS = ("pi_iso", "N_iso_direct", "N_iso_exchange")


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _load_float_array(path: Path, *, expected_length: int, field_name: str) -> list[float]:
    text = path.read_text(encoding="utf-8").replace(",", " ")
    values: list[float] = []
    for token in text.split():
        try:
            value = float(token)
        except Exception as exc:
            raise ValueError(f"{field_name} contains nonnumeric token {token!r}") from exc
        if not math.isfinite(value):
            raise ValueError(f"{field_name} contains nonfinite value {value!r}")
        values.append(value)
    if len(values) != expected_length:
        raise ValueError(f"{field_name} has length {len(values)}, expected {expected_length}")
    return values


def _dataset_lengths(dataset_index: dict[str, Any]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for item in dataset_index.get("datasets", []):
        lengths[str(item["path"])] = int(item["length"])
    return lengths


def build_inline_export(
    manifest: dict[str, Any],
    dataset_index: dict[str, Any],
    *,
    array_dir: str | Path,
) -> dict[str, Any]:
    array_dir = Path(array_dir)
    lengths = _dataset_lengths(dataset_index)
    index_by_path = {str(item["path"]): item for item in dataset_index.get("datasets", [])}
    raw_export_provenance = {
        "manifest_artifact": manifest.get("artifact"),
        "manifest_path": None,
        "correlators_hdf5": None,
        "array_dir": str(array_dir),
        "execution_class": manifest.get("execution_class") or "unclassified",
        "claim_tier": manifest.get("claim_tier"),
        "public_promotion_allowed": manifest.get("public_promotion_allowed", False),
        "profile_id": manifest.get("profile_id"),
        "backend": manifest.get("backend"),
        "physics": manifest.get("physics"),
        "solvers": manifest.get("solvers"),
        "integrator": manifest.get("integrator"),
        "boundary_conditions": manifest.get("boundary_conditions"),
        "sources": manifest.get("sources"),
        "contractions": manifest.get("contractions"),
    }
    out: dict[str, Any] = {
        "artifact": "oph_hadron_backend_raw_export_inlined",
        "format_version": int(manifest.get("format_version") or 1),
        "execution_class": manifest.get("execution_class") or "unclassified",
        "claim_tier": manifest.get("claim_tier"),
        "public_promotion_allowed": manifest.get("public_promotion_allowed", False),
        "profile_id": manifest.get("profile_id"),
        "raw_export_provenance": raw_export_provenance,
        "ensembles": {},
    }
    for ensemble in manifest.get("ensembles", []):
        ensemble_id = str(ensemble["ensemble_id"])
        cfgs: dict[str, Any] = {}
        for cfg in ensemble.get("cfgs", []):
            cfg_id = str(cfg["cfg_id"])
            sources: dict[str, Any] = {}
            for source in cfg.get("sources", []):
                src_id = normalize_source_id(str(source["src_id"]))
                source_payload: dict[str, Any] = {
                    "coord": list(source.get("coord") or []),
                }
                for channel in CHANNELS:
                    dataset_path = str((source.get("datasets") or {})[channel])
                    item = index_by_path.get(dataset_path)
                    if item is None:
                        raise ValueError(f"dataset {dataset_path!r} is missing from dataset_index")
                    array_file = item.get("array_file")
                    if not array_file:
                        raise ValueError(f"dataset {dataset_path!r} is missing array_file in dataset_index")
                    source_payload[channel] = _load_float_array(
                        array_dir / str(array_file),
                        expected_length=lengths[dataset_path],
                        field_name=dataset_path,
                    )
                sources[src_id] = source_payload
            cfgs[cfg_id] = {
                "trajectory_stop": cfg.get("trajectory_stop"),
                "sources": sources,
            }
        out["ensembles"][ensemble_id] = {
            "ensemble_id": ensemble_id,
            "cfgs": cfgs,
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Assemble an inline backend export from array files.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--dataset-index", required=True)
    parser.add_argument("--array-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    export = build_inline_export(
        _load_json(args.manifest),
        _load_json(args.dataset_index),
        array_dir=args.array_dir,
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(export, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
