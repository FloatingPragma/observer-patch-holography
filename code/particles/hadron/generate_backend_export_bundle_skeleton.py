#!/usr/bin/env python3
"""Emit a fillable HDF5 backend-export skeleton for hadron production.

This is the production-side counterpart of the existing JSON skeleton. Instead
of asking an external RHMC/HMC code to hand-build a giant JSON file, it emits a
canonical HDF5 dataset tree plus a manifest that records the exact dataset paths
and per-source coordinates.
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

from particles.hadron.backend_export_bundle import build_backend_export_skeleton


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a raw backend export bundle skeleton.")
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--payload", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--profile-id", default="oph_reference_backend_v1")
    args = parser.parse_args()

    try:
        import h5py  # type: ignore
        import numpy as np  # type: ignore
    except Exception as exc:  # pragma: no cover - dependency error path
        raise RuntimeError("h5py and numpy are required to generate the backend export skeleton") from exc

    receipt = _load_json(args.receipt)
    payload = _load_json(args.payload)
    manifest, datasets = build_backend_export_skeleton(receipt, payload, profile_id=args.profile_id)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "backend_run_manifest.json"
    h5_path = out_dir / "correlators.h5"

    with h5py.File(h5_path, "w") as h5:
        for dset_path, length in datasets:
            dset = h5.create_dataset(dset_path, shape=(length,), dtype="<f8")
            dset[...] = np.full((length,), math.nan, dtype=np.float64)
            dset.attrs["status"] = "fill_with_real_backend_output"

    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {manifest_path}")
    print(f"wrote {h5_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
