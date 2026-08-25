#!/usr/bin/env python3
"""Tests for the production backend bundle skeleton generator."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "hadron" / "generate_backend_export_bundle_skeleton.py"
RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"


def test_manifest_only_backend_bundle_skeleton_is_fail_closed_until_backend_fills_it() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = pathlib.Path(tmp) / "bundle"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--receipt",
                str(RECEIPT),
                "--payload",
                str(PAYLOAD),
                "--out-dir",
                str(out_dir),
                "--manifest-only",
            ],
            check=True,
            cwd=ROOT,
        )
        manifest = json.loads((out_dir / "backend_run_manifest.json").read_text(encoding="utf-8"))
        dataset_index = json.loads((out_dir / "dataset_index.json").read_text(encoding="utf-8"))

        assert manifest["artifact"] == "oph_hadron_backend_raw_export"
        assert manifest["execution_class"] == "production_request_template"
        assert manifest["claim_tier"] == "unfilled_backend_export_template"
        assert manifest["public_promotion_allowed"] is False
        assert manifest["files"]["correlators_hdf5"] == "correlators.h5"
        assert dataset_index["artifact"] == "oph_hadron_backend_dataset_index"
        assert dataset_index["placeholder_hdf5_written"] is False
        assert len(dataset_index["datasets"]) == 36
        assert all(item["fill_status"] == "requires_real_backend_output" for item in dataset_index["datasets"])
        assert not (out_dir / "correlators.h5").exists()
