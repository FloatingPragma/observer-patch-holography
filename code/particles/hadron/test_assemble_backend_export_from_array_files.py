#!/usr/bin/env python3
"""Tests for assembling backend exports from plain array files."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "particles" / "hadron" / "generate_backend_export_bundle_skeleton.py"
ASSEMBLER = ROOT / "particles" / "hadron" / "assemble_backend_export_from_array_files.py"
RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"


def _write_array_file(path: pathlib.Path, length: int, offset: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = [f"{1.0 + offset + 0.001 * idx:.12g}" for idx in range(length)]
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def test_assemble_backend_export_from_array_files_roundtrips_manifest_contract() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        bundle_dir = pathlib.Path(tmp) / "bundle"
        output_path = pathlib.Path(tmp) / "backend_export_inlined.json"
        subprocess.run(
            [
                sys.executable,
                str(GENERATOR),
                "--receipt",
                str(RECEIPT),
                "--payload",
                str(PAYLOAD),
                "--out-dir",
                str(bundle_dir),
                "--manifest-only",
            ],
            check=True,
            cwd=ROOT,
        )
        dataset_index = json.loads((bundle_dir / "dataset_index.json").read_text(encoding="utf-8"))
        for offset, item in enumerate(dataset_index["datasets"]):
            _write_array_file(bundle_dir / item["array_file"], int(item["length"]), offset)

        subprocess.run(
            [
                sys.executable,
                str(ASSEMBLER),
                "--manifest",
                str(bundle_dir / "backend_run_manifest.json"),
                "--dataset-index",
                str(bundle_dir / "dataset_index.json"),
                "--array-dir",
                str(bundle_dir),
                "--output",
                str(output_path),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(output_path.read_text(encoding="utf-8"))
        assert payload["artifact"] == "oph_hadron_backend_raw_export_inlined"
        assert payload["execution_class"] == "production_request_template"
        assert payload["public_promotion_allowed"] is False
        first_ensemble = payload["ensembles"]["qcd_2p1_seed_n0"]
        first_source = first_ensemble["cfgs"]["qcd_2p1_seed_n0__cfg0"]["sources"]["src0"]
        assert len(first_source["pi_iso"]) == 291
        assert len(first_source["N_iso_direct"]) == 291
        assert len(first_source["N_iso_exchange"]) == 291
