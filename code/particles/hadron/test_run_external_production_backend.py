#!/usr/bin/env python3
"""Tests for the external production backend runner."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "particles" / "hadron" / "generate_backend_export_bundle_skeleton.py"
RUNNER = ROOT / "particles" / "hadron" / "run_external_production_backend.py"
RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
SEQUENCE_POPULATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"


def _make_request(bundle_dir: pathlib.Path) -> None:
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


def _fill_manifest_provenance(bundle_dir: pathlib.Path) -> None:
    manifest_path = bundle_dir / "backend_run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["execution_class"] = "production"
    manifest["claim_tier"] = "production_backend_export_bundle"
    manifest["public_promotion_allowed"] = True
    manifest["backend"] = {
        "family": "rhmc_hmc",
        "name": "pytest-production-backend",
        "version": "1.0",
        "git_commit": "pytestabcdef",
        "run_id": "pytest-run",
        "build_id": "pytest-build",
        "machine": "pytest-node",
    }
    manifest["solvers"] = {
        "rhmc_strange": {
            "rational_coefficients": [{"alpha": 1.0, "beta": 2.0}],
        }
    }
    manifest["integrator"] = {"id": "pytest_omelyan", "n_steps": 16}
    manifest["boundary_conditions"] = {
        "gauge": {"x": "periodic", "y": "periodic", "z": "periodic", "t": "periodic"},
        "fermion": {"x": "periodic", "y": "periodic", "z": "periodic", "t": "anti_periodic"},
    }
    manifest["sources"] = {"construction": "local_point_from_payload"}
    manifest["contractions"] = {
        "pion": {"export_channel": "pi_iso"},
        "nucleon": {"export_channels": ["N_iso_direct", "N_iso_exchange"]},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_array_files(bundle_dir: pathlib.Path) -> None:
    dataset_index = json.loads((bundle_dir / "dataset_index.json").read_text(encoding="utf-8"))
    for offset, item in enumerate(dataset_index["datasets"]):
        path = bundle_dir / item["array_file"]
        path.parent.mkdir(parents=True, exist_ok=True)
        length = int(item["length"])
        dset_path = str(item["path"])
        if dset_path.endswith("/pi_iso"):
            values = [f"{(1.0 + 0.01 * offset) * math.exp(-0.02 * idx):.12g}" for idx in range(length)]
        elif dset_path.endswith("/N_iso_exchange"):
            values = [f"{(0.1 + 0.001 * offset) * math.exp(-0.05 * idx):.12g}" for idx in range(length)]
        elif dset_path.endswith("/N_iso_direct"):
            values = [
                f"{(1.2 + 0.01 * offset) * math.exp(-0.04 * idx) + (0.1 + 0.001 * offset) * math.exp(-0.05 * idx):.12g}"
                for idx in range(length)
            ]
        else:
            raise AssertionError(f"unexpected dataset path {dset_path}")
        path.write_text("\n".join(values) + "\n", encoding="utf-8")


def test_external_runner_refuses_without_backend_command_or_arrays() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        bundle_dir = pathlib.Path(tmp) / "bundle"
        _make_request(bundle_dir)
        result = subprocess.run(
            [
                sys.executable,
                str(RUNNER),
                "--request-dir",
                str(bundle_dir),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "Refusing to synthesize production correlators" in (result.stderr + result.stdout)


def test_external_runner_validates_existing_array_files_and_runs_writeback() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = pathlib.Path(tmp)
        bundle_dir = tmpdir / "bundle"
        _make_request(bundle_dir)
        _fill_manifest_provenance(bundle_dir)
        _write_array_files(bundle_dir)
        backend_export = tmpdir / "backend_export.json"
        dump = tmpdir / "dump.json"
        manifest = tmpdir / "manifest.json"
        evaluation = tmpdir / "evaluation.json"
        closure = tmpdir / "closure.json"
        readiness = tmpdir / "readiness.json"
        receipt = tmpdir / "receipt.json"
        payload = tmpdir / "payload.json"
        sequence_population = tmpdir / "sequence_population.json"
        receipt.write_text(RECEIPT.read_text(encoding="utf-8"), encoding="utf-8")
        payload.write_text(PAYLOAD.read_text(encoding="utf-8"), encoding="utf-8")
        sequence_population.write_text(SEQUENCE_POPULATION.read_text(encoding="utf-8"), encoding="utf-8")

        subprocess.run(
            [
                sys.executable,
                str(RUNNER),
                "--request-dir",
                str(bundle_dir),
                "--skip-backend-command",
                "--backend-export-output",
                str(backend_export),
                "--sequence-population",
                str(sequence_population),
                "--receipt",
                str(receipt),
                "--payload",
                str(payload),
                "--dump-output",
                str(dump),
                "--manifest-output",
                str(manifest),
                "--evaluation-output",
                str(evaluation),
                "--closure-output",
                str(closure),
                "--readiness-output",
                str(readiness),
                "--schedule-provenance",
                "pytest_existing_array_files",
            ],
            check=True,
            cwd=ROOT,
        )
        assert backend_export.exists()
        readiness_payload = json.loads(readiness.read_text(encoding="utf-8"))
        assert readiness_payload["publication_bundle_ready"] is False
        assert "custody receipt" in readiness_payload["smallest_backend_residual_object"]
