#!/usr/bin/env python3
"""Tests for raw hadron backend export bundle support."""

from __future__ import annotations

import json
import pathlib
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import h5py
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover - optional dependency path
    h5py = None
    np = None


ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.backend_export_bundle import (
    build_backend_export_skeleton,
    load_backend_input_artifact,
)
from particles.hadron.production_execution_support import build_backend_manifest, build_production_dump


class BackendExportBundleTest(unittest.TestCase):
    def setUp(self) -> None:
        self.receipt = {
            "artifact": "runtime_schedule_receipt_N_therm_and_N_sep",
            "execution_contract": {
                "ensemble_schedule": [
                    {
                        "ensemble_id": "qcd_2p1_seed_n0",
                        "cfg_ids": ["qcd_2p1_seed_n0__cfg0", "qcd_2p1_seed_n0__cfg1"],
                        "trajectory_stop_by_cfg": {
                            "qcd_2p1_seed_n0__cfg0": 2048,
                            "qcd_2p1_seed_n0__cfg1": 2560,
                        },
                    }
                ]
            },
        }
        self.payload = {
            "ensemble_payloads": [
                {
                    "ensemble_id": "qcd_2p1_seed_n0",
                    "family_index": 0,
                    "beta": 6.0,
                    "L": 8,
                    "T": 4,
                    "aLambda_msbar3": 0.1,
                    "am_l": 0.01,
                    "am_s": 0.02,
                    "source_descriptors_by_cfg": {
                        "qcd_2p1_seed_n0__cfg0": [
                            {"src_id": "src0", "coords": [0, 0, 0, 0]},
                            {"src_id": "src1", "coords": [4, 4, 4, 2]},
                        ],
                        "qcd_2p1_seed_n0__cfg1": [
                            {"src_id": "src0", "coords": [0, 0, 0, 0]},
                            {"src_id": "src1", "coords": [4, 4, 4, 2]},
                        ],
                    },
                }
            ]
        }

    def test_build_backend_export_skeleton_uses_receipt_schedule_and_payload_coords(self) -> None:
        manifest, datasets = build_backend_export_skeleton(self.receipt, self.payload)
        self.assertEqual(manifest["artifact"], "oph_hadron_backend_raw_export")
        self.assertEqual(manifest["profile_id"], "oph_reference_backend_v1")
        self.assertEqual(manifest["ensembles"][0]["cfgs"][0]["trajectory_stop"], 2048)
        self.assertEqual(
            manifest["ensembles"][0]["cfgs"][0]["sources"][1]["coord"],
            [4, 4, 4, 2],
        )
        self.assertEqual(len(datasets), 12)
        self.assertEqual(datasets[0], ("/ensembles/qcd_2p1_seed_n0/cfgs/qcd_2p1_seed_n0__cfg0/sources/src0/pi_iso", 4))

    def test_build_backend_export_skeleton_rejects_bad_source_coord(self) -> None:
        bad_payload = json.loads(json.dumps(self.payload))
        bad_payload["ensemble_payloads"][0]["source_descriptors_by_cfg"]["qcd_2p1_seed_n0__cfg0"][0]["coords"] = [0, 0, 0]
        with self.assertRaisesRegex(ValueError, "must contain exactly 4 coordinates"):
            build_backend_export_skeleton(self.receipt, bad_payload)

    def _make_bundle(self, root: Path) -> Path:
        bundle_dir = root / "bundle"
        bundle_dir.mkdir(parents=True, exist_ok=True)
        manifest, datasets = build_backend_export_skeleton(self.receipt, self.payload)
        manifest["backend"] = {
            "family": "rhmc_hmc",
            "name": "ref_backend",
            "version": "1.0",
            "git_commit": "abc123",
            "run_id": "run-001",
            "build_id": "build-001",
            "machine": "test-node",
        }
        (bundle_dir / "backend_run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with h5py.File(bundle_dir / "correlators.h5", "w") as h5:
            for offset, (dset_path, length) in enumerate(datasets):
                h5.create_dataset(dset_path, data=np.linspace(1 + offset, length + offset, length, dtype=np.float64))
        return bundle_dir

    @unittest.skipUnless(h5py is not None and np is not None, "h5py and numpy are required for HDF5 roundtrip test")
    def test_raw_bundle_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle_dir = self._make_bundle(Path(tmp))
            backend_input = load_backend_input_artifact(bundle_dir)
            self.assertEqual(backend_input["artifact"], "oph_hadron_backend_raw_export_inlined")
            dump = build_production_dump(self.receipt, self.payload, backend_input)
            manifest = build_backend_manifest(
                self.receipt,
                self.payload,
                backend_input,
                backend_input_path=str(bundle_dir),
            )
            src0_pi = dump["ensembles"]["qcd_2p1_seed_n0"]["cfgs"]["qcd_2p1_seed_n0__cfg0"]["sources"]["src0"]["pi_iso"]
            self.assertEqual(src0_pi, [1.0, 2.0, 3.0, 4.0])
            self.assertEqual(manifest["backend_name"], "ref_backend")
            self.assertEqual(manifest["backend_run_id"], "run-001")
            self.assertIn("raw_export_provenance", manifest)


if __name__ == "__main__":
    unittest.main()
