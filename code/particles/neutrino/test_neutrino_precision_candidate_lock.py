#!/usr/bin/env python3
"""Tests for the immutable neutrino precision candidate lock."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib
import sys
import tempfile


HERE = pathlib.Path(__file__).resolve().parent
SCRIPT = HERE / "build_neutrino_precision_candidate_lock.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("oph_neutrino_candidate_lock", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_lock_is_fail_closed_hashed_and_immutable() -> None:
    module = _load_module()
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_candidate_lock_") as tmpdir:
        output = pathlib.Path(tmpdir) / "v1.0.0"
        manifest = module.build_bundle(output, "2026-07-11T04:20:00Z")
        assert manifest["prediction_promotion_allowed"] is False
        assert manifest["source_closure_passes"] is False
        assert manifest["candidate_rejected_by_nufit61"] is True

        source_dag = json.loads((output / "source_dag.json").read_text(encoding="utf-8"))
        predictions = json.loads((output / "predictions.json").read_text(encoding="utf-8"))
        adjudication = json.loads((output / "adjudication.json").read_text(encoding="utf-8"))
        assert source_dag["acyclic"] is True
        assert source_dag["target_audit"]["historical_target_exposure"] is True
        assert source_dag["source_closure_gate"]["passes"] is False
        assert predictions["absolute_scale"]["absolute_mass_prediction_allowed"] is False
        assert adjudication["baseline"]["candidate_rejected_by_declared_3sigma_gate"] is True
        assert adjudication["baseline"]["convention_rescue_found"] is False

        for filename, expected_hash in manifest["files"].items():
            assert _sha256(output / filename) == expected_hash

        try:
            module.build_bundle(output, "2026-07-11T04:20:00Z")
        except FileExistsError:
            pass
        else:
            raise AssertionError("an existing lock directory must never be overwritten")
