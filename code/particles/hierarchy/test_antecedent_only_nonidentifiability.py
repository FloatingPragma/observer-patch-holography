#!/usr/bin/env python3
"""Tests for the antecedent-only hierarchy/naturality no-go receipt (#518)."""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parent
MODULE_PATH = ROOT / "computations" / "antecedent_only_nonidentifiability.py"
RECEIPT_PATH = (
    ROOT / "certificates" / "antecedent_only_nonidentifiability_receipt.json"
)

spec = importlib.util.spec_from_file_location(
    "antecedent_only_nonidentifiability",
    MODULE_PATH,
)
module = importlib.util.module_from_spec(spec)
sys.modules["antecedent_only_nonidentifiability"] = module
spec.loader.exec_module(module)


def test_tracked_receipt_matches_deterministic_recomputation() -> None:
    tracked = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))

    assert tracked == module.build_receipt()
    assert tracked["target_artifacts_consumed"] == []
    assert tracked["pass"] is True


def test_capacity_countermodels_share_antecedents_and_have_distinct_outputs() -> None:
    block = module.build_receipt()["capacity_nonidentifiability"]
    models = block["countermodels"]

    assert len({model["antecedent_fingerprint"] for model in models}) == 1
    assert {model["fixed_log_capacity"] for model in models} == {"1", "3"}
    assert all(model["strict_contraction"] is True for model in models)
    assert all(model["self_map"] is True for model in models)
    assert all(model["unique_fixed_point"] is True for model in models)
    assert block["result"] == (
        "NO_UNIQUE_CAPACITY_PRODUCER_FROM_DECLARED_ANTECEDENTS"
    )


def test_capacity_completion_rejects_invalid_defining_antecedents() -> None:
    with pytest.raises(ValueError, match="interior"):
        module.capacity_completion(
            completion_id="bad-center",
            fixed_log_capacity=Fraction(4),
        )
    with pytest.raises(ValueError, match="averaging weight"):
        module.capacity_completion(
            completion_id="bad-weight",
            fixed_log_capacity=Fraction(2),
            averaging_weight=Fraction(0),
        )


def test_capacity_family_is_not_the_bridge_target_in_disguise() -> None:
    receipt = module.build_receipt()
    source = receipt["source_packet"]

    assert receipt["target_artifacts_consumed"] == []
    assert "fixed_log_capacity" not in source
    assert "N_CRC" not in source
    assert "B_EW(P,N)=0" in source["excluded_target_conditions"]
    assert receipt["promotion_boundary"]["not_promotable"][0] == (
        "N_CRC^EW as a physical capacity"
    )


def test_naturality_defects_are_evaluated_from_maps() -> None:
    block = module.build_receipt()["naturality_nonidentifiability"]
    models = block["countermodels"]

    assert [
        module.evaluate_naturality_completion(model)["epsilon_H"]
        for model in models
    ] == [0, 1, 1]
    assert models[1]["evaluated_defects"] == {
        "epsilon_n": 1,
        "epsilon_h": 0,
        "epsilon_H": 1,
    }
    assert models[2]["evaluated_defects"] == {
        "epsilon_n": 0,
        "epsilon_h": 1,
        "epsilon_H": 1,
    }


def test_naturality_evaluator_ignores_supplied_defect_booleans() -> None:
    completion = module.build_receipt()["naturality_nonidentifiability"][
        "countermodels"
    ][1]
    completion["evaluated_defects"] = {
        "epsilon_n": 0,
        "epsilon_h": 0,
        "epsilon_H": 0,
    }

    assert module.evaluate_naturality_completion(completion) == {
        "epsilon_n": 1,
        "epsilon_h": 0,
        "epsilon_H": 1,
    }


def test_naturality_countermodels_share_scalar_packet() -> None:
    receipt = module.build_receipt()
    models = receipt["naturality_nonidentifiability"]["countermodels"]
    source_hash = receipt["source_packet_sha256"]

    assert {model["antecedent_fingerprint"] for model in models} == {
        source_hash
    }
    assert receipt["source_packet"]["naturality_type_signatures"] == {
        "rho_sH": "Q_s -> Q_H",
        "n_s": "Q_s -> Q_s",
        "n_H": "Q_H -> Q_H",
        "h_s": "Q_s -> O_s",
        "chi_sH": "O_s -> O_H",
        "h_H": "Q_H -> O_H",
    }


def test_cli_check_accepts_tracked_receipt() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--check",
            "--output",
            str(RECEIPT_PATH),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
