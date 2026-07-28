"""Regression and adversarial tests for the Workstream A action bundle."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))
sys.path.insert(0, str(ROOT / "checkers"))

import sm_eft_action_producer as producer  # noqa: E402
import check_sm_eft_action as checker  # noqa: E402


@pytest.fixture(scope="module")
def bundle() -> dict:
    return producer.build_bundle()


def test_bundle_is_deterministic(bundle: dict) -> None:
    again = producer.build_bundle()
    assert again == bundle


def test_emitted_file_matches_builder(bundle: dict) -> None:
    stored = json.loads(producer.OUT_PATH.read_text(encoding="utf-8"))
    assert stored == bundle


def test_checker_passes_on_emitted_bundle(bundle: dict) -> None:
    checker.check(bundle)


def test_census_and_anomalies(bundle: dict) -> None:
    census = bundle["field_census"]
    assert census["weyl_state_total"] == 45
    assert len(census["fermions"]) == 15
    assert len(census["ghosts"]) == 3


def test_controls_recorded_refusals(bundle: dict) -> None:
    for name, verdict in bundle["controls"].items():
        assert verdict["expected_failure"] is True, name
        assert verdict["failed"] is True, name


def test_numeric_yukawa_injection_fails_checker(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["yukawa_packet"]["matrices"]["Yu"][0][0] = {"value": "0.99"}
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_missing_generation_fails_checker(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["field_census"]["fermions"] = [
        f for f in doctored["field_census"]["fermions"] if f["generation"] != 3
    ]
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_gut_normalization_fails_checker(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["conventions"]["gprime_convention"] = "GUT_g1"
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_mass_core_fails_checker(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["mass_core"] = {"mW": "80.4"}
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_target_value_fails_numeric_scan(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["vev_types"]["v_chart"]["hint"] = 246
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_tampered_digest_fails_checker(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["subject_digest"] = "0" * 64
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_vev_types_stay_distinct(bundle: dict) -> None:
    doctored = json.loads(json.dumps(bundle))
    doctored["vev_types"]["v_F"]["type"] = doctored["vev_types"]["v_chart"]["type"]
    with pytest.raises(SystemExit):
        checker.check(doctored)


def test_producer_cli_round_trip(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "producers" / "sm_eft_action_producer.py")],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert result.returncode == 0
    assert "subject_digest" in result.stdout
