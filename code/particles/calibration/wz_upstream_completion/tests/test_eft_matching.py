"""Regression and adversarial tests for the Workstream B matching bundle."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))
sys.path.insert(0, str(ROOT / "checkers"))

import eft_matching_producer as producer  # noqa: E402
import check_eft_matching as checker  # noqa: E402


@pytest.fixture(scope="module")
def bundle() -> dict:
    return producer.build_bundle()


def test_bundle_is_deterministic(bundle: dict) -> None:
    assert producer.build_bundle() == bundle


def test_emitted_file_matches_builder(bundle: dict) -> None:
    stored = json.loads(producer.OUT_PATH.read_text(encoding="utf-8"))
    assert stored == bundle


def test_checker_passes() -> None:
    checker.check()


def test_census_derived_coefficients(bundle: dict) -> None:
    coefficients = bundle["gauge_betas"]["coefficients"]
    assert Fraction(coefficients["b1"]) == Fraction(41, 6)
    assert Fraction(coefficients["b2"]) == Fraction(-19, 6)
    assert Fraction(coefficients["b3"]) == Fraction(-7)


def test_binding_to_action_digest(bundle: dict) -> None:
    action = json.loads(producer.ACTION_PATH.read_text(encoding="utf-8"))
    assert bundle["action_subject_digest"] == action["subject_digest"]


def test_controls_recorded_refusals(bundle: dict) -> None:
    for name, verdict in bundle["controls"].items():
        assert verdict["expected_failure"] is True, name
        assert verdict["failed"] is True, name


def test_rederivation_rejects_a_doctored_census() -> None:
    action = json.loads(producer.ACTION_PATH.read_text(encoding="utf-8"))
    census = json.loads(json.dumps(action["field_census"]))
    census["fermions"] = [f for f in census["fermions"] if f["generation"] != 3]
    derived = checker.rederive_betas(census)
    assert derived["b1"] != Fraction(41, 6)
    assert derived["b3"] != Fraction(-7)


def test_scalar_weight_matters() -> None:
    action = json.loads(producer.ACTION_PATH.read_text(encoding="utf-8"))
    census = json.loads(json.dumps(action["field_census"]))
    census["scalars"] = []
    derived = checker.rederive_betas(census)
    assert derived["b2"] == Fraction(-19, 6) - Fraction(1, 6)


def test_output_vector_is_fully_symbolic(bundle: dict) -> None:
    for component in bundle["matching"]["output_vector"]["components"]:
        assert set(component) == {"symbol"}


def test_decoupling_is_recorded_empty_with_reason(bundle: dict) -> None:
    interval = bundle["matching"]["intervals"][0]
    assert interval["decoupling_maps"]["maps"] == []
    assert "recorded statement" in interval["decoupling_maps"]["reason"]
