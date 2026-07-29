"""Fail-closed tests for the sampled W/Z boundary diagnostic."""

from __future__ import annotations

import copy
import json
import sys
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "checkers"))
sys.path.insert(0, str(ROOT / "producers"))

import check_wz_pole_diagnostic as checker  # noqa: E402
import wz_pole_receipts as producer  # noqa: E402

ARTIFACT = json.loads(
    (ROOT / "outputs" / "wz_pole_receipts.json").read_text(encoding="utf-8")
)


def rehash(payload: dict) -> None:
    payload.pop("artifact_sha256", None)
    payload["artifact_sha256"] = "sha256:" + producer.canonical_sha256(payload)


def test_committed_diagnostic_contract_passes_without_promotion() -> None:
    verdict = checker.check(ARTIFACT)
    assert verdict["status"] == checker.PASS_STATUS, verdict["problems"]
    assert verdict["root_count_certified"] is False
    assert verdict["acceptance_row_discharged"] is False
    assert verdict["promotion_allowed"] is False


def test_every_certification_flag_is_fixed_false() -> None:
    assert ARTIFACT["acceptance_row_discharged"] is False
    assert all(value is False for value in ARTIFACT["promotion"].values())
    assert all(value is False for value in ARTIFACT["certification"].values())
    for receipt in ARTIFACT["receipts"].values():
        for field in (
            "root_count_certified",
            "simple_root_certified",
            "analytic_sheet_certified",
            "precision_enclosure_certified",
            "laurent_data_certified",
            "physical_current_pole_certified",
        ):
            assert receipt[field] is False


def test_boundary_sampler_visits_all_four_edges_once() -> None:
    points = producer.boundary_sample_points(0 + 0j, 2 + 2j, 4)
    assert len(points) == 16
    assert Counter(edge for edge, _ in points) == {
        "bottom": 4,
        "right": 4,
        "top": 4,
        "left": 4,
    }
    assert len({point for _, point in points}) == 16
    assert points[0] == ("bottom", 0 + 0j)
    assert points[4] == ("right", 2 + 0j)
    assert points[8] == ("top", 2 + 2j)
    assert points[12] == ("left", 0 + 2j)


def test_dimensional_prefactor_corrections_replay_exactly() -> None:
    vector = json.loads(producer.VECTOR_PATH.read_text(encoding="utf-8"))
    names = {"W": "WpWm", "Z": "ZZ", "AZ": "AZ"}
    recomputed = {}
    independently_replayed = {}
    fixture = producer.fixture_substitutions()
    for label, block in names.items():
        expression = producer.block_transverse(vector, block)
        recomputed[label] = str(
            producer.dimensional_prefactor_finite_correction(expression)
        )

        substituted = sp.together(expression.subs(fixture))
        pole_coefficients = {}
        for symbol in substituted.free_symbols:
            if "__LB__" not in str(symbol):
                continue
            head, args = producer.loop_symbol_parts(symbol)
            evaluated = [sp.simplify(arg.subs(fixture)) for arg in args]
            if head == "A0":
                coefficient = evaluated[0]
            elif head in {"A0p", "B0"}:
                coefficient = sp.Integer(1)
            elif head.startswith("C"):
                coefficient = sp.Integer(0)
            else:
                raise AssertionError(f"unknown loop head {head}")
            pole_coefficients[symbol] = coefficient
        residue_as_function_of_d = sp.together(
            substituted.xreplace(pole_coefficients)
        )
        independent = sp.factor(
            -2
            * sp.diff(residue_as_function_of_d, producer.NS["d"]).subs(
                producer.NS["d"],
                4,
            )
        )
        independently_replayed[label] = str(independent)

    assert recomputed == checker.EXPECTED_CORRECTIONS
    assert independently_replayed == checker.EXPECTED_CORRECTIONS


def test_checker_rejects_promoted_certification_even_with_valid_self_digest() -> None:
    mutated = copy.deepcopy(ARTIFACT)
    mutated["certification"]["root_count"] = True
    rehash(mutated)
    verdict = checker.check(mutated)
    assert verdict["status"] == "FAIL"
    assert any("root_count" in problem for problem in verdict["problems"])


def test_checker_rejects_source_digest_tamper_even_with_valid_self_digest() -> None:
    mutated = copy.deepcopy(ARTIFACT)
    mutated["input"]["sha256"] = "0" * 64
    rehash(mutated)
    verdict = checker.check(mutated)
    assert verdict["status"] == "FAIL"
    assert any("source vector" in problem for problem in verdict["problems"])


def test_checker_rejects_exact_correction_tamper_even_with_valid_self_digest() -> None:
    mutated = copy.deepcopy(ARTIFACT)
    mutated["conventions"]["dimensional_prefactor_finite_correction"]["W"] = "0"
    for key, receipt in mutated["receipts"].items():
        if key.startswith("W_"):
            receipt["dimensional_prefactor_finite_correction"] = "0"
    rehash(mutated)
    verdict = checker.check(mutated)
    assert verdict["status"] == "FAIL"
    assert any("finite correction" in problem for problem in verdict["problems"])


def test_checker_rejects_sample_summary_tamper_even_with_valid_self_digest() -> None:
    mutated = copy.deepcopy(ARTIFACT)
    mutated["receipts"]["Z_256"]["sampled_boundary_winding"] = 2
    rehash(mutated)
    verdict = checker.check(mutated)
    assert verdict["status"] == "FAIL"
    assert any("sampled winding" in problem for problem in verdict["problems"])
