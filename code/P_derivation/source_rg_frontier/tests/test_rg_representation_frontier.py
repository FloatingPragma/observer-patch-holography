#!/usr/bin/env python3
"""Regression and adversarial tests for the issue #32 source frontier."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import pytest


PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))

import build_rg_representation_frontier as producer  # noqa: E402
import check_rg_representation_frontier as checker  # noqa: E402


OUTPUT = PACKAGE / "outputs" / "rg_representation_frontier.json"


def load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def rehash(payload: dict) -> dict:
    payload = deepcopy(payload)
    payload.pop("subject_digest", None)
    payload["subject_digest"] = checker.canonical_sha256(payload)
    return payload


def must_fail(payload: dict, code: str) -> None:
    with pytest.raises(SystemExit, match=rf"^{code}:"):
        checker.check(rehash(payload))


def test_clean_rebuild_is_byte_exact_and_independently_resolves() -> None:
    subprocess.run(
        [sys.executable, str(PACKAGE / "build_rg_representation_frontier.py"), "--check-byte-exact"],
        cwd=PACKAGE,
        check=True,
    )
    checker.check(load())


def test_exact_parametric_law_and_declared_evaluation() -> None:
    payload = load()
    assert payload["parametric_one_loop_law"]["coefficients"] == {
        "b_Y": "(20/9) N_g + (1/6) N_H",
        "b_2": "-22/3 + (4/3) N_g + (1/6) N_H",
        "b_3": "-11 + (4/3) N_g",
    }
    row = payload["conditional_evaluations"][0]
    assert row["coefficients"] == {"b_Y": "41/6", "b_2": "-19/6", "b_3": "-7"}
    assert row["promotion_allowed"] is False
    assert payload["promotion_allowed"] is False


def test_zero_gauge_index_boundary_does_not_overclaim_full_decoupling() -> None:
    boundary = load()["invisible_sector_boundary"]
    assert boundary["exact_result"]["delta_b_Y"] == "0"
    assert boundary["exact_result"]["delta_b_2"] == "0"
    assert boundary["exact_result"]["delta_b_3"] == "0"
    assert boundary["full_WZ_decoupling_proved"] is False
    assert boundary["missing_for_full_WZ_decoupling"]


def test_full_wz_decoupling_self_attestation_is_rejected() -> None:
    payload = load()
    payload["invisible_sector_boundary"]["full_WZ_decoupling_proved"] = True
    must_fail(payload, "STERILE_OVERPROMOTION")


def test_census_fiber_remains_explicit_and_unselected() -> None:
    payload = load()
    rows = payload["nonidentifiability_witnesses"]["census_nonuniqueness"]["witnesses"]
    assert [(row["N_g"], row["N_H"]) for row in rows] == [
        (3, 0),
        (3, 1),
        (3, 2),
        (4, 1),
        (5, 1),
    ]
    assert len({tuple(row["coefficients"].values()) for row in rows}) == 5
    assert all(row["physical_selection_status"] == "not_source_selected" for row in rows)


def test_external_validation_packet_in_source_ancestry_is_rejected() -> None:
    payload = load()
    row = payload["source_inputs"]["artifacts"][0]
    row["path"] = "code/particles/calibration/wz_upstream_completion/outputs/eft_matching_1.json"
    must_fail(payload, "SOURCE_PATH")


def test_source_pin_mutation_is_rejected() -> None:
    payload = load()
    payload["source_inputs"]["artifacts"][0]["byte_sha256"] = "0" * 64
    must_fail(payload, "PIN_HASH")


def test_beta_coefficient_mutation_is_rejected() -> None:
    payload = load()
    payload["conditional_evaluations"][0]["coefficients"]["b_2"] = "-3"
    must_fail(payload, "CONDITIONAL_BETA")


def test_family_count_promotion_is_rejected() -> None:
    payload = load()
    payload["conditional_evaluations"][0]["status"] = "source_selected"
    must_fail(payload, "CONDITIONAL_STATUS")


def test_scalar_count_promotion_is_rejected() -> None:
    payload = load()
    payload["representation_indices"]["conditional_scalar_doublet"]["status"] = "source_selected"
    must_fail(payload, "SCALAR_STATUS")


def test_hypercharge_normalization_mutation_is_rejected() -> None:
    payload = load()
    payload["parametric_one_loop_law"]["hypercharge_normalization"] = "GUT normalized"
    must_fail(payload, "HYPERCHARGE_NORMALIZATION")


def test_sign_convention_formula_mutation_is_rejected() -> None:
    payload = load()
    payload["parametric_one_loop_law"]["coefficients"]["b_3"] = "11 - (4/3) N_g"
    must_fail(payload, "BETA_LAW")


def test_matching_self_attestation_is_rejected() -> None:
    payload = load()
    payload["matching_objects"]["threshold_locations"] = {
        "status": "emitted",
        "reason": "caller asserted",
    }
    must_fail(payload, "SCHEMA")


def test_acceptance_overclaim_is_rejected() -> None:
    payload = load()
    payload["acceptance_map"][0]["status"] = "complete"
    must_fail(payload, "ACCEPTANCE_OVERCLAIM")


def test_scheme_witness_collapse_is_rejected() -> None:
    payload = load()
    schemes = payload["nonidentifiability_witnesses"]["coordinate_nonuniqueness"][
        "finite_scheme_redefinitions"
    ]
    schemes[1]["c"] = "0"
    must_fail(payload, "SCHEME_WITNESSES")


def test_lean_binding_mutation_is_rejected() -> None:
    payload = load()
    payload["formal_certificate"]["byte_sha256"] = "f" * 64
    must_fail(payload, "LEAN_HASH")


def test_producer_formula_has_no_target_argument() -> None:
    assert producer.beta_coefficients(3, 1) == {
        "b_Y": "41/6",
        "b_2": "-19/6",
        "b_3": "-7",
    }
    with pytest.raises(producer.FrontierError, match=r"NEGATIVE_MULTIPLICITY"):
        producer.beta_coefficients(-1, 1)
