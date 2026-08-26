#!/usr/bin/env python3
"""Guard the D10 current-carrier tau2 obstruction artifact."""

from __future__ import annotations

import json
import pathlib
from fractions import Fraction

import derive_d10_ew_tau2_current_carrier_obstruction as lane


ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"


def _inputs() -> tuple[dict, dict, dict, dict]:
    source_pair = json.loads(lane.DEFAULT_SOURCE_PAIR.read_text(encoding="utf-8"))
    population = json.loads(lane.DEFAULT_POPULATION.read_text(encoding="utf-8"))
    fiberwise = json.loads(lane.DEFAULT_FIBERWISE_TREE_LAW.read_text(encoding="utf-8"))
    references = json.loads(lane.REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    return source_pair, population, fiberwise, references


def test_d10_tau2_current_carrier_obstruction_is_emitted() -> None:
    source_pair, population, fiberwise, references = _inputs()
    rebuilt = lane.build_artifact(source_pair, population, fiberwise, references)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload == rebuilt

    assert payload["artifact"] == "oph_d10_ew_tau2_current_carrier_obstruction"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["proof_status"] == "exact_interval_excludes_single_tau2_central_WZ_pair_on_current_carrier"
    assert payload["next_single_residual_object"] == "delta_n_tree_exact"
    assert payload["diagnostic_only"] is True
    assert payload["exactness_convention"]["status"] == "conditional_exact_diagnostic"
    assert "finite-decimal" in payload["exactness_convention"]["upstream_decimal_interpretation"]
    assert all(
        len(record["payload_sha256"]) == 64
        for name, record in payload["input_provenance"].items()
        if name != "digest_scope"
    )
    direction = payload["direction_obstruction"]
    assert direction["single_tau2_possible_at_first_order"] is False
    assert direction["global_no_go_inferred_from_germ"] is False
    assert direction["germ_coefficient_W"] > 0.0
    assert direction["germ_coefficient_Z"] > 0.0
    assert (direction["tau2_required_for_W_first_order"] > 0.0) != (direction["tau2_required_for_Z_first_order"] > 0.0)
    nonlinear = payload["closed_form_nonlinear_point_test"]
    assert nonlinear["obstruction_established"] is True
    assert nonlinear["simultaneous_exact_central_pair_possible"] is False
    assert nonlinear["MZ_residual_gev"] < -0.019
    assert "display-only" in nonlinear["floating_point_values_role"]
    exact = nonlinear["exact_interval_certificate"]
    assert exact["arithmetic"] == "exact_rational_interval_arithmetic"
    assert exact["residual_squared_sign"] == "strictly_negative"
    assert exact["zero_excluded"] is True
    upper = exact["MZ_squared_residual_interval_gev2"]["upper"]
    assert Fraction(int(upper["numerator"]), int(upper["denominator"])) < 0
    assert nonlinear["nonlinear_routes_left_open"]
    distance = payload["reference_distance"]
    assert distance["W_offset_sigma"] > 0.0
    assert distance["Z_offset_sigma"] < 0.0


def test_wide_valid_pi_enclosure_prevents_false_no_go() -> None:
    source_pair, population, fiberwise, references = _inputs()
    payload = lane.build_artifact(
        source_pair,
        population,
        fiberwise,
        references,
        pi_bounds=(Fraction(3), Fraction(4)),
    )
    assert payload["status"] == "obstruction_not_established_interval_contains_zero"
    assert payload["proof_status"] == "exact_interval_does_not_exclude_single_tau2_central_WZ_pair"
    assert payload["next_single_residual_object"] is None
    nonlinear = payload["closed_form_nonlinear_point_test"]
    assert nonlinear["simultaneous_exact_central_pair_possible"] is None
    exact = nonlinear["exact_interval_certificate"]
    assert exact["residual_squared_sign"] == "undetermined_interval_contains_zero"
    assert exact["zero_excluded"] is False


def test_custom_input_paths_do_not_masquerade_as_canonical(tmp_path: pathlib.Path) -> None:
    source_pair, population, fiberwise, references = _inputs()
    custom_source = tmp_path / "source.json"
    custom_population = tmp_path / "population.json"
    custom_fiberwise = tmp_path / "fiberwise.json"
    custom_references = tmp_path / "references.json"
    payload = lane.build_artifact(
        source_pair,
        population,
        fiberwise,
        references,
        source_pair_path=custom_source,
        population_path=custom_population,
        fiberwise_tree_law_path=custom_fiberwise,
        reference_path=custom_references,
    )
    provenance = payload["input_provenance"]
    assert provenance["source_pair"]["declared_path"] == str(custom_source.resolve())
    assert provenance["population"]["declared_path"] == str(custom_population.resolve())
    assert provenance["fiberwise_tree_law"]["declared_path"] == str(custom_fiberwise.resolve())
    assert provenance["reference_entries"]["declared_path"] == str(custom_references.resolve())
