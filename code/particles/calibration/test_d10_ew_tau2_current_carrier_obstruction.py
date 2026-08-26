#!/usr/bin/env python3
"""Guard the D10 current-carrier tau2 obstruction artifact."""

from __future__ import annotations

import json
import pathlib
from copy import deepcopy
from fractions import Fraction

import pytest

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
    assert provenance["source_pair"]["declared_path"] == custom_source.resolve().as_posix()
    assert provenance["population"]["declared_path"] == custom_population.resolve().as_posix()
    assert provenance["fiberwise_tree_law"]["declared_path"] == custom_fiberwise.resolve().as_posix()
    assert provenance["reference_entries"]["declared_path"] == custom_references.resolve().as_posix()


@pytest.mark.parametrize(
    ("field", "bad_value"),
    [
        ("artifact", "wrong_fiber_artifact"),
        ("status", "retracted"),
        ("proof_status", "unproved"),
        ("coordinate_symbol", "wrong_coordinate"),
        ("fiber_population_functional_formula", "0"),
        ("tauY_formula", "0"),
        ("n_EW_formula", "1"),
        ("MW_formula", "0"),
        ("MZ_formula", "0"),
    ],
)
def test_corrupt_fiber_contract_is_rejected(field: str, bad_value: object) -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted = deepcopy(fiberwise)
    corrupted[field] = bad_value
    with pytest.raises(ValueError, match=f"fiberwise_tree_law.{field}"):
        lane.build_artifact(source_pair, population, corrupted, references)


def test_corrupt_fiber_scalars_and_anchor_are_rejected() -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted_scalar = deepcopy(fiberwise)
    corrupted_scalar["carrier_basis_scalar"]["alpha2_mz"] = 999
    with pytest.raises(ValueError, match="carrier_basis_scalar.alpha2_mz"):
        lane.build_artifact(source_pair, population, corrupted_scalar, references)

    corrupted_anchor = deepcopy(fiberwise)
    corrupted_anchor["anchor_point"]["tau2_tree_exact"] = 0.75
    with pytest.raises(ValueError, match="anchor_point.tau2_tree_exact"):
        lane.build_artifact(source_pair, population, corrupted_anchor, references)


@pytest.mark.parametrize(
    ("field", "bad_value"),
    [
        ("eta_EW", 0.75),
        ("sigma_EW", 0.75),
        ("tau_2", 0.75),
        ("tau_Y", 0.75),
    ],
)
def test_incompatible_selected_population_point_is_rejected(
    field: str,
    bad_value: object,
) -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted = deepcopy(population)
    corrupted["selected_population_point"][field] = bad_value
    with pytest.raises(ValueError, match=f"selected_population_point.{field}"):
        lane.build_artifact(source_pair, corrupted, fiberwise, references)


def test_wrong_population_and_source_receipts_are_rejected() -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted_population = deepcopy(population)
    corrupted_population["artifact"] = "wrong_population_artifact"
    with pytest.raises(ValueError, match="population.artifact"):
        lane.build_artifact(source_pair, corrupted_population, fiberwise, references)

    corrupted_population_status = deepcopy(population)
    corrupted_population_status["status"] = "open"
    with pytest.raises(ValueError, match="population.status"):
        lane.build_artifact(
            source_pair,
            corrupted_population_status,
            fiberwise,
            references,
        )

    corrupted_source = deepcopy(source_pair)
    corrupted_source["artifact"] = "wrong_source_artifact"
    with pytest.raises(ValueError, match="source_pair.artifact"):
        lane.build_artifact(corrupted_source, population, fiberwise, references)

    corrupted_source_status = deepcopy(source_pair)
    corrupted_source_status["status"] = "open"
    with pytest.raises(ValueError, match="source_pair.status"):
        lane.build_artifact(
            corrupted_source_status,
            population,
            fiberwise,
            references,
        )


def test_incompatible_source_and_population_scalars_are_rejected() -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted_source = deepcopy(source_pair)
    corrupted_source["source_slots"]["alpha2_mz"] = 999
    with pytest.raises(ValueError, match="source_pair.source_slots.alpha2_mz"):
        lane.build_artifact(corrupted_source, population, fiberwise, references)

    corrupted_basis = deepcopy(population)
    corrupted_basis["selected_population_basis_point"]["n_EW"] = 0.5
    with pytest.raises(ValueError, match="selected_population_basis_point.n_EW"):
        lane.build_artifact(source_pair, corrupted_basis, fiberwise, references)


def test_corrupt_selector_and_mass_map_formulas_are_rejected() -> None:
    source_pair, population, fiberwise, references = _inputs()
    corrupted_source_formula = deepcopy(source_pair)
    corrupted_source_formula["population_basis"]["n_EW_formula"] = "1"
    with pytest.raises(ValueError, match="source_pair.population_basis.n_EW_formula"):
        lane.build_artifact(
            corrupted_source_formula,
            population,
            fiberwise,
            references,
        )

    corrupted_selector = deepcopy(population)
    corrupted_selector["population_functional_status"] = "open"
    with pytest.raises(ValueError, match="population.population_functional_status"):
        lane.build_artifact(source_pair, corrupted_selector, fiberwise, references)

    corrupted_mass_map = deepcopy(population)
    corrupted_mass_map["population_atomic_quartet"]["mZ_formula"] = "0"
    with pytest.raises(ValueError, match="population_atomic_quartet.mZ_formula"):
        lane.build_artifact(source_pair, corrupted_mass_map, fiberwise, references)


def test_canonical_receipts_record_fail_closed_input_binding() -> None:
    source_pair, population, fiberwise, references = _inputs()
    payload = lane.build_artifact(source_pair, population, fiberwise, references)
    validation = payload["input_contract_validation"]
    assert validation["status"] == "PASS"
    assert validation["binding"] == "fail_closed_formula_scalar_anchor_validation"
    assert validation["canonical_eta_source"] == "fiberwise_tree_law.eta_source"
    assert set(validation["validated_formula_fields"]) == set(
        lane.EXPECTED_FIBERWISE_FORMULAS
    )
    contracts = validation["validated_formula_contracts"]
    assert set(contracts["source_pair"]) == set(lane.EXPECTED_SOURCE_PAIR_FIELDS)
    assert set(contracts["population"]) == set(lane.EXPECTED_POPULATION_FIELDS)
    assert set(contracts["population_basis"]) == set(
        lane.EXPECTED_POPULATION_BASIS_FORMULAS
    )
    assert set(contracts["population_mass_map"]) == set(
        lane.EXPECTED_POPULATION_MASS_FORMULAS
    )
