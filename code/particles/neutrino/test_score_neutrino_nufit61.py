#!/usr/bin/env python3
"""Tests for the NuFIT 6.1 correlated-profile scorer."""

from __future__ import annotations

import importlib.util
import copy
import json
import math
import re
from fractions import Fraction
import pathlib
import sys

import pytest


HERE = pathlib.Path(__file__).resolve().parent
SCRIPT = HERE / "score_neutrino_nufit61.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("oph_nufit61_scorer", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_bilinear_interpolation_and_out_of_grid_failure() -> None:
    module = _load_module()
    grid = module.RectilinearGrid.from_rows(
        [
            (0.0, 0.0, 0.0),
            (0.0, 2.0, 4.0),
            (2.0, 0.0, 2.0),
            (2.0, 2.0, 6.0),
        ],
        "synthetic",
    )
    assert math.isclose(grid.interpolate(1.0, 1.0)["delta_chi2"], 3.0, abs_tol=1.0e-15)
    try:
        grid.interpolate(-0.1, 1.0)
    except ValueError as exc:
        assert "extrapolation is forbidden" in str(exc)
    else:
        raise AssertionError("out-of-grid interpolation must fail")


def test_candidate_coordinate_transform_and_phase_wrap() -> None:
    module = _load_module()
    candidate = {
        "pmns_observables": {
            "theta12_deg": 30.0,
            "theta13_deg": 10.0,
            "theta23_deg": 45.0,
            "delta_deg": 305.0,
        },
        "dimensionless_ratio_dm21_over_dm32": 0.03,
    }
    coordinates = module._candidate_coordinates(candidate)
    assert math.isclose(coordinates["sin2_theta12"], 0.25, abs_tol=1.0e-15)
    assert math.isclose(coordinates["sin2_theta23"], 0.5, abs_tol=1.0e-15)
    assert math.isclose(coordinates["delta_cp_deg_wrapped"], -55.0, abs_tol=1.0e-15)


def test_template_kernel_blocks_prediction_promotion() -> None:
    module = _load_module()
    boundary = module._source_boundary(
        {"status": "template", "proof_status": "conjugacy_riesz_candidate"}
    )
    assert boundary["source_only_prediction_eligible"] is False
    assert boundary["prospective_evidence_eligible"] is False
    assert boundary["historical_target_exposure"] is True


def test_source_emitted_kernel_is_necessary_but_history_stays_retrospective() -> None:
    module = _load_module()
    boundary = module._source_boundary(
        {"status": "source_only_frozen", "proof_status": "closed_source_emitted"}
    )
    assert boundary["source_only_prediction_eligible"] is True
    assert boundary["prospective_evidence_eligible"] is False


def _valley(module, minimum=11.6):
    xs = (math.log10(0.0005), math.log10(0.001))
    ys = (1.0, 1.5004, 1.5005, 2.0)
    heights = (minimum + 10000, minimum, minimum + 10000, minimum + 10000)
    return module.RectilinearGrid.from_rows(
        [(x, y, z) for x in xs for y, z in zip(ys, heights)], "narrow valley")


def test_unsampled_narrow_valley_cannot_be_falsely_rejected() -> None:
    module = _load_module()
    grid = _valley(module)
    sampled = module._ratio_profile(grid, 1.0, 1001)
    bracket = module._ratio_profile_bracket(grid, 1.0)
    threshold = module.THREE_SIGMA_TWO_DOF_DELTA_CHI2
    # Exact independent oracle: the interpolant is independent of x and has
    # its minimum at the y=1.5004 vertex, which lies on the continuous curve.
    assert sampled["delta_chi2"] > threshold
    assert bracket["lower_bound"] <= 11.6 <= bracket["upper_bound"] < threshold
    assert bracket["upper_bound"] - bracket["lower_bound"] < 1.01e-8


def test_production_score_uses_the_bracket_lower_endpoint(monkeypatch) -> None:
    module = _load_module()
    angular = module.RectilinearGrid.from_rows(
        [(x, y, 0.0) for x in (0.0, 1.0) for y in (-180.0, 180.0)], "angular")
    monkeypatch.setattr(module, "_verify_table", lambda *args: {"synthetic": True})
    monkeypatch.setattr(module, "_read_grids", lambda *args: {
        "T13/T12": angular, "T23/DCP": angular, "DMS/DMA": _valley(module)})
    result = module.score_table(pathlib.Path("unused"), {}, "synthetic", {
        "sin2_theta13": .02, "sin2_theta12": .3, "sin2_theta23": .5,
        "delta_cp_deg_wrapped": 0.0, "ratio_dm21_over_dm32": 1.0}, 1001)
    curve = result["profiles"]["DMS/DMA_ratio_profile"]
    assert curve["delta_chi2"] > module.THREE_SIGMA_TWO_DOF_DELTA_CHI2
    assert result["joint_fixed_candidate_delta_chi2_lower_bound"] == curve["certified_profile_bracket"]["lower_bound"]
    assert result["passes_published_3sigma_2d_compatibility"]
    assert result["published_3sigma_2d_comparison_status"] == "COMPATIBLE_PROFILES"


def test_threshold_straddling_interval_is_unresolved(monkeypatch) -> None:
    module = _load_module()
    threshold = module.THREE_SIGMA_TWO_DOF_DELTA_CHI2
    angular = module.RectilinearGrid.from_rows(
        [(x, y, 0.0) for x in (0.0, 1.0) for y in (-180.0, 180.0)], "angular")
    monkeypatch.setattr(module, "_verify_table", lambda *args: {})
    monkeypatch.setattr(module, "_read_grids", lambda *args: {
        "T13/T12": angular, "T23/DCP": angular, "DMS/DMA": _valley(module)})
    monkeypatch.setattr(module, "_ratio_profile_bracket", lambda *args: {
        "lower_bound": threshold-1e-10, "upper_bound": threshold+1e-10})
    result = module.score_table(pathlib.Path("unused"), {}, "synthetic", {
        "sin2_theta13": .02, "sin2_theta12": .3, "sin2_theta23": .5,
        "delta_cp_deg_wrapped": 0.0, "ratio_dm21_over_dm32": 1.0}, 1001)
    assert result["published_3sigma_2d_comparison_status"] == "UNRESOLVED_INTERVAL"
    assert result["passes_published_3sigma_2d_compatibility"] is None


def test_cli_does_not_turn_unresolved_treatments_into_rejection(monkeypatch, tmp_path) -> None:
    module = _load_module()
    monkeypatch.setattr(module, "score_table", lambda *args: {
        "profiles": {"T23/DCP": {"delta_chi2": 0.0}},
        "joint_fixed_candidate_delta_chi2_lower_bound": 0.0,
        "passes_published_3sigma_2d_compatibility": None,
        "published_3sigma_2d_comparison_status": "UNRESOLVED_INTERVAL"})
    output = tmp_path / "unresolved.json"
    monkeypatch.setattr(sys, "argv", ["score", "--tb-off-no", "unused", "--tb-yes-no", "unused",
                                     "--output", str(output)])
    assert module.main() == 0
    payload = json.loads(output.read_bytes())
    assert payload["decision"]["current_weighted_cycle_candidate_rejected_by_declared_gate"] is False


def test_actual_profile_exclusion_remains_available() -> None:
    module = _load_module()
    bracket = module._ratio_profile_bracket(_valley(module, 20.0), 1.0)
    assert bracket["lower_bound"] > module.THREE_SIGMA_TWO_DOF_DELTA_CHI2
    assert bracket["lower_bound"] <= 20.0 <= bracket["upper_bound"]


def test_interior_bilinear_logarithmic_minimum_has_exact_sign_oracle() -> None:
    module = _load_module()
    # On r=1, x=log10(y/2000). Thus x*(y-2000)>=0 everywhere,
    # with exact minimum zero at the interior point y=2000, x=0.
    grid = module.RectilinearGrid.from_rows(
        [(x, y, x * (y-2000)) for x in (-1.0, 1.0) for y in (1000.0, 3000.0)], "interior")
    result = module._ratio_profile_bracket(grid, 1.0)
    assert result["lower_bound"] <= 0 <= result["upper_bound"]
    assert result["upper_bound"] - result["lower_bound"] < 1.01e-8


@pytest.mark.parametrize("power", [-100, -4, 0, 3, 100])
def test_directed_log10_encloses_exact_power_identity(power) -> None:
    module = _load_module()
    lo, hi = module._log10_interval(Fraction(10) ** power)
    assert lo <= power <= hi
    assert hi-lo < Fraction(1, 10**65)


def test_directed_float_serialization_never_rounds_inwards() -> None:
    module = _load_module()
    for value in (Fraction(1, 10), Fraction(-1, 10), Fraction(1, 3), Fraction(10**20, 3)):
        assert Fraction(module._outward_float(value, False)) <= value
        assert Fraction(module._outward_float(value, True)) >= value


def test_angular_point_bounds_use_exact_weights_and_outward_rounding() -> None:
    module = _load_module()
    grid = module.RectilinearGrid.from_rows(
        [(x, y, x + 2*y + x*y) for x in (0., 2.) for y in (0., 3.)], "exact angular")
    x, y = .1, .3
    result = module._bilinear_point_bracket(grid, x, y)
    exact = Fraction(x) + 2*Fraction(y) + Fraction(x)*Fraction(y)
    assert Fraction(result["lower_bound"]) <= exact <= Fraction(result["upper_bound"])


def test_curve_clipping_preserves_valid_minimum_and_rejects_empty_domain() -> None:
    module = _load_module()
    # Allowed x=log10(y/2000) lies in [0,1] only when y>=2000.
    grid = module.RectilinearGrid.from_rows(
        [(x, y, y) for x in (0.0, 1.0) for y in (1000.0, 3000.0)], "clipped")
    result = module._ratio_profile_bracket(grid, 1.0)
    assert result["lower_bound"] <= 2000 <= result["upper_bound"]
    assert result["upper_bound"] - result["lower_bound"] < 1.01e-8
    absent = module.RectilinearGrid.from_rows(
        [(x, y, 0.0) for x in (1.0, 2.0) for y in (1.0, 2.0)], "absent")
    with pytest.raises(ValueError, match="does not cross"):
        module._ratio_profile_bracket(absent, 1.0)


@pytest.mark.parametrize("ratio", [0, -1, True, float("nan"), float("inf"), "1", 1j])
def test_profile_ratio_invalid_inputs_fail_closed(ratio) -> None:
    module = _load_module()
    with pytest.raises(ValueError, match="finite and positive"):
        module._ratio_profile_bracket(_valley(module), ratio)


def test_resource_limit_cannot_issue_a_false_certificate() -> None:
    module = _load_module()
    with pytest.raises(ValueError, match="did not reach"):
        module._ratio_profile_bracket(_valley(module), 1.0, max_subdivisions=0)


def test_nonfinite_grid_rejected_before_interval_operations() -> None:
    module = _load_module()
    grid = module.RectilinearGrid.from_rows(
        [(x, y, float("nan")) for x in (0.0, 1.0) for y in (1.0, 2.0)], "invalid")
    with pytest.raises(ValueError, match="finite real"):
        module._ratio_profile_bracket(grid, 1.0)


def test_historical_projection_records_different_bytes_without_target_change() -> None:
    module = _load_module()
    coordinates = {"ratio": .03}
    prior = {"artifact": "oph_neutrino_nufit61_retrospective_profile_score",
             "candidate": {"coordinates": coordinates, "sha256": "old"},
             "nufit_release": {"source_manifest_sha256": "manifest"},
             "source_boundary": {"family_transport_kernel_sha256": "kernel"}}
    result = module._prior_score_projection(prior, coordinates, "new", "manifest", "kernel")
    assert result["scored_coordinate_projection_exactly_equal"]
    assert result["candidate_byte_identity"] is False
    assert result["original_candidate_sha256"] == "old"
    for section, key, value in (
        ("candidate", "coordinates", {"ratio": .04}),
        ("nufit_release", "source_manifest_sha256", "different"),
        ("source_boundary", "family_transport_kernel_sha256", "different"),
    ):
        corrupted = copy.deepcopy(prior)
        corrupted[section][key] = value
        with pytest.raises(ValueError, match="changed"):
            module._prior_score_projection(corrupted, coordinates, "new", "manifest", "kernel")


def test_content_digest_ignores_the_generation_timestamp(tmp_path) -> None:
    """Re-stamping a receipt must not change the provenance digest."""
    module = _load_module()
    committed = module.DEFAULT_CANDIDATE
    payload = json.loads(committed.read_text(encoding="utf-8"))
    assert "generated_utc" in payload
    stamped = []
    for index, stamp in enumerate(("2026-01-01T00:00:00Z", "2027-12-31T23:59:59Z")):
        path = tmp_path / f"candidate_{index}.json"
        path.write_text(
            json.dumps(dict(payload, generated_utc=stamp), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        stamped.append(path)
    assert module._sha256(stamped[0]) != module._sha256(stamped[1])
    assert module._content_sha256(stamped[0]) == module._content_sha256(stamped[1])
    # And the committed receipt's own content digest is the stable one.
    assert module._content_sha256(committed) == module._content_sha256(stamped[0])


def test_published_score_receipt_keeps_its_generation_timestamp() -> None:
    module = _load_module()
    payload = json.loads(module.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", payload["generated_utc"])


def test_prior_projection_reports_content_identity_across_a_restamp() -> None:
    module = _load_module()
    coordinates = {"ratio": .03}
    prior = {"artifact": "oph_neutrino_nufit61_retrospective_profile_score",
             "candidate": {"coordinates": coordinates, "sha256": "old-bytes",
                           "content_sha256": "same-content"},
             "nufit_release": {"source_manifest_sha256": "manifest"},
             "source_boundary": {"family_transport_kernel_sha256": "kernel"}}
    result = module._prior_score_projection(
        prior, coordinates, "new-bytes", "manifest", "kernel", "same-content")
    assert result["candidate_byte_identity"] is False
    assert result["candidate_content_identity"] is True
    assert result["original_candidate_content_sha256"] == "same-content"
    legacy = copy.deepcopy(prior)
    legacy["candidate"].pop("content_sha256")
    legacy_result = module._prior_score_projection(
        legacy, coordinates, "new-bytes", "manifest", "kernel", "same-content")
    assert legacy_result["candidate_content_identity"] is False
