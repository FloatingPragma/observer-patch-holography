"""Independent checks of inversion, scheme comparisons and claim boundaries."""
import copy
import importlib.util
import math

import mpmath as mp
import pytest

from . import derive


def test_receipt_and_conditional_comparison_reproduce():
    result = derive.build()
    assert result == derive.read("results.json")
    assert result["outside_printed_intervals"] == 12
    assert result["lepton_replay_max_relative_error"] < 1e-13
    assert all(not group["intersection_nonempty"] for group in result["common_factor"].values())


@pytest.mark.parametrize("factor", [.1, .7, 1, 3, 20])
def test_exponent_inverse_solves_forward_equation_and_has_right_direction(factor):
    initial, exponent = 3.7, 6
    target = initial/factor
    solved = derive.frozen_transport_exponent(exponent, initial, target)
    assert math.isclose(initial*6**(-(solved-exponent)), target, rel_tol=1e-14)
    if factor > 1:
        assert solved > exponent
        assert initial*6**(-(2*exponent-solved-exponent)) > initial
    if factor < 1:
        assert solved < exponent


def test_source_law_matches_existing_register_producer_on_same_coordinates():
    path = derive.HERE.parents[1]/"code/particles/leptons/derive_charged_mcpr_completion_conditional.py"
    spec = importlib.util.spec_from_file_location("charged_source_replay", path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    packet = derive.read("inputs.json")["boundary_packet"]
    values = derive.forward_coefficients(packet)
    with mp.workdps(80):
        p, u = mp.mpf(packet["P_decimal"]), mp.mpf(packet["alpha_U_decimal"])
        geometry = source.icosahedral_incidence_solve()
        fixed = source.response_fixed_point(p, u, source.register_dimensions(geometry))
        triple = source.charged_triple(p, u, fixed, geometry)
        for key, coefficient in zip(("e", "mu", "tau"), triple["masses_over_v"]):
            assert math.isclose(values["masses_GeV"][key], float(coefficient*packet["v_GeV"]), rel_tol=1e-14)


def test_reference_mutation_changes_residuals_without_changing_forward_coefficients():
    inputs = derive.read("inputs.json")
    reference = derive.read("pdg2026_reference.json")
    baseline = derive.build(inputs, reference)
    changed = copy.deepcopy(reference)
    changed["masses"]["e"]["central"] = "0.52"
    mutated = derive.build(inputs, changed)
    assert mutated["forward_leptons"] == baseline["forward_leptons"]
    assert mutated["rows"][0]["value_GeV"] == baseline["rows"][0]["value_GeV"]
    assert mutated["rows"][0]["residual_ppm"] != baseline["rows"][0]["residual_ppm"]


def test_common_factor_rectangle_has_both_feasible_and_disjoint_cases():
    refs = {"a": {"central": "2", "error_minus": ".1", "error_plus": ".1", "unit": "GeV"},
            "b": {"central": "4", "error_minus": ".1", "error_plus": ".1", "unit": "GeV"}}
    assert derive.factor_intersection({"a": 1, "b": 2}, refs)["intersection_nonempty"]
    assert not derive.factor_intersection({"a": 1, "b": 1}, refs)["intersection_nonempty"]
    with pytest.raises(ValueError):
        derive.factor_intersection({"a": 0, "b": 1}, refs)


def test_ideal_circulant_does_not_pass_joint_lepton_ratio_check():
    ratio = derive.build()["lepton_ratios"]["mu_over_e"]
    assert not ratio["within_printed_rectangle"]
    assert not ratio["exact_point_inside_rectangle"]
    assert 9 < ratio["exact_Koide_Brannen_residual_ppm"] < 11
    assert .07 < ratio["residual_ppm"] < .08


def test_neutral_diagnostic_uses_finite_bound_hypotheses():
    from model.particle_finite_response.certificate import phase_bound
    neutral = derive.build()["neutral"]
    bound = phase_bound()
    assert math.isclose(neutral["conditional_unshifted_phase_bound"],
                        bound["normal_adjacent_squared_gap_ratio_upper"], rel_tol=1e-14)
    assert neutral["conditional_unshifted_phase_bound"] < neutral["normal_gap_ratio_reference_rectangle"][0]
    assert neutral["minimal_action_squared_gaps_eV2"] == [0, 0]
    assert not neutral["arbitrary_neutral_models_excluded"]


def test_scope_cannot_promote_reference_controls_or_runtime_separation():
    result = derive.build()
    assert result["forward_leptons"]["target_informed_architecture"]
    assert result["scope"]["physical_spectrum_identity_established"] is False
    assert result["scope"]["particle_abundance_distribution_computed"] is False
    assert result["scope"]["joint_statistical_significance_computed"] is False
    assert result["scope"]["theory_error_bound"] is None
    assert all(not c["prospective_prediction"] for c in result["vev_substitution_controls"])
    assert result["massless_gauge"]["photon"]["reference_kind"] == "upper_limit"
    assert result["massless_gauge"]["gluon"]["reference_kind"] == "theoretical_mass_parameter"
