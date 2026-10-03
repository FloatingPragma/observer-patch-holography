"""Tests retain disagreements and prohibit target/provenance laundering."""
import importlib.util
import json
from pathlib import Path

import pytest

PATH = Path(__file__).with_name("compare_conditional_quark_masses.py")
SPEC = importlib.util.spec_from_file_location("conditional_comparison", PATH)
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


def inputs():
    return [json.loads(p.read_text()) for p in (comparison.PACKET, comparison.PDG, comparison.FLAG)]


def test_receipt_and_all_comparisons_retained():
    receipt = comparison.build()
    assert receipt == json.loads(comparison.OUTPUT.read_text())
    assert tuple(receipt["branches"]) == comparison.BRANCHES
    for branch in receipt["branches"].values():
        assert len(branch["pdg_individual_masses"]) == 6
        assert len(branch["pdg_light_quantities"]) == 3
        assert len(branch["flag_light_ratios"]) == 4
    full = receipt["branches"]["full_rscc"]
    assert all(row["inside_quoted_marginal_range"] for row in full["pdg_individual_masses"])
    ratios = {row["reference"]["quantity"]: row for row in full["pdg_light_quantities"]}
    assert not ratios["ms_over_mud"]["inside_quoted_marginal_range"]
    assert receipt["joint_likelihood"] is None
    assert receipt["theory_uncertainty"] is None
    top = receipt["criticality_top_running_coordinate_comparisons"]
    assert len(top) == 3
    assert top["coupled_two_loop_top_only"]["inside_quoted_marginal_range"]
    assert not top["one_loop"]["inside_quoted_marginal_range"]


@pytest.mark.parametrize("guard", ["formula_discovery_target_informed", "quark_reference_values_used_at_runtime", "common_scale_yukawa_matrix"])
def test_provenance_mutations_rejected(guard):
    packet, pdg, flag = inputs()
    packet["guards"][guard] = not packet["guards"][guard]
    with pytest.raises(ValueError):
        comparison.compare(packet, pdg, flag)


@pytest.mark.parametrize("value", ["NaN", "Infinity", "-1", "0"])
def test_invalid_mass_rejected(value):
    packet, pdg, flag = inputs()
    packet["full_rscc"]["native_mass_chart_GeV"]["u"] = value
    with pytest.raises(ValueError):
        comparison.compare(packet, pdg, flag)


def test_new_disagreement_is_reported_without_branch_selection():
    packet, pdg, flag = inputs()
    packet["full_rscc"]["native_mass_chart_GeV"]["t"] = "190"
    output = comparison.compare(packet, pdg, flag)
    assert output["primary_branch"] == "full_rscc"
    top = output["branches"]["full_rscc"]["pdg_individual_masses"][-1]
    assert not top["inside_quoted_marginal_range"]


def test_asymmetric_range_is_not_a_gaussian_error():
    row = {"central": "27.33", "minus": "0.14", "plus": "0.18", "quoted_confidence": "90%"}
    output = comparison.quoted_range(comparison.D("27.50"), row)
    assert output["inside_quoted_marginal_range"]
    assert output["quoted_upper"] == "27.51"
    assert "sigma" not in output


@pytest.mark.parametrize("mutation", ["control_selection", "producer_chart", "reference_chart", "promotion", "missing_disagreement", "missing_flavor_theory", "criticality_chart"])
def test_comparison_laundering_rejected(mutation):
    packet, pdg, flag = inputs()
    if mutation == "control_selection":
        packet["guards"]["controls_used_for_model_selection"] = True
    elif mutation == "producer_chart":
        packet["model"]["charts"]["t"]["scheme"] = "MSbar"
    elif mutation == "reference_chart":
        pdg["masses"][-1]["scheme"] = "MSbar"
    elif mutation == "promotion":
        pdg["interpretation"]["new_independent_evidence_claimed"] = True
    elif mutation == "missing_disagreement":
        pdg["light_quantities"].pop(1)
    elif mutation == "missing_flavor_theory":
        flag["averages"].pop()
    else:
        pdg["criticality_top_comparison"]["scheme"] = "pole_extraction"
    with pytest.raises(ValueError):
        comparison.compare(packet, pdg, flag)
