"""False-green controls for the independent conditional six-mass verifier."""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path

import pytest

import verify_conditional_quark_mass_replay as check


@pytest.fixture
def receipt() -> dict:
    return json.loads(check.DEFAULT_RECEIPT.read_text())


def replace_at(payload: dict, path: tuple[str, ...], value: object) -> None:
    destination = payload
    for key in path[:-1]:
        destination = destination[key]
    destination[path[-1]] = value


def test_committed_receipt_passes_independent_equations(receipt):
    result = check.verify_payload(receipt)
    assert result["status"] == "PASS"
    assert result["rows_verified"] == 12


@pytest.mark.parametrize("path,value", [
    (("schema",), "oph.source_only_mass_prediction.v1"),
    (("claim_class",), "independent_blind_prediction"),
    (("status",), "physical_mass_derivation_complete"),
    (("inputs", "alpha_inverse", "value"), "136.99483517741294"),
    (("inputs", "alpha_inverse", "role"), "source_derived"),
    (("inputs", "E_star_GeV", "value"), "1.221e19"),
    (("inputs", "E_star_GeV", "role"), "derived_absolute_unit"),
    (("inputs", "P", "definition"), "legacy_flavor_dependent_fixed_point"),
    (("inputs", "P", "value"), "1.6309720958588974"),
    (("inputs", "structural_counts", "oriented_repair_slots"), 25),
    (("source_packet", "P"), "1.6309720958588974"),
    (("source_packet", "alpha_U"), "0.041124247441816685"),
    (("source_packet", "v_over_E_star"), "2.019811407857633e-17"),
    (("source_packet", "mz_over_E_star"), "7.501767385088419e-18"),
    (("source_packet", "mu_U_over_E_star"), "0.002"),
    (("source_packet", "alpha1_mz"), "0.016"),
    (("source_packet", "alpha2_mz"), "0.034"),
    (("source_packet", "alpha3_mz"), "0.1184"),
    (("source_packet", "alpha_U"), "nan"),
    (("source_packet", "alpha_U"), "1e999"),
    (("source_packet", "alpha_U"), 0.04112433619563048),
    (("model", "primary"), "lower_order_control"),
    (("model", "control"), "selected_by_smallest_mass_error"),
    (("model", "charts", "u", "scale"), "MZ_common_scale"),
    (("model", "charts", "c", "scheme"), "physical_dimensionless_Yukawa"),
    (("model", "charts", "t", "scheme"), "MSbar"),
    (("model", "coefficient_menu", "mean_u_quadratic"), "1/28"),
    (("model", "equations", "mean_up"), "3*P+79*w/15+w^2/28"),
    (("model", "precision_boundary"), "printed digits are physical precision"),
    (("model", "gauge_solver", "su3_cutoff"), 20),
    (("hypotheses", "H8_chart_attachment"), "physical common-scale masses proved"),
    (("full_rscc", "coefficients", "a_u"), "5.6431"),
    (("full_rscc", "coefficients", "delta_g"), "0"),
    (("full_rscc", "coefficients", "B"), "0.81"),
    (("full_rscc", "mass_over_v", "t"), "0.7"),
    (("full_rscc", "mass_over_v", "u"), "-0.00000875"),
    (("full_rscc", "mass_over_E_star", "s"), "7.61e-21"),
    (("full_rscc", "native_mass_chart_GeV", "c"), "1.273"),
    (("full_rscc", "light_ratios", "u_over_d"), "0.5"),
    (("lower_order_control", "coefficients", "delta_g"), "0.001"),
    (("lower_order_control", "native_mass_chart_GeV", "b"), "4.19"),
])
def test_altered_scientific_input_or_result_fails(receipt, path, value):
    replace_at(receipt, path, value)
    with pytest.raises(ValueError):
        check.verify_payload(receipt)


@pytest.mark.parametrize("name", sorted(check.GUARDS))
def test_every_claim_and_provenance_guard_is_enforced(receipt, name):
    receipt["guards"][name] = not receipt["guards"][name]
    with pytest.raises(ValueError, match="guard"):
        check.verify_payload(receipt)


def test_integer_zero_cannot_impersonate_false_guard(receipt):
    receipt["guards"]["independent_prediction"] = 0
    with pytest.raises(ValueError, match="guard"):
        check.verify_payload(receipt)


@pytest.mark.parametrize("location", [(), ("inputs",), ("source_packet",),
                                      ("model",), ("full_rscc",)])
def test_hidden_target_or_selector_field_is_rejected(receipt, location):
    target = receipt
    for key in location:
        target = target[key]
    target["observed_quark_mass_target"] = "173.0"
    with pytest.raises(ValueError):
        check.verify_payload(receipt)


def test_control_cannot_replace_primary_law(receipt):
    receipt["full_rscc"] = copy.deepcopy(receipt["lower_order_control"])
    with pytest.raises(ValueError):
        check.verify_payload(receipt)


def test_ratio_preserving_sector_rescale_still_fails(receipt):
    row = receipt["full_rscc"]
    for name in ("u", "c", "t"):
        for field in ("mass_over_v", "mass_over_E_star", "native_mass_chart_GeV"):
            row[field][name] = str(float(row[field][name]) * 1.001)
    with pytest.raises(ValueError, match="geometric mean"):
        check.verify_payload(receipt)


def test_geometric_mean_preserving_shape_change_still_fails(receipt):
    row = receipt["full_rscc"]
    for name, factor in (("u", 1.001), ("c", 1 / 1.001)):
        for field in ("mass_over_v", "mass_over_E_star", "native_mass_chart_GeV"):
            row[field][name] = str(float(row[field][name]) * factor)
    with pytest.raises(ValueError, match="log ratio"):
        check.verify_payload(receipt)


def test_observed_target_path_cannot_enter_numeric_ancestry(receipt):
    receipt["model"]["numeric_input_paths"].append("code/particles/quark_references.json")
    with pytest.raises(ValueError, match="ancestry"):
        check.verify_payload(receipt)


def test_observed_target_source_pin_is_rejected(receipt):
    receipt["source_pins"]["code/particles/quark_references.json"] = "0" * 64
    with pytest.raises(ValueError, match="ancestry"):
        check.verify_payload(receipt)


def test_actual_source_digest_is_checked(receipt):
    receipt["source_pins"][check.PRODUCER] = "0" * 64
    with pytest.raises(ValueError, match="source pin drift"):
        check.verify_payload(receipt)


@pytest.mark.parametrize("mutation", [
    lambda source: source + "\nfrom quark_reference_data import measured_masses\n",
    lambda source: source.replace("ALPHA_PATH.read_text()", "Path('quark_targets.json').read_text()"),
    lambda source: source + "\ndef hidden(path):\n    return json.loads(path.read_text())\n",
    lambda source: source.replace("GAUGE_PATH)\n", "ALPHA_PATH)\n"),
])
def test_target_reads_and_dynamic_import_changes_are_rejected(mutation):
    original = (check.REPO / check.PRODUCER).read_text()
    modified = mutation(original)
    assert original != modified
    with pytest.raises(ValueError):
        check.check_producer_inputs(modified)


def test_constructor_cannot_restore_stage5_flavor_ancestry():
    producer = (check.REPO / check.PRODUCER).read_text()
    gauge = (check.REPO / check.GAUGE).read_text()
    modified = producer.replace("return {}", "return super()._derive_stage5_integer_vectors()")
    assert modified != producer
    with pytest.raises(ValueError, match="Stage-5"):
        check.check_gauge_call_graph(modified, gauge)


def test_reachable_gauge_method_cannot_consume_quark_template():
    producer = (check.REPO / check.PRODUCER).read_text()
    gauge = (check.REPO / check.GAUGE).read_text()
    modified = gauge.replace("v_ev = self.v_from_transmutation(alpha_u, p)",
                             "v_ev = self.diagonal_quark_masses(p)")
    assert modified != gauge
    with pytest.raises(ValueError, match="nongauge method"):
        check.check_gauge_call_graph(producer, modified)


def test_verifier_imports_neither_producer_nor_gauge_solver():
    tree = ast.parse(Path(check.__file__).read_text())
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    assert not any("conditional_quark_mass_replay" in name or "paper_math" in name
                   for name in imports)


@pytest.mark.parametrize("key,value", [("family_dimension",6),("line_complement_dimension",3),("scalar_fraction","1/6"),("centered_fraction","3/5"),("common_exposure","4/5")])
def test_certified_source_projection_corruption_fails(receipt,key,value):
    receipt["model"]["source_response"]["projection"][key]=value
    with pytest.raises(ValueError):check.verify_payload(receipt)


def test_source_response_cannot_be_promoted_to_selected_flavor(receipt):
    receipt["model"]["source_response"]["physical_boundary"]="All physical readouts selected by repair dynamics"
    with pytest.raises(ValueError):check.verify_payload(receipt)
