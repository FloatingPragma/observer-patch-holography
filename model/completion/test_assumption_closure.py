"""Independent coverage and adversarial disclosure tests for the scoped inventory."""
import copy
import json

import pytest

from model.completion import assumption_closure as ac


@pytest.fixture(scope="module")
def inventory():
    return json.loads(ac.PATH.read_text())


def row(rows, identifier):
    return next(item for item in rows if item["id"] == identifier)


def rejects(inventory, change):
    mutated = copy.deepcopy(inventory)
    change(mutated)
    with pytest.raises(ValueError):
        ac.validate(mutated)


def test_retained_inventory(inventory):
    counts = ac.validate(inventory)
    assert counts["canonical_rows"] == 84
    assert counts["selected_claims"] == 16
    assert counts["joint_assumptions"] == 10


def test_exact_canonical_identity_status_and_consumption(inventory):
    canonical = json.loads((ac.ROOT / "tracking/premise_register.json").read_text())["rows"]
    assert [r["id"] for r in inventory["canonical_rows"]] == [r["id"] for r in canonical]
    assert "PR-83" not in {r["id"] for r in canonical}
    for original, recorded in zip(canonical, inventory["canonical_rows"]):
        assert original["discharge"] == recorded["canonical_contract"]["discharge"]
    claims = {r["claim_id"]: r for r in inventory["claim_manifest"]}
    actual = {(r["id"], u["claim_id"], u["role"]) for r in inventory["canonical_rows"] for u in r["uses"]}
    expected = {(premise, c["claim_id"], role) for c in claims.values() for role in ("consumed", "open", "boundary") for premise in c["premise_dependencies"][role]}
    assert actual == expected
    native = claims["OPH-MODEL-NATIVE-OBSERVER-EXECUTION"]
    assert native["premise_dependencies"]["classification"] == "explicit_non_consumer"
    assert len(native["assumptions"]) == 48
    assert not any(claim == native["claim_id"] for _, claim, _ in actual)


def test_named_claim_inputs_have_exact_coverage(inventory):
    expected = {(c["claim_id"], name) for c in inventory["claim_manifest"] for name in c["assumptions"]}
    actual = [(r["claim_id"], r["statement"]) for r in inventory["claim_inputs"]]
    assert len(actual) == len(set(actual))
    assert set(actual) == expected


def test_joint_graph_consumers_are_exact(inventory):
    graph = json.loads((ac.ROOT / "model/particle_joint_reconstruction/dependencies.json").read_text())
    expected = {(name, n["id"]) for n in graph["nodes"] for name in n["conditional_assumptions"]}
    actual = {(r["id"], n) for r in inventory["joint_inputs"] for n in r["consumer_nodes"]}
    assert actual == expected
    # An empty gravity-node list does not erase the separately consumed gravity inputs.
    assert any(r["source_domain"] == "gravity_branch" for r in inventory["claim_inputs"])


def test_every_proof_has_actual_interface_and_finite_hypotheses(inventory):
    for proof in inventory["proved_inputs"]:
        assert set(proof["interface"]) == set(ac.FIELDS)
        assert all(proof["interface"].values())
        assert not proof["canonical_premises_discharged"]
        assert not proof["physical_map_supplied"]
        for input_id in proof["hypothesis_inputs"]:
            inp = row(inventory["service_inputs"], input_id)
            assert inp["classification"] not in {"open_physical_map", "open_implementation"}
            assert inp["source_domain"] == proof["source_domain"]
    assert row(inventory["proved_inputs"], "finite_population_absorption")["hypothesis_inputs"] == ["population_preparation"]
    assert "population_federation" in row(inventory["proved_inputs"], "finite_population_absorption")["downstream_obligations"]
    conjugation = row(inventory["proved_inputs"], "finite_conjugation_evenness")
    assert conjugation["classification"] == "imported_theorem_with_matched_hypotheses"
    assert {"finite_state", "permutation", "kernel_equivariance", "even_preparation", "odd_observable"} <= conjugation["matched_hypotheses"].keys()


def test_declared_interpretations_keep_open_obligations(inventory):
    thermal = row(inventory["service_inputs"], "thermal_density_statistics")
    assert thermal["classification"] == "declared_architectural_primitive"
    assert thermal["physical_interpretation"]
    assert thermal["obligation"]["classification"] == "open_physical_map"
    for item in inventory["service_inputs"]:
        assert item["obligation"]["smallest_witness"] and item["obligation"]["counterexample"]
    for batch in inventory["batch_admissions"]:
        assert not batch["physical_map_supplied"]
    assert {"independent_run_sets", "supplied_marginal_errors", "linear_influence_coefficients"} <= set(row(inventory["batch_admissions"], "electroweak_covariance")["inherited_inputs"])


def test_finite_numbers_and_scopes_are_separate(inventory):
    weak = row(inventory["proved_inputs"], "complete_lepton_instrument")
    assert "diameter 2" in weak["matched_hypotheses"]["generator"]
    assert "dimensionless" in weak["matched_hypotheses"]["time"]
    cube = inventory["domains"]["charged_cube"]
    assert "c=21/2" in cube["algebra"]
    assert "c=1/2" not in cube["algebra"]
    population = row(inventory["proved_inputs"], "finite_population_absorption")
    assert "nonincreasing" in population["statement"]
    assert "every initial distribution" in population["matched_hypotheses"]["uniform_tail"]
    assert "independent uniform" in population["matched_hypotheses"]["occupied_limit"]


@pytest.mark.parametrize("section", ["canonical_rows", "joint_inputs", "claim_inputs", "imported_inputs", "service_inputs", "proved_inputs", "batch_admissions", "source_join_obligations", "source_manifest", "theorem_premise_references"])
def test_missing_input_or_witness_rejected(inventory, section):
    rejects(inventory, lambda d: d[section].pop())


@pytest.mark.parametrize("section", ["canonical_rows", "joint_inputs", "claim_inputs", "service_inputs"])
def test_duplicate_not_counted_as_coverage(inventory, section):
    rejects(inventory, lambda d: d[section].__setitem__(0, copy.deepcopy(d[section][1])))


@pytest.mark.parametrize("category", [None, "", "unclassified", "proved_on_this_source", "imported_theorem_with_matched_hypotheses"])
def test_unclassified_or_promoted_service_rejected(inventory, category):
    rejects(inventory, lambda d: d["service_inputs"][0].__setitem__("classification", category))


def test_false_canonical_discharge_rejected(inventory):
    rejects(inventory, lambda d: d["canonical_rows"][0]["canonical_contract"]["discharge"].__setitem__("status", "proved"))


def test_calibration_cannot_become_independent_prediction(inventory):
    rejects(inventory, lambda d: row(d["service_inputs"], "lepton_calibration_ancestry").__setitem__("independent_prediction", True))


def test_missing_open_interpretation_obligation_rejected(inventory):
    rejects(inventory, lambda d: row(d["service_inputs"], "thermal_density_statistics").pop("obligation"))


def test_false_source_join_rejected(inventory):
    rejects(inventory, lambda d: d["source_join_obligations"][0].__setitem__("join_established", True))


def test_source_relabeling_rejected_even_with_same_numbers(inventory):
    rejects(inventory, lambda d: row(d["proved_inputs"], "fixed_charged_minmax").__setitem__("source_domain", "primary_wilson_gauss"))


@pytest.mark.parametrize("field", ac.FIELDS)
def test_each_source_interface_field_required(inventory, field):
    rejects(inventory, lambda d: d["proved_inputs"][0]["interface"].pop(field))


def test_arbitrary_reference_not_replaced_by_marginal(inventory):
    rejects(inventory, lambda d: d["proved_inputs"][0]["interface"].__setitem__("reference", "source marginal only"))


def test_receipt_hash_is_not_proof(inventory):
    rejects(inventory, lambda d: d["proved_inputs"][0]["proof"].__setitem__("path", "model/lepton_cycle/results.json"))


def test_imported_theorem_needs_hypothesis_match(inventory):
    rejects(inventory, lambda d: row(d["proved_inputs"], "finite_conjugation_evenness")["matched_hypotheses"].pop("kernel_equivariance"))


def test_downstream_map_not_added_as_finite_hypothesis(inventory):
    rejects(inventory, lambda d: row(d["proved_inputs"], "finite_population_absorption")["hypothesis_inputs"].append("population_federation"))


def test_stale_source_hash_rejected(inventory):
    rejects(inventory, lambda d: d["source_manifest"][0].__setitem__("sha256", "0" * 64))


def test_wrong_revision_and_missing_owner_rejected(inventory):
    rejects(inventory, lambda d: d["audit"].__setitem__("baseline_revision", "0" * 40))
    rejects(inventory, lambda d: d["service_inputs"][0].__setitem__("owner", "CP-99"))


def test_join_cannot_hide_in_batch_admission(inventory):
    rejects(inventory, lambda d: d["batch_admissions"][0].__setitem__("physical_map_supplied", True))


def test_path_escape_rejected():
    for path in ("../outside", "/etc/passwd", "missing-scientific-input"):
        with pytest.raises(ValueError):
            ac.local_file(path)


def test_every_import_edge_retains_application_boundary(inventory):
    expected = {(c["claim_id"], name) for c in inventory["claim_manifest"] for name in c["imported_results"]}
    actual = [(r["claim_id"], r["statement"]) for r in inventory["imported_inputs"]]
    assert len(actual) == len(set(actual)) and set(actual) == expected
    for item in inventory["imported_inputs"]:
        assert not item["hypotheses_matched_here"]
        assert not item["source_specific_discharge"]
        if "CODATA" in item["statement"] or "PDG" in item["statement"]:
            assert item["classification"] == "empirical_calibration"
            assert not item["independent_prediction"]
    rejects(inventory, lambda d: d["imported_inputs"][0].__setitem__("hypotheses_matched_here", True))
