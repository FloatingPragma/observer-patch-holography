from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest


HERE = Path(__file__).resolve().parents[1]
REPO_ROOT = HERE.parents[3]
OUTPUT = HERE / "outputs" / "source_parent_inventory.json"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_module("wz594_frontier_builder", HERE / "build_source_parent_inventory.py")
checker = load_module("wz594_frontier_checker", HERE / "check_source_parent_inventory.py")


def load_inventory() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def canonical_sha256(value) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def rehash_inventory(payload: dict) -> None:
    payload["finite_parent_bundle_digest"] = canonical_sha256(
        checker._bundle_seed(payload["positive_parent_bindings"])
    )
    payload["conditional_context_digest"] = canonical_sha256(
        checker._bundle_seed(payload["conditional_context"])
    )
    body = {key: value for key, value in payload.items() if key != "inventory_digest"}
    payload["inventory_digest"] = canonical_sha256(body)


def copy_input_closure(tmp_path: Path, inventory: dict) -> Path:
    root = tmp_path / "sealed"
    paths = set(inventory["target_firewall"]["resolved_source_paths"])
    paths.add(builder.POLICY_REL)
    for relative in sorted(paths):
        source = REPO_ROOT / relative
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return root


def update_pin_for_path(inventory: dict, root: Path, relative: str) -> None:
    path = root / relative
    raw = path.read_bytes()
    parsed = json.loads(raw)
    pin_groups = (
        inventory["positive_parent_bindings"]
        + inventory["conditional_context"]
        + inventory["resolved_boundaries"]
    )
    for binding in pin_groups:
        for pin in binding["files"]:
            if pin["path"] == relative:
                pin["bytes"] = len(raw)
                pin["byte_sha256"] = hashlib.sha256(raw).hexdigest()
                pin["canonical_json_sha256"] = canonical_sha256(parsed)
                rehash_inventory(inventory)
                return
    raise AssertionError(f"pin not found: {relative}")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def assert_rejected(payload: dict, match: str) -> None:
    with pytest.raises(checker.FrontierVerificationError, match=match):
        checker.verify_inventory(
            payload,
            repo_root=REPO_ROOT,
            run_native_verifiers=False,
        )


def test_committed_inventory_passes_independent_resolution() -> None:
    result = checker.verify_inventory(
        load_inventory(),
        repo_root=REPO_ROOT,
        run_native_verifiers=False,
    )
    assert result["status"] == "PASS"
    assert result["positive_parent_count"] == 4
    assert result["conditional_context_count"] == 7
    assert result["promotion_allowed"] is False


def test_build_is_byte_exact_and_deterministic(tmp_path: Path) -> None:
    first = builder.build_inventory(REPO_ROOT)
    second = builder.build_inventory(REPO_ROOT)
    assert first == second
    out = tmp_path / "frontier.json"
    builder.write_inventory(out, first)
    assert out.read_bytes() == OUTPUT.read_bytes()


def test_only_allowlisted_inputs_are_sufficient_for_replay(tmp_path: Path) -> None:
    committed = load_inventory()
    sealed_root = copy_input_closure(tmp_path, committed)
    rebuilt = builder.build_inventory(sealed_root)
    result = checker.verify_inventory(
        rebuilt,
        repo_root=sealed_root,
        run_native_verifiers=False,
    )
    assert result["status"] == "PASS"
    assert rebuilt["inventory_digest"] == committed["inventory_digest"]


def test_unknown_inventory_field_fails_closed() -> None:
    mutated = load_inventory()
    mutated["producer_says_pass"] = True
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation")


def test_path_traversal_fails_closed() -> None:
    mutated = load_inventory()
    mutated["positive_parent_bindings"][0]["files"][0]["path"] = (
        "code/a5_closure/manifests/../../particles/calibration/"
        "wz_pdg_2026_target_fixture.json"
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "file set or order drifted|path traversal")


def test_parent_byte_drift_fails_closed(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/echosahedral_federation_reference.json"
    path = sealed_root / relative
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(checker.FrontierVerificationError, match="byte count mismatch"):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_rehashed_source_firewall_contamination_fails_closed(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/receipts/"
        "echosahedral_federation_reference.receipt.json"
    )
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["source_firewall"]["forbidden_dependency_hits"] = ["experimental_target"]
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(checker.FrontierVerificationError, match="#565 source firewall"):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_rehashed_structured_target_injection_fails_closed(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/manifests/"
        "charged_response_pole_residue_artifact.json"
    )
    path = sealed_root / relative
    artifact = json.loads(path.read_text(encoding="utf-8"))
    artifact["experimental_target"] = {"mounted": True}
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    artifact["artifact_sha256"] = "sha256:" + canonical_sha256(body)
    write_json(path, artifact)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(checker.FrontierVerificationError, match="structured target content"):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_family_candidate_cannot_be_promoted_to_positive_parent() -> None:
    mutated = load_inventory()
    mutated["conditional_context"][0]["status"] = "verified_finite_parent"
    rehash_inventory(mutated)
    assert_rejected(mutated, "wrong binding status")


def test_declared_parent_scope_cannot_be_rewritten() -> None:
    mutated = load_inventory()
    mutated["conditional_context"][0]["usable_exports"] = [
        "physical family and chiral gauge action",
    ]
    mutated["conditional_context"][0]["excluded_promotions"] = [
        "no excluded physical promotions"
    ]
    rehash_inventory(mutated)
    assert_rejected(mutated, "declared export scope drifted")


def test_family_open_receipts_are_load_bearing(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/manifests/family_band_attachment_reference.json"
    )
    path = sealed_root / relative
    certificate = json.loads(path.read_text(encoding="utf-8"))
    certificate["named_interface"]["open_receipts"] = []
    body = {key: value for key, value in certificate.items() if key != "manifest_sha256"}
    certificate["manifest_sha256"] = "sha256:" + canonical_sha256(body)
    write_json(path, certificate)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(checker.FrontierVerificationError, match="physical-family boundary"):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_family_tensor_identity_limitation_is_load_bearing(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/matter_attachment_receipt.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["gap_inheritance_certificate"]["structure"] = (
        "physical family-sensitive chiral gauge operator"
    )
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="tensor-identity limitation",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_spin_packet_cannot_be_promoted_onto_local_domain(
    tmp_path: Path,
) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/matter_attachment_receipt.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["spin_layer"]["spin_to_local_domain_bridge_certified"] = True
    receipt["spin_layer"]["same_source_domain_certified"] = True
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="silently attached",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_ambiguous_matter_domain_clause_is_rejected(
    tmp_path: Path,
) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/matter_attachment_receipt.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["clause_verdicts"]["same_source_domain_binding"] = True
    receipt["clause_verdicts"].pop(
        "local_stage2_same_source_domain_binding",
        None,
    )
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="conflate the spin and local-domain packets",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_declared_matter_tensor_cannot_be_marked_source_selected(
    tmp_path: Path,
) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/matter_attachment_receipt.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["matter_operator_certificate"]["source_selected"] = True
    receipt["gap_inheritance_certificate"][
        "matter_action_source_selected"
    ] = True
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="relabeled as source-selected",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_finite_classical_spectral_control_is_load_bearing(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/receipts/"
        "flux_defect_criterion_reference.receipt.json"
    )
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["exact_support_classical_realification"][
        "declared_adjacency_spectral_family_recoverable"
    ] = False
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="exact-support vector-spring realification",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_local_and_exact_flux_domains_cannot_be_identified(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/receipts/"
        "flux_defect_criterion_reference.receipt.json"
    )
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["local_domain_classical_spectral_context"][
        "exact_flux_identity_bridge"
    ] = True
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="separate local-domain vector-spring context",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_exact_flux_classical_metric_is_load_bearing(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/receipts/"
        "flux_defect_criterion_reference.receipt.json"
    )
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["exact_support_classical_realification"]["per_class"][0][
        "classical_hessian_certificate"
    ]["coordinate_metric"] = [[1, 0], [0, 1]]
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="positive-metric edge-Hessian proof",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


@pytest.mark.parametrize(
    ("role", "field", "forged_value", "match"),
    [
        (
            "scalar_yukawa_source_frontier",
            "promotion_allowed",
            True,
            "#630 frontier promoted itself",
        ),
        (
            "source_clock_frontier",
            "physical_promotion_allowed",
            True,
            "#633 frontier fabricated",
        ),
    ],
)
def test_partial_frontier_promotions_fail_closed(
    tmp_path: Path,
    role: str,
    field: str,
    forged_value,
    match: str,
) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    binding = next(item for item in inventory["conditional_context"] if item["role"] == role)
    relative = binding["files"][0]["path"]
    path = sealed_root / relative
    frontier = json.loads(path.read_text(encoding="utf-8"))
    frontier[field] = forged_value
    write_json(path, frontier)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(checker.FrontierVerificationError, match=match):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_local_carrier_physical_promotion_fails_closed(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    binding = next(
        item
        for item in inventory["conditional_context"]
        if item["role"] == "local_ew_order_unit_frontier"
    )
    relative = binding["files"][0]["path"]
    path = sealed_root / relative
    frontier = json.loads(path.read_text(encoding="utf-8"))
    first_key = next(iter(frontier["promotion"]))
    frontier["promotion"][first_key] = True
    body = {key: value for key, value in frontier.items() if key != "manifest_sha256"}
    frontier["manifest_sha256"] = "sha256:" + canonical_sha256(body)
    write_json(path, frontier)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="#631 finite line isomorphism was given a physical promotion",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_v_chart_cannot_be_relabelled_v_f() -> None:
    mutated = load_inventory()
    mutated["coordinate_bridge"]["source_coordinate"] = "v_F"
    mutated["coordinate_bridge"]["relabel_allowed"] = True
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|relabelled")


def test_physical_units_or_poles_cannot_be_emitted() -> None:
    mutated = load_inventory()
    mutated["unit_scope"]["dimensionful_values_present"] = True
    mutated["unit_scope"]["emitted_observables"] = [{"m_W_GeV": "80"}]
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation")


def test_physical_unit_verdict_cannot_be_reclassified() -> None:
    mutated = load_inventory()
    mutated["unit_scope"]["physical_unit_verdict"] = "OPEN"
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|physical-unit row drifted")


def test_resolved_finite_boundary_cannot_be_promoted_to_continuum() -> None:
    mutated = load_inventory()
    boundary = next(
        item
        for item in mutated["resolved_boundaries"]
        if item["provenance_issue"] == 634
    )
    boundary["classification"] = "positive_physical_source"
    boundary["effect"] = "continuum Lorentzian spacetime attained"
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|classification drifted")


def test_resolved_boundary_receipt_verdict_is_load_bearing(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/clock_unit_verdict.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["verdict"] = "PHYSICAL_UNITS_EVALUABLE"
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="receipt verdict drift",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_local_classical_receipt_cannot_claim_exact_flux_identity(
    tmp_path: Path,
) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = (
        "code/a5_closure/manifests/classical_realization_receipt.json"
    )
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["spectral_interface_identity"][
        "rer_exact_flux_12_42_vertex_identity_bridge"
    ] = True
    receipt["spectral_interface_identity"][
        "separate_from_rer_exact_flux_certificate"
    ] = False
    write_json(path, receipt)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="conflated with the exact flux support",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_resolved_boundary_target_injection_fails_closed(tmp_path: Path) -> None:
    inventory = load_inventory()
    sealed_root = copy_input_closure(tmp_path, inventory)
    relative = "code/a5_closure/manifests/classical_realization_receipt.json"
    path = sealed_root / relative
    receipt = json.loads(path.read_text(encoding="utf-8"))
    receipt["experimental_target"] = {"mounted": True}
    write_json(path, receipt)
    with pytest.raises(
        builder.FrontierBuildError,
        match="forbidden structured target content",
    ):
        builder.build_inventory(sealed_root)
    update_pin_for_path(inventory, sealed_root, relative)
    with pytest.raises(
        checker.FrontierVerificationError,
        match="structured target content",
    ):
        checker.verify_inventory(
            inventory,
            repo_root=sealed_root,
            run_native_verifiers=False,
        )


def test_forged_promotion_flag_fails_closed() -> None:
    mutated = load_inventory()
    mutated["promotion_allowed"] = True
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation")


def test_consumer_cannot_be_marked_frozen_early() -> None:
    mutated = load_inventory()
    mutated["consumer_contract"]["frozen_algorithm_substitution_ready"] = True
    mutated["consumer_contract"]["common_subject_digest_ready"] = True
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|marked ready")


def test_consumer_schema_cannot_be_promoted_or_reclassified() -> None:
    mutated = load_inventory()
    mutated["consumer_contract"]["schemas"][0]["status"] = (
        "nonpromoting_specification_schema"
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "consumer schema status drift")


def test_target_path_cannot_enter_resolved_input_closure() -> None:
    mutated = load_inventory()
    mutated["target_firewall"]["resolved_source_paths"].append(
        "code/particles/calibration/wz_pdg_2026_target_fixture.json"
    )
    mutated["target_firewall"]["resolved_source_paths"].sort()
    rehash_inventory(mutated)
    assert_rejected(
        mutated,
        "outside allowlisted roots|forbidden target file|closure drifted",
    )


def test_source_dag_cycle_fails_closed() -> None:
    mutated = load_inventory()
    mutated["source_dag"]["edges"].append(
        {
            "from": "oph_native_dimensionless_packet",
            "to": "physical_family_and_matter_pole_attachment",
        }
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "dependency edges drifted|source DAG is cyclic")


def test_field_census_does_not_depend_on_common_carrier() -> None:
    mutated = load_inventory()
    mutated["source_dag"]["edges"].append(
        {
            "from": "common_screen_electroweak_carrier",
            "to": "source_complete_field_census",
        }
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "dependency edges drifted|incorrectly depends")


def test_scientific_interface_provenance_is_explicit() -> None:
    inventory = load_inventory()
    interfaces = {
        row["gate_id"]: row
        for row in inventory["required_interfaces"]
    }
    provenance = {
        gate_id: row["provenance_issues"]
        for gate_id, row in interfaces.items()
    }
    assert provenance["finite_to_lorentzian_quantum_eft_transfer"] == [635]
    assert provenance["source_scalar_action"] == [636]
    assert provenance["full_yukawa_operator_and_coefficients"] == [637]
    assert provenance["source_to_fj_coordinate_map"] == [638]
    assert provenance["scalar_yukawa_fj_integration"] == [630]
    rg_prerequisites = set(
        interfaces["target_clean_rg_threshold_matching"]["prerequisites"]
    )
    assert {
        "finite_local_domain_boundary",
        "physical_family_and_matter_pole_attachment",
        "source_scalar_action",
        "full_yukawa_operator_and_coefficients",
        "common_screen_electroweak_carrier",
        "source_complete_field_census",
    }.issubset(rg_prerequisites)
    assert "scalar_yukawa_fj_integration" not in rg_prerequisites
    assert "source_to_fj_coordinate_map" not in rg_prerequisites


def test_scientific_boundary_cannot_be_overpromoted() -> None:
    mutated = load_inventory()
    mutated["scientific_boundary_map"][0]["classification"] = (
        "positive_physical_source"
    )
    mutated["scientific_boundary_map"][0]["missing_interfaces"] = []
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|overpromoted")


def test_scientific_boundary_summary_cannot_overclaim() -> None:
    mutated = load_inventory()
    mutated["scientific_boundary_map"][0]["summary"] = (
        "the receipt supplies a physical Lorentzian quantum field theory"
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "claim scope drifted")


def test_boundary_row_nine_cannot_drop_a_replay_interface() -> None:
    mutated = load_inventory()
    row = next(
        item
        for item in mutated["scientific_boundary_map"]
        if item["boundary_index"] == 9
    )
    row["missing_interfaces"].remove("runtime_and_human_target_firewall")
    rehash_inventory(mutated)
    assert_rejected(
        mutated,
        "scientific boundary row 9 classification, missing interfaces, or claim scope drifted",
    )


def test_inventory_digest_is_load_bearing() -> None:
    mutated = copy.deepcopy(load_inventory())
    mutated["inventory_digest"] = "0" * 64
    assert_rejected(mutated, "self-digest mismatch")


def test_checker_does_not_import_the_producer() -> None:
    source = (HERE / "check_source_parent_inventory.py").read_text(encoding="utf-8")
    assert "import build_source_parent_inventory" not in source
    assert "from build_source_parent_inventory" not in source
