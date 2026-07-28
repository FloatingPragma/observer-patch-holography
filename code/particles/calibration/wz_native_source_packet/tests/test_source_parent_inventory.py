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
    for binding in inventory["positive_parent_bindings"] + inventory["conditional_context"]:
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
    assert result["conditional_context_count"] == 6
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
            "from": "external_qft_pole_consumer",
            "to": "oph_native_dimensionless_packet",
        }
    )
    rehash_inventory(mutated)
    assert_rejected(mutated, "source DAG is cyclic")


def test_acceptance_rows_cannot_self_close() -> None:
    mutated = load_inventory()
    mutated["acceptance_map"][0]["status"] = "complete"
    mutated["acceptance_map"][0]["blocking_gates"] = []
    rehash_inventory(mutated)
    assert_rejected(mutated, "schema validation|overpromoted")


def test_acceptance_row_nine_cannot_drop_a_native_action_gate() -> None:
    mutated = load_inventory()
    row = next(
        item
        for item in mutated["acceptance_map"]
        if item["acceptance_index"] == 9
    )
    row["blocking_gates"].remove("full_yukawa_operator_and_coefficients")
    rehash_inventory(mutated)
    assert_rejected(mutated, "acceptance row 9 lost")


def test_inventory_digest_is_load_bearing() -> None:
    mutated = copy.deepcopy(load_inventory())
    mutated["inventory_digest"] = "0" * 64
    assert_rejected(mutated, "self-digest mismatch")


def test_checker_does_not_import_the_producer() -> None:
    source = (HERE / "check_source_parent_inventory.py").read_text(encoding="utf-8")
    assert "import build_source_parent_inventory" not in source
    assert "from build_source_parent_inventory" not in source
