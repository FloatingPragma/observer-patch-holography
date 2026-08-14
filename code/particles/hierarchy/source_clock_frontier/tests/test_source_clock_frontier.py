from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import shutil

import pytest


HERE = Path(__file__).resolve().parents[1]
REPO_ROOT = HERE.parents[3]
OUTPUT = HERE / "outputs" / "source_clock_frontier.json"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_module("clock633_frontier_builder", HERE / "build_source_clock_frontier.py")
checker = load_module("clock633_frontier_checker", HERE / "check_source_clock_frontier.py")


def load_frontier() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def canonical_sha256(value) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def rehash(payload: dict) -> None:
    body = {key: value for key, value in payload.items() if key != "frontier_digest"}
    payload["frontier_digest"] = canonical_sha256(body)


def assert_rejected(payload: dict, match: str) -> None:
    with pytest.raises(checker.FrontierVerificationError, match=match):
        checker.verify_frontier(payload, repo_root=REPO_ROOT)


def copy_input_closure(tmp_path: Path, frontier: dict) -> Path:
    root = tmp_path / "sealed"
    paths = set(frontier["target_firewall"]["source_input_paths"])
    paths.update(frontier["target_firewall"]["diagnostic_paths"])
    for relative in sorted(paths):
        source = REPO_ROOT / relative
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return root


def test_committed_frontier_passes_independent_resolution() -> None:
    result = checker.verify_frontier(load_frontier(), repo_root=REPO_ROOT)
    assert result["status"] == "PASS"
    assert load_frontier()["issue"] == 633
    assert result["hard_issue_dependencies"] == [634]
    assert result["alternative_routes_allowed"] is True
    assert result["optional_candidate_route_count"] == 1
    assert result["optional_candidate_component_count"] == 5
    assert result["physical_source_payload_count"] == 0
    assert result["physical_promotion_allowed"] is False


def test_builder_is_deterministic_and_byte_exact(tmp_path: Path) -> None:
    first = builder.build_frontier(REPO_ROOT)
    second = builder.build_frontier(REPO_ROOT)
    assert first == second
    out = tmp_path / "frontier.json"
    builder.write_frontier(out, first)
    assert out.read_bytes() == OUTPUT.read_bytes()


def test_sealed_input_closure_is_sufficient(tmp_path: Path) -> None:
    committed = load_frontier()
    sealed = copy_input_closure(tmp_path, committed)
    rebuilt = builder.build_frontier(sealed)
    result = checker.verify_frontier(rebuilt, repo_root=sealed)
    assert result["status"] == "PASS"
    assert rebuilt["frontier_digest"] == committed["frontier_digest"]


def test_interval_inversion_is_exact_and_order_reversing() -> None:
    receipt = load_frontier()["interval_inversion_receipt"]
    eps_lo = Fraction(receipt["epsilon_interval"]["lower"])
    eps_hi = Fraction(receipt["epsilon_interval"]["upper"])
    e_j_lo = Fraction(receipt["energy_interval_J"]["lower"])
    e_j_hi = Fraction(receipt["energy_interval_J"]["upper"])
    numerator = checker.H_SI * checker.NU_CS_HZ
    assert Fraction(0) < eps_lo <= eps_hi
    assert e_j_lo == numerator / eps_hi
    assert e_j_hi == numerator / eps_lo
    assert e_j_lo <= e_j_hi


def test_unknown_top_level_field_fails_closed() -> None:
    mutated = load_frontier()
    mutated["producer_says_physical"] = True
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_component_status_cannot_be_relabelled() -> None:
    mutated = load_frontier()
    mutated["candidate_routes"][0]["component_contracts"][0][
        "current_status"
    ] = "absent"
    rehash(mutated)
    assert_rejected(mutated, "candidate routes do not match")


def test_hard_dependency_set_is_exactly_issue_634() -> None:
    frontier = load_frontier()
    assert frontier["dependency_semantics"]["hard_issue_dependencies"] == [634]
    assert frontier["dependency_semantics"]["optional_route_owner_issues"] == [
        32,
        34,
        317,
        318,
        425,
        522,
        545,
        546,
        569,
        633,
    ]
    assert frontier["dependency_semantics"]["downstream_only_issues"] == [334]
    mutated = load_frontier()
    mutated["dependency_semantics"]["hard_issue_dependencies"].append(32)
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_issue_334_cannot_move_upstream() -> None:
    mutated = load_frontier()
    mutated["dependency_semantics"]["hard_issue_dependencies"] = [334]
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_alternative_routes_cannot_be_disabled() -> None:
    mutated = load_frontier()
    mutated["dependency_semantics"]["alternative_routes_allowed"] = False
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_cesium_candidate_cannot_become_required() -> None:
    mutated = load_frontier()
    mutated["candidate_routes"][0]["required_for_issue_633"] = True
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_optional_cesium_component_cannot_feed_generic_clock_gap() -> None:
    mutated = load_frontier()
    mutated["provenance_dag"]["edges"].append(
        {
            "from": "R_U",
            "to": "dimensionless_clock_gap",
        }
    )
    rehash(mutated)
    assert_rejected(mutated, "optional candidate component feeds the generic clock gap")


def test_route_neutral_gate_owner_cannot_be_reassigned() -> None:
    mutated = load_frontier()
    mutated["open_gates"][0]["owner_issues"] = [32]
    rehash(mutated)
    assert_rejected(mutated, "open gates do not match|route-neutral gate ownership")


def test_source_payload_injection_fails_closed() -> None:
    mutated = load_frontier()
    mutated["source_payloads"] = [
        {
            "component_id": "R_U",
            "source_map": "target_residual_fit",
            "positive_interval": ["1", "2"],
        }
    ]
    mutated["target_firewall"]["source_payload_count"] = 1
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_diagnostic_cannot_become_source_ancestor() -> None:
    mutated = load_frontier()
    mutated["provenance_dag"]["edges"].append(
        {
            "from": "legacy_gravity_checksum_skeleton",
            "to": "dimensionless_clock_gap",
        }
    )
    rehash(mutated)
    assert_rejected(mutated, "diagnostic artifact has a path")


def test_diagnostic_ancestry_flag_cannot_be_enabled() -> None:
    mutated = load_frontier()
    mutated["diagnostic_reference_pins"][0]["source_ancestry_allowed"] = True
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_newton_g_composition_is_downstream_only() -> None:
    frontier = load_frontier()
    assert frontier["downstream_consumers"] == [
        {
            "issue": 334,
            "role": "newton_g_composition",
            "required_input": "source_energy_interval",
            "possible_output": "source_g_si_interval",
            "status": "downstream_optional_context_not_source_clock_evidence",
        }
    ]
    mutated = load_frontier()
    mutated["provenance_dag"]["edges"].append(
        {
            "from": "newton_g_composition",
            "to": "source_energy_interval",
        }
    )
    rehash(mutated)
    assert_rejected(mutated, "cyclic|source-energy ancestor")


def test_physical_output_flags_cannot_be_enabled() -> None:
    for key in (
        "dimensionless_clock_gap_emitted",
        "source_energy_interval_emitted",
        "source_g_si_interval_emitted",
        "physical_promotion_allowed",
    ):
        mutated = load_frontier()
        mutated[key] = True
        rehash(mutated)
        assert_rejected(mutated, "schema validation")


def test_exact_si_definition_drift_fails_closed() -> None:
    mutated = load_frontier()
    mutated["exact_unit_chart"]["nu_Cs_Hz"] = "9192631771/1"
    rehash(mutated)
    assert_rejected(mutated, "exact SI definition drifted")


def test_interval_endpoint_tamper_fails_closed() -> None:
    mutated = load_frontier()
    mutated["interval_inversion_receipt"]["energy_interval_J"]["lower"] = "1/1"
    rehash(mutated)
    assert_rejected(mutated, "joule interval inversion mismatch")


def test_policy_pin_tamper_fails_closed() -> None:
    mutated = load_frontier()
    mutated["policy_pin"]["byte_sha256"] = "0" * 64
    rehash(mutated)
    assert_rejected(mutated, "policy byte digest mismatch")


def test_diagnostic_byte_pin_tamper_fails_closed() -> None:
    mutated = load_frontier()
    mutated["diagnostic_reference_pins"][1]["byte_sha256"] = "f" * 64
    rehash(mutated)
    assert_rejected(mutated, "diagnostic byte digest mismatch")


def test_diagnostic_path_traversal_fails_closed() -> None:
    mutated = load_frontier()
    mutated["diagnostic_reference_pins"][0]["path"] = (
        "code/particles/hierarchy/../runs/clock/cs133_feshbach_scalarization.json"
    )
    rehash(mutated)
    assert_rejected(mutated, "diagnostic path mismatch|path traversal")


def test_acceptance_cannot_be_marked_complete() -> None:
    mutated = load_frontier()
    mutated["acceptance_progress"][0]["status"] = "complete"
    rehash(mutated)
    assert_rejected(mutated, "schema validation")


def test_frontier_digest_tamper_fails_closed() -> None:
    mutated = copy.deepcopy(load_frontier())
    mutated["frontier_digest"] = "a" * 64
    assert_rejected(mutated, "frontier digest mismatch")
