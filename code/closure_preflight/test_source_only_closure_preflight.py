from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Callable

import pytest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUTPUT = HERE / "outputs/source_only_closure_preflight.json"
RC_OUTPUT = HERE / "outputs/rc_load_source_only_certificate.json"
CONTRACT = HERE / "data/source_only_closure_contract_v3.json"
PRODUCER = HERE / "source_only_closure_preflight.py"
RC_PRODUCER = HERE / "rc_load_source_only.py"
VERIFIER = HERE / "verify_source_only_closure_preflight_independent.py"

sys.path.insert(0, str(HERE))
import source_only_closure_preflight as producer  # noqa: E402


def canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def reseal(value: dict) -> bytes:
    value.pop("certificate_sha256", None)
    value["certificate_sha256"] = "sha256:" + hashlib.sha256(
        canonical_json_bytes(value)
    ).hexdigest()
    return canonical_json_bytes(value)


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def row(value: dict, candidate_id: str) -> dict:
    return next(
        item for item in value["candidate_rows"] if item["candidate_id"] == candidate_id
    )


def test_generated_outputs_are_current_and_independently_valid() -> None:
    generated = run(str(PRODUCER), "--verify")
    assert generated.returncode == 0, generated.stdout + generated.stderr
    rc = run(str(RC_PRODUCER), "--verify")
    assert rc.returncode == 0, rc.stdout + rc.stderr
    independent = run(str(VERIFIER))
    assert independent.returncode == 0, independent.stdout + independent.stderr


def test_inventory_crosswalk_and_open_boundary_are_exact() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    required = certificate["qualification"]["required_source_only_gates"]
    completion = certificate["preflight_completion"]["gates"]
    assert certificate["schema"] == "oph.source_only_closure_preflight.v3"
    assert len(certificate["candidate_rows"]) == 30
    assert certificate["inventory_crosswalk"]["authoritative_menu_row_count"] == 17
    assert certificate["inventory_crosswalk"]["family_counts"] == {
        "CAP-B": 1,
        "CAP-K": 4,
        "CAP-L": 5,
        "CAP-P": 6,
        "coupled": 1,
    }
    assert len([item for item in certificate["candidate_rows"] if item["family"] == "CAP-P"]) == 6
    assert certificate["inventory_crosswalk"]["contract_candidate_order_used"] is False
    assert certificate["inventory_crosswalk"]["hierarchy_packet_count"] == 5
    assert certificate["status"] == producer.derive_progress_status(
        certificate["candidate_rows"], required, completion
    )
    assert certificate["scope"]["issue_closure_claimed"] is False
    assert certificate["scope"]["future_source_laws_exhausted"] is False
    assert certificate["scope"]["new_observational_comparison_opened"] is False
    assert certificate["qualification"]["qualifying_candidate_count"] == 0
    assert certificate["preflight_completion"]["all_complete"] is False


def test_status_rule_requires_row_and_global_completion() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    required = contract["source_only_required_gates"]
    completion_names = contract["preflight_completion_gates"]
    all_attained = {
        "row_kind": "candidate_map",
        "solver_target_payload_read": False,
        "exposure_class": "sealed_prospective",
        "exact_fixed_point_progress": True,
        "source_only_gates": {
            name: {
                "attained": True,
                "state": "ATTAINED",
                "evidence": [{"synthetic": True}],
            }
            for name in required
        },
    }
    incomplete = {
        name: {"attained": name == "decision_rule_frozen"}
        for name in completion_names
    }
    assert producer.derive_progress_status(
        [all_attained], required, incomplete
    ) == "OPEN_GATE_COMPLETE_CANDIDATE_BLOCKED_BY_PREFLIGHT_COMPLETENESS"
    complete = {name: {"attained": True} for name in completion_names}
    assert producer.derive_progress_status(
        [all_attained], required, complete
    ) == "SOURCE_ONLY_CANDIDATE_PRESENT"
    exposed = copy.deepcopy(all_attained)
    exposed["solver_target_payload_read"] = True
    assert producer.derive_progress_status(
        [exposed], required, complete
    ) == "OPEN_REGISTERED_CANDIDATES_FAIL_SOURCE_ONLY_GATES"
    no_fixed_point = copy.deepcopy(all_attained)
    no_fixed_point["source_only_gates"][required[0]]["attained"] = False
    no_fixed_point["exact_fixed_point_progress"] = False
    assert producer.derive_progress_status(
        [no_fixed_point], required, complete
    ) == "OPEN_NO_FIXED_POINT_EVIDENCE"


def test_source_only_rule_excludes_downstream_attachment() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    required = set(contract["source_only_required_gates"])
    downstream = set(contract["explicitly_downstream_gates"])
    assert required.isdisjoint(downstream)
    assert contract["scope_boundary"]["physical_or_laboratory_attachment_required"] is False
    assert contract["decision_rule"]["archive_enumeration_is_not_a_source_return_map"] is True
    assert contract["decision_rule"]["fixed_cutoff_controls_are_not_cosmic_fixed_points"] is True


def test_row_boundaries_are_fail_closed() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    archive = row(certificate, "n.capacity_map_archive_190")
    direct = row(certificate, "n.direct_correctable_record")
    coupled = row(certificate, "n.menu.coupled.cp1_cp2_cp3")
    assert archive["source_only_gates"]["source_return_map_complete"]["attained"] is False
    assert archive["details"]["source_return_map_constructed"] is False
    assert archive["solver_target_payload_read"] is True
    assert direct["row_kind"] == "fixed_cutoff_control"
    assert direct["details"]["fixed_cutoff_control"] is True
    assert direct["exact_fixed_point_progress"] is False
    assert certificate["progress"]["direct_D_equals_M0_control_counted_as_fixed_point_progress"] is False
    assert coupled["solver_target_payload_read"] is True
    assert coupled["exposure_class"] == "target_exposed_theorem_coupled"
    assert row(certificate, "hierarchy.p_public_endpoint")["solver_target_payload_read"] is True
    assert row(certificate, "hierarchy.p_source_audit_witness")["solver_target_payload_read"] is False
    assert row(certificate, "hierarchy.pn_product_contraction_adapter")["solver_target_payload_read"] is True
    assert row(certificate, "hierarchy.n_global_repair_tick_adapter")["solver_target_payload_read"] is True
    assert row(certificate, "hierarchy.n_ew_bridge_defined_capacity")["solver_target_payload_read"] is True


def test_formulas_and_rc_payload_are_exact() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rc = json.loads(RC_OUTPUT.read_text(encoding="utf-8"))
    assert rc["equation"] == producer.RC_EQUATION
    assert row(certificate, "n.rc_load")["details"]["equation"] == producer.RC_EQUATION
    assert row(certificate, "n.common_load_baseline")["details"]["formula"] == producer.COMMON_LOAD_FORMULA
    for branch_id, formula in producer.RESERVE_FORMULAS.items():
        assert row(certificate, f"n.reserve_{branch_id}")["details"]["formula"] == formula
    assert rc["scope"]["observational_comparison_imported"] is False
    assert rc["scope"]["physical_or_laboratory_attachment_evaluated"] is False
    assert rc["scope"]["workspace_root_checkout_required"] is False
    assert rc["verdict"] == "PASS"
    text = RC_PRODUCER.read_text(encoding="utf-8").lower()
    for forbidden in (
        "n_planck",
        "n_bao",
        "z_planck",
        "z_bao",
        "3.312593218645471",
        "3.2632655911878014",
    ):
        assert forbidden not in text


def mutate_status(value: dict) -> None:
    value["status"] = "SOURCE_ONLY_CANDIDATE_PRESENT"


def mutate_remove_candidate(value: dict) -> None:
    value["candidate_rows"].pop()


def mutate_gate(value: dict) -> None:
    gate = row(value, "p.thomson_structured_running")["source_only_gates"][
        "same_quantity_constructed"
    ]
    gate["attained"] = True
    gate["state"] = "ATTAINED"


def mutate_parent_pin(value: dict) -> None:
    value["parent_pins"][0]["sha256"] = "sha256:" + "0" * 64


def mutate_rc_root(value: dict) -> None:
    row(value, "n.rc_load")["details"]["fixed_point_X"]["lo"] = "1"


def mutate_archive_count(value: dict) -> None:
    row(value, "n.capacity_map_archive_190")["details"]["total_rows"] = 189


def mutate_attachment_gate(value: dict) -> None:
    value["scope"]["physical_or_laboratory_attachment_used_as_gate"] = True


def mutate_future_exhaustion(value: dict) -> None:
    value["scope"]["future_source_laws_exhausted"] = True


def mutate_new_comparison(value: dict) -> None:
    value["scope"]["new_observational_comparison_opened"] = True


def mutate_qualifier_list(value: dict) -> None:
    value["qualification"]["qualifying_candidate_ids"] = ["n.rc_load"]
    value["qualification"]["qualifying_candidate_count"] = 1


def mutate_evidence_value(value: dict) -> None:
    row(value, "p.thomson_structured_running")["source_only_gates"][
        "existence_certified"
    ]["evidence"][0]["observed"] = False


def mutate_wildcard_evidence(value: dict) -> None:
    row(value, "p.thomson_structured_running")["source_only_gates"][
        "existence_certified"
    ]["evidence"][0]["json_pointer"] = "/modes/*/banach/existence"


def mutate_false_gate_state(value: dict) -> None:
    row(value, "p.thomson_structured_running")["source_only_gates"][
        "same_quantity_constructed"
    ]["state"] = "ATTAINED"


def mutate_mathematical_status(value: dict) -> None:
    row(value, "n.rc_load")["mathematical_status"] = "UNCONDITIONAL"


def mutate_fixed_point_progress(value: dict) -> None:
    row(value, "n.direct_correctable_record")["exact_fixed_point_progress"] = True


def mutate_solver_target_payload_read(value: dict) -> None:
    row(value, "n.rc_load")["solver_target_payload_read"] = True


def mutate_exposure_class(value: dict) -> None:
    row(value, "n.rc_load")["exposure_class"] = "sealed_prospective"


def mutate_rc_equation_target_injection(value: dict) -> None:
    row(value, "n.rc_load")["details"]["equation"] += " using N_observed"


def mutate_common_load_formula(value: dict) -> None:
    row(value, "n.common_load_baseline")["details"]["formula"] = "N0 = pi * exp(6*pi/P)"


def mutate_reserve_formula(value: dict) -> None:
    row(value, "n.reserve_finite_presence")["details"]["formula"] = "N = N0"


def mutate_cap_p_count(value: dict) -> None:
    cap_p = next(
        item for item in value["candidate_rows"] if item["family"] == "CAP-P"
    )
    value["candidate_rows"].remove(cap_p)
    value["inventory_crosswalk"]["family_counts"]["CAP-P"] = 5


def mutate_archive_source_return(value: dict) -> None:
    gate = row(value, "n.capacity_map_archive_190")["source_only_gates"][
        "source_return_map_complete"
    ]
    gate["attained"] = True
    gate["state"] = "ATTAINED"


def mutate_completion_gate(value: dict) -> None:
    gate = value["preflight_completion"]["gates"]["input_catalogue_frozen"]
    gate["attained"] = True
    gate["state"] = "ATTAINED"


def mutate_hierarchy_disposition(value: dict) -> None:
    row(value, "hierarchy.pn_product_contraction_adapter")["disposition"] = "SOURCE_ONLY"


MUTATIONS: tuple[Callable[[dict], None], ...] = (
    mutate_status,
    mutate_remove_candidate,
    mutate_gate,
    mutate_parent_pin,
    mutate_rc_root,
    mutate_archive_count,
    mutate_attachment_gate,
    mutate_future_exhaustion,
    mutate_new_comparison,
    mutate_qualifier_list,
    mutate_evidence_value,
    mutate_wildcard_evidence,
    mutate_false_gate_state,
    mutate_mathematical_status,
    mutate_fixed_point_progress,
    mutate_solver_target_payload_read,
    mutate_exposure_class,
    mutate_rc_equation_target_injection,
    mutate_common_load_formula,
    mutate_reserve_formula,
    mutate_cap_p_count,
    mutate_archive_source_return,
    mutate_completion_gate,
    mutate_hierarchy_disposition,
)


@pytest.mark.parametrize("mutation", MUTATIONS, ids=lambda fn: fn.__name__)
def test_independent_verifier_rejects_resealed_mutations(
    tmp_path: Path, mutation: Callable[[dict], None]
) -> None:
    value = json.loads(OUTPUT.read_text(encoding="utf-8"))
    mutation(value)
    mutant = tmp_path / "mutant.json"
    mutant.write_bytes(reseal(value))
    result = run(str(VERIFIER), "--certificate", str(mutant))
    assert result.returncode != 0, result.stdout + result.stderr
