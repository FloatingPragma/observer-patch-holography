from __future__ import annotations

import copy
import json

import pytest

import direct_n_closure_verdict as verdict_module
from capacity_indexed_source_family import canonical_json_bytes, tagged_sha256


def _write_payload(tmp_path, name: str, payload: dict):
    path = tmp_path / name
    path.write_bytes(canonical_json_bytes(payload))
    return path


def test_direct_n_verdict_consumes_the_exact_negative_source_result():
    verdict = verdict_module.build_verdict()

    assert verdict["status"] == (
        "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
    )
    assert verdict["issues"] == [505, 551]
    assert verdict["fixed_cutoff_result"] == {
        "D": 24,
        "M0": 24,
        "status": "SOURCE_DERIVED_FIXED_CUTOFF_PHYSICAL_PACKET",
        "cosmic_value_selected": False,
    }
    assert verdict["capacity_indexed_result"]["zero_sets_differ"] is True
    assert verdict["capacity_indexed_result"][
        "unique_source_zero_entailed"
    ] is False
    assert verdict["capacity_indexed_result"][
        "strange_loop_identity_rejected"
    ] is False
    assert verdict["source_controls"]["all_rung_lean_theorem"].endswith(
        "boundedCompletionClass_doesNotForceUniqueZero"
    )


def test_every_required_control_and_fiber_status_is_explicit():
    verdict = verdict_module.build_verdict()
    assert verdict["fiber_classifier"] == {
        "ambiguous": "AMBIGUOUS",
        "empty": "EMPTY",
        "incomplete": "INCOMPLETE",
        "singleton": "SINGLETON",
    }
    assert verdict["fiber_closure_branches"]["empty"][
        "existentially_closed"
    ] is False
    assert verdict["fiber_closure_branches"]["empty"][
        "robustly_closed"
    ] is False
    assert verdict["fiber_closure_branches"]["ambiguous"][
        "existentially_closed"
    ] is True
    assert verdict["fiber_closure_branches"]["ambiguous"][
        "robustly_closed"
    ] is False
    assert verdict["fiber_closure_branches"]["singleton"][
        "existentially_closed"
    ] is True
    assert verdict["fiber_closure_branches"]["singleton"][
        "robustly_closed"
    ] is True
    controls = verdict["source_controls"]
    assert controls["constructor_reads_desired_capacity"] is False
    assert all(
        value is True
        for key, value in controls.items()
        if key
        not in {
            "constructor_reads_desired_capacity",
            "all_rung_lean_theorem",
        }
    )
    assert verdict["comparison_boundary"]["direct_numeric_N_emitted"] is False
    assert verdict["comparison_boundary"][
        "cosmological_comparison_permitted"
    ] is False
    assert verdict["comparison_boundary"][
        "horizon_record_attachment_evaluable"
    ] is False


def test_parent_pin_or_verdict_mutation_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    family = json.loads(verdict_module.CERTIFICATE_PATH.read_text())
    family["status"] = "SOURCE_CLASS_IDENTIFIED"
    mutated = _write_payload(tmp_path, "family.json", family)
    monkeypatch.setattr(verdict_module, "CERTIFICATE_PATH", mutated)
    with pytest.raises(ValueError, match="not attained"):
        verdict_module.build_verdict()


def test_independent_projection_mismatch_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    receipt = json.loads(
        verdict_module.INDEPENDENT_RECEIPT_PATH.read_text()
    )
    receipt["projection_sha256"] = "sha256:" + "0" * 64
    mutated = _write_payload(tmp_path, "independent.json", receipt)
    monkeypatch.setattr(verdict_module, "INDEPENDENT_RECEIPT_PATH", mutated)
    with pytest.raises(ValueError, match="different projection"):
        verdict_module.build_verdict()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda family: family.__setitem__("zero_sets_differ", False),
            "zero_sets_differ drift",
        ),
        (
            lambda family: family["branch_receipts"][
                "reversible_identity"
            ].__setitem__("extension", []),
            "extension receipts are empty",
        ),
        (
            lambda family: family["branch_receipts"][
                "copy_collapse_erasure"
            ]["sewing"][0].__setitem__("status", "FAIL"),
            "sewing receipt did not pass",
        ),
    ],
)
def test_family_scientific_mutations_fail_closed(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
    mutation,
    message,
):
    family = json.loads(verdict_module.CERTIFICATE_PATH.read_text())
    mutation(family)
    family.pop("certificate_sha256", None)
    family["certificate_sha256"] = tagged_sha256(
        canonical_json_bytes(family)
    )
    mutated = _write_payload(tmp_path, "family-mutation.json", family)
    monkeypatch.setattr(verdict_module, "CERTIFICATE_PATH", mutated)
    with pytest.raises(ValueError, match=message):
        verdict_module.build_verdict()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda receipt: receipt.__setitem__("target_clean", False),
            "not target-clean",
        ),
        (
            lambda receipt: receipt["scope"].__setitem__(
                "producer_implementation_independent", False
            ),
            "scope drift",
        ),
        (
            lambda receipt: receipt.__setitem__("branch_ids_replayed", []),
            "branch coverage is incomplete",
        ),
    ],
)
def test_independent_replay_mutations_fail_closed(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
    mutation,
    message,
):
    receipt = json.loads(
        verdict_module.INDEPENDENT_RECEIPT_PATH.read_text()
    )
    mutation(receipt)
    mutated = _write_payload(tmp_path, "independent-mutation.json", receipt)
    monkeypatch.setattr(verdict_module, "INDEPENDENT_RECEIPT_PATH", mutated)
    with pytest.raises(ValueError, match=message):
        verdict_module.build_verdict()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda custody: custody.__setitem__("commit", "0" * 40),
            "custody commit is invalid",
        ),
        (
            lambda custody: custody.__setitem__("artifacts", []),
            "custody artifacts are empty",
        ),
        (
            lambda custody: custody["artifacts"][0].__setitem__(
                "sha256", "1" * 64
            ),
            "custody projection hash mismatch",
        ),
    ],
)
def test_custody_mutations_fail_closed(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
    mutation,
    message,
):
    custody = json.loads(verdict_module.INDEPENDENT_CUSTODY_PATH.read_text())
    mutation(custody)
    mutated = _write_payload(tmp_path, "custody-mutation.json", custody)
    monkeypatch.setattr(verdict_module, "INDEPENDENT_CUSTODY_PATH", mutated)
    with pytest.raises(ValueError, match=message):
        verdict_module.build_verdict()


def test_fixed_fiber_classification_mutation_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
):
    fixed = json.loads(verdict_module.FIXED_CERTIFICATE_PATH.read_text())
    fixed["controls"]["terminal_fibers"]["ambiguous"][
        "terminal_fiber_capacity_set"
    ] = [24]
    mutated = _write_payload(tmp_path, "fixed-mutation.json", fixed)
    monkeypatch.setattr(verdict_module, "FIXED_CERTIFICATE_PATH", mutated)
    with pytest.raises(ValueError, match="ambiguous fiber classifier drift"):
        verdict_module.build_verdict()


def test_runtime_verdict_is_canonical_and_byte_exact():
    expected = canonical_json_bytes(verdict_module.build_verdict())
    assert verdict_module.OUTPUT_PATH.read_bytes() == expected
