from __future__ import annotations

import json

import pytest

import horizon_record_attachment_verdict as verdict_module


def test_verdict_consumes_incomplete_direct_source_antecedent() -> None:
    verdict = verdict_module.build_verdict()
    assert verdict["status"] == "NOT_EVALUABLE_NO_HORIZON_RECORD_ATTACHMENT"
    assert verdict["issue"] == 589
    assert verdict["consumed_verdicts"]["direct_n_closure"] == (
        "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
    )
    boundary = verdict["comparison_boundary"]
    assert boundary["cosmological_payload_read"] is False
    assert boundary["lambda_comparison_permitted"] is False
    assert boundary["forecast_entry_permitted"] is False
    assert len(verdict["parent_pins"]) == 2


def test_runtime_verdict_is_byte_exact() -> None:
    committed = verdict_module.OUTPUT_PATH.read_bytes()
    assert committed == verdict_module.canonical_json_bytes(
        verdict_module.build_verdict()
    )


def test_direct_n_status_drift_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tampered = json.loads(
        verdict_module.DIRECT_N_VERDICT_PATH.read_text(encoding="utf-8")
    )
    tampered["status"] = "SOURCE_POSITIVE_UNIQUE_ZERO"
    tampered_path = tmp_path / "direct_n_closure_verdict.json"
    tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(verdict_module, "DIRECT_N_VERDICT_PATH", tampered_path)
    with pytest.raises(ValueError, match="status drift"):
        verdict_module.build_verdict()


def test_comparison_leak_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tampered = json.loads(
        verdict_module.DIRECT_N_VERDICT_PATH.read_text(encoding="utf-8")
    )
    tampered["comparison_boundary"]["cosmological_comparison_permitted"] = True
    tampered.pop("verdict_sha256", None)
    tampered["verdict_sha256"] = verdict_module.tagged_sha256(
        verdict_module.canonical_json_bytes(tampered)
    )
    tampered_path = tmp_path / "direct_n_closure_verdict.json"
    tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(verdict_module, "DIRECT_N_VERDICT_PATH", tampered_path)
    with pytest.raises(ValueError, match="permitted a comparison"):
        verdict_module.build_verdict()


def test_universal_membership_promotion_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tampered = json.loads(
        verdict_module.DIRECT_N_VERDICT_PATH.read_text(encoding="utf-8")
    )
    tampered["bounded_generation_register_result"][
        "universal_all_rung_membership_proved"
    ] = True
    tampered.pop("verdict_sha256", None)
    tampered["verdict_sha256"] = verdict_module.tagged_sha256(
        verdict_module.canonical_json_bytes(tampered)
    )
    tampered_path = tmp_path / "direct_n_closure_verdict.json"
    tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(verdict_module, "DIRECT_N_VERDICT_PATH", tampered_path)
    with pytest.raises(ValueError, match="promoted universal membership"):
        verdict_module.build_verdict()
