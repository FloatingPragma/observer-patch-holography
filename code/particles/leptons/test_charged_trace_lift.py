#!/usr/bin/env python3
"""Guard the fail-closed charged sector-isolated trace-lift theorem gate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
LEPTON_DIR = ROOT / "particles" / "leptons"
sys.path.insert(0, str(LEPTON_DIR))

from derive_charged_trace_lift import (  # noqa: E402
    CHARGED_BUDGET_JSON,
    D10_REDUCTION_JSON,
    DETERMINANT_FRONTIER_JSON,
    EXACT_READOUT_JSON,
    NORMALIZATION_NO_GO_JSON,
    SAME_LABEL_READBACK_JSON,
    build_artifact,
    classify_source_certificate,
)


SCRIPT = LEPTON_DIR / "derive_charged_trace_lift.py"


def _load(path: pathlib.Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _inputs() -> tuple[dict, dict, dict, dict, dict, dict]:
    return (
        _load(EXACT_READOUT_JSON),
        _load(SAME_LABEL_READBACK_JSON),
        _load(DETERMINANT_FRONTIER_JSON),
        _load(D10_REDUCTION_JSON),
        _load(NORMALIZATION_NO_GO_JSON),
        _load(CHARGED_BUDGET_JSON),
    )


def _complete_source_certificate() -> dict:
    return {
        "artifact": "oph_charged_trace_lift_source_certificate",
        "theorem_grade": True,
        "source_only": True,
        "no_target_leak": True,
        "ancestry": [
            "oph_d10_source_packet",
            "oph_charged_central_projector_receipt",
            "oph_charged_determinant_reference_receipt",
        ],
        "multiplicity_vector_M_ch": {"psi12": 1, "psi23": 1, "psi31": 1},
        "source_closed_stage_indexed_q": True,
        "physical_label_map": {
            "psi12": "electron",
            "psi23": "muon",
            "psi31": "tau",
        },
        "charged_central_projector": {"artifact": "oph_charged_central_projector_receipt"},
        "factorization_lemma": {
            "status": "certified",
            "leakage_bound": {"interval": ["0", "0"]},
        },
        "uncentered_lift_constant": {
            "source_object_name": "oph_charged_determinant_reference_receipt",
            "reference_stage": "r_0",
            "value": "-3",
            "interval": ["-3", "-3"],
        },
        "determinant_scalar_interval": ["-9", "-9"],
        "P_interval": ["1.6310", "1.6311"],
        "mass_space_affine_anchor_log_gev_interval": ["-3.1", "-3.0"],
        "mass_space_anchor_source_object": "oph_d10_electroweak_scale_conversion_receipt",
        "attachment_identity_residual": {
            "certified": True,
            "value": "0",
            "interval": ["0", "0"],
        },
        "anchor_bridge_verdict": "open",
    }


def test_current_corpus_emits_strengthened_no_go_and_ratio_checksum(tmp_path: pathlib.Path) -> None:
    output = tmp_path / "charged_trace_lift_theorem.json"
    missing_source = tmp_path / "absent_source_certificate.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--source-certificate",
            str(missing_source),
            "--output",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )
    payload = _load(output)

    assert payload["artifact"] == "oph_charged_trace_lift_theorem"
    assert payload["claim_label"] == "no_go_confirmed_new_source_needed"
    assert payload["factorization_lemma"]["leakage_bound"]["certified_zero"] is False
    assert payload["uncentered_lift_constant"]["source_object_name"] is None
    assert payload["uncentered_lift_constant"]["value"] is None
    residual = payload["attachment_identity_residual"]
    assert residual["interval"] == [None, None]
    assert residual["certified_zero"] is False
    assert residual["computable"] is False
    assert payload["promotion"]["conditional_on_P_allowed"] is False
    assert payload["conditional_mass_rows"] == []

    regression = payload["exact_ratio_regression"]
    assert regression["status"] == "passed"
    assert regression["theorem_scope"] == "current_family_only"
    assert regression["ratios"]["mu_over_e"]["witness"] == 206.76828270846718
    assert regression["max_relative_residual_abs"] <= 5.0e-15


def test_complete_source_certificate_is_the_only_certified_path() -> None:
    source = _complete_source_certificate()
    payload = build_artifact(*_inputs(), source)

    assert payload["claim_label"] == "trace_lift_certified"
    assert payload["factorization_lemma"]["leakage_bound"]["interval"] == ["0", "0"]
    assert payload["attachment_identity_residual"]["interval"] == ["0", "0"]
    assert payload["attachment_identity_residual"]["certified_zero"] is True
    assert payload["promotion"]["conditional_on_P_allowed"] is True
    assert payload["promotion"]["public_promotion_allowed"] is False
    assert len(payload["conditional_mass_rows"]) == 3
    assert payload["conditional_mass_rows"][0]["P_interval"] == ["1.6310", "1.6311"]
    assert payload["conditional_mass_rows"][0]["mass_interval"][0] < (
        payload["conditional_mass_rows"][0]["mass_interval"][1]
    )


def test_absent_lift_constant_fails_closed() -> None:
    source = _complete_source_certificate()
    source["uncentered_lift_constant"] = {}
    label, checks = classify_source_certificate(source)

    assert label != "trace_lift_certified"
    assert checks["lift_constant_source_named"] is False
    assert checks["lift_constant_numeric"] is False
    payload = build_artifact(*_inputs(), source)
    assert payload["promotion"]["conditional_on_P_allowed"] is False
    assert payload["attachment_identity_residual"]["certified_zero"] is False


def test_interval_containing_zero_is_not_a_zero_identity_certificate() -> None:
    source = _complete_source_certificate()
    source["attachment_identity_residual"]["interval"] = ["-1e-12", "1e-12"]
    label, checks = classify_source_certificate(source)

    assert label != "trace_lift_certified"
    assert checks["residual_singleton_zero"] is False


def test_source_open_q_or_nonzero_leakage_fails_closed() -> None:
    source = _complete_source_certificate()
    source["source_closed_stage_indexed_q"] = False
    source["factorization_lemma"]["leakage_bound"]["interval"] = ["0", "1e-9"]
    label, checks = classify_source_certificate(source)

    assert label != "trace_lift_certified"
    assert checks["source_closed_stage_indexed_q"] is False
    assert checks["zero_leakage"] is False
