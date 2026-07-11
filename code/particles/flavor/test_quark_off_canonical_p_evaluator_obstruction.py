#!/usr/bin/env python3
"""Guard the off-canonical quark P-evaluator obstruction certificate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
PUBLIC_SOURCE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_public_source_payload.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_off_canonical_p_evaluator_obstruction.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_off_canonical_p_evaluator_obstruction.json"


def test_quark_off_canonical_p_evaluator_obstruction_records_hard_blockers() -> None:
    subprocess.run([sys.executable, str(PUBLIC_SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_off_canonical_p_evaluator_obstruction"
    assert payload["proof_status"] == "off_canonical_p_evaluator_underdetermined_current_corpus"
    assert payload["theorem_grade_closure"] is False
    assert payload["selected_public_exact_surface"]["closed"] is False
    assert payload["selected_public_pure_B_payload"]["closed"] is True
    assert payload["lane_closure_verdict"]["closure_kind"] == "hard_no_go_current_corpus"
    assert payload["lane_closure_verdict"]["issue_212_acceptance_met"] is False
    assert payload["lane_closure_verdict"]["closed_theorem_grade_surface"] == (
        "selected_public_physical_quark_frame_class_only"
    )

    edge = payload["edge_statistics_sigma_lift_obstruction"]
    assert edge["missing_object"] == "oph_edge_statistics_sigma_lift"
    assert edge["bridge_status"] == "candidate_only"
    assert abs(edge["residuals"]["c_d_used_minus_fit"]) > 0.1
    assert abs(edge["residuals"]["sigma_d_candidate_minus_active"]) > 0.05
    canonical_point = payload["formal_countermodel_witness"]["canonical_point"]
    assert canonical_point["sigma_u_star"] == 5.579692209267639
    assert canonical_point["sigma_d_star"] == 3.300314452061615
    assert edge["active_closed_sigmas"]["sigma_u_total_log_per_side"] != canonical_point["sigma_u_star"]
    assert edge["active_closed_sigmas"]["sigma_d_total_log_per_side"] != canonical_point["sigma_d_star"]

    odd = payload["odd_response_scale_obstruction"]
    assert odd["missing_object"] == "oph_off_canonical_odd_response_kappa_value_law"
    assert odd["lift_parameterization_kind"] == "single_kappa_signed_projector_ray"
    assert len(odd["admissible_kappa_witnesses"]) == 3
    assert odd["admissible_kappa_witnesses"][0]["kappa"] == 0.0

    countermodels = payload["formal_countermodel_witness"]
    assert countermodels["nonidentifiability_result"].startswith("Because these countermodels agree")
    assert {item["name"] for item in countermodels["even_sigma_countermodels"]} == {
        "constant_selected_class_extension",
        "edge_candidate_response_extension",
        "canonical_vanishing_free_perturbation_family",
    }
    assert {item["name"] for item in countermodels["odd_kappa_countermodels"]} == {
        "zero_odd_extension",
        "canonical_vanishing_odd_response_family",
        "unit_representative_extension",
    }
    edge_candidate = next(
        item for item in countermodels["even_sigma_countermodels"] if item["name"] == "edge_candidate_response_extension"
    )
    assert abs(edge_candidate["canonical_residuals"]["sigma_d_candidate_minus_selected"]) > 0.1
    assert "arbitrary_P_public_quark_frame_classification_not_derived" in payload["closure_blockers"]
