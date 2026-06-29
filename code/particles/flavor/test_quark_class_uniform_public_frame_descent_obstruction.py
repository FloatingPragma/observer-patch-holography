#!/usr/bin/env python3
"""Guard the global public quark-frame descent obstruction."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_class_uniform_public_frame_descent_obstruction.py"
OFF_CANONICAL_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_off_canonical_p_evaluator_obstruction.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_class_uniform_public_frame_descent_obstruction.json"


def test_quark_class_uniform_public_frame_descent_obstruction_closes_issue_199_as_no_go() -> None:
    subprocess.run([sys.executable, str(OFF_CANONICAL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_class_uniform_public_frame_descent_obstruction"
    assert payload["github_issue"] == 199
    assert payload["issue_199_acceptance_met_as_obstruction"] is True
    assert payload["theorem_grade_global_descent"] is False
    assert payload["selected_class_descent"]["closed"] is True
    assert payload["selected_class_descent"]["source_only_sigma_emitted"] is False
    assert payload["selected_class_descent"]["theorem_scope"] == "selected_public_physical_quark_frame_class_only"
    assert payload["lane_closure_verdict"]["closure_kind"] == "hard_no_go_current_corpus"
    assert "oph_arbitrary_P_public_quark_frame_transport_classification" in payload["missing_global_objects"]
    assert payload["formal_countermodel_witness"]["free_even_sigma_family"]["name"] == (
        "canonical_vanishing_free_perturbation_family"
    )
    assert "promote_selected_fiber_constancy_to_all_public_quark_frame_classes" in payload["forbidden_promotions"]
