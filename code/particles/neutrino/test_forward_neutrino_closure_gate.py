#!/usr/bin/env python3
"""The forward neutrino bundle must inherit every upstream failure."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "export_forward_neutrino_closure_bundle.py"
RUNS = ROOT / "particles" / "runs" / "neutrino"


def test_live_bundle_fails_closed_on_template_and_profile_rejection() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_closure_gate_") as tmpdir:
        output = pathlib.Path(tmpdir) / "bundle.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--weighted-cycle",
                str(RUNS / "neutrino_weighted_cycle_repair.json"),
                "--bridge-rigidity",
                str(RUNS / "neutrino_bridge_rigidity_theorem.json"),
                "--absolute-attachment",
                str(RUNS / "neutrino_absolute_attachment_theorem.json"),
                "--profile-score",
                str(RUNS / "nufit61_weighted_cycle_retrospective_score.json"),
                "--output",
                str(output),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(output.read_text(encoding="utf-8"))
        assert payload["prediction_promotion_allowed"] is False
        assert payload["public_surface_candidate_allowed"] is False
        assert payload["transitive_promotion_gate"]["weighted_cycle_prediction_promotion_allowed"] is False
        assert payload["transitive_promotion_gate"]["bridge_prediction_promotion_allowed"] is False
        assert payload["transitive_promotion_gate"]["absolute_attachment_prediction_promotion_allowed"] is False
        assert payload["transitive_promotion_gate"]["correlated_profile_rejected"] is True
        assert payload["legacy_absolute_fields_are_compare_only"] is True
        assert payload["pmns_status"] == "target_informed_template_candidate_rejected_by_nufit61_profile"
