#!/usr/bin/env python3
"""Validate the compare-only quark transport-frame diagnostic orbit."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_transport_frame_diagnostic_orbit.py"


def test_transport_frame_diagnostic_orbit_is_derived_and_non_promotable() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_transport_frame_") as tmpdir:
        out = pathlib.Path(tmpdir) / "transport_frame.json"
        subprocess.run([sys.executable, str(SCRIPT), "--output", str(out)], check=True, cwd=ROOT)
        payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["proof_status"] == "compare_only_not_sector_attached"
    assert payload["self_overlap"]["symbol"] == "F0^dagger F1"
    assert payload["line_lift_labels"] == ["f1", "f2", "f3"]
    assert payload["missing_sector_attachment"]["cannot_emit"] == [
        "sigma_id",
        "canonical_token",
        "U_u_left",
        "U_d_left",
        "V_CKM",
    ]
    assert payload["debug_log_shell_loss"]["transport_frame_self_overlap"] < payload["debug_log_shell_loss"]["current_same_sheet"]
    assert payload["debug_log_shell_loss"]["improvement_factor_vs_current_same_sheet"] > 20.0
    assert payload["residual_gauge_quotient"]["entrywise_moduli_invariant"] is True
    assert payload["residual_gauge_quotient"]["ckm_moduli_invariant"] is True
    assert payload["smallest_exact_blocker"]["id"] == "sector_attached_transport_frame_to_sigma_ud_lift_equivalence"
    assert math.isclose(payload["ckm_invariants"]["theta_12"], 0.05303513965374766, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(payload["ckm_invariants"]["theta_23"], 0.03505328791223491, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(payload["ckm_invariants"]["theta_13"], 0.004481306693226519, rel_tol=0.0, abs_tol=1.0e-15)
