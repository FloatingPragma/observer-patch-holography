#!/usr/bin/env python3
"""Validate the exact current-family quark PDG reconstruction theorem artifact."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
AFFINE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_affine_anchor_theorem.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_pdg_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_pdg_theorem.json"


def test_quark_current_family_exact_pdg_theorem() -> None:
    subprocess.run([sys.executable, str(AFFINE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_exact_pdg_theorem"
    assert payload["proof_status"] == (
        "closed_target_anchored_mixed_convention_coordinate_reconstruction"
    )
    assert payload["theorem_scope"] == "current_family_only"
    assert payload["public_promotion_allowed"] is False

    for residual in payload["exact_fit_residuals_gev"].values():
        assert math.isclose(float(residual), 0.0, rel_tol=0.0, abs_tol=1e-12)

    frontier = payload["next_target_free_bridge"]["remaining_public_frontier"]
    assert frontier == [
        "quark_source_spread_pair_action_breaking_theorem",
        "quark_rg_covariant_scheme_readout_or_invariant",
        "quark_common_scale_dimensionless_yukawa_certificate",
    ]
    assert "target-anchored current-family witness" in payload["next_target_free_bridge"]["note"]
    assert payload["target_anchored"] is True
    assert payload["source_only_prediction"] is False
    assert payload["single_running_quark_sextet_claim_allowed"] is False
