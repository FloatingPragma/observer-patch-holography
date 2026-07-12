#!/usr/bin/env python3
"""Validate the exact current-family quark affine-anchor theorem artifact."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
EXACT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_affine_anchor_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_affine_anchor_theorem.json"


def test_quark_current_family_affine_anchor_theorem() -> None:
    subprocess.run([sys.executable, str(EXACT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_affine_anchor_theorem"
    assert payload["proof_status"] == "closed_current_family_affine_anchor"
    assert payload["theorem_scope"] == "current_family_only"
    assert payload["public_promotion_allowed"] is False
    assert payload["target_anchored"] is True
    assert payload["source_only_prediction"] is False
    assert payload["single_running_quark_sextet_claim_allowed"] is False

    anchor = float(payload["current_family_affine_anchor"]["value"])
    split = float(payload["current_family_sector_split"]["value"])
    g_mean = float(payload["current_family_shared_sector_mean"]["value"])
    sector_means = payload["current_family_sector_means"]

    assert math.isclose(g_mean, math.exp(anchor), rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(float(sector_means["reconstructed_g_u"]), math.exp(anchor + split), rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(float(sector_means["reconstructed_g_d"]), math.exp(anchor - split), rel_tol=0.0, abs_tol=1e-15)
    for sector in ("u_sector", "d_sector"):
        for residual in payload["exact_fit_residuals"][sector]:
            assert math.isclose(float(residual), 0.0, rel_tol=0.0, abs_tol=1e-12)
