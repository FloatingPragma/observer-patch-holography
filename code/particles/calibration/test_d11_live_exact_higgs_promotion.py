#!/usr/bin/env python3
"""Validate the exact D11 Higgs-only promotion theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_higgs_promotion.py"


def test_d11_live_exact_higgs_promotion_hits_exact_higgs_without_top_inverse_readback() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_d11_exact_higgs_") as tmpdir:
        out = pathlib.Path(tmpdir) / "d11_live_exact_higgs_promotion.json"
        subprocess.run([sys.executable, str(SCRIPT), "--output", str(out)], check=True, cwd=ROOT)
        payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d11_live_exact_higgs_promotion"
    assert payload["proof_status"] == "closed_target_anchored_live_exact_higgs_promotion"
    assert payload["theorem_id"] == "D11LiveForwardExactHiggsPromotion"
    assert payload["mass_readout"]["mH_gev"] == pytest.approx(125.1995304097179, abs=1.0e-12)
    assert payload["mass_readout"]["exact_residual_gev"] == pytest.approx(0.0, abs=1.0e-12)
    assert payload["exact_higgs_seed"]["value"] == pytest.approx(payload["exactifier"]["pi_H_exact"], abs=1.0e-15)
    assert "exact_d11_top_promotion_on_this_surface" in payload["strictly_not_claimed"]
    assert "full_higgs_top_inverse_slice_promotion" in payload["strictly_not_claimed"]
