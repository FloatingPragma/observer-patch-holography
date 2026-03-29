#!/usr/bin/env python3
"""Guard the smaller fixed-eta affine germ beneath the open D10 tree identity."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fixed_eta_post_transport_affine_germ.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fixed_eta_post_transport_affine_germ.json"


def test_d10_fixed_eta_affine_germ_is_emitted() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_fixed_eta_post_transport_affine_germ"
    assert payload["object_id"] == "EWFixedEtaPostTransportAffineGerm_D10"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["strictly_smaller_than"] == "EWSinglePostTransportTreeIdentity_D10"
    assert payload["diagnostic_only"] is True
    assert payload["next_residual_object"] == "EWSinglePostTransportTreeIdentity_D10"
    anchor = payload["selected_anchor_point"]
    assert abs(anchor["sigma_EW"] + anchor["eta_EW"]) < 1.0e-15
