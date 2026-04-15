#!/usr/bin/env python3
"""Validate the D10/gravity shared edge-entropy bridge artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FORWARD_CERTIFICATE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_forward_transmutation_certificate.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_gravity_shared_edge_entropy_bridge.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_gravity_shared_edge_entropy_bridge.json"


def test_d10_gravity_shared_edge_entropy_bridge_closes_product_branch_identity() -> None:
    subprocess.run([sys.executable, str(FORWARD_CERTIFICATE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_gravity_shared_edge_entropy_bridge"
    assert payload["status"] == "closed"
    assert payload["object_id"] == "D10GravitySharedEdgeEntropyBridge"
    assert payload["branch_data"]["realized_product_group_branch"] == "(SU(3) x SU(2) x U(1)) / Z6"
    assert abs(payload["branch_data"]["P"] - 1.63094) < 1.0e-12
    assert abs(payload["branch_data"]["ellbar_shared"] - 0.407735) < 1.0e-12
    assert (
        payload["theorem"]["formulas"]["shared_entropy_identity"]
        == "ellbar_shared = ellbar_SU3(t3_run) + ellbar_SU2(t2_run)"
    )
    assert payload["downstream_consequences"]["gravity_nat"] == "G_nat = a_cell / P"
    assert payload["downstream_consequences"]["gravity_si"] == "G_SI = c^3 * a_cell / (hbar * P)"
