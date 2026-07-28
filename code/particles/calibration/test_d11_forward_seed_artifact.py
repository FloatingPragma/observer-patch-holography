#!/usr/bin/env python3
"""Guard the compact D11 forward-seed artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
MISSING_CERT = ROOT / "particles" / "runs" / "calibration" / "_missing_d11_forward_seed_promotion_certificate.json"


def test_d11_forward_seed_exports_single_scalar_candidate() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT), "--promotion-certificate", str(MISSING_CERT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_forward_seed"
    assert payload["forward_seed_object"] == "sigma_D11_HT"
    assert payload["seed_status"] == "declared_surface_fixed_ray_candidate"
    assert payload["predictive_promotion_allowed"] is False
    assert payload["source_surface_promotable"] is False
    assert payload["public_surface_candidate_allowed"] is False
    assert payload["comparison_surface_allowed"] is True
    assert payload["closure_state"] == "fixed_ray_algebra_only__source_promotion_open"
    assert payload["smallest_predictive_missing_object"] == "source_emitted_higgs_yukawa_fj_packet"
    assert payload["source_artifact"] == "oph_d11_declared_calibration_surface"
    assert payload["source_seed_law"]["sigma_D11_HT_formula"] == "alpha_u * cos(2*theta_W0) / sqrt(pi)"


def test_closed_status_without_explicit_source_gate_cannot_promote() -> None:
    forged = {
        "certificate_id": "forward_seed_promotion_certificate",
        "status": "closed",
        "predictive_promotion_allowed": True,
    }
    with tempfile.TemporaryDirectory(prefix="oph_d11_source_gate_") as tmpdir:
        certificate = pathlib.Path(tmpdir) / "forged_closed_certificate.json"
        output = pathlib.Path(tmpdir) / "forward_seed.json"
        certificate.write_text(json.dumps(forged), encoding="utf-8")
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--promotion-certificate",
                str(certificate),
                "--output",
                str(output),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(output.read_text(encoding="utf-8"))

    assert payload["predictive_promotion_allowed"] is False
    assert payload["source_surface_promotable"] is False
    assert payload["public_surface_candidate_allowed"] is False
