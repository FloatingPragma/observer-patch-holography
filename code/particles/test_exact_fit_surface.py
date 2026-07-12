#!/usr/bin/env python3
"""Validate the exact-fits-only diagnostic surface."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_exact_fit_surface.py"


def test_exact_fit_surface_contains_only_exact_hits() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_exact_fit_surface_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        md = tmp / "EXACT_FITS_ONLY.md"
        js = tmp / "exact_fits_only.json"
        forward = tmp / "exact_fits_only_current.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--markdown-out",
                str(md),
                "--json-out",
                str(js),
                "--forward-out",
                str(forward),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(js.read_text(encoding="utf-8"))
        assert payload["artifact"] == "oph_exact_fits_only_surface"
        ids = {entry["id"] for entry in payload["entries"]}
        assert ids == {
            "higgs_top_reference_exact_adapter",
            "charged_current_family_exact_witness",
            "quark_current_family_exact_witness",
            "quark_selected_class_target_audit",
            "neutrino_two_parameter_exact_adapter",
            "neutrino_atmospheric_only_exact_adapter",
            "neutrino_solar_only_exact_adapter",
        }

        d11 = next(entry for entry in payload["entries"] if entry["id"] == "higgs_top_reference_exact_adapter")
        charged = next(entry for entry in payload["entries"] if entry["id"] == "charged_current_family_exact_witness")
        quark = next(entry for entry in payload["entries"] if entry["id"] == "quark_current_family_exact_witness")
        quark_selected_audit = next(
            entry for entry in payload["entries"] if entry["id"] == "quark_selected_class_target_audit"
        )
        neutrino_exact = next(entry for entry in payload["entries"] if entry["id"] == "neutrino_two_parameter_exact_adapter")
        assert d11["max_abs_residual"] == pytest.approx(0.0, abs=1.0e-12)
        assert charged["max_abs_residual"] == pytest.approx(0.0, abs=1.0e-12)
        assert quark["max_abs_residual"] == pytest.approx(0.0, abs=1.0e-10)
        assert quark_selected_audit["max_abs_residual"] == pytest.approx(0.0, abs=1.0e-10)
        assert quark_selected_audit["promotable"] is False
        assert quark_selected_audit["fit_kind"] == (
            "selected_class_target_anchored_mixed_convention_mass_texture_audit"
        )
        assert neutrino_exact["max_abs_residual"] == pytest.approx(0.0, abs=1.0e-18)

        markdown = md.read_text(encoding="utf-8")
        assert "Electroweak Frozen-Target Exact Pair" not in markdown
        assert "Higgs/Top Reference Exact Adapter" in markdown
        assert "Charged Current-Family Exact Witness" in markdown
        assert "Quark Selected-Class Target Audit" in markdown
        assert "Quark Current-Family Mixed-Convention Target Audit" in markdown
        assert "Neutrino Two-Parameter Exact Adapter" in markdown
        assert "Neutrino Atmospheric Only Exact Adapter" in markdown
        assert "Neutrino Solar Only Exact Adapter" in markdown
