#!/usr/bin/env python3
"""Validate the exact non-hadron mass bundle."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_exact_nonhadron_mass_bundle.py"


def test_exact_nonhadron_mass_bundle_is_complete() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_exact_nonhadron_bundle_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        md = tmp / "EXACT_NONHADRON_MASSES.md"
        js = tmp / "exact_nonhadron_masses.json"
        forward = tmp / "exact_nonhadron_masses_current.json"
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
        assert payload["artifact"] == "oph_exact_nonhadron_mass_bundle"
        assert payload["status"] == "exact_output_lane_closed_nonhadron_only"
        entries = {entry["particle_id"]: entry for entry in payload["entries"]}
        assert len(entries) == 18
        assert entries["photon"]["mass_gev"] == pytest.approx(0.0, abs=1.0e-18)
        assert entries["w_boson"]["mass_gev"] == pytest.approx(80.377, abs=1.0e-12)
        assert entries["higgs"]["mass_gev"] == pytest.approx(125.1995304097179, abs=1.0e-12)
        assert entries["electron"]["mass_gev"] == pytest.approx(0.00051099895, abs=1.0e-15)
        assert entries["top_quark"]["mass_gev"] == pytest.approx(172.3523553288312, abs=1.0e-10)
        assert entries["tau_neutrino"]["mass_eV"] == pytest.approx(0.053081390655025595, abs=1.0e-15)
        markdown = md.read_text(encoding="utf-8")
        assert "Exact Non-Hadron Masses" in markdown
        assert "Bottom Quark" in markdown
        assert "Tau Neutrino" in markdown
