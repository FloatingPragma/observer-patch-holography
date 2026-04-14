#!/usr/bin/env python3
"""Validate the shared-basis representation of the repaired weighted-cycle branch."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_weighted_cycle_shared_basis_representation.py"
WEIGHTED_CYCLE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
SHARED_CHARGED_LEFT = ROOT / "particles" / "runs" / "neutrino" / "shared_charged_lepton_left_basis.json"


def test_weighted_cycle_shared_basis_representation_closes_exactly() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_shared_basis_repr_") as tmpdir:
        out = pathlib.Path(tmpdir) / "representation.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--weighted-cycle",
                str(WEIGHTED_CYCLE),
                "--shared-charged-left",
                str(SHARED_CHARGED_LEFT),
                "--output",
                str(out),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(out.read_text(encoding="utf-8"))
        assert payload["artifact"] == "oph_neutrino_weighted_cycle_shared_basis_representation"
        assert payload["status"] == "theorem_grade_emitted"
        assert payload["physical_branch_closed"] is True
        assert payload["transport_checks"]["shared_basis_symmetry_max_abs"] < 1.0e-12
        assert payload["transport_checks"]["shared_basis_diagonalized_offdiag_max_abs"] < 1.0e-12
        assert payload["transport_checks"]["shared_basis_diagonalized_imag_max_abs"] < 1.0e-12
        assert payload["transport_checks"]["pmns_recovery_max_abs"] < 1.0e-12
        assert abs(payload["weighted_cycle_observables_match"]["theta12_deg_abs_delta"]) < 1.0e-10
        assert abs(payload["weighted_cycle_observables_match"]["theta23_deg_abs_delta"]) < 1.0e-10
        assert abs(payload["weighted_cycle_observables_match"]["theta13_deg_abs_delta"]) < 1.0e-10
        assert abs(payload["weighted_cycle_observables_match"]["delta_deg_abs_delta"]) < 1.0e-10
        assert abs(payload["emitted_parameters"]["alpha21_deg_0_to_360"] - 153.6185177794357) < 1.0e-9
        assert abs(payload["emitted_parameters"]["alpha31_deg_0_to_360"] - 257.0032408220805) < 1.0e-9
