#!/usr/bin/env python3
"""Tests for the fail-closed carrier-mode / quantum-particle gate."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_carrier_mode_acceptance.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_carrier_mode_acceptance", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_group_only_and_classical_equation_only_do_not_pass_particle_gate() -> None:
    module = _load_module()
    group_only = {"abstract_symmetry_group_reconstructed": True}
    group_gate = module.evaluate_acceptance(group_only)
    assert group_gate["classical_carrier_gate"]["passed"] is False
    assert group_gate["quantum_particle_gate"]["passed"] is False
    assert group_gate["particle_promotion_allowed"] is False

    classical_only = {name: True for name in module.CLASSICAL_REQUIRED}
    classical_gate = module.evaluate_acceptance(classical_only)
    assert classical_gate["classical_carrier_gate"]["passed"] is True
    assert classical_gate["quantum_particle_gate"]["passed"] is False
    assert classical_gate["particle_promotion_allowed"] is False

    complete_conditional_receipt = {
        name: True for name in (*module.CLASSICAL_REQUIRED, *module.QUANTUM_REQUIRED)
    }
    complete_gate = module.evaluate_acceptance(complete_conditional_receipt)
    assert complete_gate["quantum_particle_gate"]["passed"] is True
    assert complete_gate["particle_promotion_allowed"] is True


def test_default_receipt_records_modes_but_no_quantum_particles() -> None:
    module = _load_module()
    payload = module.build_payload()
    assert payload["github_issue"] == 536
    assert payload["abstract_symmetry_group_alone_passes_quantum_gate"] is False
    carriers = {row["carrier_id"]: row for row in payload["carriers"]}
    assert set(carriers) == {"photon", "gluon", "graviton"}
    for carrier in carriers.values():
        assert carrier["hard_quadratic_mass_parameter_squared"] == 0
        assert "mass_gev" not in carrier
        assert carrier["classical_carrier_gate"]["passed"] is True
        assert "positive" in carrier["reduced_classical_hamiltonian"]
        assert carrier["quantum_particle_gate"]["passed"] is False
        assert carrier["particle_promotion_allowed"] is False
        assert carrier["abstract_symmetry_group_alone_sufficient"] is False
    assert "confining QCD" in carriers["gluon"]["phase_boundary"]
    assert "higher-derivative" in carriers["graviton"]["phase_boundary"]


def test_builder_is_deterministic(tmp_path: pathlib.Path) -> None:
    json_a, md_a = tmp_path / "a.json", tmp_path / "a.md"
    json_b, md_b = tmp_path / "b.json", tmp_path / "b.md"
    for json_out, md_out in ((json_a, md_a), (json_b, md_b)):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--json-out", str(json_out), "--markdown-out", str(md_out)],
            check=True,
            cwd=ROOT,
        )
    assert json_a.read_bytes() == json_b.read_bytes()
    assert md_a.read_bytes() == md_b.read_bytes()
    payload = json.loads(json_a.read_text(encoding="utf-8"))
    assert all(row["particle_promotion_allowed"] is False for row in payload["carriers"])
