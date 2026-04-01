#!/usr/bin/env python3
"""Validate the theorem-side quark relative-sheet selector scaffold."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
REPAIR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
ORBIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sigma_ud_orbit.py"
SELECTOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_relative_sheet_selector.py"


def test_reference_singleton_orbit_emits_sigma_ref_when_uniqueness_theorem_closes() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_selector_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        repair_out = tmp / "repair.json"
        orbit_out = tmp / "orbit.json"
        selector_out = tmp / "selector.json"

        subprocess.run([sys.executable, str(REPAIR_SCRIPT), "--output", str(repair_out)], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(ORBIT_SCRIPT), "--output", str(orbit_out)], check=True, cwd=ROOT)
        subprocess.run(
            [
                sys.executable,
                str(SELECTOR_SCRIPT),
                "--repair-theorem",
                str(repair_out),
                "--orbit",
                str(orbit_out),
                "--output",
                str(selector_out),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(selector_out.read_text(encoding="utf-8"))
        assert payload["selection_rule_status"] == "theorem_grade_rule_emitted_on_current_solver_surface"
        assert payload["selection_status"] == "theorem_grade_value_emitted"
        assert payload["quark_relative_sheet_selector"]["value"]["sigma_id"] == "sigma_ref"
        assert payload["quark_relative_sheet_selector"]["value"]["canonical_token"] == "D12::same_label_left::reference_sheet"
        assert payload["quark_relative_sheet_selector"]["value"]["branch_key"] == ["D12", "sigma_ref"]
        assert payload["branch_key_after_repair"] == ["D12", "sigma_ref"]
        assert payload["debug_best_candidate"]["sigma_id"] == "sigma_ref"
        assert payload["debug_best_candidate_promotable"] is False
        assert payload["next_after_selection"]["id"] == "quark_d12_t1_value_law"


def test_theorem_witness_selects_exactly_one_sigma() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_selector_witness_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        repair_out = tmp / "repair.json"
        elements_json = tmp / "elements.json"
        orbit_out = tmp / "orbit.json"
        selector_out = tmp / "selector.json"

        elements_json.write_text(
            json.dumps(
                [
                    {
                        "sigma_id": "sig-a",
                        "canonical_token": "A",
                        "branch_key": ["D12", "sig-a"],
                        "ckm": {"theta_12": 0.01, "theta_23": 0.002, "theta_13": 0.00003, "delta_ckm": 1.0, "jarlskog": 1e-9},
                    },
                    {
                        "sigma_id": "sig-b",
                        "canonical_token": "B",
                        "branch_key": ["D12", "sig-b"],
                        "ckm": {"theta_12": 0.02, "theta_23": 0.003, "theta_13": 0.00004, "delta_ckm": 1.1, "jarlskog": 2e-9},
                        "selection_proof": {"theorem_grade_select": True},
                    },
                ],
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )

        subprocess.run([sys.executable, str(REPAIR_SCRIPT), "--output", str(repair_out)], check=True, cwd=ROOT)
        subprocess.run(
            [sys.executable, str(ORBIT_SCRIPT), "--elements-json", str(elements_json), "--output", str(orbit_out)],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(SELECTOR_SCRIPT),
                "--repair-theorem",
                str(repair_out),
                "--orbit",
                str(orbit_out),
                "--output",
                str(selector_out),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(selector_out.read_text(encoding="utf-8"))
        assert payload["selection_rule_status"] == "theorem_grade_rule_emitted_on_current_solver_surface"
        assert payload["selection_status"] == "theorem_grade_value_emitted"
        assert payload["quark_relative_sheet_selector"]["value"]["sigma_id"] == "sig-b"


def test_debug_ranking_is_never_promoted() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_selector_debug_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        repair_out = tmp / "repair.json"
        elements_json = tmp / "elements.json"
        orbit_out = tmp / "orbit.json"
        selector_out = tmp / "selector.json"

        elements_json.write_text(
            json.dumps(
                [
                    {
                        "sigma_id": "sig-close",
                        "canonical_token": "close",
                        "branch_key": ["D12", "sig-close"],
                        "ckm": {"theta_12": 0.2256, "theta_23": 0.0438, "theta_13": 0.00347, "delta_ckm": 1.0, "jarlskog": 1e-9},
                    }
                ],
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )

        subprocess.run([sys.executable, str(REPAIR_SCRIPT), "--output", str(repair_out)], check=True, cwd=ROOT)
        subprocess.run(
            [sys.executable, str(ORBIT_SCRIPT), "--elements-json", str(elements_json), "--output", str(orbit_out)],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(SELECTOR_SCRIPT),
                "--repair-theorem",
                str(repair_out),
                "--orbit",
                str(orbit_out),
                "--output",
                str(selector_out),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(selector_out.read_text(encoding="utf-8"))
        assert payload["selection_rule_status"] == "open_target_free_rule_unemitted"
        assert payload["debug_best_candidate_promotable"] is False
        assert payload["quark_relative_sheet_selector"]["value"] is None


def test_closed_singleton_branch_key_guard_rejects_mismatched_repair_theorem() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_quark_selector_guard_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        repair_out = tmp / "repair.json"
        orbit_out = tmp / "orbit.json"
        selector_out = tmp / "selector.json"

        subprocess.run([sys.executable, str(REPAIR_SCRIPT), "--output", str(repair_out)], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(ORBIT_SCRIPT), "--output", str(orbit_out)], check=True, cwd=ROOT)

        broken_repair = json.loads(repair_out.read_text(encoding="utf-8"))
        broken_repair["current_d12_sheet"]["branch_key"] = ["D12", None]
        repair_out.write_text(json.dumps(broken_repair, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        result = subprocess.run(
            [
                sys.executable,
                str(SELECTOR_SCRIPT),
                "--repair-theorem",
                str(repair_out),
                "--orbit",
                str(orbit_out),
                "--output",
                str(selector_out),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "closed singleton quark selector surface has inconsistent branch keys" in result.stderr
