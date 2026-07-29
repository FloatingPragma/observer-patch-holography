"""Tests for the converted-chart receipts and the ST/Nielsen replay."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "checkers"))

CONVERTED = json.loads((ROOT / "outputs" / "fj_converted_equivalence.json").read_text(encoding="utf-8"))
VERDICT = json.loads((ROOT / "outputs" / "ward_st_nielsen_check.json").read_text(encoding="utf-8"))


def test_all_equivalence_receipts_pass() -> None:
    assert CONVERTED["problems"] == []
    for name, receipt in CONVERTED["equivalence_receipts"].items():
        assert receipt["passed"], name
    assert set(CONVERTED["equivalence_receipts"]) == {"ZZ", "WpWm", "hh", "G0G0", "GpGm"}


def test_chart_markers_distinct() -> None:
    markers = CONVERTED["chart"]["markers"]
    assert set(markers) == {"h", "hbar", "loop_order"}
    assert "gauge_fixing_v" in CONVERTED["chart"]


def test_st_replay_is_scoped_to_the_completed_projections() -> None:
    assert VERDICT["status"] == "PARTIAL_REPLAY_PASS__GHOST_ST_RESIDUAL_OPEN"
    assert VERDICT["replay_checks_passed"] is True
    assert VERDICT["acceptance_complete"] is False
    assert VERDICT["problems"] == []


def test_st_replay_covers_declared_projections() -> None:
    replays = VERDICT["replays"]
    for key in ("per_diagram_resummation", "photon_ward_full_expression",
                "charge_universality_nielsen", "custodial", "fj_equivalence",
                "ct_reachability"):
        assert key in replays, key
    assert replays["photon_ward_full_expression"] == "0"
    assert "gauge_condition_correlator_protection" in VERDICT["residual"]


def test_replay_checker_rejects_mutation(tmp_path, monkeypatch) -> None:
    import check_ward_st_nielsen as checker
    mutated = json.loads((ROOT / "outputs" / "fj_direct_vector_blocks.json").read_text(encoding="utf-8"))
    mutated["blocks"]["AA"]["diagrams"][0]["p_pole"] = "17*g2**2*p2"
    path = tmp_path / "vector.json"
    path.write_text(json.dumps(mutated), encoding="utf-8")
    monkeypatch.setattr(checker, "VECTOR_PATH", path)
    verdict = checker.check()
    assert verdict["status"] == "FAIL"
