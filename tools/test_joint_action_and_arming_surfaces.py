"""Cross-surface gates for the joint-action, dispersion-arming,
source-clock, and export-validator rungs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

ACTION_CLAIM = "OPH-GEOMETRY-COMMON-WORLD-JOINT-ACTION"
ARM_CLAIM = "OPH-EM-DISPERSION-ARMING-INTERFACE"
CLK_CLAIM = "OPH-QFT-SOURCE-CLOCK-CANDIDATE"
VAL_CLAIM = "OPH-QUANTUM-INS03-EXPORT-VALIDATOR"

ACTION_LEAN = "Lean/Geometry/CommonWorldJointAction.lean"
ARM_LEAN = "Lean/Screen/DispersionArmingInterface.lean"
CLK_LEAN = "Lean/QFT/SourceClockCandidate.lean"
VAL_PY = "code/phase_instrument_export/ins03_export_validator.py"

ARMING_ROWS = tuple(f"PR-7{i}" for i in range(7))


def _registry() -> list[dict]:
    return yaml.safe_load((ROOT / "claims/claim_registry.yaml").read_text(
        encoding="utf-8"))["claims"]


def _claim(claim_id: str) -> dict:
    matches = [row for row in _registry() if row["claim_id"] == claim_id]
    assert len(matches) == 1, (claim_id, len(matches))
    return matches[0]


def _collapsed(relative_path: str) -> str:
    return " ".join((ROOT / relative_path).read_text(encoding="utf-8").split())


def test_claims_registered_with_gates() -> None:
    gates = {ACTION_CLAIM: [740], ARM_CLAIM: [733, 742],
             CLK_CLAIM: [736, 730], VAL_CLAIM: [730, 737]}
    for claim_id, expected in gates.items():
        assert _claim(claim_id)["gates"] == expected, claim_id


def test_action_claim_carries_direct_sum_boundary() -> None:
    statement = _claim(ACTION_CLAIM)["statement"]
    for token in ("unique stationary path", "formal direct sum",
                  "exactly decoupled", "action -2",
                  "supplies no shared carrier", "relative coefficient one"):
        assert token in statement, token


def test_arming_claim_and_register_rows() -> None:
    statement = _claim(ARM_CLAIM)["statement"]
    for token in ("PR-70 through PR-76", "schema-level checklist",
                  "stipulated mock", "not a universal no-go",
                  "nothing here arms or scores"):
        assert token in statement, token
    deps = _claim(ARM_CLAIM)["premise_dependencies"]
    assert deps["open"] == list(ARMING_ROWS)
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    rows = {r["id"]: r for r in register["rows"]}
    assert len(register["rows"]) == 76
    for pid in ARMING_ROWS:
        row = rows[pid]
        assert row["disposition"] in ("remove", "axiomatize"), pid
        assert 733 in row["consuming_lanes"] and 742 in row["consuming_lanes"]
        assert ARM_LEAN in row["evidence"], pid


def test_clock_claim_carries_exact_invariants() -> None:
    statement = _claim(CLK_CLAIM)["statement"]
    for token in ("197/1754", "94/1754", "103/1754", "1754",
                  "multiplicity", "does not replace or reduce PR-15"):
        assert token in statement, token
    assert _claim(CLK_CLAIM)["premise_dependencies"]["open"] == ["PR-15"]


def test_validator_claim_and_package() -> None:
    statement = _claim(VAL_CLAIM)["statement"]
    for token in ("0 <= E_i <= I", "renamed-context rejection",
                  "STATIC_COMMITTED_FIXTURE_CONFORMANT",
                  "PRODUCER_AUTHENTICATION_UNIMPLEMENTED"):
        assert token in statement, token
    base = ROOT / "code/phase_instrument_export"
    for name in ("ins03_export_validator.py", "test_ins03_export_validator.py",
                 "sample_conforming_export.json", "VALIDATOR_CONTRACT.md"):
        assert (base / name).exists(), name
    text = (base / "ins03_export_validator.py").read_text(encoding="utf-8")
    assert "PARTIAL_COMMITTED_CONTEXT_SET" in text


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        ACTION_LEAN: ("clockStationary_iff_uniform",
                      "worldline_is_stationary_path",
                      "jointAction_one_principle",
                      "jointAction_sector_decoupling",
                      "stationarity_not_minimality",
                      "jointAction_not_forcing_calibration"),
        ARM_LEAN: ("DispersionArmingInterface", "ArmedComparison",
                   "prop_fields_carry_no_separation"),
        CLK_LEAN: ("SourceClockCandidate", "197", "1754",
                   "identification_field_free"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)


def test_owner_papers_carry_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("CommonWorldJointAction", "direct sum",
                  "zero-coupling", "relative sector coefficient"):
        assert token in observers, token
    screen = _collapsed(
        "paper/screen_microphysics_and_observer_synchronization.tex")
    for token in ("schema-level arming checklist", "DispersionArmingInterface",
                  "natural-number digest labels"):
        assert token in screen, token


def test_ledger_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert ACTION_LEAN in rows["OL-N1"]["evidence"]
    assert "exactly decoupled sum" in rows["OL-N1"]["notes"]
    assert "Viable carrier-map" in rows["OL-N1"]["notes"]
    assert rows["OL-N1"]["status"] == "owed"
    assert ARM_LEAN in rows["OL-F4"]["evidence"]
    for pid in ARMING_ROWS:
        assert pid in rows["OL-F4"]["open_premises"], pid
    assert CLK_LEAN in rows["OL-H8"]["evidence"]
    assert VAL_PY in rows["OL-C5"]["evidence"]
