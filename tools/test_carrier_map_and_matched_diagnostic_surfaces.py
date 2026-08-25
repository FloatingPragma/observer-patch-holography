"""Cross-surface gates for the carrier-map candidate, matched dark-sector
diagnostic, INS-03 v2 replay gate, and fiber-rate comparison rungs."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

MAP_CLAIM = "OPH-GEOMETRY-SCREEN-CARRIER-MAP-CANDIDATE"
DM_CLAIM = "OPH-DM-MATCHED-OBSERVABLE-DIAGNOSTIC"
REPLAY_CLAIM = "OPH-QUANTUM-INS03-V2-REPLAY-GATE"
YM_CLAIM = "OPH-YM-FIBER-RATE-COMPARISON"

MAP_LEAN = "Lean/Geometry/ScreenCarrierMapCandidate.lean"
DM_PY = "code/cosmology/rar_deep_regime/matched_observable_diagnostic.py"
DM_RECEIPT = "code/cosmology/rar_deep_regime/runtime/matched_observable_receipt.json"
REPLAY_PY = "code/phase_instrument_export/v2/transcript_replay_verifier.py"
YM_LEAN = "Lean/Screen/KogutSusskindFiberRateComparison.lean"


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
    gates = {MAP_CLAIM: [740], DM_CLAIM: [742, 751],
             REPLAY_CLAIM: [730, 737], YM_CLAIM: [743, 744]}
    for claim_id, expected in gates.items():
        assert _claim(claim_id)["gates"] == expected, claim_id


def test_carrier_map_claim_stays_a_candidate() -> None:
    statement = _claim(MAP_CLAIM)["statement"]
    for token in ("candidate throughout",
                  "no clause identifies a port with a physical direction",
                  "sixty-word certificate",
                  "declared selection among valid equivariant embeddings",
                  "not consumed and not discharged"):
        assert token in statement, token
    assert _claim(MAP_CLAIM)["premise_dependencies"]["open"] == [
        "PR-15", "PR-52", "PR-53"]


def test_matched_diagnostic_claim_is_direction_neutral() -> None:
    statement = _claim(DM_CLAIM)["statement"]
    for token in ("v_A^2 = v_obs^2 - v_bar^2",
                  "containing zero",
                  "proxy mismatch",
                  "no confirmation of OPH is claimed",
                  "not the preregistered joint likelihood"):
        assert token in statement, token


def test_matched_receipt_carries_the_quoted_numbers() -> None:
    receipt = json.loads((ROOT / DM_RECEIPT).read_text(encoding="utf-8"))
    blob = json.dumps(receipt)
    for token in ("8.75", "8.41", "115"):
        assert token in blob, token
    statement = _claim(DM_CLAIM)["statement"]
    for token in ("8.42e-11", "8.75e-11", "[-0.080, +0.046]", "-0.017"):
        assert token in statement, token


def test_replay_claim_and_package() -> None:
    statement = _claim(REPLAY_CLAIM)["statement"]
    for token in ("TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED",
                  "PRODUCER_AUTHENTICATION_UNIMPLEMENTED",
                  "fails closed",
                  "owner freezes it before any producer transcript"):
        assert token in statement, token
    base = ROOT / "code/phase_instrument_export/v2"
    for name in ("PRIMITIVE_GENERATION_SEMANTICS.md",
                 "transcript_replay_verifier.py",
                 "test_transcript_replay_verifier.py",
                 "sample_synthetic_transcript.json",
                 "AUTHENTICATED_BINDING_SPEC.md"):
        assert (base / name).exists(), name


def test_ym_claim_keeps_finite_scope() -> None:
    statement = _claim(YM_CLAIM)["statement"]
    for token in ("minimum fiber rate",
                  "uniformly in every fiber ratio",
                  "no mass gap, no Clay-problem step",
                  "exact open gap"):
        assert token in statement, token
    receipt = json.loads((ROOT / (
        "code/yang_mills/receipts/"
        "kogut_susskind_fiber_rate_instances.json")).read_text(
        encoding="utf-8"))
    assert receipt, "instance receipt must be nonempty"


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        MAP_LEAN: ("ScreenCarrierMapCandidate", "candidate_gram_bridge",
                   "candidate_dot_table"),
        YM_LEAN: ("OPH.KogutSusskindFiberRate", "two_mul_le_fiberRate",
                  "fiberRate_eq_two_mul_iff"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)


def test_owner_papers_carry_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("ScreenCarrierMapCandidate",
                  "declared rather than forced"):
        assert token in observers, token
    dark = _collapsed("cosmology/oph_dark_matter_paper.tex")
    for token in ("matched diagnostic removes the proxy mismatch",
                  "contains zero", "Consistency on seen data selects nothing"):
        assert token in dark, token
    ym = _collapsed("extra/yang_mills_gap_clay_problem.tex")
    for token in ("KogutSusskindFiberRateComparison",
                  "minimum fiber rate", "exact open step"):
        assert token in ym, token


def test_ledger_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert MAP_LEAN in rows["OL-N1"]["evidence"]
    assert "declared, not forced" in rows["OL-N1"]["notes"]
    assert rows["OL-N1"]["status"] == "owed"
    assert DM_PY in rows["OL-I3"]["evidence"]
    assert "neither confirms nor selects" in rows["OL-I3"]["notes"]
    assert rows["OL-I3"]["status"] == "owed"
    assert REPLAY_PY in rows["OL-C5"]["evidence"]
    assert rows["OL-C5"]["status"] == "partial"


def test_premise_rows_stay_open() -> None:
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    rows = {r["id"]: r for r in register["rows"]}
    assert MAP_LEAN in rows["PR-53"]["evidence"]
    assert "stays open and unconsumed" in rows["PR-53"]["notes"]
    for pid in ("PR-64", "PR-65"):
        assert REPLAY_PY in rows[pid]["evidence"], pid
        assert "stays open" in rows[pid]["notes"], pid
