"""Cross-surface gates for the carrier-map candidate, matched dark-sector
diagnostic, INS-03 v2 replay gate, and fiber-rate comparison rungs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
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
                  "positive refinement factor",
                  "declared selection among valid equivariant embeddings",
                  "not consumed and not discharged"):
        assert token in statement, token
    assert _claim(MAP_CLAIM)["premise_dependencies"]["open"] == [
        "PR-15", "PR-52", "PR-53"]


def test_matched_diagnostic_claim_is_direction_neutral() -> None:
    statement = _claim(DM_CLAIM)["statement"]
    for token in ("97-galaxy common set",
                  "containing zero",
                  "consistent with a mixed-proxy explanation",
                  "do not causally identify",
                  "not a calibrated likelihood or OPH confirmation"):
        assert token in statement, token


def test_matched_receipt_carries_the_quoted_numbers() -> None:
    receipt = json.loads((ROOT / DM_RECEIPT).read_text(encoding="utf-8"))
    primary = receipt["primary_fraction_0p1"]
    assert primary["channel_a_rar"]["a0_m_s2"] == pytest.approx(8.750402e-11)
    assert primary["paired_common_set"]["channel_a_common_a0_m_s2"] == pytest.approx(
        8.533413e-11)
    assert primary["channel_b_matched"]["a0_unweighted_log_mean_m_s2"] == pytest.approx(
        7.688961e-11)
    assert primary["paired_common_set"]["n_common_galaxies"] == 97
    assert primary["channel_b_matched"][
        "n_excluded_bulge_luminosity_ambiguous"] == 18
    statement = _claim(DM_CLAIM)["statement"]
    for token in ("8.750e-11", "8.533e-11", "7.689e-11",
                  "[-0.1166, +0.0264]", "-0.0453"):
        assert token in statement, token


def test_replay_claim_and_package() -> None:
    statement = _claim(REPLAY_CLAIM)["statement"]
    for token in ("TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED",
                  "PRODUCER_AUTHENTICATION_UNIMPLEMENTED",
                  "fails closed",
                  "must be owner-frozen before a producer transcript",
                  "larger remaining conjunction"):
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
                  "any mass-gap or Clay-problem step",
                  "committed non-product L=2,3 orbit laws"):
        assert token in statement, token
    receipt = json.loads((ROOT / (
        "code/yang_mills/receipts/"
        "kogut_susskind_fiber_rate_instances.json")).read_text(
        encoding="utf-8"))
    assert receipt, "instance receipt must be nonempty"


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        MAP_LEAN: ("ScreenCarrierMapCandidate", "candidate_gram_bridge",
                   "candidate_dot_table", "evalVec_ne_zero",
                   "baryCarrier_refine_meshRay"),
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
    for token in ("matched diagnostic removes the mixed observable",
                  "contains zero", "Consistency on seen data selects nothing",
                  "does not causally isolate"):
        assert token in dark, token
    ym = _collapsed("extra/yang_mills_gap_clay_problem.tex")
    for token in ("KogutSusskindFiberRateComparison",
                  "minimum fiber rate", "named residual steps"):
        assert token in ym, token


def test_ledger_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert MAP_LEAN in rows["OL-N1"]["evidence"]
    assert "exact equivariant carrier embedding candidate" in rows["OL-N1"]["notes"]
    assert "physical identifications remain declared" in rows["OL-N1"]["notes"]
    assert rows["OL-N1"]["status"] == "owed"
    assert DM_PY in rows["OL-I3"]["evidence"]
    assert "selects no model" in rows["OL-I3"]["notes"]
    assert "does not causally isolate" in rows["OL-I3"]["notes"]
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
