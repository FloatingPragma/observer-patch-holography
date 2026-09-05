"""Cross-surface gates for the V3.33 third wave: the transported-charge force
law, the seam-step speed limit, the source-clock rate along worldlines, and
the golden-sector characters."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

FORCE_CLAIM = "OPH-GEOMETRY-TRANSPORTED-CHARGE-FORCE-LAW"
SPEED_CLAIM = "OPH-GEOMETRY-SEAM-STEP-SPEED-LIMIT"
CLOCK_CLAIM = "OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES"
GOLD_CLAIM = "OPH-EM-GOLDEN-SECTOR-CHARACTERS"

FORCE_LEAN = "Lean/Geometry/TransportedChargeForceLaw.lean"
SPEED_LEAN = "Lean/Geometry/SeamStepSpeedLimit.lean"
CLOCK_LEAN = "Lean/Geometry/SourceClockRateAlongWorldlines.lean"
GOLD_LEAN = "Lean/Screen/GoldenSectorCharacters.lean"
GOLD_BRIDGE_LEAN = "Lean/Screen/A5PortSixAxesBridge.lean"
GOLD_PSL_BRIDGE_LEAN = "Lean/Screen/PSL2F5SixAxesBridge.lean"


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
    assert _claim(FORCE_CLAIM)["gates"] == [740]
    assert _claim(SPEED_CLAIM)["gates"] == [733, 740, 736]
    assert _claim(CLOCK_CLAIM)["gates"] == [736, 739]
    assert _claim(GOLD_CLAIM)["gates"] == [733, 728]
    for claim_id in (FORCE_CLAIM, SPEED_CLAIM, CLOCK_CLAIM, GOLD_CLAIM):
        assert _claim(claim_id)["premise_dependencies"]["consumed"] == [], claim_id
    for claim_id in (FORCE_CLAIM, SPEED_CLAIM, CLOCK_CLAIM):
        assert "discharges none" in _claim(claim_id)["statement"], claim_id
    assert "discharges neither" in _claim(GOLD_CLAIM)["statement"]


def test_force_claim_keeps_declared_class() -> None:
    statement = _claim(FORCE_CLAIM)["statement"]
    for token in ("closed two-step variation", "q h E_(j+1)(e) with zero clock cost",
                  "-8 - q h E_(j+1)(e)", "zero circulation",
                  "every difference forgets the unit", "declared enrichment"):
        assert token in statement, token


def test_speed_claim_keeps_identification_declared() -> None:
    statement = _claim(SPEED_CLAIM)["statement"]
    for token in ("declared identification", "2 phi^2 = 3 + sqrt 5",
                  "at least two rests per crossing", "2/3 < h",
                  "approached and never attained", "thresholds are independent",
                  "nothing selects it"):
        assert token in statement, token


def test_clock_claim_keeps_rules_declared() -> None:
    statement = _claim(CLOCK_CLAIM)["statement"]
    for token in ("E c (tau - sqrt(tau^2 - 4))", "zero exactly on resting windows",
                  "half seams are not seams", "sqrt 5 / 3 at tau = 3",
                  "both readings are declared",
                  "Neither rule is selected by the source or the join"):
        assert token in statement, token


def test_golden_claim_separates_character_and_later_irreducibility_theorems() -> None:
    statement = _claim(GOLD_CLAIM)["statement"]
    for token in ("(1/20) Z[phi]", "entrywise Galois conjugate", "phi or 1 - phi",
                  "twelve each", "This character module proves the projector and character identities, not irreducibility",
                  "The later dedicated GoldenSectorIrreducibility module proves",
                  "The subsequent dedicated GoldenSectorComplexIrreducibility module proves complex scalar-extension irreducibility",
                  "exactly the repository's explicit six-axis PSL(2,F5) model",
                  "canonical SL(2,ZMod 5) to PSL(2,ZMod 5) center quotient",
                  "its image is exactly A5SixAxes.L60",
                  "does not identify PSL2F5 with abstract A5",
                  "nor transport the Golden sectors as typed PSL representations",
                  "Identification with the abstract A5 or icosahedral character table remains an inference outside Lean"):
        assert token in statement, token
    scope_if_false = _claim(GOLD_CLAIM)["scope_if_false"]
    for token in ("Failure of the separate PSL2F5 bridge removes only the abstract-cover interface",
                  "the projector, character, and later irreducibility results remain intact"):
        assert token in scope_if_false, token


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        FORCE_LEAN: ("transportedAction_closed_variation", "exchange_action_difference",
                     "exchange_stationary_iff", "roundTrip_action_difference",
                     "rest_not_localMin_zero_field", "static_interaction_difference",
                     "difference_forgets_unit"),
        SPEED_LEAN: ("window_timelike_incompatible", "threshold_gap",
                     "block_timelike_iff", "block_two_rests",
                     "two_rest_admissible_interval", "block_speed_sup_one",
                     "thresholds_independent", "incompatibility_from_identification"),
        CLOCK_LEAN: ("properLength_generated_counts", "accrual_difference",
                     "accruals_agree_iff_resting", "refine_not_generated",
                     "indexAccrual_not_refinementInvariant",
                     "index_accrual_compatible_iff_resting", "dilationFactor_lt_one",
                     "dilationFactor_tendsto_one", "lengthRule_return_eq_inv_dilation"),
        GOLD_LEAN: ("plus_add_minus", "plus_mul_minus", "normal_plus", "minus_eq_conj_plus",
                    "chi_plus_values", "chi_conj", "chi_differ_iff", "five_class_counts",
                    "character_norms_golden", "faceAct_comm_plus"),
        GOLD_BRIDGE_LEAN: ("quotient_respects_antipode", "axisRelabel",
                           "rowEquiv", "bridged_axis_eq_six_axis",
                           "every_six_axis_row_realized", "quotient_action_faithful"),
        GOLD_PSL_BRIDGE_LEAN: ("slToPsl_ker_center", "p1EquivSix_affine",
                              "p1EquivSix_infinity", "slProjectiveAction_ker",
                              "pslProjectiveAction_injective", "tClass_action",
                              "sClass_action", "center_eq_plus_minus_one",
                              "center_card_two", "pslToSix_range",
                              "psl_equiv_six_axis_group"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path


def test_psl_bridge_is_registered_and_exported() -> None:
    lake = (ROOT / "Lean/lakefile.lean").read_text(encoding="utf-8")
    umbrella = (ROOT / "Lean/Screen/OPHScreen.lean").read_text(encoding="utf-8")
    assert "`PSL2F5SixAxesBridge" in lake
    assert "import PSL2F5SixAxesBridge" in umbrella


def test_owner_paper_carries_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("TransportedChargeForceLaw", "SeamStepSpeedLimit",
                  "two rests per crossing", "SourceClockRateAlongWorldlines",
                  "GoldenSectorCharacters", "PSL2F5SixAxesBridge",
                  "canonical center quotient", "pointwise port bridge is not appended"):
        assert token in observers, token


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert FORCE_LEAN in rows["OL-N1"]["evidence"]
    assert SPEED_LEAN in rows["OL-N1"]["evidence"]
    assert SPEED_LEAN in rows["OL-H8"]["evidence"]
    assert CLOCK_LEAN in rows["OL-H8"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert FORCE_LEAN in prows["PR-54"]["evidence"]
    assert SPEED_LEAN in prows["PR-15"]["evidence"]
    assert CLOCK_LEAN in prows["PR-15"]["evidence"]
    assert GOLD_LEAN in prows["PR-53"]["evidence"]
    assert GOLD_LEAN in rows["OL-F2"]["evidence"]
    assert prows["PR-15"]["disposition"] == "import"
