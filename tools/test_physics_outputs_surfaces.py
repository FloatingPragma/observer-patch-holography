"""Cross-surface gates for the fixed-capacity w-law, integer-k comb
instrument, and SI calibration import rungs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

WLAW_CLAIM = "OPH-DE-FIXED-CAPACITY-WLAW"
COMB_CLAIM = "OPH-GW-INTEGER-K-COMB-TEMPLATE"
CAL_CLAIM = "OPH-CAL-SI-ANCHOR-IMPORT"

WLAW_LEAN = "Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean"
COMB_LEAN = "Lean/Geometry/IntegerKCombInvariance.lean"
CAL_LEAN = "Lean/Thermodynamics/PhysicalCalibrationImport.lean"


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
    gates = {WLAW_CLAIM: [736, 742], COMB_CLAIM: [729], CAL_CLAIM: [736, 732]}
    for claim_id, expected in gates.items():
        claim = _claim(claim_id)
        assert claim["gates"] == expected, (claim_id, claim["gates"])


def test_wlaw_claim_carries_law_and_pending_status() -> None:
    statement = _claim(WLAW_CLAIM)["statement"]
    for token in ("(w0,wa)=(-1,0)", "w = -1 + (1/3) d ln N/d ln a",
                  "w >= -1", "not anchorable", "no OPH confirmation",
                  "seen diagnostics only"):
        assert token in statement, token


def test_comb_claim_carries_template_and_pending_status() -> None:
    statement = _claim(COMB_CLAIM)["statement"]
    for token in ("ln k/(8 pi)", "2 sqrt(1 - chi^2)/(1 + sqrt(1 - chi^2))",
                  "64 pi^2", "d_before = k d_after", "-ln k",
                  "M_det = (1+z) M_source", "registered_pending_freeze",
                  "Posterior samples alone",
                  "no event likelihood has been evaluated"):
        assert token in statement, token


def test_cal_claim_carries_anchor_exactness_and_import_boundary() -> None:
    statement = _claim(CAL_CLAIM)["statement"]
    for token in ("9192631770", "299792458", "6.62607015e-34",
                  "exactly one inhabitant", "m/(2 pi tau)",
                  "import, never derived"):
        assert token in statement, token
    deps = _claim(CAL_CLAIM)["premise_dependencies"]
    assert deps["classification"] == "explicit_edges"
    assert deps["open"] == ["PR-15"]


def test_frozen_register_rows_pending_with_owner_slots() -> None:
    register = json.loads((ROOT / "claims/frozen_prediction_register.json"
                           ).read_text(encoding="utf-8"))
    rows = {r["id"]: r for r in register["rows"]}
    for fz, issue in (("FZ-13", 742), ("FZ-14", 729)):
        row = rows[fz]
        assert row["status"] == "registered_pending_freeze", fz
        assert row["owning_issue"] == issue, fz
        assert row["frozen_utc"] is None and row["content_sha256"] is None, fz
        assert "owner" in row["kill_band"], fz
    assert "DESI DR1 and DR2" in rows["FZ-13"]["comparison_protocol"] or \
        "DESI DR1 and DR2" in rows["FZ-13"]["content"] or \
        "seen data" in rows["FZ-13"]["comparison_protocol"]
    assert "GW150914" in rows["FZ-14"]["comparison_protocol"]


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        WLAW_LEAN: ("w_eq_neg_one_of_constN", "w_eq_drift",
                    "w_ge_neg_one_of_monotone",
                    "deriv_N_neg_of_w_lt_neg_one", "cpl_forced_of_constN",
                    "decreasingCapacity_w_lt_neg_one"),
        COMB_LEAN: ("offsetSubtracted_ratio", "ratio_template_independent",
                    "pow_bound_32_lower", "log_ladder_32_lower",
                    "after_eq_before_div", "signed_entropy_change",
                    "positive_entropy_loss"),
        CAL_LEAN: ("SIAnchors", "ClockCalibration", "labFrequency",
                   "tick_declaration_not_forced", "9192631770"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)


def test_comb_code_package_complete_and_receipt_committed() -> None:
    base = ROOT / "code/gravitation/ringdown_comb"
    for name in ("integer_k_comb_template.py",
                 "verify_integer_k_comb_independent.py",
                 "test_integer_k_comb.py",
                 "REGISTRATION_CONTRACT_DRAFT.md",
                 "runtime/integer_k_comb_template_receipt.json"):
        assert (base / name).exists(), name
    receipt = json.loads((base / "runtime/integer_k_comb_template_receipt.json"
                          ).read_text(encoding="utf-8"))
    text = json.dumps(receipt)
    assert "imported continuation premise" in text
    assert receipt["schema"].endswith(".v3")
    assert receipt["frame_contract"]["detector_mass"] == (
        "M_det = (1+z)*M_source"
    )
    assert "posterior samples alone" in receipt["comparison_contract_boundary"].lower()
    assert "not a normalized cross-k" in receipt["imported_continuation_law"]["kms_scope"]
    assert receipt["imported_continuation_law"]["linewidth_scope"].startswith(
        "Gamma/Delta_E_k"
    )
    contract = (base / "REGISTRATION_CONTRACT_DRAFT.md").read_text(
        encoding="utf-8")
    for token in ("GW150914", "GW250114", "M_det=(1+z)M_source",
                  "Posterior samples alone", "source derivation open"):
        assert token in contract, token


def test_owner_papers_carry_the_results() -> None:
    desitter = _collapsed(
        "extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex")
    for token in ("The fixed-capacity equation of state", "(w_0,w_a)=(-1,0)",
                  "w(a)\\ge-1", "pending the owner's freeze"):
        assert token in desitter, token
    bh = _collapsed("cosmology/oph_black_hole_information_ledger.tex")
    for token in ("integer-\\(k\\) ringdown transition comb",
                  "x_k=\\ln k/(8\\pi)", "imported continuation premise",
                  "M_{\\rm det}=(1+z)M_{\\rm source}",
                  "not a cross-\\(k\\) transition probability"):
        assert token in bh, token
    technical = _collapsed("paper/tex_fragments/PAPER.tex")
    for token in ("signed black-hole entropy change is \\(-\\ln k\\)",
                  "M_{\\rm det}=(1+z)M_{\\rm source}",
                  "not a derived hierarchy of transition probabilities"):
        assert token in technical, token
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("cesium hyperfine", "one declared tick duration",
                  "an import and never a derivation"):
        assert token in observers, token


def test_ledger_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert WLAW_LEAN in rows["OL-H2"]["evidence"]
    assert "FZ-13" in rows["OL-H2"]["notes"]
    assert "frozen_targets" not in rows["OL-H2"]
    assert rows["OL-B4"]["frozen_targets"] == ["FZ-06", "FZ-14"]
    assert rows["OL-B4"]["status"] == "owed"
    assert "Posterior samples alone" in rows["OL-B4"]["notes"]
    assert "d_before=k d_after" in rows["OL-B4"]["notes"]
    assert CAL_LEAN in rows["OL-H8"]["evidence"]
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    pr15 = [r for r in register["rows"] if r["id"] == "PR-15"][0]
    assert CAL_LEAN in pr15["evidence"]
    assert pr15["disposition"] == "import"
