"""Cross-surface gates for the V3.37 fifth wave: the certain-state selection
of the Lueders instrument, the sourced field consistency criterion, and the
complex irreducibility of the golden pieces."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

INSTR_CLAIM = "OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES"
GAUSS_CLAIM = "OPH-EM-SOURCED-FIELD-CONSISTENCY"
CPLX_CLAIM = "OPH-EM-GOLDEN-SECTOR-COMPLEX-IRREDUCIBILITY"

INSTR_LEAN = "Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean"
GAUSS_LEAN = "Lean/Screen/SeamChargeContinuity.lean"
CPLX_LEAN = "Lean/Screen/GoldenSectorComplexIrreducibility.lean"
HILB_CLAIM = "OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT"
HILB_LEAN = "Lean/Screen/FieldSectorEnergyInnerProduct.lean"
CURL_CLAIM = "OPH-EM-CURL-SECTOR-EIGENBASIS"
CURL_LEAN = "Lean/Screen/CurlSectorEigenbasis.lean"


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
    gates = {INSTR_CLAIM: [730], GAUSS_CLAIM: [740, 733], CPLX_CLAIM: [733, 728],
             HILB_CLAIM: [730, 733], CURL_CLAIM: [733, 730]}
    for claim_id, expected in gates.items():
        claim = _claim(claim_id)
        assert claim["gates"] == expected, claim_id
        assert claim["premise_dependencies"]["consumed"] == [], claim_id
        assert claim["claim_class"] == "conditional_implication", claim_id


def test_claims_keep_declared_scope() -> None:
    instr = _claim(INSTR_CLAIM)["statement"]
    for token in ("K = K P", "Lueders outcome map X -> P X P", "conversely",
                  "strictly weaker than certain-state invariance",
                  "declared operational hypothesis", "not source-selected",
                  "discharges none of PR-02, PR-03, PR-04, PR-52, PR-64, or PR-65"):
        assert token in instr, token
    gauss = _claim(GAUSS_CLAIM)["statement"]
    for token in ("if and only if the initial load is neutral",
                  "solvable exactly when the charge vanishes",
                  "neutral partner or a neutralising background",
                  "sign of the committed load relative to a physical charge",
                  "discharges none of PR-15, PR-53, PR-54"):
        assert token in gauss, token
    cplx = _claim(CPLX_CLAIM)["statement"]
    for token in ("extension of scalars to the complex numbers",
                  "complex Burnside span has dimension nine",
                  "traces 1 - phi and phi", "real type",
                  "discharges neither PR-52 nor PR-53"):
        assert token in cplx, token
    hilb = _claim(HILB_CLAIM)["statement"]
    for token in ("static gradient-amplitude subspace",
                  "gradient velocity direction remains non-null",
                  "does not itself supply the carrier's nineteen-vector basis",
                  "private-algebra Stone surface"):
        assert token in hilb, token
    curl = _claim(CURL_CLAIM)["statement"]
    for token in ("nineteen nonzero curl modes", "positive definite",
                  "coefficient state is recovered",
                  "not identified with a physical photon Hilbert space"):
        assert token in curl, token


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        INSTR_LEAN: ("kraus_mul_compl_eq_zero", "rank_one_decomposition_collinear",
                     "kraus_component_fixing_certain_states_eq_compress",
                     "instrument_component_fixing_certain_states_is_lueders",
                     "hasKrausFormWithEffect_and_fixesCertainStates_iff_eq_luedersOutcomeMap",
                     "instrument_with_projective_effects_fixing_certain_states_is_lueders",
                     "effect_does_not_select_but_certain_states_do",
                     "phaseInstrument_eq_lueders_iff",
                     "committed_instruments_certain_state_invariance_selects_lueders",
                     "classical_subkernel_fixing_certain_laws_is_restriction"),
        GAUSS_LEAN: ("hoppingCurrent_continuity", "family_continuity",
                     "gauss_propagates_under_sourced_ampere", "family_gauss_propagation",
                     "sourced_maxwell_solvable_iff", "single_hopping_charge_obstruction",
                     "family_obstruction", "crossing_continuity_table"),
        CPLX_LEAN: ("complex_irreducible", "finrank_SC", "no_unit_intertwiner",
                    "chi_differ_iff"),
        HILB_LEAN: ("modeForm_posDef_iff", "stepMatrix_energy_isometry",
                    "invariant_inner_product_unique_up_to_positive_scale",
                    "staggered_energy_fixes_scale", "mode_unitary_group_explicit_generator",
                    "orthogonal_family_hilbert_reading", "assembledInner_radical",
                    "assembledInner_bilinear_radical",
                    "degenerate_block_invariant_forms_not_unique"),
        CURL_LEAN: ("curl_eigen", "curl_orth", "curl_linearIndependent", "curl_span",
                    "curl_committed_members", "curl_grad_isCompl", "fullLam_admissible_iff",
                    "curl_sector_energy_hilbert_packet", "carrier_flow_full_curl"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path
        assert "by native_decide" not in text, relative_path
        assert "sorry" not in text.replace("sorry-free", ""), relative_path


def test_lean_modules_registered() -> None:
    assert "import EventAlgebra.InstrumentSelectionByCertainStates" in (
        ROOT / "Lean/EventAlgebra.lean").read_text(encoding="utf-8")
    lakefile = (ROOT / "Lean/lakefile.lean").read_text(encoding="utf-8")
    assert "`SeamChargeContinuity" in lakefile
    assert "`GoldenSectorComplexIrreducibility" in lakefile
    assert "`FieldSectorEnergyInnerProduct" in lakefile
    assert "`CurlSectorEigenbasis" in lakefile


def test_owner_papers_carry_the_results() -> None:
    assert "InstrumentSelectionByCertainStates" in _collapsed(
        "extra/machine_checked_finite_event_algebras.tex")
    assert "SeamChargeContinuity" in _collapsed(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex")
    assert "GoldenSectorComplexIrreducibility" in _collapsed(
        "paper/observers_are_all_you_need.tex")
    assert "FieldSectorEnergyInnerProduct" in _collapsed(
        "paper/observers_are_all_you_need.tex")
    assert "CurlSectorEigenbasis" in _collapsed(
        "paper/observers_are_all_you_need.tex")


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert INSTR_LEAN in rows["OL-C1"]["evidence"]
    assert GAUSS_LEAN in rows["OL-N1"]["evidence"]
    assert CPLX_LEAN in rows["OL-F2"]["evidence"]
    assert HILB_LEAN in rows["OL-C2"]["evidence"]
    assert CURL_LEAN in rows["OL-C2"]["evidence"]
    assert rows["OL-C2"]["status"] == "partial"
    assert rows["OL-N1"]["status"] == "owed"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert INSTR_LEAN in prows["PR-64"]["evidence"]
    assert GAUSS_LEAN in prows["PR-54"]["evidence"]
    assert CPLX_LEAN in prows["PR-53"]["evidence"]
    assert HILB_LEAN in prows["PR-15"]["evidence"]
    assert CURL_LEAN in prows["PR-15"]["evidence"]
    assert prows["PR-64"]["disposition"] != "discharged"
