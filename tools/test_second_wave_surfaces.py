"""Cross-surface gates for the V3.32 second wave: worldline-to-hop transport,
the collar temperature reading, and the carrier mode equivariance."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

TRANS_CLAIM = "OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT"
COLLAR_CLAIM = "OPH-THERMO-COLLAR-TEMPERATURE-READING"
EQUI_CLAIM = "OPH-EM-CARRIER-MODE-EQUIVARIANCE"

TRANS_LEAN = "Lean/Geometry/WorldlineHopTransport.lean"
COLLAR_LEAN = "Lean/Thermodynamics/CollarTemperatureReading.lean"
EQUI_LEAN = "Lean/Screen/CarrierModeEquivariance.lean"

GRAV = "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex"


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
    assert _claim(TRANS_CLAIM)["gates"] == [740]
    assert _claim(COLLAR_CLAIM)["gates"] == [729, 736]
    assert _claim(EQUI_CLAIM)["gates"] == [733, 736, 728]
    for claim_id in (TRANS_CLAIM, COLLAR_CLAIM, EQUI_CLAIM):
        assert _claim(claim_id)["premise_dependencies"]["consumed"] == [], claim_id


def test_transport_claim_keeps_declared_scope() -> None:
    statement = _claim(TRANS_CLAIM)["statement"]
    for token in ("ray difference of its endpoints", "the ray map is injective",
                  "the agreement is endpoint-local",
                  "projection onto the stepped seam divided by 4",
                  "timelike exactly when tau^2 > 4",
                  "the unit is unconstrained by the field sector",
                  "discharges none"):
        assert token in statement, token
    assert _claim(TRANS_CLAIM)["premise_dependencies"]["open"] == [
        "PR-15", "PR-52", "PR-53", "PR-54"]


def test_collar_claim_is_an_exact_negative_with_dictionary() -> None:
    statement = _claim(COLLAR_CLAIM)["statement"]
    for token in ("pins no inverse temperature",
                  "realized branch of the declared maximum-entropy constraint",
                  "Gibbs state of K / beta at beta",
                  "exactly when K = beta H + constant",
                  "trade against the energy scale",
                  "are declared", "discharges neither"):
        assert token in statement, token
    assert _claim(COLLAR_CLAIM)["premise_dependencies"]["open"] == ["PR-15", "PR-52"]


def test_equivariance_claim_keeps_irrep_reading_as_inference() -> None:
    statement = _claim(EQUI_CLAIM)["statement"]
    for token in ("all 36000 triples", "an inference outside the theorems",
                  "is an observation", "character norms 60, 60, 60, 120",
                  "none discharged"):
        assert token in statement, token
    assert _claim(EQUI_CLAIM)["premise_dependencies"]["open"] == ["PR-15", "PR-52", "PR-53"]


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        TRANS_LEAN: ("seamVectorZ_eq_ray_difference", "ray_port_eq_spatialZ",
                     "candidateRayZ_injective", "transport_load_endpoints",
                     "transport_load_charge_fixed", "hoppingCurrent_eq_projection",
                     "transported_field_equations", "seam_step_timelike_iff",
                     "unit_unconstrained_by_field_sector"),
        COLLAR_LEAN: ("collarRef_eq_gibbs_one", "collar_beta_free",
                      "collar_beta_not_pinned", "collar_modular_hamiltonian_split",
                      "collar_kms_one", "collar_kms_beta", "collar_lab_not_forced"),
        EQUI_LEAN: ("incidence_equivariant", "localMaxwellOperator_pull",
                    "projector_images_invariant", "ampereEvolutionScaled_pull",
                    "fieldEnergyScaled_pull", "seamAct_comp", "facePerm_comp",
                    "order_counts", "projector_characters", "character_norms"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)
        assert "#print axioms" in text, relative_path


def test_owner_papers_carry_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("WorldlineHopTransport", "timelike exactly when",
                  "CarrierModeEquivariance", "split of the golden sector by eigenvalue is an observation"):
        assert token in observers, token
    gravity = _collapsed(GRAV)
    for token in ("CollarTemperatureReading",
                  "does not supply that inverse temperature"):
        assert token in gravity, token


def test_ledger_and_premise_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert TRANS_LEAN in rows["OL-N1"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"
    assert COLLAR_LEAN in rows["OL-B1"]["evidence"]
    assert rows["OL-B1"]["status"] == "partial"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    for pid in ("PR-52", "PR-53", "PR-54"):
        assert TRANS_LEAN in prows[pid]["evidence"], pid
    assert COLLAR_LEAN in prows["PR-15"]["evidence"]
    assert EQUI_LEAN in prows["PR-53"]["evidence"]
    assert EQUI_LEAN in rows["OL-F2"]["evidence"]
    assert prows["PR-15"]["disposition"] == "import"
