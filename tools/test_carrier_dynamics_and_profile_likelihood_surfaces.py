"""Cross-surface gates for the carrier-dynamics, profile-likelihood,
enlarged-target, and quotient-descent rungs."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]

DYN_CLAIM = "OPH-GEOMETRY-CARRIER-DYNAMICS-COMPATIBILITY"
LIK_CLAIM = "OPH-DM-ROTATION-CURVE-PROFILE-LIKELIHOOD"
ENL_CLAIM = "OPH-QFT-ENLARGED-TARGET-RECORD-EMBEDDING"

QUO_CLAIM = "OPH-YM-GAUGE-ORBIT-QUOTIENT-DESCENT"

DYN_LEAN = "Lean/Geometry/CarrierDynamicsCompatibility.lean"
ENL_LEAN = "Lean/QFT/EnlargedTargetRecordEmbedding.lean"
QUO_LEAN = "Lean/Screen/GaugeOrbitQuotientGap.lean"
LIK_DIR = "code/cosmology/rar_deep_regime/joint_likelihood"
LIK_RECEIPT = f"{LIK_DIR}/runtime/joint_likelihood_receipt.json"


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
    gates = {DYN_CLAIM: [740], LIK_CLAIM: [742, 751], ENL_CLAIM: [730, 743],
             QUO_CLAIM: [743, 744]}
    for claim_id, expected in gates.items():
        assert _claim(claim_id)["gates"] == expected, claim_id


def test_quotient_claim_keeps_invariant_subspace_scope() -> None:
    statement = _claim(QUO_CLAIM)["statement"]
    for token in ("descent-and-restriction",
                  "gauge-invariant subspace",
                  "orbit-space carrier formulation the named residual",
                  "moved-link ratios are pinned to one",
                  "no mass gap and no Clay-problem step"):
        assert token in statement, token


def test_dynamics_claim_stays_candidate_level() -> None:
    statement = _claim(DYN_CLAIM)["statement"]
    for token in ("compatibility of the candidate",
                  "uniformly in the step, so no step and no calibration",
                  "word-level matrix products are not assembled",
                  "so the coupling is declared",
                  "no frame identification is forced"):
        assert token in statement, token
    assert _claim(DYN_CLAIM)["premise_dependencies"]["open"] == [
        "PR-15", "PR-52", "PR-53"]


def test_likelihood_claim_is_conditional_and_direction_neutral() -> None:
    statement = _claim(LIK_CLAIM)["statement"]
    for token in ("no combined two-channel likelihood formed",
                  "7.9e-11 at rho = 0 to 5.8e-11",
                  "3.6 and 16.4",
                  "understates the observed scatter",
                  "consistency is shared with the standard null",
                  "The source does not fix a0"):
        assert token in statement, token


def test_likelihood_receipt_matches_quoted_numbers() -> None:
    receipt = json.loads((ROOT / LIK_RECEIPT).read_text(encoding="utf-8"))
    assert receipt["seen_data_postdiction"] is True
    assert receipt["physical_claim"] is False
    assert receipt["source_derived_output"] is False
    deep = receipt["subset_results"]["deep_f_0p1"]
    assert deep["n_points"] == 960 and deep["n_galaxies"] == 115
    by_rho = {row["rho"]: row for row in deep["per_rho"]}
    assert abs(by_rho[0.0]["a0_ml_m_s2"] - 7.94e-11) < 0.02e-11
    assert abs(by_rho[0.6]["a0_ml_m_s2"] - 5.82e-11) < 0.02e-11
    for row in deep["per_rho"]:
        paired = row["paired_btfr"]
        low, high = paired["log10_ratio_95pct"]
        assert low < 0.0 < high
        assert paired["verdict"] == "CONSISTENT_INTERVAL_CONTAINS_ZERO"


def test_enlarged_claim_keeps_finite_matrix_scope() -> None:
    statement = _claim(ENL_CLAIM)["statement"]
    for token in ("no AQFT vocabulary attached",
                  "single clause whose target algebra changes",
                  "proved of non-product form",
                  "both required",
                  "197/1754"):
        assert token in statement, token
    assert _claim(ENL_CLAIM)["premise_dependencies"]["open"] == [
        "PR-15", "PR-52"]


def test_lean_modules_carry_headline_declarations() -> None:
    expectations = {
        DYN_LEAN: ("DynamicsTransport", "ampere_intertwine",
                   "committed_perm_dynamics_compatibility",
                   "interaction_normalization_not_forced",
                   "transport_selection_not_forced"),
        ENL_LEAN: ("recordEmbedding", "enlargedDensity_ne_product",
                   "fixed_target_obstructed_enlarged_target_inhabited",
                   "no_product_state_bimodule_map_into_enlargedTarget"),
        QUO_LEAN: ("gaugeGroup", "generator_translate",
                   "dirichlet_bound_restricts",
                   "quotient_kogutSusskind_floor",
                   "instance_quotient_gap",
                   "distinct incident sites"),
    }
    for relative_path, tokens in expectations.items():
        text = (ROOT / relative_path).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, (relative_path, token)


def test_likelihood_package_files_exist() -> None:
    base = ROOT / LIK_DIR
    for name in ("joint_rar_likelihood.py",
                 "verify_joint_likelihood_independent.py",
                 "test_joint_rar_likelihood.py",
                 "runtime/joint_likelihood_receipt.json",
                 "LIKELIHOOD_CONVENTIONS.md"):
        assert (base / name).exists(), name
    conventions = _collapsed(f"{LIK_DIR}/LIKELIHOOD_CONVENTIONS.md")
    assert "no combined two-channel likelihood over both data channels" in \
        conventions
    assert "point-count dof convention" in conventions


def test_owner_papers_carry_the_results() -> None:
    observers = _collapsed("paper/observers_are_all_you_need.tex")
    for token in ("CarrierDynamicsCompatibility",
                  "declared rather than forced",
                  "EnlargedTargetRecordEmbedding".replace(
                      "EnlargedTargetRecordEmbedding",
                      "enlarged-target path is constructed")):
        assert token in observers, token
    dark = _collapsed("cosmology/oph_dark_matter_paper.tex")
    for token in ("Profile likelihood",
                  "dominates the systematic budget",
                  "understates the observed scatter",
                  "measurement-side interval"):
        assert token in dark, token
    ym = _collapsed("extra/yang_mills_gap_clay_problem.tex")
    for token in ("GaugeOrbitQuotientGap",
                  "descent-and-restriction",
                  "named residual steps"):
        assert token in ym, token


def test_ledger_rows_cite_without_promotion() -> None:
    ledger = json.loads((ROOT / "tracking/observation_ledger.json").read_text(
        encoding="utf-8"))
    rows = {row["id"]: row for row in ledger["rows"]}
    assert DYN_LEAN in rows["OL-N1"]["evidence"]
    assert rows["OL-N1"]["status"] == "owed"
    assert f"{LIK_DIR}/joint_rar_likelihood.py" in rows["OL-I3"]["evidence"]
    assert "dominates the systematic budget" in rows["OL-I3"]["notes"]
    assert rows["OL-I3"]["status"] == "owed"
    assert ENL_LEAN in rows["OL-C6"]["evidence"]
    assert rows["OL-C6"]["status"] == "partial"
    register = json.loads((ROOT / "tracking/premise_register.json").read_text(
        encoding="utf-8"))
    prows = {r["id"]: r for r in register["rows"]}
    assert DYN_LEAN in prows["PR-53"]["evidence"]
