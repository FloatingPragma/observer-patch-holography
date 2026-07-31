"""Focused tests for the generated OPH null-model scorecard."""

from __future__ import annotations

import sys
from decimal import Decimal
from pathlib import Path

import pytest


TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import check_null_models as nulls  # noqa: E402


def _synthetic_constant_rows(hit_keys: set[str]) -> list[dict]:
    rows = []
    for c1, c2 in nulls.declared_constant_grid():
        numeric_key = f"{c1.equivalence_key}|{c2.equivalence_key}"
        canonical = c1.key == "phi" and c2.key == "pi"
        hit = canonical or numeric_key in hit_keys
        rows.append(
            {
                "canonical_pair": canonical,
                "numeric_pair_key": numeric_key,
                "threshold_results": {
                    threshold: {
                        "certified_inside": hit,
                        "certified_outside": not hit,
                    }
                    for threshold in ("1e-4", "1e-5", "2.5e-6")
                },
            }
        )
    return rows


def test_declared_grid_is_explicit_and_alias_aware() -> None:
    grid = nulls.declared_constant_grid()
    assert len(grid) == 48
    numeric_keys = {
        f"{c1.equivalence_key}|{c2.equivalence_key}" for c1, c2 in grid
    }
    assert len(numeric_keys) == 42
    assert sum(c1.key == "phi" and c2.key == "pi" for c1, c2 in grid) == 1
    assert {count for _, count in nulls.DELTAHEDRAL_MENU} == {
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        12,
    }


def test_w_f6_expression_grammar_is_exact_bounded_and_nonrecursive() -> None:
    atoms = nulls.expression_grammar_atoms(nulls.ROOT)
    grammar = nulls.declared_expression_grammar(nulls.ROOT)
    assert [atom.expression for atom in atoms] == [
        "P",
        "phi",
        "pi",
        "e",
        "sqrt(P)",
        "sqrt(pi)",
        "ln(2)",
        "P-phi",
        "sqrt(5)",
        "1",
        "2",
        "3",
        "4",
        "5",
    ]
    assert len(grammar) == 602
    assert sum(row.binary_operation_count == 0 for row in grammar) == 14
    assert sum(row.binary_operation_count == 1 for row in grammar) == 588
    assert all(row.binary_operation_count <= 1 for row in grammar)
    assert next(row.value for row in atoms if row.expression == "P") == Decimal(
        "1.63097209585889737696451390350695562847912625"
    )


def test_w_f6_random_targets_are_deterministic_and_open_interval() -> None:
    first = nulls.deterministic_expression_null_targets()
    second = nulls.deterministic_expression_null_targets()
    assert first == second
    assert len(first) == 2000
    assert len(set(first)) == 2000
    assert all(Decimal("0.6") < target < Decimal("0.8") for target in first)
    assert first[0] == Decimal(
        "0.665825579273475419769942260095341168835147982463240623474121093750"
    )


def test_w_f6_expression_null_rates_and_corrections_fail_closed() -> None:
    report = nulls.build_expression_grammar_null_model(nulls.ROOT)
    assert report["source_claimed_280_expression_rates_reproduced"] is False
    assert report["grammar"]["syntactic_expression_count"] == 602
    assert report["grammar"]["expressions_inside_random_target_interval"] == 34
    assert report["random_target_sampler"]["target_count"] == 2000
    assert {
        tolerance: row["random_target_hit_count"]
        for tolerance, row in report["thresholds"].items()
    } == {
        "0.20%": 773,
        "0.10%": 413,
        "0.05%": 209,
    }
    assert report["promotion_allowed"] is False
    assert report["evidential_weight_granted"] is False
    assert (
        report["claim_entry_gate"]["current_candidate_corrections_promotable"]
        is False
    )

    rows = {
        row["label"]: row
        for row in report["candidate_correction_diagnostics"]
    }
    bottom = rows["m_b multiplicative correction"]
    assert bottom["nearest_expression"] == "ln(2)"
    assert bottom["threshold_hits"] == {
        "0.20%": True,
        "0.10%": True,
        "0.05%": True,
    }
    strange = rows["m_s multiplicative correction"]
    assert strange["nearest_expression"] == "(2/3)"
    assert strange["threshold_hits"] == {
        "0.20%": True,
        "0.10%": False,
        "0.05%": False,
    }
    down = rows["m_d multiplicative correction"]
    assert not any(down["threshold_hits"].values())
    bridge = rows["m_s/m_d bridge factor"]
    assert bridge["inside_random_calibration_interval"] is False
    assert not any(bridge["threshold_hits"].values())
    assert all(row["promotion_allowed"] is False for row in rows.values())


def test_constant_interpretation_rule_uses_distinct_alternatives() -> None:
    rows = _synthetic_constant_rows(set())
    summary = nulls.summarize_constant_scan(rows)
    assert summary["declared_pair_count"] == 48
    assert summary["numerically_distinct_pair_count"] == 42
    assert summary["unique_alternative_pair_count"] == 41
    assert summary["duplicate_alias_pairs"] == 6
    assert summary["thresholds"]["2.5e-6"]["certified_alternative_hits"] == 0
    assert summary["interpretation_rule"]["triggered"] is False
    assert "NO_POSITIVE_WEIGHT" in summary["interpretation_rule"]["current_interpretation"]

    hit_keys = {
        "81/50|2",
        "33/20|e",
        "sqrt(e)|3",
    }
    triggered = nulls.summarize_constant_scan(
        _synthetic_constant_rows(hit_keys)
    )
    assert triggered["thresholds"]["2.5e-6"]["certified_alternative_hits"] == 3
    assert triggered["interpretation_rule"]["triggered"] is True
    assert "NO_EVIDENTIAL_WEIGHT" in triggered["interpretation_rule"]["current_interpretation"]


def test_canonical_pixel_root_has_interval_banach_certificate() -> None:
    pytest.importorskip("mpmath")
    source = (
        nulls.ROOT
        / "code/P_derivation/interval_contraction_certificate.py"
    )
    module = nulls._load_module(source, "_oph_test_null_model_interval")
    row = nulls._certify_constant_pair(
        module,
        nulls.C1_CHOICES[0],
        nulls.C2_CHOICES[0],
    )
    assert row["certificate"]["status"] == "INTERVAL_BANACH_ROOT_CERTIFIED"
    assert row["certificate"]["contraction"] is True
    assert row["certificate"]["self_map"] is True
    assert row["alpha_inv_point"].startswith("137.03566013694657")
    assert row["threshold_results"]["2.5e-6"]["certified_inside"] is True


def test_selector_menu_separates_axiom_domain_from_counterfactuals() -> None:
    report = nulls.build_selector_ablation(nulls.ROOT)
    summary = report["summary"]
    assert summary["declared_probe_menu_size"] == 8
    assert summary["configurations_with_executable_port_response_producer"] == 1
    assert summary["configurations_with_undeclared_port_response_producer"] == 7
    assert summary["sm_lie_type_available_count"] == 1
    assert summary["sm_lie_type_uniquely_selected_count"] == 0
    icosahedron = next(
        row for row in report["rows"] if row["configuration"] == "icosahedron"
    )
    assert (
        icosahedron["status"]
        == "AXIOM_DOMAIN_FORCED__COUNTERFACTUAL_MENU_UNTESTED"
    )
    assert icosahedron["physical_source_binding"] is False
    assert len(icosahedron["compatible_compact_lie_types"]) == 3
    scorecard = (nulls.ROOT / nulls.DEFAULT_SCORECARD).read_text(encoding="utf-8")
    assert "forced in the A1/A2 domain" in scorecard
    assert "physical source" in scorecard
    assert "binding and laboratory gauge-current identification remain open" in scorecard
    assert "does not contradict the A1/A2 abstract theorem" in scorecard
    assert report["matter_nonuniqueness"] == {
        "rank15_projector_verified": True,
        "inequivalent_completions": [
            "rank15_exterior_packet",
            "rank15_plus_sterile_singlet",
        ],
        "sterile_singlet_excluded": False,
        "completion_unique": False,
    }


def test_rscc_ablation_is_recomputed_and_beats_full_model() -> None:
    pytest.importorskip("mpmath")
    report = nulls.build_rscc_ablation(nulls.ROOT)
    assert report["status"] == "NEGATIVE_CONTROL_BEATS_FULL_RSCC"
    assert report["ablation_beats_full_maximum_error"] is True
    assert report["ablation_beats_full_residual_sum"] is True
    assert (
        report["zero_w2_zero_delta_g_ablation"][
            "max_abs_relative_error_percent"
        ]
        < report["full_rscc"]["max_abs_relative_error_percent"]
    )
    assert (
        report["zero_w2_zero_delta_g_ablation"]["raw_diagonal_residual_sum"]
        < report["full_rscc"]["raw_diagonal_residual_sum"]
    )


def test_clebsch_audit_enumerates_six_and_all_are_flag_rejected() -> None:
    report = nulls.build_quark_clebsch_audit(nulls.ROOT)
    assert report["menu"] == {
        "factor_set": ["1", "1/3", "3"],
        "assignment_slots": ["b_over_tau", "s_over_mu", "d_over_e"],
        "permutation_count": 6,
        "exhaustive": True,
    }
    assignments = {
        tuple(row["assignment"][slot] for slot in report["menu"]["assignment_slots"])
        for row in report["rows"]
    }
    assert len(assignments) == 6
    winners = [
        row
        for row in report["rows"]
        if row["retrospective_unique_least_discrepant"]
    ]
    assert len(winners) == 1
    assert winners[0]["assignment"] == {
        "b_over_tau": "1",
        "s_over_mu": "1/3",
        "d_over_e": "3",
    }
    metric = report["retrospective_target_informed_metric"]
    assert metric["current_assumed_order_uniquely_least_discrepant"] is True
    assert metric["target_informed"] is True
    assert metric["preregistered"] is False
    assert metric["physical_generation_order_selected"] is False
    assert metric["all_six_permutations_rejected"] is True
    assert "GENERATION_REGISTER_ORDER" in metric["remaining_open_premises"]
    assert all(
        row["conservative_flag_rejected_for_all_nf_rows"]
        for row in report["rows"]
    )
    pairing = report["pairing_and_weight_set_provenance"]
    assert pairing["constrains_candidate_channel_pairing"] is True
    assert pairing["equates_independent_yukawa_coefficients"] is False
    assert "CONDITIONAL_DECLARED_ALGEBRA" in pairing["F1_F2_status"]


def test_clebsch_rg_identity_is_scale_qualified_and_menu_wide() -> None:
    report = nulls.build_quark_clebsch_audit(nulls.ROOT)
    assert all(row["rg_ratio_identity"]["passes"] for row in report["rows"])
    adopted = next(row for row in report["rows"] if row["adopted_assignment"])
    assert adopted["rg_ratio_identity"]["factor_ratio_s_over_d"] == "1/9"
    identity = report["scale_qualified_rg_identity"]
    assert identity["passes_for_all_six_permutations"] is True
    assert Decimal(identity["adopted_register_scale_relative_gap"]) < Decimal(
        "1e-13"
    )
    assert identity["literal_low_scale_identity_exact"] is False
    assert Decimal(identity["literal_low_scale_relative_gap"]) > Decimal("1e-8")
    assert "generation-dependent" in identity["qualification"]
    assert "modified lane" in identity["qualification"]


def test_flag_comparison_uses_correct_inputs_and_conservative_gate() -> None:
    report = nulls.build_quark_clebsch_audit(nulls.ROOT)
    flag = report["flag_2024_compare_only"]
    assert flag["input_covariance_available"] is False
    assert flag["oph_theory_uncertainty_supplied"] is False
    assert flag["prediction_solve_or_physical_selection_input"] is False
    assert (
        flag["uncertainty_policy"][
            "conservative_rejection_gate_uses_rho_plus_one"
        ]
        is True
    )
    assert flag["rejection_gate_triggered_for_all_declared_nf_rows"] is True
    assert flag["prediction_preexisted_audit"] is True
    assert flag["significance_gate_preregistered"] is False
    assert flag["comparison_is_retrospective"] is True

    by_nf = {row["nf"]: row for row in flag["rows"]}
    nf_211 = by_nf["2+1+1"]
    assert nf_211["transcribed_inputs"]["ms_over_mud"][
        "published_notation"
    ] == "27.227(81)"
    assert nf_211["transcribed_inputs"]["mu_over_md"][
        "published_notation"
    ] == "0.465(24)"
    assert Decimal(nf_211["derived_ms_over_md"]) == Decimal("19.9437775")
    assert abs(
        Decimal(nf_211["independent_propagation"]["gap_sigma"])
        - Decimal("9.1260946")
    ) < Decimal("1e-7")
    assert abs(
        Decimal(
            nf_211[
                "maximally_positive_correlation_propagation"
            ]["gap_sigma"]
        )
        - Decimal("7.8498373")
    ) < Decimal("1e-7")

    nf_21 = by_nf["2+1"]
    assert Decimal(nf_21["derived_ms_over_md"]) == Decimal("20.35935")
    assert abs(
        Decimal(nf_21["independent_propagation"]["gap_sigma"])
        - Decimal("9.4981596")
    ) < Decimal("1e-7")
    assert abs(
        Decimal(
            nf_21[
                "maximally_positive_correlation_propagation"
            ]["gap_sigma"]
        )
        - Decimal("7.4799285")
    ) < Decimal("1e-7")


def test_koide_invariant_is_conditional_not_physically_promoted() -> None:
    report = nulls.build_quark_clebsch_audit(nulls.ROOT)
    koide = report["koide_and_sqrt_mass_invariant"]
    assert koide["mcpr_is_exact_two_thirds"] is False
    assert koide["runtime_charged_reference_consumed"] is False
    assert koide["finite_gns_public_koide_promotion_allowed"] is False
    assert "DECLARED_ARCHITECTURE" in koide["classification"]
    assert "not independently derived" in koide["stage5_audit_interpretation"]
    assert "historically target-informed" in koide["classification"]
    assert "blind/source-derived" in koide["classification"]


def test_rscc_win_must_be_disclosed_before_w3a() -> None:
    report = {"rscc_ablation": {"ablation_beats_full_model": True}}
    with pytest.raises(nulls.NullModelError, match="header does not disclose"):
        nulls.enforce_rscc_front_page_disclosure(
            report,
            "# scorecard\n\n## W3a — hidden warning\n",
        )
    nulls.enforce_rscc_front_page_disclosure(
        report,
        (
            "# scorecard\n\n"
            + nulls.RSCC_DISCLOSURE_MARKER
            + "\n\n## W3a — visible warning\n"
        ),
    )


def test_scorecard_check_detects_drift_and_supports_explicit_path(
    tmp_path: Path,
) -> None:
    report = {"rscc_ablation": {"ablation_beats_full_model": True}}
    expected = (
        "# scorecard\n\n"
        + nulls.RSCC_DISCLOSURE_MARKER
        + "\n\n## W3a\n"
    )
    path = tmp_path / "tracking" / "null_model_scorecard.md"
    path.parent.mkdir()
    path.write_text(expected, encoding="utf-8")
    nulls.check_scorecard(path, expected, report)
    path.write_text(expected + "tampered\n", encoding="utf-8")
    with pytest.raises(nulls.NullModelError, match="has drifted"):
        nulls.check_scorecard(path, expected, report)

    args = nulls.parse_args(
        ["--root", str(tmp_path), "--output", "tracking/custom.md", "--check"]
    )
    assert args.root == tmp_path
    assert args.output == Path("tracking/custom.md")
    assert args.check is True
