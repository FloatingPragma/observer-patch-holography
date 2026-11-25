"""Focused tests for the generated OPH null-model scorecard."""

from __future__ import annotations

import sys
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


def test_scorecard_check_detects_drift_and_supports_explicit_path(
    tmp_path: Path,
) -> None:
    report = {}
    expected = "# scorecard\n\n## W3a\n"
    path = tmp_path / "scorecards" / "null_model_scorecard.md"
    path.parent.mkdir()
    path.write_text(expected, encoding="utf-8")
    nulls.check_scorecard(path, expected, report)
    path.write_text(expected + "tampered\n", encoding="utf-8")
    with pytest.raises(nulls.NullModelError, match="has drifted"):
        nulls.check_scorecard(path, expected, report)

    args = nulls.parse_args(
        ["--root", str(tmp_path), "--output", "scorecards/custom.md", "--check"]
    )
    assert args.root == tmp_path
    assert args.output == Path("scorecards/custom.md")
    assert args.check is True
