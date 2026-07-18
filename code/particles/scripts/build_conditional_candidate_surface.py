#!/usr/bin/env python3
"""Build the conditional candidate output surface.

The strict source-only promotion policy reports `n/a` in every public
OPH-value column while the sector gates stay open. The numeric conditional
candidates exist in the run artifacts. This surface renders those
values in one place together with their recorded claim labels and open
gates. Every row stays `promotable=False`; the surface adds no new claims
and must not feed back into any promotion decision.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
D11_SPLIT_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json"
D10_VALUE_LAW_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
HIERARCHY_WITNESS_JSON = ROOT / "particles" / "hierarchy" / "computations" / "hierarchy_numeric_witness.json"
DEFAULT_MD_OUT = ROOT / "particles" / "CONDITIONAL_CANDIDATES.md"
DEFAULT_JSON_OUT = ROOT / "particles" / "conditional_candidates.json"
DEFAULT_FORWARD_OUT = ROOT / "particles" / "runs" / "status" / "conditional_candidates_current.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing required artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _require(payload: dict[str, Any], key: str, context: str) -> Any:
    if key not in payload:
        raise SystemExit(f"missing required field `{key}` in {context}")
    return payload[key]


def _rel(path: pathlib.Path) -> str:
    return path.relative_to(ROOT.parent).as_posix()


def build_gev_chart_candidates() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    split = _load_json(D11_SPLIT_PAIR_JSON)
    pair = _require(split, "exact_split_pair", _rel(D11_SPLIT_PAIR_JSON))
    rows.append(
        {
            "observable": "m_H",
            "value": _require(pair, "mH_gev", "exact_split_pair"),
            "unit": "GeV",
            "claim_label": "conditional_declared_surface_higgs_top_split_candidate",
            "chart": "declared_d10_d11_running_matching_threshold_surface",
            "open_gates": [
                "d10_target_free_repair_closure",
                "qt1_qt5_quotient_path_source_entailment",
                "ds1_ds5_d11_split_character_and_rigidity",
                "physical_complex_pole_conversion",
            ],
            "source_artifact": _rel(D11_SPLIT_PAIR_JSON),
            "promotable": False,
        }
    )
    rows.append(
        {
            "observable": "m_t",
            "value": _require(pair, "mt_pole_gev", "exact_split_pair"),
            "unit": "GeV",
            "claim_label": "conditional_declared_surface_higgs_top_split_candidate",
            "chart": "declared_d10_d11_running_matching_threshold_surface",
            "open_gates": [
                "d10_target_free_repair_closure",
                "qt1_qt5_quotient_path_source_entailment",
                "ds1_ds5_d11_split_character_and_rigidity",
                "physical_complex_pole_conversion",
                "top_threshold_control",
            ],
            "source_artifact": _rel(D11_SPLIT_PAIR_JSON),
            "promotable": False,
        }
    )

    value_law = _load_json(D10_VALUE_LAW_JSON)
    if _require(value_law, "status", _rel(D10_VALUE_LAW_JSON)) != "candidate_only":
        raise SystemExit(
            "d10 target-free repair value law left `candidate_only`; "
            "re-audit this surface before rebuilding"
        )
    quintet = _require(value_law, "coherent_emitted_quintet", _rel(D10_VALUE_LAW_JSON))
    for observable, key in (("M_W", "MW_pole"), ("M_Z", "MZ_pole")):
        rows.append(
            {
                "observable": observable,
                "value": _require(quintet, key, "coherent_emitted_quintet"),
                "unit": "GeV",
                "claim_label": "conditional_d10_target_free_value_law_chart_coordinate",
                "chart": "archived_d10_value_law_candidate",
                "open_gates": [
                    "qt1_qt5_quotient_path_source_entailment",
                    "physical_complex_pole_attachment",
                ],
                "source_artifact": _rel(D10_VALUE_LAW_JSON),
                "promotable": False,
            }
        )

    return rows


def build_scale_free_chart_ratios() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    value_law = _load_json(D10_VALUE_LAW_JSON)
    quintet = _require(value_law, "coherent_emitted_quintet", _rel(D10_VALUE_LAW_JSON))
    mw = _require(quintet, "MW_pole", "coherent_emitted_quintet")
    mz = _require(quintet, "MZ_pole", "coherent_emitted_quintet")
    for observable, value in (
        ("M_W / M_Z", mw / mz),
        ("sin2_theta_W_eff", _require(quintet, "sin2w_eff", "coherent_emitted_quintet")),
        ("alpha_em_eff_inv", _require(quintet, "alpha_em_eff_inv", "coherent_emitted_quintet")),
    ):
        rows.append(
            {
                "observable": observable,
                "value": value,
                "unit": "dimensionless",
                "claim_label": "conditional_d10_target_free_value_law_scale_free_coordinate",
                "chart": "archived_d10_value_law_candidate",
                "open_gates": [
                    "qt1_qt5_quotient_path_source_entailment",
                    "physical_complex_pole_attachment",
                ],
                "source_artifact": _rel(D10_VALUE_LAW_JSON),
                "promotable": False,
            }
        )

    split = _load_json(D11_SPLIT_PAIR_JSON)
    pair = _require(split, "exact_split_pair", _rel(D11_SPLIT_PAIR_JSON))
    rows.append(
        {
            "observable": "m_H / m_t",
            "value": _require(pair, "mH_gev", "exact_split_pair")
            / _require(pair, "mt_pole_gev", "exact_split_pair"),
            "unit": "dimensionless",
            "claim_label": "conditional_declared_surface_higgs_top_ratio",
            "chart": "declared_d10_d11_running_matching_threshold_surface",
            "open_gates": [
                "d10_target_free_repair_closure",
                "ds1_ds5_d11_split_character_and_rigidity",
                "physical_complex_pole_conversion",
            ],
            "source_artifact": _rel(D11_SPLIT_PAIR_JSON),
            "promotable": False,
        }
    )
    return rows


def build_dimensionless_candidates() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    # The source-root-plus-public-width coordinate is a diagnostic, not a fixed
    # point of a single declared map, so it is excluded from this candidate set.

    witness = _load_json(HIERARCHY_WITNESS_JSON)
    endpoint_branch = _require(witness, "public_endpoint_branch", _rel(HIERARCHY_WITNESS_JSON))
    rows.append(
        {
            "observable": "v_over_E_star",
            "value": _require(endpoint_branch, "v_over_E_star", "public_endpoint_branch"),
            "unit": "dimensionless",
            "claim_label": "conditional_public_endpoint_hierarchy_ratio",
            "chart": "closed_form_hierarchy_law_on_public_endpoint_packet",
            "open_gates": [
                "e_star_physical_normalization",
                "clock_scale_certificate",
                "frozen_running_threshold_convention",
            ],
            "source_artifact": _rel(HIERARCHY_WITNESS_JSON),
            "promotable": False,
        }
    )
    return rows


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Conditional Candidate Outputs",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        "This surface renders the numeric conditional candidates that the strict",
        "source-only promotion policy withholds from the public OPH-value columns.",
        "Each row carries its recorded claim label and its open gates. No row is",
        "promotable, and this surface must not feed back into promotion decisions.",
        "The public claim boundary stays in `RESULTS_STATUS.md`; target-anchored",
        "exact fits stay in `EXACT_FITS_ONLY.md`.",
        "",
        "The electroweak row is the target-free D10 value law. It is the OPH",
        "derivation for the pair: no measured mass enters, and the value follows as",
        "an exact implication of the five quotient-transport premises. The",
        "diagnostic base-running chart and the measured-reference inverse adapter",
        "stay in the technical audit README and off this surface.",
        "",
        "## Conditional GeV Chart Candidates",
        "",
        "| Observable | Value | Chart | Open gates | Source |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in payload["gev_chart_candidates"]:
        lines.append(
            f"| `{row['observable']}` | `{row['value']} {row['unit']}` "
            f"| `{row['chart']}` | `{', '.join(row['open_gates'])}` "
            f"| `{row['source_artifact']}` |"
        )
    lines += [
        "",
        "## Conditional Scale-Free Chart Ratios",
        "",
        "These ratios cancel the chart scale. The value-law triple shares one",
        "emitted quintet, so the three coordinates stand or fall together. The",
        "Higgs/top ratio shares the declared-surface split pair.",
        "",
        "| Observable | Value | Chart | Open gates | Source |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in payload["scale_free_chart_ratios"]:
        lines.append(
            f"| `{row['observable']}` | `{row['value']}` "
            f"| `{row['chart']}` | `{', '.join(row['open_gates'])}` "
            f"| `{row['source_artifact']}` |"
        )
    lines += [
        "",
        "## Conditional Dimensionless Candidates",
        "",
        "| Observable | Value | Chart | Open gates | Source |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in payload["dimensionless_candidates"]:
        lines.append(
            f"| `{row['observable']}` | `{row['value']}` "
            f"| `{row['chart']}` | `{', '.join(row['open_gates'])}` "
            f"| `{row['source_artifact']}` |"
        )
    lines.append("")
    return "\n".join(lines)


def build_payload() -> dict[str, Any]:
    return {
        "artifact": "oph_conditional_candidate_surface",
        "generated_utc": _timestamp(),
        "surface_policy": "conditional_candidates_rendered_promotion_withheld",
        "must_not_feed_back": True,
        "gev_chart_candidates": build_gev_chart_candidates(),
        "scale_free_chart_ratios": build_scale_free_chart_ratios(),
        "dimensionless_candidates": build_dimensionless_candidates(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--forward-out", default=str(DEFAULT_FORWARD_OUT))
    args = parser.parse_args()

    payload = build_payload()
    markdown = render_markdown(payload)
    for path, text in (
        (pathlib.Path(args.markdown_out), markdown),
        (pathlib.Path(args.json_out), json.dumps(payload, indent=2, sort_keys=True) + "\n"),
        (pathlib.Path(args.forward_out), json.dumps(payload, indent=2, sort_keys=True) + "\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
