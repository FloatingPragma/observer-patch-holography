#!/usr/bin/env python3
"""Guard the Majorana-phase rows on the particle status surface."""

from __future__ import annotations

import importlib.util
import json
import pathlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_results_status_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_results_status_table", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_majorana_phase_surface_rows_emit_on_live_theorem() -> None:
    module = _load_module()
    rows = module.build_majorana_phase_surface_rows(module.build_surface_state(with_hadrons=False))
    by_id = {row["observable_id"]: row for row in rows}
    assert {"alpha21_majorana", "alpha31_majorana"} <= set(by_id)
    assert by_id["alpha21_majorana"]["status"] == "theorem_grade"
    assert by_id["alpha21_majorana"]["prediction_display"] == "153.618518 deg"
    assert by_id["alpha31_majorana"]["prediction_display"] == "257.003241 deg"


def test_majorana_phase_surface_rows_stay_absent_without_public_promotion(tmp_path: pathlib.Path) -> None:
    module = _load_module()
    theorem_path = tmp_path / "majorana_theorem.json"
    theorem_path.write_text(
        json.dumps(
            {
                "artifact": "oph_neutrino_physical_majorana_phase_theorem",
                "status": "candidate_only",
                "public_surface_candidate_allowed": False,
                "public_promotion_blocker": "synthetic blocker for candidate-only branch coverage",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    original = module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM
    module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM = theorem_path
    try:
        rows = module.build_majorana_phase_surface_rows(module.build_surface_state(with_hadrons=False))
    finally:
        module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM = original

    by_id = {row["observable_id"]: row for row in rows}
    assert by_id["alpha21_majorana"]["status"] == "still_absent"
    assert by_id["alpha21_majorana"]["prediction_display"] == "n/a"
    assert by_id["alpha31_majorana"]["prediction_display"] == "n/a"
    assert "synthetic blocker" in by_id["alpha21_majorana"]["note"]


def test_majorana_phase_surface_rows_require_theorem_grade_status(tmp_path: pathlib.Path) -> None:
    module = _load_module()
    theorem_path = tmp_path / "majorana_theorem.json"
    theorem_path.write_text(
        json.dumps(
            {
                "artifact": "oph_neutrino_physical_majorana_phase_theorem",
                "status": "candidate_only",
                "public_surface_candidate_allowed": True,
                "public_promotion_blocker": "synthetic inconsistent status gate",
                "candidate_parameters": {
                    "alpha21_deg_0_to_360": 153.6185177794357,
                    "alpha31_deg_0_to_360": 257.00324082207993,
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    original = module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM
    module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM = theorem_path
    try:
        rows = module.build_majorana_phase_surface_rows(module.build_surface_state(with_hadrons=False))
    finally:
        module.NEUTRINO_PHYSICAL_MAJORANA_PHASE_THEOREM = original

    by_id = {row["observable_id"]: row for row in rows}
    assert by_id["alpha21_majorana"]["status"] == "still_absent"
    assert by_id["alpha21_majorana"]["prediction_display"] == "n/a"
    assert "synthetic inconsistent status gate" in by_id["alpha21_majorana"]["note"]


def test_render_markdown_includes_majorana_phase_section() -> None:
    module = _load_module()
    surface_state = module.build_surface_state(with_hadrons=False)
    premise_boundaries = module.build_premise_boundaries()
    reference_payload = {
        "source": {
            "label": "Particle Data Group",
            "edition": "2025",
            "api_info_url": "https://pdg.lbl.gov/api",
        }
    }
    markdown = module.render_markdown(
        rows=[
            {
                "particle": "nu_e",
                "group": "Leptons",
                "status": "continuation",
                "prediction_display_gev": "n/a",
                "reference_display": "not directly measured",
                "delta_display": "n/a",
                "note": "stub",
            }
        ],
        comparison_rows=module.build_neutrino_oscillation_comparison_rows(surface_state),
        majorana_rows=module.build_majorana_phase_surface_rows(surface_state),
        generated_utc="2026-03-30T00:00:00Z",
        P=1.63094,
        log_dim_H=1.0e122,
        loops=4,
        with_hadrons=False,
        hadron_profile="serious",
        reference_payload=reference_payload,
        surface_state=surface_state,
        premise_boundaries=premise_boundaries,
    )
    assert "## Majorana Phase Surface" in markdown
    assert "| alpha21^(Maj) | theorem_grade | 153.618518 deg |" in markdown
