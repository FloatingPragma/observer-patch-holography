#!/usr/bin/env python3
"""Guard the explicit neutrino oscillation comparison rows on the status surface."""

from __future__ import annotations

import importlib.util
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


def test_neutrino_oscillation_comparison_rows_are_emitted_on_repaired_branch() -> None:
    module = _load_module()
    rows = module.build_neutrino_oscillation_comparison_rows(module.build_surface_state(with_hadrons=False))
    observable_ids = {row["observable_id"] for row in rows}
    assert {
        "theta12_deg",
        "theta23_deg",
        "theta13_deg",
        "delta_cp_deg",
        "delta_m21_sq_over_delta_m32_sq",
        "delta_m21_sq_eV2",
        "delta_m32_sq_eV2",
    } <= observable_ids


def test_neutrino_oscillation_comparison_rows_preserve_compare_only_boundary() -> None:
    module = _load_module()
    rows = module.build_neutrino_oscillation_comparison_rows(module.build_surface_state(with_hadrons=False))
    by_id = {row["observable_id"]: row for row in rows}
    assert by_id["theta12_deg"]["status"] == "weighted_cycle_dimensionless"
    assert by_id["delta_m21_sq_eV2"]["status"] == "compare_only"
    assert by_id["delta_m32_sq_eV2"]["status"] == "compare_only"
    assert "two-parameter positive-segment neutrino adapter" in by_id["delta_m32_sq_eV2"]["note"]


def test_render_markdown_includes_neutrino_oscillation_section() -> None:
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
    assert "## Neutrino Oscillation Comparison" in markdown
    assert "| theta12 | weighted_cycle_dimensionless |" in markdown
