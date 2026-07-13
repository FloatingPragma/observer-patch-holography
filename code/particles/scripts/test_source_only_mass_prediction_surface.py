"""Tests for the per-family source-only mass prediction surface."""

from __future__ import annotations

import build_source_only_mass_prediction_surface as lane


def test_surface_covers_every_family():
    surface = lane.build()
    names = [family["family"] for family in surface["families"]]
    assert names == [
        "charged leptons",
        "electroweak bosons",
        "Higgs boson and top quark",
        "quarks (u, d, s, c, b)",
        "neutrinos",
        "massless carriers (photon, gluon)",
        "hadrons",
        "operational scale E_star",
    ]
    assert surface["promotion_allowed"] is False


def test_lepton_block_carries_both_lanes_and_open_routes():
    surface = lane.build()
    leptons = surface["families"][0]
    lanes = {row["lane"] for row in leptons["rows"]}
    assert any("MCPR" in name for name in lanes)
    assert any("kappa interval" in name for name in lanes)
    mcpr_row = next(row for row in leptons["rows"] if "MCPR" in row["lane"])
    assert mcpr_row["tier"] == "T2"
    assert "W5_ORB" in mcpr_row["blocking_objects"]
    assert len(mcpr_row["masses_MeV_display"]) == 3


def test_boson_block_carries_the_discrete_two_law_pair():
    surface = lane.build()
    bosons = surface["families"][1]
    assert len(bosons["rows"]) == 3
    width = bosons["discrete_ambiguity_width_GeV"]
    assert 0.001 < width["MW"] < 0.02
    zero = next(r for r in bosons["rows"] if "zero-selector" in r["lane"])
    assert abs(zero["MW_GeV_display"] - 80.3301) < 0.001
    assert abs(zero["MZ_GeV_display"] - 91.1191) < 0.001


def test_markdown_renders_with_tier_ladder_and_all_families():
    surface = lane.build()
    text = lane.render_markdown(surface)
    assert "# Source-Only Mass Prediction Surface" in text
    assert "| T4 |" in text
    for family in surface["families"]:
        assert f"## {family['family']}" in text
    assert "Global blocking objects" in text


def test_every_input_artifact_exists():
    surface = lane.build()
    assert all(entry["exists"] for entry in surface["inputs"].values())
