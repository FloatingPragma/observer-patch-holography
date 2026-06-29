#!/usr/bin/env python3
"""Tests for the final public-frontier quark theorem package."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_public_physical_sigma_datum_descent import build_artifact as build_public_sigma_descent
from derive_quark_public_strengthened_physical_sigma_lift_frontier import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_public_strengthened_physical_sigma_lift_frontier_records_final_routes() -> None:
    public_sigma_theorem = build_public_sigma_descent(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )
    payload = build_artifact(
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("quark_absolute_readout_algebraic_collapse.json"),
        _load("quark_exact_pdg_end_to_end_theorem.json"),
        _load("quark_exact_yukawa_end_to_end_theorem.json"),
        _load("overlap_edge_line_lift.json"),
        _load("generation_bundle_branch_generator.json"),
        public_sigma_theorem,
    )

    assert payload["proof_status"] == "blocked_target_derived_sigma_datum_descent"
    assert payload["public_promotion_allowed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "quark_sigma_source_datum_no_target_leak_required"
    )
    assert payload["resolved_by_theorem_artifact"] == "oph_quark_public_physical_sigma_datum_descent"
    final = payload["final_public_theorem_candidate"]
    assert final["id"] == "selected_bridge_fiber_sigma_descent_not_source_selection"
    assert final["induces_global_contract"]["id"] == "strengthened_quark_physical_sigma_ud_lift"
    sigma = final["must_emit"]["physical_sigma_datum"]
    assert abs(float(sigma["sigma_u"]) - 5.573928426395543) < 1.0e-12
    assert abs(float(sigma["sigma_d"]) - 3.296264198808688) < 1.0e-12
    exact_yukawas = payload["algebraic_consequence_after_closure"]["forced_exact_yukawas"]
    assert exact_yukawas["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert exact_yukawas["forward_certified"] is True
    assert payload["alternate_upstream_route"]["id"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["alternate_upstream_route"]["status"] == "upstream_alternative_route_currently_deprioritized"
    assert payload["alternate_upstream_route"]["smaller_exact_missing_clause"] == (
        "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split"
    )
    obstruction = payload["alternate_upstream_route"]["current_attached_data_obstruction"]
    assert abs(float(obstruction["commutator_operator_norm"]) - 0.04861550547372144) < 1.0e-12
    assert abs(float(obstruction["projector_defect_operator_norm"]) - 0.06363734112184061) < 1.0e-12
