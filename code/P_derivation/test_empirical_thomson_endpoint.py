#!/usr/bin/env python3
"""Tests for the empirical Thomson endpoint and its payload."""

from __future__ import annotations

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from empirical_thomson_endpoint import evaluate  # noqa: E402

PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"
SCHEMA = ROOT / "particles" / "hadron" / "empirical_ee_hadronic_spectral_measure.schema.json"


def _payload() -> dict:
    return json.loads(PAYLOAD.read_text(encoding="utf-8"))


def test_payload_satisfies_schema_required_fields() -> None:
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    payload = _payload()
    for key in schema["required"]:
        assert key in payload, key
    assert payload["artifact"] == "oph_empirical_ee_hadronic_spectral_measure"
    assert payload["row_class"] == "oph_plus_empirical_hadron_closure"
    guards = payload["guards"]
    assert guards["source_only"] is False
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["external_cross_section_data_used"] is True
    grid = payload["energy_grid"]["values"]
    r_vals = payload["R_values"]["values"]
    assert len(grid) == len(r_vals) and len(grid) > 10


def test_payload_integral_is_sane_and_cross_checked() -> None:
    integ = _payload()["integral"]
    value, unc = integ["value"], integ["uncertainty"]
    assert 0.024 < value < 0.030
    assert 0.0 < unc < 0.15 * value
    xc = integ["external_cross_check"]
    assert xc["role"] == "compare_only_not_an_input"
    # payload agrees with the external compilation within two payload sigmas
    assert abs(value - xc["value"]) < 2.0 * unc


def test_endpoint_guards_and_solve_path() -> None:
    report = evaluate()
    assert report["row_class"] == "oph_plus_empirical_hadron_closure"
    guards = report["guards"]
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["measured_alpha_in_solve_path"] is False
    lo, hi = (float(x) for x in report["endpoint"]["alpha_inv_interval"])
    central = float(report["endpoint"]["alpha_inv_central"])
    assert lo < central < hi
    p_lo, p_hi = (float(x) for x in report["endpoint"]["P_interval"])
    assert 1.60 < p_lo < p_hi < 1.66


def test_endpoint_reports_anchor_gap_consistently() -> None:
    report = evaluate()
    cmp_blk = report["compare_only"]
    gap_lo, gap_hi = (float(x) for x in cmp_blk["gap_interval_inv_alpha"])
    a_lo, a_hi = (float(x) for x in
                  cmp_blk["same_scheme_anchor_gap_interval_inv_alpha"])
    assert gap_lo < gap_hi and a_lo < a_hi
    # the anchor gap and the endpoint gap describe the same discrepancy with
    # opposite signs, up to the second-order insertion correction
    assert abs((-0.5 * (gap_lo + gap_hi)) - 0.5 * (a_lo + a_hi)) < 0.05
    inside = cmp_blk["codata_inside_endpoint_interval"]
    codata = float(cmp_blk["codata_alpha_inv"])
    lo, hi = (float(x) for x in report["endpoint"]["alpha_inv_interval"])
    assert inside == (lo <= codata <= hi)


def test_endpoint_requires_spectral_measure_export() -> None:
    report = evaluate()
    export = report["inputs"]["hadronic_spectral_measure_export"]
    assert export["artifact"] == "oph_empirical_ward_projected_hadronic_spectral_measure"
    assert float(export["requadrature_abs_difference"]) <= float(
        export["requadrature_tolerance"])
    assert export["positivity_status"] == (
        "verified_nonnegative_on_exported_grids_and_atoms")


def test_transport_split_carries_hadronic_spectral_object() -> None:
    report = evaluate()
    split = report["transport_split"]
    assert split["a0_anchor_inv_alpha"] == report["inputs"]["source_anchor_inv_alpha_MZ"]
    hadronic = split["hadronic_spectral_object"]
    assert hadronic["status"] == "declared_empirical_not_source_emitted"
    # the spectral object's timelike moment is the solve-path hadronic input
    assert abs(float(hadronic["delta_alpha_had_5_MZ_timelike"])
               - float(report["inputs"]["delta_alpha_had_5_MZ"])) < 2.0e-05
    packet = float(hadronic["inverse_alpha_packet_spacelike"])
    assert 2.5 < packet < 4.5


def test_anchor_bridge_block_references_scheme_analysis() -> None:
    report = evaluate()
    bridge = report["compare_only"]["anchor_bridge"]
    assert bridge["route"] == "route_A_empirical_class_supplied_by_this_lane"
    assert "#425" in bridge["source_only_branch"]
