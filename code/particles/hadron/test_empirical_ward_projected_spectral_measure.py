#!/usr/bin/env python3
"""Tests for the empirical Ward-projected spectral-measure export."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "derive_empirical_ward_projected_spectral_measure.py"
SCHEMA = HERE / "empirical_ward_projected_spectral_measure.schema.json"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "empirical_ward_projected_spectral_measure.json"
UPSTREAM = ROOT / "particles" / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"

M_Z2 = 91.1876 * 91.1876
ALPHA0 = 1.0 / 137.035999177
KNT19_VALUE = 0.02761
KNT19_UNCERTAINTY = 0.00011


def _export() -> dict:
    if not OUTPUT.exists():
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_export_regenerates_and_satisfies_schema() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = _export()
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    try:
        import jsonschema
    except ImportError:
        for key in schema["required"]:
            assert key in payload, key
    else:
        jsonschema.validate(instance=payload, schema=schema)
    assert payload["artifact"] == "oph_empirical_ward_projected_hadronic_spectral_measure"
    assert payload["row_class"] == "oph_plus_empirical_hadron_closure"


def test_guards_keep_empirical_class_non_promotable() -> None:
    guards = _export()["guards"]
    assert guards["source_only"] is False
    assert guards["empirical_hadron_closure"] is True
    assert guards["external_cross_section_data_used"] is True
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["usable_for_public_final_values"] is True
    assert guards["surrogate_hadron_artifact"] is False
    assert guards["satisfies_production_constructive_next_artifact"] is False


def test_measure_is_positive() -> None:
    payload = _export()
    for segment in payload["continuum_density"]:
        assert all(rv >= 0.0 for rv in segment["R_values"]), segment["segment_id"]
        assert len(segment["grid_sqrt_s"]) == len(segment["R_values"])
    for atom in payload["spectral_atoms"]:
        assert atom["weight"] > 0.0, atom["atom_id"]
    assert payload["rho_had_or_measure"]["positivity_status"] == (
        "verified_nonnegative_on_exported_grids_and_atoms")


def test_requadrature_consistency_gate_holds() -> None:
    payload = _export()
    upstream = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    consistency = payload["consistency"]
    assert consistency["within_tolerance"] is True
    assert consistency["abs_difference"] <= consistency["tolerance"]
    assert abs(consistency["upstream_integral_value"]
               - upstream["integral"]["value"]) < 1e-12
    # the gate must sit far below the physics budget
    assert consistency["tolerance"] < 0.1 * upstream["integral"]["uncertainty"]


def test_transport_moments_reintegrate_from_exported_objects() -> None:
    """Third, test-side reconstruction of the full timelike moment from the
    exported objects only: midpoint rule on the grids (the derive script uses
    trapezoid), analytic atoms, and the perturbative tail rebuilt from the
    exported parametric form at a different point count. No import of the
    ingest module."""
    payload = _export()

    def kernel(s: float) -> float:
        return M_Z2 / (s * (M_Z2 - s))

    total = 0.0
    for segment in payload["continuum_density"]:
        grid = segment["grid_sqrt_s"]
        r_vals = segment["R_values"]
        for i in range(len(grid) - 1):
            rs = 0.5 * (grid[i] + grid[i + 1])
            r_mid = 0.5 * (r_vals[i] + r_vals[i + 1])
            total += r_mid * kernel(rs * rs) * 2.0 * rs * (grid[i + 1] - grid[i])
    for atom in payload["spectral_atoms"]:
        total += atom["weight"] * kernel(atom["s"])

    tail_spec = payload["pqcd_tail"]
    assert tail_spec["alpha_s_loop_order"] == 1
    lam2 = tail_spec["lambda_qcd5_gev"] ** 2

    def r_pqcd(s: float) -> float:
        a = 4.0 * math.pi / (23.0 / 3.0 * math.log(s / lam2)) / math.pi
        return (11.0 / 3.0) * (1.0 + a + 1.409 * a * a - 12.767 * a ** 3)

    s_lo = tail_spec["sqrt_s_start"] ** 2
    s_hi = tail_spec["numeric_end_sqrt_s"] ** 2
    r_mz = r_pqcd(M_Z2)
    n = 30001
    for i in range(n):
        s = s_lo + (s_hi - s_lo) * (i + 0.5) / n
        total += (r_pqcd(s) - r_mz) * kernel(s) * (s_hi - s_lo) / n
    total += r_mz * (math.log(s_hi / abs(M_Z2 - s_hi)) - math.log(s_lo / abs(M_Z2 - s_lo)))
    a_end = 12.0 / (23.0 * math.log(s_hi / lam2))
    r_end = (11.0 / 3.0) * (1.0 + a_end)
    total += -r_end * math.log(1.0 / (1.0 - M_Z2 / s_hi))

    reconstructed = ALPHA0 / (3.0 * math.pi) * total
    reported = payload["transport_moments"]["timelike_on_shell_mz"]["value"]
    # cross-rule agreement far below the physics uncertainty (~7.5e-4)
    assert abs(reconstructed - reported) < 5.0e-06, (reconstructed, reported)


def test_external_cross_check_within_budget() -> None:
    """The compilation-class moment agrees with the data-driven KNT19 row
    within two combined standard deviations, compare-only."""
    moments = _export()["transport_moments"]["timelike_on_shell_mz"]
    combined = math.hypot(moments["uncertainty"], KNT19_UNCERTAINTY)
    assert abs(moments["value"] - KNT19_VALUE) <= 2.0 * combined


def test_both_kernels_present_and_packet_alpha_free() -> None:
    moments = _export()["transport_moments"]
    timelike = moments["timelike_on_shell_mz"]["value"]
    spacelike = moments["spacelike_mz"]["value"]
    packet = moments["inverse_alpha_packet_spacelike"]["value"]
    assert 0.02 < timelike < 0.035
    assert 0.02 < spacelike < 0.035
    # packet = spacelike moment with the alpha prefactor removed
    assert abs(packet - spacelike / ALPHA0) < 1e-9
    assert 2.5 < packet < 4.5
