import json

import numpy as np
from scipy.linalg import expm

from . import geometry, measure, phase, preparation


def test_all_native_columns_and_unique_bulk_centre():
    route = list(range(20))
    h, ports = geometry.microscopic(route)
    assert len(h) == len(route) + 2 and set(ports) == set(route)
    sign = (-1.0) ** np.arange(len(h))
    target = 2 * (h - 36 * np.eye(len(h)))
    # Only the two winding-end onsite costs break the bulk anticommutation.
    error = target + 2 * sign[:, None] * (h - 36 * np.eye(len(h))) * sign
    assert set(map(tuple, np.argwhere(error))) == {(0, 0), (len(h) - 1, len(h) - 1)}
    np.testing.assert_array_equal(h @ sign, 72 * sign)


def test_counted_spectra_and_changed_stock():
    result = json.loads((measure.HERE / "results.json").read_text())
    stocks = result["impulse"]["stocks"]
    assert stocks["alternating"]["lattice"]["passed"]
    assert not stocks["alternating"]["fitted_hopping"]["passed"]
    assert not stocks["triple_blank_16"]["lattice"]["passed"]
    assert all(row["passed"] for row in result["packets"])
    rec = json.loads((measure.HERE / "phase.json").read_text())
    assert rec["passed"] and rec["max_hopping_error"] > 0.008
    survey = json.loads((measure.HERE / "geometry.json").read_text())
    assert len(survey["flights"]) == len(survey["belts"]) == 6
    assert all(row["passed"] for row in survey["flights"] + survey["belts"])


def test_native_pulses_compile_all_clock_and_packet_states():
    terms = measure.source.contacts()
    for theta in (0.227, -2.7, np.pi):
        np.testing.assert_allclose(
            preparation.phase_gates(theta),
            np.diag([np.exp(1j * theta), 1, 1]),
            atol=3e-13,
        )
        expected = expm(-1j * theta * preparation.contact(terms["blank"], 4, 0, 1))
        np.testing.assert_allclose(
            preparation.endpoint_pulse(theta), expected, atol=3e-13
        )
    ports, energy, vectors = measure.generator("01")
    r, b = phase.clocks(energy, vectors)
    loc = np.zeros(len(r))
    loc[[72, 74]] = np.array([1, -1]) / np.sqrt(2)
    probes = [b, loc, *[phase.marked_source(j) for j in (69, 64, 57, 49, 36, 21)]]
    targets = [
        (r + alpha * p) / np.sqrt(2) for p in probes for alpha in (1, -1, 1j, -1j)
    ]
    targets.extend(
        measure.packet(ports, q, sign)
        for q in (0.2, 0.4, 0.6, 1, 1.6, -0.2, -0.4, -0.6, -1, -1.6)
        for sign in (-1, 1)
    )
    for target in targets:
        np.testing.assert_allclose(
            preparation.compile_state(target, terms), target, atol=1e-12
        )


def test_unreferenced_density_discriminates_selected_and_sharp_sources():
    rec = json.loads((measure.HERE / "density.json").read_text())
    tested = rec["lines"][:6]
    assert [
        row["sampled"]["resolved"] and row["sampled"]["spatial_pass"] for row in tested
    ] == [True, True, True, True, True, False]
    assert all(abs(row["sampled"]["lattice_error"]) < 0.012 for row in tested[:5])
    assert not rec["single_band"]["sampled"]["resolved"]
    assert all(row["passed"] for row in rec["packets"])
    assert all(row["sampled_single_line_error"] > 0.015 for row in rec["sharp"])
    assert all(row["sampled_full_law_error"] < 0.015 for row in rec["sharp"])
    audit = json.loads((measure.HERE / "review.json").read_text())
    assert all(v < 0.008 for v in audit["max_long_errors"].values())
    assert audit["max_all_errors"]["continuum"] > 0.4
    assert audit["max_all_errors"]["native"] < 0.001


def test_all_recoil_modes_against_known_and_fitted_hopping():
    r = json.loads((measure.HERE / "transfer.json").read_text())["recoil"]
    assert len(r["rows"]) == 72
    assert r["all_held_spatial_max"] < 0.0005
    assert r["errors"]["native"]["n_held_above_gate"] == 0
    assert r["errors"]["nearest"]["n_held_above_gate"] == 31
    assert r["errors"]["next_nearest"]["n_held_above_gate"] == 0


def test_reciprocal_hopping_moments_and_charged_gateway():
    census = json.loads((measure.HERE / "harmonics.json").read_text())
    assert census["all_spatial_certified"] and all(
        v["held_pass"] for v in census["meter"].values()
    )
    for mode in census["harmonics"]:
        lo, hi = mode["sampling_interval"]
        a, b = mode["contact_forecast"]
        assert lo < mode["exact"] < hi < 0 and a < mode["exact"] < b
        assert min(hi, b) > max(lo, a)
        assert mode["foil_interval"][0] < 0 < mode["foil_interval"][1]
    read = json.loads((measure.HERE / "transfer.json").read_text())["gateway"]
    modes = {mode["mode"]: mode for mode in read["modes"]}
    assert all(r["largest_phase_error"] > 2 for r in modes["free"]["phase_rows"])
    assert all(
        r["largest_phase_error"] < 0.02 for r in modes["check_1024"]["phase_rows"]
    )
    assert modes["dephased"]["minimum_interferometric_visibility"] == 0
    charge = json.loads((measure.HERE / "conjugate.json").read_text())
    assert charge["native_generator_max_integer_error"] == 0
    assert charge["max_positive_phase_error"] < 0.02
