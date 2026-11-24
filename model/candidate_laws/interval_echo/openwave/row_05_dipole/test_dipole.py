"""Independent checks of stationary subspaces, native energies and recorded pointer phases."""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell, federation

HERE = Path(__file__).parent


def read(name):
    return json.loads((HERE / f"{name}.json").read_text())


def test_integer_nullity_and_no_exterior_transition():
    board = cell()
    h = sp.eye(12)
    for a, b in board.seams:
        h[a, b] = h[b, a] = 1
    data = read("results")
    for r in data["census"]:
        form = h.col_join(sp.eye(12)[r["gates"], :])
        assert 12 - form.rank() == r["dimension"]
        u = sp.Matrix.hstack(*form.nullspace())
        assert h * u == sp.zeros(12, u.cols)
        assert u[r["gates"], :] == sp.zeros(3, u.cols)
    for n in (20, 80, 320):
        sources = [
            a for r in data["census"] if r["owners"] == n for a in r["addresses"]
        ]
        assert sorted(sources) == list(range(n))


def test_directional_source_current_and_tilt_counts():
    d = read("results")
    for panel in d["current"].values():
        counts = np.array(panel["counts"])
        assert np.all(counts.sum(axis=1) == d["shots"])
        observed = (
            -(counts[:, 0] - counts[:, 1])
            @ np.array(d["axial_areas"])
            / (2 * d["shots"])
        )
        np.testing.assert_allclose(observed, panel["moment"], atol=1e-10)
    assert abs(d["current"]["real"]["moment"]) < 0.001
    assert min(d["current"]["plus"]["moment"], -d["current"]["minus"]["moment"]) > 0.137
    z = complex(*d["rotation"]["amplitude"])
    assert abs(z**3 - 1) < 2e-10
    v = read("vector")
    assert sum(r["singular_values"][2] > 0.01 for r in v["layouts"]) > 0
    assert max(
        max(abs(np.array(r["measured"]) - r["projected_forecast"])) for r in v["rows"]
    ) < min(v["simultaneous_radius"])
    t = {r["field"]: r for r in v["rows"] if r["time"] == 128}
    assert t[1 / 32]["measured"][1] > 0.068 and t[-1 / 32]["measured"][1] < -0.068
    assert abs(t[0]["measured"][1]) < 0.001


def test_native_flux_calibration_and_clock_counts():
    data = read("calibration")
    counts = np.array(data["clock_counts_z_y"])
    assert np.all(counts.sum(axis=2) == data["shots"])
    z, y = ((counts[:, :, 0] - counts[:, :, 1]) / data["shots"]).T
    theta = np.unwrap(np.angle(z + 1j * y))
    np.testing.assert_allclose(theta - theta[0], data["clock_read"], atol=1e-10)
    assert max(abs(theta - theta[0] - np.arange(len(theta)))) < 0.001
    for field, trial in data["fields"].items():
        x, y = (np.array(trial["counts"]) @ [1, -1]) / data["shots"]
        b = np.angle(x + 1j * y) / (2 * data["charge"] * data["oriented_area"])
        assert abs(b - trial["measured"]) < 1e-10
        assert abs(b - float(field)) < 0.0005
    assert max(abs(r["forecast"] - r["phase"]) for r in data["held"]) < 0.004
    assert sum(data["vacancy_address_counts"]) == data["shots"]


def test_exact_contact_residual_from_local_term_and_flat_recoil():
    h = law.local_term_parts()
    blank = (h[0] + 8 * h[1])[: law.N_SRC, : law.N_SRC]

    def cost(a, b):
        i = law.SRC.index((a, b))
        return blank[i, i]

    delta = cost(0, 2) - cost(1, 2) - cost(0, 1) + cost(1, 1)
    d = read("exterior")
    assert abs(delta - 15.5) < 1e-14
    for r in d["exterior"]:
        assert abs(r["control_connected_energy"] - delta) < 1e-14
        assert r["max_residual"] < 1e-9
        assert all(
            sum(c) == d["shots"]
            for s in r["receiver"]
            for c in s["directions_and_absent"]
        )
        if r["owners"] == 20:
            assert r["native_columns_checked"] == federation(0).n - 12
    for p in d["kick"]:
        for r in p["rows"]:
            a, b, sharp = r["sources"]
            assert a["return_counts"] == b["return_counts"] == [d["shots"], 0]
            if r["time"] == 16:
                assert sharp["return_counts"][1] > d["shots"] / 4


def test_slow_precession_frequency_from_opposite_field_counts():
    data = read("vector")
    odd = {}
    for t in (128, 256):
        same = {r["field"]: r["measured"][1] for r in data["rows"] if r["time"] == t}
        odd[t] = (same[1 / 32] - same[-1 / 32]) / 2
    radius = data["simultaneous_radius"][1]
    low = np.arccos((odd[256] + radius) / (2 * (odd[128] - radius))) / 128
    high = np.arccos((odd[256] - radius) / (2 * (odd[128] + radius))) / 128
    assert 0.01093 < low < 1 / (32 * np.sqrt(8)) < high < 0.01117
