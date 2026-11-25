"""Independent reconstruction from the saved readouts, geometry and original counts."""

import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.special import eval_legendre

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "records.json").read_text())
CLOCK = DATA["oscillator"]


def test_independent_symbol_and_clock_readouts():
    vertices = np.array(CLOCK["native_vertices"])
    distances = np.linalg.norm(vertices[:, None] - vertices, axis=2)
    length = distances[distances > 0.01].min()
    edges = [
        (i, j)
        for i, j in combinations(range(12), 2)
        if abs(distances[i, j] - length) < 1e-7
    ]
    assert edges == list(map(tuple, CLOCK["edges"])) and len(edges) == 30
    assert len(set(CLOCK["rail_hubs"])) == 2 and CLOCK["polarization_route"] > 0
    assert all(len(v) == 30 and min(v) > 0 for v in CLOCK["route_lengths"])
    assert np.all(np.abs(np.diff(CLOCK["marker"])) == 1)
    for mode in CLOCK["modes"]:
        k = np.asarray(mode["wavevector"])
        phases = np.array([(vertices[b] - vertices[a]) @ k for a, b in edges])
        sign = 1 if mode["mode"] == "wrong_sign" else -1
        # Positive expression avoids loss near the zero-gradient mode.
        b = (1 + sign * np.exp(1j * phases)) / np.sqrt(60)
        omega = np.linalg.norm(b)
        assert abs(omega**2 - mode["stiffness"]) < 3e-9
        eps = CLOCK["signed_halfwidth"]
        first = next(r for r in mode["reads"] if r["t"] == 1 and r["input"] == "vertex")
        counts = np.array(first["counts"])[0, 0]
        pilot = (counts[0] - counts[1]) / CLOCK["shots"]
        err = mode["refinements"][-1][2] + eps + 1e-8
        low, high = np.arccos(np.clip([pilot + err, pilot - err], -1, 1))
        assert low <= omega + 1e-7 and omega <= high + 1e-7
        for r in mode["reads"]:
            t = r["t"]
            coeff = (1 if r["input"] == "edge_plus" else -1) * b[
                mode["edge"]
            ].conjugate()
            a = (
                np.cos(t * omega)
                if r["input"] == "vertex"
                else t * np.sinc(t * omega / np.pi) * coeff
            )
            observed = np.array(r["counts"])
            assert np.all(observed.sum(axis=2) == CLOCK["shots"])
            z = (observed[:, :, 0] - observed[:, :, 1]) / CLOCK["shots"]
            assert (
                np.max(abs(z - [a.real, a.imag]))
                < eps + t * mode["refinements"][-1][2] + 2e-8
            )
        for idx, row in enumerate(mode["coherent_polarization"]):
            assert abs(sum(row) - 1) < 2e-8 and row[1 - idx % 2] < 2e-8


def test_refinement_and_sign_control():
    banks = {v["mode"]: v for v in CLOCK["modes"]}
    for mode in ("port", "oblique", "wrong_sign"):
        refs = np.array(banks[mode]["refinements"])
        assert np.array_equal(refs[:, 0], [8, 32, 128])
        assert np.all(refs[:-1, 1] > 15 * refs[1:, 1])
        assert np.all(refs[:, 1] <= refs[:, 2] + 2e-9)
    assert banks["zero"]["stiffness"] == 0 and banks["wrong_sign"]["stiffness"] == 2
    assert banks["zero"]["reads"][-3]["response"][0] > 0.99999
    assert banks["wrong_sign"]["reads"][-3]["response"][0] < -0.4


def test_last_encounters_and_independent_blind_prediction():
    v = DATA["visibility"]
    prev = json.loads((HERE.parent / "angular_transfer" / "records.json").read_text())
    n = np.array(prev["directions"])
    targets = {0, *(end for _, end in prev["endpoints"]["arrivals"])}
    assert {(src, t) for src, t, _, _ in v["messages"]} == {
        (src, t) for src in v["owners"] for t in targets
    }
    assert all(
        length > 0 and outcomes == [[0, 0], [0, 1]]
        for _, _, length, outcomes in v["messages"]
    )
    score = eval_legendre(6, n[:12] @ n.T).mean(axis=0)
    c = np.r_[np.ones(12) / 12, -np.ones(20) / 20]
    raw = {}
    with np.load(HERE / "directions.npz", allow_pickle=False) as f:
        for key, weights in [
            ("layered", [3, 9, 1, 3]),
            ("fresh_fair", [4, 4, 4, 4]),
            ("reused_fair", [8, 0, 0, 8]),
        ]:
            counts = f[key].astype(float)
            assert counts.shape == (32, 4, 34)
            assert np.all(counts.sum(axis=(1, 2)) == v["angular_starts"])
            r = v["banks"][key]
            w = np.array(weights) / 16
            last = w * np.array([1, 2 / 3, 2 / 3, 4 / 9])
            assert (
                np.max(
                    abs(
                        np.asarray(r["surviving_last"])
                        - [last[0], last[2], last[1] + last[3]]
                    )
                )
                < 3e-9
            )
            assert abs(r["survival"] - sum(last)) < 1e-8
            two = w @ np.array([1, 1 / 15, 1 / 15, 7 / 225])
            assert abs(r["raw_two"] - two) < 1e-8
            raw[key] = (
                c
                @ (counts[:, :, :32].sum(axis=1) / v["angular_starts"])
                @ score
                / (c @ score)
            )
            assert abs(raw[key] - r["counted_six"]) < 2e-9
            assert abs(raw[key] - w[0]) < v["gain_halfwidth"]
            assert abs(r["accepted_six"] - w[0] / sum(last)) < 1e-8
    p = np.asarray(v["calibration"]) / v["calibration_starts"]
    assert (
        abs((1 - p[0]) * (1 - p[1]) - raw["layered"])
        < v["gain_halfwidth"] + 2 * v["calibration_halfwidth"]
    )
    assert raw["reused_fair"] - raw["fresh_fair"] > 0.24
