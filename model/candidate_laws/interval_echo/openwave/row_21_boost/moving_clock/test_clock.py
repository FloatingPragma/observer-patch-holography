import json
from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from . import apparatus as lab, current, preparation, readback


def test_entire_local_sequence_against_physical_one_particle_matrix():
    n = 5
    rng = np.random.default_rng(447)
    psi = rng.normal(size=2 * n) + 1j * rng.normal(size=2 * n)
    local_h = law.seam_term(8)
    ids = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    direct = psi.copy()
    for phase, t in enumerate(lab.HOLDS):
        full = np.zeros((2 * n, 2 * n))
        for i in range(n):
            edge = [(2 * i + phase) % (2 * n), (2 * i + phase + 1) % (2 * n)]
            full[np.ix_(edge, edge)] += local_h[np.ix_(ids, ids)]
        direct = expm(-1j * t * full) @ direct
    np.testing.assert_allclose(
        lab.step(psi.reshape(n, 2), lab.gates(1, 8)).ravel(), direct, atol=2e-13
    )
    sample = psi / np.linalg.norm(psi)
    for s in (1, 0.8):
        assert preparation.build(sample, s)[1] < 2e-13


def test_nonrounded_phase_and_current_with_opposite_arms():
    for s, h in product((1.0, 0.8), (-5, 0, 5)):
        n = 57
        layers = lab.gates(s, 8.0)
        exact = lab.exact(n, h, s)
        z = lab.arms(n, h, layers, s, 8.0, plane=True)
        out = lab.step(z, layers)
        rail = int(h < 0)
        w = lab.step(z, layers)
        for i, arm in enumerate(z):
            middle = arm @ layers[0].T
            observed = (abs(w[i, :, 0]) ** 2).sum() - (abs(middle[:, 0]) ** 2).sum()
            np.testing.assert_allclose(observed, exact["group_velocity"], atol=1e-12)
            bins = np.array(current.route_meter(arm, rail, s))
            im, re = bins[:, 1] - bins[:, 0]
            np.testing.assert_allclose(
                np.angle(re + 1j * im), (-1) ** i * exact["k"], atol=9e-13
            )
        q = (
            np.conj(out[0, 0, rail])
            * out[1, 0, rail]
            / (np.conj(z[0, 0, rail]) * z[1, 0, rail])
        )
        np.testing.assert_allclose(q, np.exp(2j * exact["band_energy"]), atol=2e-13)


def test_complete_radar_return_and_all_front_refusals():
    n = 61
    for d in (-11, 11):
        x = int(d < 0)
        layers = lab.gates(1.0, 8.0)
        u = lab.go(lab.sharp(n, 0, x), abs(d), layers)
        np.testing.assert_allclose(
            abs(u[d, x]) ** 2, np.cos(0.2) ** (2 * abs(d)), atol=1e-13
        )
        assert np.max(abs(lab.go(lab.sharp(n, 0, x), abs(d) - 1, layers)[d])) == 0
        rng = np.random.default_rng(47)
        report = lab.radar(n, 0, d, layers, 2**16, rng)
        assert sum(report["forward_then_return"]) == sum(report["parked_turn"]) == 2**16
        assert report["parked_turn"][-1] == report["early_clicks"] == 0
        assert report["forward_then_return"][-1] > 10_000


def test_deposited_counts_decisions():
    data = json.loads((readback.HERE / "records.json").read_text())
    with np.load(readback.HERE / "records.npz") as b:
        out = readback.report(data, b)
    assert len(out) == 28
    assert all(
        x["packet"]["resolved"]
        and x["plane"]["resolved"]
        and not x["mixture"]["resolved"]
        for x in out
    )
    assert all(not x["packet"]["passes"] for x in out if x["h"])
    proof = current.read(json.loads((readback.HERE / "current.json").read_text()))
    for row in proof:
        if row["key"] in ("2_native_21", "2_native_-21"):
            assert row["ideal_interval"][1] < row["lorentz_interval"][0]
        if row["key"] in ("2_native_5", "2_native_-5"):
            assert (
                row["ideal_interval"][0] < row["lorentz_interval"][1]
                and row["lorentz_interval"][0] < row["ideal_interval"][1]
            )
    controls = json.loads((readback.HERE / "preparation.json").read_text())
    assert sum(r["states"] for r in controls) == 336
    assert max(r["maximum_column_error"] for r in controls) < 1e-11
