"""Independent DESI data, distance, propagation and train/holdout checks."""

import hashlib
import json

import camb
import numpy as np
from scipy.integrate import quad

from .measure import C, HERE, PREFIX, DistanceData, evaluate, physical_background


def test_primary_data_and_official_parameter_recovery():
    expected = {
        "mean": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585",
        "cov": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509",
    }
    for kind, digest in expected.items():
        assert (
            hashlib.sha256((HERE / f"{PREFIX}{kind}.txt").read_bytes()).hexdigest()
            == digest
        )
    d = DistanceData()
    assert sum(d.train) == 7 and sum(d.hold) == 6
    assert list(d.kind[-2:]) == ["DH_over_rs", "DM_over_rs"]
    assert np.count_nonzero(d.B) == 0 and d.cov[11, 12] < 0
    r = json.loads((HERE / "results.json").read_text())
    for stock in r["stocks"].values():
        released = stock["runs"]["0.02237"]["all_H0_and_free_ruler"]
        np.testing.assert_allclose(released["Omega_m"], 0.2975, atol=0.0007)
        np.testing.assert_allclose(
            released["H0"] / 100 * released["effective_ruler_Mpc"], 101.54, atol=0.025
        )


def test_physical_ruler_distance_integral_and_each_desi_observable():
    d = DistanceData()
    b, q, h0 = 0.02237, 1037 / 202, 70.0
    base, pred = physical_background(d, b, q, h0)
    cp = camb.CAMBparams()
    cp.set_cosmology(H0=h0, ombh2=b, omch2=b * q, mnu=0.06)
    bg = camb.get_background(cp)
    rd = bg.sound_horizon(bg.get_derived_params()["zdrag"])
    np.testing.assert_allclose(rd, base["rd_Mpc"], rtol=3e-6)
    independently = []
    for z, kind in zip(d.z, d.kind, strict=True):
        dm = quad(lambda zp: C / bg.hubble_parameter(zp), 0, z, epsabs=1e-7)[0]
        dh = C / bg.hubble_parameter(z)
        v = {
            "DM_over_rs": dm,
            "DH_over_rs": dh,
            "DV_over_rs": (z * dm * dm * dh) ** (1 / 3),
        }[kind]
        independently.append(v / rd)
    np.testing.assert_allclose(pred, independently, rtol=5e-6)
    res = d.obs - pred
    np.testing.assert_allclose(
        d.score(pred)["full"], res @ np.linalg.solve(d.cov, res), rtol=1e-12
    )


def test_held_points_do_not_set_training_H0_and_replay_score():
    d = DistanceData()
    stock = {"ratio": 1037 / 202}
    previous = json.loads((HERE / "results.json").read_text())["stocks"]["80"]["runs"][
        "0.02237"
    ]["low_z_H0"]
    # The released high-z vector gets a large perturbation for this control.
    d.obs[d.hold] += np.linspace(-1.0, 1.0, d.hold.sum())
    trial = evaluate(d, 0.02237, stock)["low_z_H0"]
    np.testing.assert_allclose(trial["H0"], previous["H0"], atol=1e-5)
    np.testing.assert_allclose(
        trial["H0_delta_chi2_1_interval"],
        previous["H0_delta_chi2_1_interval"],
        atol=1e-5,
    )
    assert (
        abs(
            trial["scores"]["held_given_train"] - previous["scores"]["held_given_train"]
        )
        > 1.0
    )
    untouched = DistanceData()
    score = untouched.score(trial["predicted"])
    np.testing.assert_allclose(score["held_given_train"], 1.00775, rtol=1e-4)
    np.testing.assert_allclose(
        score["full"], score["train"] + score["held_given_train"], rtol=1e-10
    )
