"""Independent checks on the Planck covariance, band windows and train/holdout split."""

import json
from fractions import Fraction

import numpy as np
from scipy.linalg import cholesky

from .boltzmann import RATIOS, compute
from .measure import HERE, Likelihood, norm


def test_official_check_and_covariance_factorization():
    lk = Likelihood()
    assert abs(lk.check - lk.reference_check) < 5e-8
    assert lk.n == 155 and lk.low[0] == 100 and lk.high[-1] == 1494
    rng = np.random.default_rng(54321)
    y = lk.y + cholesky(lk.cov, lower=True) @ rng.normal(size=3 * lk.n)
    scores = lk.check_polarization(y)
    np.testing.assert_allclose(
        norm(cholesky(lk.cov, lower=True), lk.y - y),
        scores["TT"] + scores["conditional_P"],
        rtol=1e-12,
    )


def test_band_windows_by_independent_loop_and_holdout_blind_calibration():
    lk = Likelihood()
    ells = np.arange(2801)
    D = 1200 + 320 * np.sin(ells / 60)
    bands = lk.bands(D, 1.0, 0.0)
    with np.load(HERE / "plik_lite.npz") as raw:
        for j, (lo, hi) in enumerate(zip(lk.low, lk.high, strict=True)):
            L = np.arange(lo, hi + 1)
            direct = sum(D[L] * raw["weight"][L - 30] * 2 * np.pi / (L * (L + 1)))
            np.testing.assert_allclose(bands[j], direct, rtol=1e-12)
    before = lk.calibrate(D, True, 0.2)
    lk.y[lk.n :] += 10
    np.testing.assert_array_equal(before, lk.calibrate(D, True, 0.2))


def test_fresh_federation_camb_and_held_scores():
    lk = Likelihood()
    ratio = float(Fraction(3, 7))
    assert RATIOS["federation_matching"] == ratio
    theory, zeq = compute(0.02237, ratio)
    result = lk.evaluate(theory, True)
    stored = json.loads((HERE / "results.json").read_text())["models"]["0.02237"][
        "federation_matching"
    ]
    np.testing.assert_allclose(zeq, stored["z_equality"], rtol=1e-6)
    for metric in (
        "TT",
        "joint_P",
        "conditional_TE",
        "conditional_EE",
        "conditional_P",
    ):
        np.testing.assert_allclose(
            result["scores"][metric],
            stored["amplitude_scale_tilt"]["scores"][metric],
            rtol=3e-5,
        )


def test_published_chain_marginals_and_weighted_ratio_replay():
    from .chains import measure

    primary = json.loads((HERE / "chain_provenance.json").read_text())["channels"]
    replay = measure()
    stored = json.loads((HERE / "chain_results.json").read_text())
    for ch in ("TE", "EE"):
        assert replay[ch] == stored[ch]
        assert replay[ch]["distinct_rows"] == sum(primary[ch]["chain_lengths"])
        refs = [primary[ch]["source_marginals"][p] for p in ("omegabh2", "omegach2")]
        for key, ref in (("density_means", "mean"), ("density_sds", "sd")):
            np.testing.assert_allclose(
                replay[ch][key], [m[ref] for m in refs], rtol=1e-4
            )
        with np.load(HERE / "chain_density_columns.npz") as raw:
            w, b, c = np.vstack([raw[f"{ch}_{i}"] for i in range(1, 5)]).T
        # Independent inverse weighted CDF on sorted, unexpanded rows.
        order = np.argsort(c / b)
        empirical = np.cumsum(w[order]) / w.sum()
        probs = [0.025, 0.16, 0.5, 0.84, 0.975]
        quantiles = (c / b)[order][np.searchsorted(empirical, probs)]
        np.testing.assert_allclose(
            quantiles, replay[ch]["ratio_quantiles_025_16_50_84_975"], atol=4e-4
        )
