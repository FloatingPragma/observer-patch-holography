"""The joint modular price at two homes and its all-rate exterior distance contribution."""

from fractions import Fraction

import numpy as np
from scipy.stats import beta

from . import peer_acquisition, sector_home


def reference(words, subset, full, mats):
    heads = (words[:, [0, 3]] > 0).astype(int)
    cell = 2 * heads[:, 0] + heads[:, 1]
    ordinary = np.bincount(cell, minlength=4)
    odds = Fraction(int(ordinary[0] * ordinary[3]), int(ordinary[1] * ordinary[2]))
    ell = np.log(float(odds))
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(words))
    tables = peer_acquisition.returns.ticket_tables(mats)
    end, _ = peer_acquisition.acquire(
        np.random.default_rng(884991), 2**19, subset, native >= 0, tables
    )
    ids = native[end]
    counts = np.bincount(cell[ids], minlength=4)
    lo, hi = (
        beta.ppf(0.001 / 8, counts, len(ids) - counts + 1),
        beta.isf(0.001 / 8, counts + 1, len(ids) - counts),
    )
    log_band = np.log(
        [lo[0] * lo[3] / (hi[1] * hi[2]), hi[0] * hi[3] / (lo[1] * lo[2])]
    )
    assert log_band[0] < ell < log_band[1] < 0
    # The equal-mixture arms were found in a search for fixed home
    # marginals, mean exchange heats and total physical board energy.
    arms = np.array([[98, 19], [16, 100]])
    assert np.all(np.bincount(ids, minlength=len(words))[arms] > 0)
    np.testing.assert_array_equal(
        np.sort(words[arms[0], :4], axis=0), np.sort(words[arms[1], :4], axis=0)
    )
    heat = sector_home.source.thermal.apparatus()[-1][0].diagonal(axis1=1, axis2=2).T
    total = sector_home.source.balls.energy.returned(words)[-1]
    for read in (heat, total):
        np.testing.assert_array_equal(read[arms[0]].sum(0), read[arms[1]].sum(0))
    assert cell[arms].tolist() == [[3, 0], [2, 1]]
    return (
        arms,
        ell,
        {
            "reference_counts": ordinary.tolist(),
            "ordinary_joint_odds": str(odds),
            "native_reference_counts": counts.tolist(),
            "log_odds_999": log_band.tolist(),
            "control_word_ids": arms.tolist(),
            "copies_found": np.bincount(ids, minlength=len(words))[arms].tolist(),
            "matched_heat_times_four": heat[arms[0]].mean(0).tolist(),
            "matched_total_energy": float(total[arms[0]].mean() / 4),
            "matched_charge_gain": -ell / 2,
            "charge_band": (-log_band[::-1] / 2).tolist(),
        },
    )


def mul(x, y):
    ends = np.multiply.outer(x, y)
    return np.array([ends.min() - 1e-12, ends.max() + 1e-12])


def all_rates(
    words, delta, subset, mats, means, price_error, standard, standard_error, ell
):
    occ = (words[:, :4] > 0).astype(int)
    x, y = occ[:, 0], occ[:, 3]
    X, Y = Fraction(int(x.sum()), len(x)), Fraction(int(y.sum()), len(y))
    assert (X, Y) == (Fraction(45, 56), Fraction(11, 14))
    joint = x * y
    centered = 56 * joint - 44 * x - 45 * y
    names = ("8", "20", "32", "44", "joint", "centered")
    stats = sector_home.rate_bounds(
        np.column_stack((occ, joint, centered)), delta, subset, mats, names
    )
    bins = {"full_geometry": [], "connected_geometry": []}
    cert = []
    assert standard > standard_error
    scale = 1 / np.array([standard + standard_error, standard - standard_error])
    for i, row in enumerate(stats["8"]):
        p = row["rate"]
        v = np.array([stats[name][i]["per_rate"] for name in names])
        covariance = v[-1] / 56 - mul(p, mul(v[0], v[3]))[::-1]
        covariance += np.array([-1e-10, 1e-10])
        cert.append(
            {
                "rate": p,
                "covariance_per_rate": covariance.tolist(),
                "charge_per_rate": mul([-ell, -ell], covariance).tolist(),
            }
        )
        full = sum(
            mul([m - price_error, m + price_error], z)
            for m, z in zip(means, v[:5], strict=True)
        )
        conn = mul([means[-1] - price_error, means[-1] + price_error], covariance)
        for name, value in (("full_geometry", full), ("connected_geometry", conn)):
            bins[name].append({"rate": p, "per_rate": mul(value, scale).tolist()})
    return {
        "bins": bins,
        "connected_certificate": cert,
        "minimum_covariance_gain_per_rate": min(
            z["covariance_per_rate"][0] for z in cert
        ),
        "minimum_modular_gain_per_rate": min(z["charge_per_rate"][0] for z in cert),
    }
