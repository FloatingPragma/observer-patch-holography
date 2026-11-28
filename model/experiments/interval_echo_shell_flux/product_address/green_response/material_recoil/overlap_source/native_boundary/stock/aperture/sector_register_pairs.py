"""Joint laws of two distinct, independently reported members of the same register."""

import json

import numpy as np

from . import sector_fronts_collective as collective
from . import sector_passage_chart as chart


def type_weights(category, counts, gamma, size):
    # Protected words have fixed quotas. The free counts vary across registers.
    k = len(counts)
    mass = (
        np.r_[
            np.asarray(counts, float), np.full(sum(category < 0), size // len(category))
        ]
        / size
    )
    M = np.outer(mass, mass) - np.diag(mass)
    M[:k, :k] += np.asarray(gamma, float) / size
    return mass, M


def distinct_law(hword, category, counts, gamma, size, error):
    h = collective.sector_fronts_counts.group_readings(hword, category)
    mass, M = type_weights(category, counts, gamma, size)
    mu = mass @ h
    delta = h.T @ M @ h / (size - 1)
    np.testing.assert_allclose(delta.sum(1), 0, atol=1e-12)
    # Each row of the dyadic conditional table is within d*error in L1.
    d = h.shape[1] * float(error)
    err = (2 * d + d * d) * abs(M).sum() / (size - 1) + 1e-12
    assert np.all(np.outer(mu, mu) + delta > 0)
    return mu, delta, err


def apparatus():
    words, _, ports, clear, adj, uv, loads, periods = chart.apparatus()
    cfg, ids, freq = np.unique(
        loads > 0, axis=0, return_inverse=True, return_counts=True
    )
    cfg = cfg.astype(int)
    base = freq / len(words)
    _, _, sub, _full, mats = chart.source.sources((8,))
    size, bank = get_bank()
    cat = collective.peer_population.fibers.apparatus()[3]
    counts, gamma = collective.sector_fronts_counts.moments(bank)
    flags = np.eye(len(cfg), dtype=int)[ids]
    report, grid, eps = collective.report_table(sub, mats, flags, 8)
    mu, change, err = distinct_law(report / grid, cat, counts, gamma, size, eps)
    cross = float(cfg[:, 0] @ change @ cfg[:, 3])
    previous = json.loads(
        (collective.HERE / "sector_fronts_collective.json").read_bytes()
    )["forecasts"]["one_eighth"]
    assert abs(cross - previous["per_distinct_pair_covariance"]) <= err + previous[
        "whole_per_member_error"
    ] / (size - 1)
    iid = np.outer(mu, mu)
    calibrated = json.loads((collective.HERE / "sector_calorimetry.json").read_bytes())[
        "ordinary_clock_calibration"
    ]
    bounds = chart.source.sector_regional.calibrate(calibrated)
    refs = {"ordinary": (base, np.stack((base, base), axis=1))}
    thermal = chart.source.thermal.bath()
    for key in ("retained", "exchanged"):
        p = np.bincount(ids, weights=thermal[key].diagonal())
        bounds_cfg = np.array([bounds[key][ids == k].sum(0) for k in range(len(cfg))])
        assert np.all((bounds_cfg[:, 0] < p) & (p < bounds_cfg[:, 1]))
        refs[key] = (p, bounds_cfg)
    # Ordinary delivery through the reserved rail into an empty home 20.
    empty = base * (1 - cfg[:, 1])
    empty /= empty.sum()
    fill = np.zeros_like(base)
    for i in np.flatnonzero(empty):
        v = cfg[i].copy()
        v[1] = 1
        out = np.flatnonzero((cfg == v).all(1)).item()
        fill[out] += empty[i]
    ordinary = np.outer(fill - empty, base)
    delta = {
        "cross": change,
        "full": iid + change - np.outer(base, base),
        "ordinary": ordinary,
    }
    errors = {
        "cross": err,
        "full": err + 4 * len(cfg) * float(eps) + 1e-12,
        "ordinary": 1e-13,
    }
    return (ports, clear, adj, uv, cfg, ids, periods), refs, mu, delta, errors


def get_bank():
    bank = []
    info, _ = collective.peer_population.histogram(16, draw_bank=bank)
    return info["size"], bank
