"""Exact native graph, independent receptions, and frozen record replay."""

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import runpy

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model.carrier import isomorphism
from . import apparatus, measure, source_check
from .readback import decode


def test_exact_native_orbit_and_current():
    base = Path(__file__).resolve().parents[5]
    b14 = runpy.run_path(str(base / "code/b14_jacobi/port_dual_metric_selector.py"))
    group = b14["group_rows"]()
    proj = b14["spectral_projectors"](group)
    Q5 = b14["Q5"]
    orbit = {tuple(sorted((g[0], g[1]))) for g in group}
    graph = apparatus.cell()
    mapped = isomorphism(orbit, graph.seams)
    assert mapped == [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 11]
    assert {tuple(sorted((mapped[a], mapped[b]))) for a, b in orbit} == set(graph.seams)
    rings, anti, _ = apparatus.geometry()
    adj = np.zeros((12, 12), int)
    for a, b in graph.seams:
        adj[a, b] = adj[b, a] = 1
    assert np.array_equal(adj[np.ix_(anti, anti)], adj)
    # sin(2*pi*k/5)/sin(2*pi/5) in Q(sqrt(5)).
    golden = Q5(Fraction(-1, 2), Fraction(1, 2))
    sine = [Q5(), Q5(1), golden, -golden, Q5(-1)]
    for m, sector in ((1, "three_plus"), (2, "three_minus")):
        P = proj[sector]
        for p, rf in enumerate(rings):
            sign = np.zeros(12, int)
            phase = np.zeros(12, int)
            for j, x in enumerate(rf):
                sign[x], sign[anti[x]] = 1, -1
                phase[x] = phase[anti[x]] = j * m % 5
            # Integer coefficient vectors in Z[zeta5]. A polynomial of
            # degree at most four vanishes at zeta5 iff all entries coincide.
            f = np.zeros((12, 5), int)
            f[np.arange(12), phase] = sign
            for v, eigen_product in (
                (abs(f), -abs(f)),
                (
                    f,
                    sum(
                        c * np.roll(f, k * m, axis=1)
                        for k, c in ((1, 1), (-1, 1), (2, -1), (-2, -1))
                    ),
                ),
            ):
                residual = adj @ v - eigen_product
                assert np.all(residual == residual[:, :1])

            # A pole has zero current; the sign and relative phase of
            # each remaining edge follow from physical ring addresses.
            def raw(r):
                return sum(
                    sine[(phase[b] - phase[a]) % 5] * int(sign[a] * sign[b])
                    for a, b in zip(r, r[1:] + r[:1], strict=True)
                )

            d = raw(rf)
            bp = mapped.index(p)
            for q, rg in enumerate(rings):
                assert raw(rg) / d == P[bp][mapped.index(q)] / P[bp][bp]


def records(name):
    return json.loads((measure.HERE / name).read_text())


def test_frozen_acquisitions_and_held_detectors():
    primary = records("records.json")
    repeat = records("source_counts.json")
    assert measure.acquire() == primary
    assert source_check.acquire() == repeat
    report = decode(primary)
    assert report["calibration_max_error"] < report["calibration_radius"]
    assert report["absolute_forecasts"]["isolated"][:2] == [348, 0]
    assert report["absolute_forecasts"]["unsorted"][:2] == [24, 24]
    assert all(h[-2:] == [5, True] for h in report["arrival_frequencies"].values())
    for m, sign in ((1, 1), (2, -1)):
        *_, lo, hi = report["marked_pair"][f"winding{m}"]
        assert lo < sign / np.sqrt(5) < hi and lo * hi > 0
        lo, hi = report["all_neighbor_ratios"][f"winding{m}"]
        assert lo * hi > 0
    shape = report["open_marked_ratio_and_shape_exclusions"]
    assert all(r[-2] * (1 if r[2] == 1 else -1) > 0 for r in shape)
    assert sum(r[-1] for r in shape) > 0
    result = source_check.compare(primary, repeat)
    for lev, _, m, _, source, _, neighbor, nb, nj, _ in result["detector_comparisons"]:
        if lev == -1 and source == "odd":
            assert nb == nj == 0
        else:
            assert nb > 0 and nj > 0
        if source == "odd":
            assert (neighbor > 0) == (m == 1)
    for rate, hand in result["clocks"].items():
        expected = np.array([0, 0.5, 1, 2, 3, 4, 5, 6]) * float(rate)
        assert max(abs(expected - hand)) < result["hand_angular_radius"]


def test_incomplete_panels_and_changed_calibration_rejected():
    p = records("records.json")
    bad = deepcopy(p)
    bad["circulation"][1] = bad["circulation"][0]
    with pytest.raises(ValueError):
        decode(bad)
    bad = deepcopy(p)
    bad["antipodes"][0] = 0
    with pytest.raises(ValueError):
        decode(bad)
    silent = deepcopy(p)
    silent["circulation"][0]["edges"] = [[0, 0, p["shots"]] for _ in range(30)]
    report = decode(silent)
    assert report["marked_pair"]["winding1"][-2:] == [None, None]
    assert report["all_neighbor_ratios"]["winding1"] is None
    repeat = records("source_counts.json")
    repeat["banks"][2] = repeat["banks"][1]
    with pytest.raises(ValueError):
        source_check.compare(p, repeat)


def test_matched_half_rate_pulse():
    # Each addressed native exchange on the slowed material reproduces the full
    # native splitter after doubling dwell and compensating the uniform phase
    # of its two endpoints with their private phase references.
    _, h = apparatus.local_source(1, 8)
    _, hh = apparatus.local_source(0.5, 8)
    for t in (np.pi / 5, np.pi / 2, np.pi):
        np.testing.assert_allclose(
            np.exp(8j * t) * expm(-2j * t * hh), expm(-1j * t * h), atol=3e-14
        )
