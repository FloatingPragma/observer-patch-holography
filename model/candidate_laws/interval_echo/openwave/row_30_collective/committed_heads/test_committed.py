"""Replay the observer verdict and compare native chain propagation independently."""

import json
import numpy as np
from . import measure as m, readback as reader
from scipy.linalg import expm


def test_complete_native_columns_and_committed_six_dimer_source():
    prob, _ = m.readings()
    for weights in prob["full"].values():
        np.testing.assert_allclose(weights, np.full(4, 1 / 16), atol=2e-14)
    n = 6
    meter = 2
    sensor = 5
    sources = [m.initial_word(n, j) for j in (1, 0)]
    p = m.probabilities(8, range(n), sources, meter, sensor)
    # Direct independently assembled physical orientation matrix for one
    # perturbed even-coin, odd-sign sector, without a KG dispersion input.
    nn = 2**n
    eye = np.eye(nn)
    h = np.zeros((nn, nn))
    det = []
    for w in range(nn):
        b = [(w >> i) & 1 for i in range(n)]
        for i in range(n):
            h[w, w ^ (1 << i)] = -0.5
        for i in range(1, n):
            if (b[i - 1], b[i]) != (1, 0):
                h[w, w] += 33 if i == meter else 32
        if b[sensor - 1] == b[sensor]:
            det.append(w)
    u = expm(-2j * h)
    state = eye[:, sources].astype(complex)
    for it in range(4):
        np.testing.assert_allclose(
            p[2, it], (abs(state[det]) ** 2).sum(axis=0), atol=3e-12, rtol=0
        )
        state = u @ state


def test_frozen_route_counts_and_effect_of_retaining_the_head():
    fresh = json.loads((m.HERE / "results.json").read_text())
    route, src, ref, meter, sensor = m.geometry()
    assert [list(r) for r in route] == fresh["route"] and (src, ref, meter, sensor) == (
        4,
        1,
        5,
        7,
    )
    result = reader.evaluate(fresh["banks"], fresh["shots"], fresh["sources"], m.TIMES)
    assert result["clock_max_error"] < 0.003
    for k in (0, 8):
        for meter in m.NAMES:
            p = fresh["summary"][str(k)][meter]["source_response"]
            for exact, ci in zip(
                p, result["source_response"][str(k)][meter], strict=True
            ):
                assert ci[0] <= exact <= ci[1]
        for row in result["conditional"]:
            if row["k"] == k:
                ci = [c["interval"] for c in row["rows"] if c["t"] == 4]
                assert all(reader.nonzero(c) == bool(k) for c in ci)
    controls = result["model_contrasts"]["8"]
    for contrast in ("coin minus 0.5 magnitude", "full minus 0.5 sign"):
        assert all(not reader.nonzero(c) for c in controls[contrast])
    for contrast in ("sign minus 1 magnitude", "full minus 1 magnitude"):
        assert reader.nonzero(controls[contrast][2])
    # Earlier, independent analytic uncommitted transfer uses a Majorana
    # reduction; the new counts were generated from the physical orientation matrix.
    analytic = [
        (
            m.strand.walls(len(route), 8, t)[sensor, src]
            - m.strand.walls(len(route), 8, t)[sensor, ref]
        )
        for t in m.TIMES
    ]
    np.testing.assert_allclose(
        analytic,
        fresh["summary"]["8"]["magnitude"]["source_response"],
        atol=2e-11,
        rtol=0,
    )
