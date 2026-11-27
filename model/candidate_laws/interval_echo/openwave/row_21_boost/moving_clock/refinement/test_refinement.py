import json

import numpy as np
from scipy.linalg import expm

from .. import apparatus as lab
from . import measure


def test_exact_floquet_matrix_and_dirac_flow_at_fixed_physical_time():
    z = np.diag([1.0, -1.0])
    x = np.array([[0.0, 1.0], [1.0, 0.0]])
    m = 0.52
    for p in (-1.1, 0.0, 1.1):
        target = expm(-1j * (p * z + m * x))
        prior = 1.0
        for steps in (4, 8, 16, 32):
            a = 1 / steps
            holds = (np.pi, np.pi - 2 * a * m)
            first, last = lab.gates(1.0, 8.0, holds=holds)
            k = a * p
            right = np.diag(np.exp([-1j * k, 1j * k]))
            matched = (
                np.array(
                    [
                        [last[1, 1], last[1, 0] * np.exp(-1j * k)],
                        [last[0, 1] * np.exp(1j * k), last[0, 0]],
                    ]
                )
                @ first
            )
            common = -np.exp(-32.5j * sum(holds))
            expected = np.cos(a * m) * right - 1j * np.sin(a * m) * x
            np.testing.assert_allclose(matched / common, expected, atol=2e-13)
            difference = np.linalg.norm(
                np.linalg.matrix_power(matched / common, steps) - target, ord=2
            )
            if p:
                assert difference < 0.8 * prior
            else:
                assert difference < 1e-11
            prior = difference


def test_covariant_local_action_and_reflection_at_compactification_boundary():
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    gx = x @ z
    gamma = (x, gx)
    metric = np.diag([1, -1])
    for mu in range(2):
        for nu in range(2):
            np.testing.assert_allclose(
                gamma[mu] @ gamma[nu] + gamma[nu] @ gamma[mu],
                2 * metric[mu, nu] * np.eye(2),
                atol=1e-14,
            )
    eta = 0.74
    spin = expm(-eta * z / 2)
    np.testing.assert_allclose(spin.T @ x @ spin, x)
    boost = np.array([[np.cosh(eta), -np.sinh(eta)], [-np.sinh(eta), np.cosh(eta)]])
    for i in range(2):
        np.testing.assert_allclose(
            np.linalg.solve(spin, gamma[i] @ spin), boost[i, 0] * x + boost[i, 1] * gx
        )
    board = lab.federation(2)
    routes = measure.guides(board)
    for n, setup in routes.items():
        ps = setup["ports"]
        assert len(set(ps)) == 2 * n
        assert all(
            b in board.neighbours[a] for a, b in zip(ps, ps[1:] + ps[:1], strict=True)
        )
        assert all(
            v in board.neighbours[a] and v not in ps
            for a, v in zip(ps[:2], setup["vacant"], strict=True)
        )


def test_recorded_terminal_clocks_radar_and_all_controls():
    data = json.loads((measure.HERE / "counts.json").read_text())
    report = measure.read(data)
    assert len(report["rows"]) == 18
    assert all(
        x["clock_resolved"] == {"plane": True, "packet": True} for x in report["rows"]
    )
    assert all(
        not x["mixture_resolved"] and x["early"] == 0 and x["source_error"] < 1e-11
        for x in report["rows"]
    )
    for x in report["rows"]:
        if x["h"] == 0:
            continue
        lo, _, hi = x["ratios"]["plane"]
        assert x["radar"] == 0.6 and x["echo"] > 0
        if x["exchange"] == 1:
            assert hi < 0.6 and x["parked"] == 0
            assert x["h"] * x["current"][1] > 0
        else:
            assert lo > 0.6 and x["parked"] > 0
    shrink = report["contractions"]
    assert len(shrink) == 4
    assert all(
        x["exchange"] == 1 and 0.125 < x["bounds"][0] < x["bounds"][2] < 0.375
        for x in shrink
    )
    for row in data["rows"]:
        s, n = row["exchange"], row["n"]
        if s != 1:
            continue
        k = 2 * np.pi * row["h"] / n
        energy = np.arccos(np.cos(row["delta"]) * np.cos(k))
        ideal = (energy - (row["d"] / row["dt"]) * k) / row["delta"]
        r = next(
            x
            for x in report["rows"]
            if (x["exchange"], x["n"], x["h"]) == (s, n, row["h"])
        )
        assert r["ratios"]["plane"][0] <= ideal <= r["ratios"]["plane"][-1]
        speed = np.cos(row["delta"]) * np.sin(k) / np.sin(energy)
        assert r["current"][0] <= speed <= r["current"][-1]
