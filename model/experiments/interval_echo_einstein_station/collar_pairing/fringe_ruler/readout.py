"""Fringe phase differences and simultaneous intervals, from counts alone."""

import numpy as np


def radius(observed, shots):
    # Each terminal bin has a Bernoulli marginal. Inverting the two-tail
    # Bernstein bound gives sqrt(2*h*x)+4*x for observed frequency h.
    x = np.log(2 * 400 / 0.01) / shots
    return np.sqrt(2 * observed * x) + 4 * x


def quadratures(rows, shots, systematic=0.0):
    points = np.array(rows)[:, :4] / shots
    widths = radius(points, shots) + systematic
    dc = points[:4].mean(axis=0)
    dc_err = widths[:4].mean(axis=0)
    cross = (points[0] - points[2] + 1j * (points[3] - points[1])) / 2
    xerr = (widths[0] + widths[2]) / 2
    yerr = (widths[3] + widths[1]) / 2
    for row, angle in ((4, np.pi / 4), (5, -np.pi / 4)):
        pred = dc + np.cos(angle) * cross.real - np.sin(angle) * cross.imag
        tol = (
            widths[row] + dc_err + abs(np.cos(angle)) * xerr + abs(np.sin(angle)) * yerr
        )
        assert np.all(abs(pred - points[row]) <= tol)
    return cross, np.hypot(xerr, yerr)


def source_change(ref, other):
    """Continuous change per port on the affine mixture. Every disk excludes zero."""
    z, err = ref
    v, verr = other
    assert np.all(abs(z) > err) and np.all(abs(v) > verr)
    d = np.angle(v * z.conj())
    w = np.arcsin(err / abs(z)) + np.arcsin(verr / abs(v))
    assert np.all(abs(d) + w < np.pi)
    # The linear interpolation between the two endpoint fringes cannot pass
    # through zero. Its phase change is uniquely the principal argument.
    return d - d[0], w + w[0]


def interval(value, width):
    return [float(value - width), float(value + width)]


def examine(receipt, systematic=0.0):
    n = receipt["shots"]
    shifts = {}
    for q, panels in receipt["counts"].items():
        raw = {s: quadratures(rows, n, systematic) for s, rows in panels.items()}
        scrambled, allowance = raw["scrambled"]
        assert np.all(abs(scrambled) < allowance)
        z, e = raw["product"]
        v, a = raw["formed"]
        half, he = raw["half"]
        assert np.all(abs(half - (z + v) / 2) < he + (e + a) / 2)
        shifts[q] = {}
        for s in ("half", "formed", "ordinary"):
            delta, err = source_change(raw["product"], raw[s])
            shifts[q][s] = [
                interval(x, y) for x, y in zip(delta[1:], err[1:], strict=True)
            ]
    clock = receipt["four_equal_ticks"]
    c = clock["one_tick_hits"] / n
    w = radius(c, n)
    tc = list(map(float, 8 * np.arcsin(np.sqrt([c - w, c + w]))))
    return dict(
        baselines=[[receipt["ports"][0], b] for b in receipt["ports"][1:]],
        four_tick_interval=tc,
        phase_changes_radians=shifts,
    )
