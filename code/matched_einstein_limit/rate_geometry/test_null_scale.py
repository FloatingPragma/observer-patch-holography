"""Curvature and clock stencil for the fixed-track conformal ambiguity."""

import sympy as sp


def test_ricci_from_christoffel_at_reference_event():
    t, x, y, z = coordinates = sp.symbols("t x y z")
    bx, by, bz = b = sp.symbols("b_x b_y b_z")
    origin = dict.fromkeys(coordinates, 0)
    f = t * (bx * x + by * y + bz * z)
    eta = sp.diag(-1, 1, 1, 1)
    metric = sp.exp(2 * f) * eta
    inverse = metric.inv()

    def gamma(a, m, n):
        return (
            sum(
                inverse[a, j]
                * (
                    metric[j, n].diff(coordinates[m])
                    + metric[j, m].diff(coordinates[n])
                    - metric[m, n].diff(coordinates[j])
                )
                for j in range(4)
            )
            / 2
        )

    assert all(
        gamma(a, m, n).subs(origin) == 0
        for a in range(4)
        for m in range(4)
        for n in range(4)
    )
    ricci = sp.Matrix(
        4,
        4,
        lambda m, n: sum(
            gamma(a, m, n).diff(coordinates[a]) - gamma(a, m, a).diff(coordinates[n])
            for a in range(4)
        ).subs(origin),
    )
    n = sp.symbols("n_x n_y n_z")
    k = sp.Matrix((1, *n))
    assert sp.expand((k.T @ ricci @ k)[0]) == -4 * sum(p * q for p, q in zip(b, n))
    assert metric.subs(t, 0) == eta


def test_scale_stencil_normalization_and_noise():
    w, slope, curvature, jerk, fourth, h = sp.symbols("w v a j m h", nonzero=True)
    s = sp.Symbol("s")
    scale = w + slope * s + curvature * s**2 / 2 + jerk * s**3 / 6 + fourth * s**4 / 24
    result = 2 * (scale.subs(s, h) - 2 * w + scale.subs(s, -h)) / h**2
    assert sp.simplify(result - 2 * curvature) == fourth * h**2 / 6
    assert sum(abs(c) for c in (2, -4, 2)) == 8
