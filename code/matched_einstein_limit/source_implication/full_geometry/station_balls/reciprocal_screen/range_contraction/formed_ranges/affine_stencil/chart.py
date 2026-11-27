"""Positive centered five-point weights from the five-clock spherical branch."""

from itertools import permutations
import json
from math import nextafter, inf
from pathlib import Path
from mpmath import iv

from .. import measure as formed

HERE = Path(__file__).parent


def bounds(v):
    return [nextafter(float(v.a), -inf), nextafter(float(v.b), inf)]


def determinant(a):
    n = len(a)
    result = iv.mpf(0)
    for perm in permutations(range(n)):
        term = iv.mpf(1)
        for i, j in enumerate(perm):
            term *= a[i][j]
        if sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n)) % 2:
            result -= term
        else:
            result += term
    return result


def minor(a, n):
    return [r[:n] for r in a[:n]]


def distance_table(archive, stock, names):
    routes = archive["directed_routes"]
    means = {
        (i, j): sum(
            formed.interval(archive["directed_mean_bands"][q][stock][k]) for q in names
        )
        / len(names)
        for k, (i, j) in enumerate(routes)
        if i < 5 and j < 5
    }
    g = [
        sum(iv.ln(means[i, j] / means[j, i]) for j in range(5) if j != i) / 5
        for i in range(5)
    ]
    return [
        [
            iv.mpf(0)
            if i == j
            else iv.sqrt(means[i, j] * means[j, i] * iv.exp(-g[i] - g[j]))
            for j in range(5)
        ]
        for i in range(5)
    ]


def spherical(d):
    def matrix(z):
        return [[iv.cos(z * x) for x in r] for r in d]

    lo, hi = iv.mpf("0.62"), iv.mpf("0.623")
    box = iv.mpf([lo.a, hi.b])
    m = matrix(box)
    assert all(determinant(minor(m, n)).a > 0 for n in range(1, 5))
    # Expanding the derivative one row at a time checks uniqueness in the
    # enrolled branch. The distances all lie in the minimizing arc range.
    dm = [[-v * iv.sin(box * v) for v in row] for row in d]
    assert sum(determinant(m[:i] + [dm[i]] + m[i + 1 :]) for i in range(5)).b < 0
    assert all((box * v).b < iv.pi.a for row in d for v in row)
    assert determinant(matrix(lo)).a > 0 and determinant(matrix(hi)).b < 0
    for _ in range(45):
        mid = (lo + hi) / 2
        det = determinant(matrix(mid))
        if det.a > 0:
            lo = mid
        elif det.b < 0:
            hi = mid
        else:
            break
    root = iv.mpf([lo.a, hi.b])
    assert float(root.delta) < 1e-8
    return root, matrix(root)


def stencil(d):
    z, c = spherical(d)
    gram = [
        [
            d[0][i]
            * d[0][j]
            * (c[i][j] - c[0][i] * c[0][j])
            / (iv.sin(z * d[0][i]) * iv.sin(z * d[0][j]))
            for j in range(1, 5)
        ]
        for i in range(1, 5)
    ]
    g = minor(gram, 3)
    assert all(determinant(minor(g, j)).a > 0 for j in range(1, 4))
    coordinates = []
    for j in range(3):
        rep = [[gram[i][3] if k == j else g[i][k] for k in range(3)] for i in range(3)]
        coordinates.append(determinant(rep) / determinant(g))
    raw = [-a for a in coordinates] + [iv.mpf(1)]
    weights = [a / sum(raw) for a in raw]
    assert min(w.a for w in weights) > 0
    moment = sum(w * gram[i][i] for i, w in enumerate(weights))
    return z, weights, gram, moment, coordinates


def room(matrix):
    return all(determinant(minor(matrix, j)).a > 0 for j in range(1, 4))


def analyze():
    iv.dps = 43
    archive = json.loads(formed.RECORD.read_text())
    assert archive["station_ports"][:5] == [0, 1, 38, 39, 5]
    results = {}
    for name, probes in (
        ("pooled", ("vacancy", "excess")),
        ("vacancy", ("vacancy",)),
        ("excess", ("excess",)),
    ):
        charts = {
            s: stencil(distance_table(archive, s, probes))
            for s in ("formed", "independent")
        }
        z, w, g, m, a = charts["formed"]
        old = charts["independent"][1]
        # Express the old reference weights on the three independent new
        # vectors. This avoids cancellation from squaring an almost-zero sum.
        v = [old[i] + old[3] * a[i] for i in range(3)]
        residual = iv.sqrt(
            sum(v[i] * g[i][j] * v[j] for i in range(3) for j in range(3))
        )
        data = dict(
            formed_inverse_radius=bounds(z),
            independent_inverse_radius=bounds(charts["independent"][0]),
            weights=[bounds(x) for x in w],
            reference_weights=[bounds(x) for x in old],
            second_moment=bounds(m),
            affine_residual_length=bounds(residual),
            fifth_position_prediction_error=bounds(residual / old[3]),
        )
        if name == "pooled":
            # Congruence in the first three tangent vectors encloses all
            # eigenvalues of the trace-one second moment without a coordinate fit.
            cov = [
                [sum(w[k] * g[i][k] * g[j][k] for k in range(4)) / m for j in range(3)]
                for i in range(3)
            ]
            assert room(
                [
                    [cov[i][j] - iv.mpf(".29") * g[i][j] for j in range(3)]
                    for i in range(3)
                ]
            )
            assert room(
                [
                    [iv.mpf(".40") * g[i][j] - cov[i][j] for j in range(3)]
                    for i in range(3)
                ]
            )
            data["directional_moment_eigenvalues_between"] = [0.29, 0.40]
        results[name] = data
    return results


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2))
