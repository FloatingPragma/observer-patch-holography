"""Rigorous angular area transfer at exposure four on the certified head equilibrium."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ..full_algebra import certificate as full
from ..full_algebra import native_precision as native
from ..null_optics.certificate import bounded_sources
from . import measure


def moments(board, route, rings, q, source, coeff):
    """Build exact rational centered moments from certified native kernel approximants."""
    stations = np.r_[measure.indices(rings[0]), measure.indices(rings[1])]
    raw = []
    for odd in (False, True):
        h = native.generator(board, route, q, odd)
        # The four-column subblock norm error is <1e-19 under the verified
        # fixed-point kernel bound. The Gram operator error is <3e-19.
        u = native.kernel_panels(h, stations, coeff)[:, 7, 5:, :5]
        re, im = u

        def tr(a):
            return a.transpose(0, 1, 3, 2)

        raw.append([tr(re) @ re + tr(im) @ im, tr(re) @ im - tr(im) @ re])
    raw = np.array(raw, dtype=object)
    denom = native.prior.FRAC**2
    probe = np.array(raw, dtype=float) / denom
    p = measure.probabilities(probe[:, 0] + 1j * probe[:, 1], source)

    # Fixed rational translations improve the interval conditioning of the determinant.
    def axes(marginal):
        center = (
            np.rint(1e8 * marginal @ measure.BEARINGS / p.sum()).astype("int64") * 10
        )
        return np.c_[np.full(5, measure.UNIT), measure.PIXELS - center].astype(object)

    src, dst = axes(p.sum(0)), axes(p.sum(1))
    operators = []
    for v in dst.T:
        for w in src.T:
            weights = np.outer(v, w)
            assert sum(np.abs(weights).ravel()) <= 100 * measure.UNIT**2
            re, im = (
                (weights[None, None, :, :, None, None] * raw)
                .sum(axis=(2, 3))
                .transpose(1, 0, 2, 3)
            )
            operators.append(
                [
                    (sp.Matrix(a) + sp.I * sp.Matrix(b)) / (denom * measure.UNIT**2)
                    for a, b in zip(re, im)
                ]
            )
    return operators


def cofactor(matrix, i, j):
    a, b = [k for k in range(3) if k != i]
    c, d = [k for k in range(3) if k != j]
    return (-1) ** (i + j) * (matrix[a, c] * matrix[b, d] - matrix[a, d] * matrix[b, c])


def ranges(effects):
    src = bounded_sources(sp.eye(4) / 4)
    x, eps = src["state"]
    d, derr = src["unpair"]
    ordinary = [sum(x, sp.zeros(4)) / 2 - w for w in x]
    angles = json.loads(
        (measure.HERE.parent / "null_optics/certificate.json").read_text()
    )["matched_pairing"]

    def interval(name):
        a, b = angles[name]
        return iv.mpf([int(np.floor(a)), int(np.ceil(b))])

    pauli = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    generators = [sp.kronecker_product(sp.eye(2), a) for a in pauli]

    def read(m, s):
        return sum(full.inner(a, b) for a, b in zip(m, s))

    bloch = [full.ri_band(read([g, g], x), eps) for g in generators]
    bn = iv.sqrt(sum(a**2 for a in bloch))
    rots = [
        [(-sp.I * (g * t - t * g) / 2).applyfunc(sp.expand) for t in x]
        for g in generators
    ]

    def panel(e):
        n = max(native.prior.one_norm(h) for h in e)

        def val(rho, radius):
            return full.ri_band(read(e, rho), n * radius + sp.Rational("1e-16"))

        rot = sum(a * val(r, eps) / bn for a, r in zip(bloch, rots))
        return [val(x, eps), val(ordinary, 2 * eps), val(d, derr), rot]

    entries = np.array([panel(e) for e in effects], dtype=object).T.reshape(4, 3, 3)
    m = entries[0]
    cof = np.array(
        [[cofactor(m, i, j) for j in range(3)] for i in range(3)], dtype=object
    )
    total = m[0, 0]
    det = sum(e * c for e, c in zip(m[0], cof[0]))
    area = 4 * det / total**3
    cof[0, 0] -= 3 * det / total

    def slope(e):
        return 4 * sum((cof * e).ravel()) / total**3

    rotate = slope(entries[3])
    ordinary = slope(entries[1]) + interval("ordinary_angle_interval") * rotate
    pairing = slope(entries[2]) + interval("angle_interval") * rotate
    ref = json.loads(
        (measure.HERE.parent / "full_algebra/certificate.json").read_text()
    )["signed_transfers"]

    def band(name):
        lo, hi = json.loads(ref[name], parse_float=str)
        return iv.mpf([lo, hi])

    return dict(
        signed_area=str(area),
        ordinary=str(ordinary),
        unpairing=str(pairing),
        ordinary_calibrated_ratio=str(
            pairing / ordinary * band("ordinary_entropy") / band("pairing_entropy")
        ),
    )


def run():
    iv.dps = 30
    f = measure.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = measure.screens(board, route)
    coeff = native.prior.weights()
    x = measure.front.optics.source()["state"]
    return {str(q): ranges(moments(board, route, rings, q, x, coeff)) for q in (-1, 1)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    a = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if a.write:
        (measure.HERE / "certificate.json").write_text(data)
    else:
        print(data, end="")
