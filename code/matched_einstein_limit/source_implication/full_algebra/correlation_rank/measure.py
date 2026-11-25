"""Certify the two-head correlation rank needed on the fixed-volume survey."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from .. import certificate as cert
from .. import enclosure
from ..native_precision import exact_effects, prior

HERE = Path(__file__).parent
iv = enclosure.iv
# A rational grid stores the native derivative intervals without float casts.
SCALE = 10**19


def dyadic(point):
    sign, mantissa, exponent, _ = point
    return (-1 if sign else 1) * mantissa * Fraction(2) ** exponent


def pack(bound):
    lo, hi = (dyadic(p) * SCALE for p in bound._mpi_)
    return [lo.__floor__(), hi.__ceil__()]


def native_derivatives():
    nums, unit, err = exact_effects()
    effects = cert.rational_effects(nums, unit)
    error = sp.Rational(err)
    record = json.loads((HERE.parent / "results.json").read_text())
    x = [
        prior.sym_matrix(block)
        for block in np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    ]
    x = [(block + block.H) / 2 for block in x]
    mass = sum(block.trace() for block in x)
    x = [block / mass for block in x]
    bounds = json.loads((HERE.parent / "certificate.json").read_text())
    radius = cert.ri(str(bounds["state_HS_error_upper"]))
    pauli = (
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    )
    joint = [sp.kronecker_product(a, b) / 2 for a in pauli for b in pauli]
    local = [sp.kronecker_product(a, sp.eye(2)) / 2 for a in pauli] + [
        sp.kronecker_product(sp.eye(2), b) / 2 for b in pauli
    ]
    basis = joint + local
    center, widths, partials = [], [], []
    for panel in effects:
        mean = sum(h.trace() for h in panel) / 8
        centered = [h - mean * sp.eye(4) for h in panel]
        width = iv.sqrt(cert.ri(sum(cert.inner(h, h) for h in centered))) + iv.sqrt(
            8
        ) * cert.ri(error)
        center.append(cert.ri(sum(cert.inner(block, h) for block, h in zip(x, panel))))
        widths.append(cert.ri(error) + radius * width)
        partials.append(
            [cert.ri_band(cert.inner(h, b), 2 * error) for h in panel for b in basis]
            + [cert.ri_band((panel[0] - panel[1]).trace(), 8 * error) / iv.sqrt(8)]
        )
    symmetric = iv.mpf([-1, 1])
    outer = [c + symmetric * w for c, w in zip(center, widths)]
    jac = np.array(partials, dtype=object)
    matrices = {}
    points = enclosure.geometry_jet(center)
    boxes = enclosure.geometry_jet(outer)

    def coefficients(point, box):
        gradient = [
            g + symmetric * sum(abs(h) * w for h, w in zip(row, widths))
            for g, row in zip(point.g, box.h)
        ]
        return np.array(gradient, dtype=object) @ jac

    for label, point, box in zip(("area", "volume"), points, boxes):
        result = coefficients(point, box)[:30].reshape(2, 15)[:, :9].reshape(2, 3, 3)
        matrices[label] = [
            [[pack(v) for v in row] for row in block] for block in result
        ]
    # Compile a concrete signed reader. Round every coefficient and certify its
    # complete traceless derivative against the geometry, including local terms.
    mids = [
        np.array(matrices[k], dtype=float).mean(axis=-1) / SCALE
        for k in ("area", "volume")
    ]
    root = brentq(lambda t: np.linalg.det(mids[0][0] - t * mids[1][0]), 2.8, 2.85)
    price = f"{root:.12f}"
    t = cert.ri(price)
    gain = iv.mpf([str(v) for v in bounds["gain_interval"]])
    assert gain.a > cert.ri("23.1")
    exact = coefficients(points[0] - points[1] * t, boxes[0] - boxes[1] * t) * (-gain)
    approx = np.array([float(v.mid) for v in exact])

    def fixed(value):
        return f"{value:.15f}"

    dyads, reconstructed = [], []
    for sector in range(2):
        offset = 15 * sector
        left, weights, right = np.linalg.svd(approx[offset : offset + 9].reshape(3, 3))
        block = []
        coeff = np.full((3, 3), iv.mpf(0), dtype=object)
        for i in range(2 + sector):
            weight = fixed(weights[i])
            sign = np.sign(left[np.argmax(abs(left[:, i])), i])
            a = [fixed(v * sign) for v in left[:, i]]
            b = [fixed(v * sign) for v in right[i, :]]
            block.append(dict(weight=weight, left=a, right=b))
            coeff += np.outer(
                [cert.ri(v) for v in a], [cert.ri(v) for v in b]
            ) * cert.ri(weight)
        dyads.append(block)
        reconstructed.extend(coeff.flat)
        reconstructed.extend(
            cert.ri(fixed(v)) for v in approx[offset + 9 : offset + 15]
        )
    reconstructed.append(cert.ri(fixed(approx[-1])))
    uncertainty = iv.sqrt(sum((v - w) ** 2 for v, w in zip(exact, reconstructed)))
    assert uncertainty.b < cert.ri("3.106e-7")
    upper = pack(uncertainty)[1]
    reader = dict(
        volume_price=price,
        dyads=dyads,
        one_head=[[fixed(v) for v in approx[o + 9 : o + 15]] for o in (0, 15)],
        parity_contrast=fixed(approx[-1]),
        error_upper_times_scale=upper,
    )
    return matrices, reader


def unpack(blocks):
    return np.array(
        [
            iv.mpf([lo, hi]) / SCALE
            for block in blocks
            for row in block
            for lo, hi in row
        ],
        dtype=object,
    ).reshape(2, 3, 3)


def cofactors_and_det(matrix):
    a, b, c = matrix

    # Cross products are the rows of the cofactor matrix.
    def cross(x, y):
        return [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ]

    bc, ca, ab = cross(b, c), cross(c, a), cross(a, b)
    determinant = sum(x * y for x, y in zip(a, bc))
    return bc + ca + ab, determinant


def singular_bounds(matrix):
    cofactors, determinant = cofactors_and_det(matrix)
    matrix_norm = iv.sqrt(sum(x**2 for x in matrix.flat))
    # Each two-by-two minor has modulus <= sigma_1 sigma_2.
    second = max(abs(h).a for h in cofactors) / matrix_norm.b
    if 0 in determinant:
        return second, iv.mpf(0)
    cofactor_norm = iv.sqrt(sum(x**2 for x in cofactors))
    return second, abs(determinant) / cofactor_norm


def certify(matrices):
    area, volume = (unpack(matrices[key]) for key in ("area", "volume"))
    # Every real price t has t = 4s or t = 4/s with |s| <= 1.
    threshold = cert.ri("4e-8")
    counts = []
    for reciprocal in (False, True):
        stack = [(-1, 0), (0, 0)]  # [n / 2^depth, (n+1) / 2^depth]
        leaves = 0
        depth_max = 0
        while stack:
            n, depth = stack.pop()
            s = iv.mpf([n, n + 1]) / 2**depth
            pencil = (
                (area * s - 4 * volume) if reciprocal else (area - volume * (4 * s))
            )
            singular = [singular_bounds(m) for m in pencil]
            if all(second.a > threshold for second, _ in singular) and any(
                third.a > threshold for _, third in singular
            ):
                leaves += 1
                depth_max = max(depth_max, depth)
            else:
                if depth >= 36:
                    raise ValueError(
                        f"uncertified pencil at numerator {n}, depth {depth}, reciprocal {reciprocal}"
                    )
                stack.extend(((2 * n, depth + 1), (2 * n + 1, depth + 1)))
        counts.append({"leaves": leaves, "maximum_depth": depth_max})
    # One fixed-marginal correlation direction also changes volume.
    assert volume[0, 2, 2].a > cert.ri(".03")
    # Either block admits a price at which its rank falls to exactly two.
    brackets = (("2.8", "2.85"), ("0", "1"))
    for block, (lo, hi) in enumerate(brackets):
        signs = [
            cofactors_and_det(area[block] - volume[block] * cert.ri(t))[1]
            for t in (lo, hi)
        ]
        assert signs[0].a > 0 and signs[1].b < 0
    return {
        "area_discrepancy_lower": "4e-8",
        "entropy_discrepancy_lower": "9.24e-7",
        "chart_subdivisions": counts,
        "price_brackets": brackets,
        "minimum_total_product_terms": 5,
    }


def measure():
    iv.dps = 32
    matrices, reader = native_derivatives()
    return {"scale": SCALE, **matrices, **certify(matrices), "signed_reader": reader}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = measure()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in result.items()
                if k not in ("area", "volume", "signed_reader")
            }
        )
    )
