"""Verify action-readout constraints using the exact-source radar-gradient enclosure."""

import json
from fractions import Fraction
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

HERE = Path(__file__).resolve().parent
RADAR = HERE.parent / "results.json"
CALIBRATION = HERE.parents[1] / "certificate.json"
iv = MPIntervalContext()
iv.dps = 32
ORBIT_SCALE = 10**11


def matrices(source, field):
    return [
        [[iv.mpf(v) / source["scale"] for v in row] for row in block]
        for block in source[field]
    ]


def enclosure(value, scale):
    """Return outward integer endpoints, using the exact directed binary bounds."""
    exact = []
    for negative, mantissa, exponent, _ in value._mpi_:
        exact.append(
            (-1 if negative else 1) * mantissa * Fraction(2) ** exponent * scale
        )
    return [exact[0].__floor__(), exact[1].__ceil__()]


def block_difference(blocks):
    return [a - b for even, odd in zip(*blocks) for a, b in zip(even, odd)]


def blind_tangent(area, volume, gain):
    a, v = block_difference(area), block_difference(volume)
    v2 = sum(x**2 for x in v)
    wedge2 = sum(
        (a[i] * v[j] - a[j] * v[i]) ** 2 for i in range(len(a)) for j in range(i)
    )
    if v2.a <= 0 or wedge2.a <= 0:
        raise ValueError("no certified hidden fixed-volume response")
    w_norm = iv.sqrt(wedge2 / v2)
    price = sum(x * y for x, y in zip(a, v)) / v2
    yz = a[5] - price * v[5]
    return {
        "scale": 10**9,
        "volume_difference_norm": enclosure(iv.sqrt(v2), 10**9),
        "entropy_slope_per_unit": enclosure(gain * w_norm / iv.sqrt(2), 10**9),
        "parity_times_yz_over_two_slope": enclosure(-iv.sqrt(2) * yz / w_norm, 10**9),
    }


def shape(matrix):
    norm2 = sum(x**2 for row in matrix for x in row)
    if norm2.a <= 0:
        return None
    cofactor = [
        [
            matrix[(i + 1) % 3][(j + 1) % 3] * matrix[(i + 2) % 3][(j + 2) % 3]
            - matrix[(i + 1) % 3][(j + 2) % 3] * matrix[(i + 2) % 3][(j + 1) % 3]
            for j in range(3)
        ]
        for i in range(3)
    ]
    norm = iv.sqrt(norm2)
    second = iv.sqrt(sum(x**2 for row in cofactor for x in row)) / norm2
    third = abs(sum(x * y for x, y in zip(matrix[0], cofactor[0]))) / (norm2 * norm)
    return norm, second, third


def orbit_lower(left, right):
    even, odd = shape(left), shape(right)
    if even is None or odd is None:
        return iv.mpf(0)
    gap = max(abs(even[1] - odd[1]).a / 2, abs(even[2] - odd[2]).a / 3)
    return min(even[0].a, odd[0].a) * gap


def certify_chart(area, volume, reciprocal):
    floor = iv.mpf(1) / ORBIT_SCALE
    stack = [(-1, 1, 0)]
    depths = []
    while stack:
        left, right, depth = stack.pop()
        price = iv.mpf([left, right]) / (2**depth)
        blocks = [
            [
                [price * a - v if reciprocal else a - price * v for a, v in zip(ar, vr)]
                for ar, vr in zip(ab, vb)
            ]
            for ab, vb in zip(area, volume)
        ]
        if orbit_lower(*blocks).a > floor.b:
            depths.append(depth)
        else:
            if depth >= 20:
                raise ValueError("shared shape has not been excluded")
            middle = left + right
            stack.extend(
                ((left * 2, middle, depth + 1), (middle, right * 2, depth + 1))
            )
    return {"leaves": len(depths), "max_depth": max(depths)}


def certify():
    source = json.loads(RADAR.read_text())
    calibration = json.loads(CALIBRATION.read_text(), parse_float=str)
    gain = iv.mpf(calibration["gain_interval"])
    area, volume = matrices(source, "area"), matrices(source, "volume")
    return {
        "blind_signed_curve": blind_tangent(area, volume, gain),
        "shared_shape_price_charts": [
            certify_chart(area, volume, c) for c in (False, True)
        ],
        "shared_shape_entropy_error_lower_times_1e15": enclosure(
            gain / (2 * iv.sqrt(2) * ORBIT_SCALE), 10**15
        )[0],
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
