"""Certify the ten-frame decoder on rational boxes for measured velocities."""

import json
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

from ..frame import PAIRS, inverse

HERE = Path(__file__).parent


def multiply(left, right):
    values = [a * b for a, b in product(left, right)]
    return min(values), max(values)


def square(bounds):
    a, b = bounds
    return (0 if a <= 0 <= b else min(a * a, b * b), max(a * a, b * b))


def row(velocity):
    components = ((F(1), F(1)), *velocity)
    return [
        square(components[i])
        if i == j
        else multiply((2, 2), multiply(components[i], components[j]))
        for i, j in PAIRS
    ]


def symmetric_menu(speed, width=0):
    if not 0 < speed or 2 * speed**2 >= 1 or width < 0:
        raise ValueError(
            "a positive speed with subluminal diagonal frames and nonnegative errors is required"
        )
    centers = [[0, 0, 0]]
    for i in range(3):
        for sign in (1, -1):
            v = [0, 0, 0]
            v[i] = sign * speed
            centers.append(v)
    for i, j in combinations(range(3), 2):
        v = [0, 0, 0]
        v[i] = v[j] = speed
        centers.append(v)
    return [[(F(x) - width, F(x) + width) for x in v] for v in centers]


def certify(velocities):
    if len(velocities) != 10 or any(len(v) != 3 for v in velocities):
        raise ValueError("ten three-component velocity boxes are required")
    velocities = [[(F(lo), F(hi)) for lo, hi in v] for v in velocities]
    if any(lo > hi for v in velocities for lo, hi in v):
        raise ValueError("velocity bounds are reversed")
    dilation = []
    for v in velocities:
        magnitudes = [square(b) for b in v]
        if sum(hi for lo, hi in magnitudes) >= 1:
            raise ValueError("a velocity box meets or crosses the light cone")
        dilation.append(1 - sum(lo for lo, hi in magnitudes))
    boxes = [row(v) for v in velocities]
    nominal = []
    for v in velocities:
        center = [((a + b) / 2, (a + b) / 2) for a, b in v]
        nominal.append([lo for lo, hi in row(center)])
    try:
        decoder = inverse(nominal)
    except ValueError as exc:
        raise ValueError(
            "nominal observers do not distinguish every tensor component"
        ) from exc
    error = [
        [max(abs(lo - n), abs(hi - n)) for (lo, hi), n in zip(r, nr, strict=True)]
        for r, nr in zip(boxes, nominal, strict=True)
    ]
    beta = max(
        sum(abs(a) * sum(er) for a, er in zip(r, error, strict=True)) for r in decoder
    )
    if beta >= 1:
        raise ValueError(
            "the uncertainty intervals are too wide for this inverse certificate"
        )
    gain = max(
        sum(abs(a) * z for a, z in zip(r, dilation, strict=True)) for r in decoder
    )
    return beta, gain / (1 - beta)


def examples():
    result = []
    for v in (F(1, 2), F(1, 4), F(1, 10)):
        exact = certify(symmetric_menu(v))[1]
        assert exact == 2 / v**2 - 1
        width = F(1, 2048)
        beta, bound = certify(symmetric_menu(v, width))
        result.append(
            dict(
                axial_speed=str(v),
                component_half_width=str(width),
                exact_gain=str(exact),
                inverse_residual_upper=str(beta),
                measurable_gain_upper=str(bound),
            )
        )
    return result


if __name__ == "__main__":
    print(json.dumps(examples(), indent=2))
