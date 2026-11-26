"""Finite orientation average which removes a common affine shape response."""

from fractions import Fraction
from itertools import product


def pullback(tensor, shift, signs):
    # Rotate e_i to signs[i] * e_{(i+shift) mod 3}.
    return tuple(
        tuple(
            signs[i] * signs[j] * tensor[(i + shift) % 3][(j + shift) % 3]
            for j in range(3)
        )
        for i in range(3)
    )


def orientation_average(tensor):
    frames = [
        pullback(tensor, shift, signs)
        for shift in range(3)
        for signs in product((-1, 1), repeat=3)
        if signs[0] * signs[1] * signs[2] == 1
    ]
    return tuple(
        tuple(Fraction(sum(q[i][j] for q in frames), len(frames)) for j in range(3))
        for i in range(3)
    )
