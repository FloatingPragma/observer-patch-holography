"""Recover a real Fourier quadrature on the seven published string directions."""
from fractions import Fraction as F

from matched_einstein_limit.frame import PAIRS, decode

from ...proper_balls.measure import FRAMES


def longitudinal_prediction(rest, speed, ward=(0, 0, 0, 0)):
    """Predict the positive-z reading given c B_0n + B_zn = ward[n]."""
    c = F(speed)
    w0, _, _, wz = map(F, ward)
    return ((5 - 4*c)**2 * F(rest) + (40 - 16*c)*w0 + 16*wz) / 9


def reconstruct(readings, speed, ward=(0, 0, 0, 0)):
    """Use the six other shapes and three Ward components, then predict the seventh.

    Each reading is the same real quadrature of an independently measured
    quadratic tensor response. The Fourier wave covector is (-omega,0,0,k)
    with k nonzero and speed omega/k. Inputs and outputs use exact rationals.
    """
    c = F(speed)
    if c == F(1, 4):
        raise ValueError("two transverse currents are invisible at speed 1/4")
    w0, wx, wy, wz = map(F, ward)
    predicted_z = longitudinal_prediction(readings["rest"], c, (w0, wx, wy, wz))
    # Substitute this prediction before decoding, so the measured z reading
    # remains an independent control even when it disagrees.
    seven = {**readings, "z": predicted_z}
    even = dict(zip(PAIRS, decode([F(seven[shape]) for shape in FRAMES.values()]), strict=True))
    jz = w0 - c*even[0, 0]
    currents = ((wx - even[1, 3] + jz/4) / (c - F(1, 4)),
                (wy - even[2, 3] + jz/4) / (c - F(1, 4)), jz)
    tensor = dict(even)
    for i, ji in enumerate(currents, start=1):
        tensor[0, i] = ji
        tensor[i, i] -= 5*ji/2
        for j in range(i + 1, 4):
            tensor[i, j] -= (ji + currents[j - 1])/4
    return tuple(tensor[pair] for pair in PAIRS), F(readings["z"]) - predicted_z
