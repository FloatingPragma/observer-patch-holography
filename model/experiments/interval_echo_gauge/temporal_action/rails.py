"""Two-rail phase accumulation and digit-by-digit physical X readback."""

from functools import lru_cache
from math import prod
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

WIDTH = 40
GRID = 1 << WIDTH


@lru_cache(maxsize=1)
def hardware():
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    h = h0 + 8 * h1
    words = ((2, 0), (0, 2))
    support = [law.SRC.index(w) for w in words]
    block = h[np.ix_(support, support)]
    np.testing.assert_allclose(np.delete(h[:, support], support, axis=0), 0, atol=1e-13)
    x = np.array([[0, 1], [1, 0]])
    np.testing.assert_allclose(block, 32.5 * np.eye(2) - x / 2, atol=1e-13)
    # Orient the marked private-reference port as the right rail. Its
    # squared-load generator is 0 or 32 on these two alternatives.
    for load in range(3):
        column = np.eye(law.DIM)[:, law.SRC.index((load, 0))]
        np.testing.assert_allclose(h1 @ column, load**2 * column, atol=1e-13)
    private = 8 * np.diag([0.0, 4.0])
    quarter = expm(-1j * np.pi * private / 64)
    splitter = quarter @ expm(-1j * np.pi * block / 2) @ quarter
    splitter *= np.exp(-1j * np.angle(splitter[0, 0]))
    np.testing.assert_allclose(
        splitter, np.array([[1, 1], [1, -1]]) / np.sqrt(2), atol=1e-13
    )
    b = federation(0)
    p, q = next(
        (a, c) for a, c in b.seams if a // 12 == c // 12 == 1 and b.degree[c] == 5
    )
    # A zero standard on the spare contact reads logical 1 as load two
    # on the right rail. Echo and inverse restore its sharp occupation.
    echo = law.echo()
    bright = np.tile([pair == (1, 1) for pair in law.SRC], law.DIM // law.N_SRC)
    for bit in (0, 1):
        initial = np.eye(law.DIM)[:, law.SRC.index((2 * bit, 0))]
        after = echo @ initial
        np.testing.assert_allclose(
            echo.T @ (after * (bright == bit)), initial, atol=1e-13
        )
    return int(p), int(q)


def read_reference():
    left, right = hardware()
    return min(i for i in federation(0).neighbours[right] if i != left)


@lru_cache(maxsize=None)
def rail_pulse(power, sign):
    """Native phase on each Fourier digit from ONE selected signed input bit."""
    denominators = 1 << np.arange(1, WIDTH + 1, dtype=np.int64)
    numerator = (sign * (1 << power)) % denominators
    return np.exp(2j * np.pi * (numerator / denominators))


def significant(value):
    n = abs(int(value))
    while n:
        least = n & -n
        yield least.bit_length() - 1
        n -= least


class Meter:
    def __init__(self):
        hardware()
        self.words = 0
        self.max_magnitude = 0
        self.max_budget = 0

    def read(self, terms):
        """Polynomial terms each supply a sign and one or more input words.

        Only their active digit settings control pulse gates. No sum or product
        of the operands is used as an angle or as a detector outcome.
        """
        terms = list(terms)
        budget = sum(prod(abs(int(v)) for v in term) for term in terms)
        if budget >= GRID // 2:
            raise OverflowError("Window exceeds the declared signed word range")
        self.max_budget = max(budget, self.max_budget)
        rails = np.ones(WIDTH, complex)
        for term in terms:
            if 0 in term:
                continue
            sign = 1
            powers = [0]
            for value in term:
                if value < 0:
                    sign = -sign
                powers = [
                    p + j for p in powers for j in significant(value) if p + j < WIDTH
                ]
            for p in powers:
                rails *= rail_pulse(p, sign)
        word = 0
        residual = 0.0
        for j in range(WIDTH):
            # Lower public bits set an inverse reference on this rail.
            phase = rails[j] * np.exp(-2j * np.pi * word / (1 << (j + 1)))
            bit = int(phase.real < 0)
            residual = max(residual, abs(phase - (-1) ** bit))
            word |= bit << j
        if residual > 2e-9:
            raise ArithmeticError("A nominally sharp X read is ambiguous")
        signed = word if word < GRID // 2 else word - GRID
        self.max_magnitude = max(self.max_magnitude, abs(signed))
        self.words += 1
        return signed


def fringe(value, scale, altered=False):
    """X and Y read after independent digit-selected phases of the public word."""
    phase = 1 + 0j
    kappa = 4 if altered else 8
    for p in significant(value):
        angle = np.pi / 2 * ((1 << p) / scale) * np.sign(value)
        # Positive laboratory exposure modulo a full native marked-rail turn.
        dwell = (-angle / (4 * 8)) % (2 * np.pi / (4 * 8))
        phase *= np.exp(-4j * kappa * dwell)
    # Changed private strength also affects the analyzer reference when it is
    # changed. Here only the data-pulse contact is changed; preparation and
    # the independent X/Y analyzers keep their calibrated native settings.
    return np.clip((1 + np.array([phase.real, phase.imag])) / 2, 0, 1)
