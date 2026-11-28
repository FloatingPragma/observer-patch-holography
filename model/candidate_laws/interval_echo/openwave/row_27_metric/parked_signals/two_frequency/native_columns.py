"""Verify the enrolled one-seam columns of Interval Echo by exact rational arithmetic."""

from fractions import Fraction as F


def check():
    loads = [(a, b) for a in range(3) for b in range(3)]
    labels = [(c, d) for c in range(2) for d in range(-2, 3)]
    writes = []
    for a, b in loads:
        low, odd = divmod(a + b, 2)
        writes.append(
            [
                9 * labels.index((c, a - b))
                + loads.index((low + c * odd, low + (1 - c) * odd))
                for c in (0, 1)
            ]
        )

    # Coordinates B and R denote the physical column (B, R/sqrt(2)).
    def echo(B, R):
        new_blank = [F(R[i] + R[j], 2) for i, j in writes]
        new_record = R.copy()
        for x, (i, j) in enumerate(writes):
            new_record[i] += B[x] - new_blank[x]
            new_record[j] += B[x] - new_blank[x]
        return new_blank, new_record

    def act(B, R):
        _, v = echo(B, R)
        coeffs = []
        for i, (c, d) in enumerate(labels):
            rev = labels.index((c, -d))
            for x in range(9):
                coeffs.append(
                    F(v[9 * i + x] - v[9 * rev + x], 2) + 8 * d * d * v[9 * i + x]
                )
        return echo([0] * 9, coeffs)

    checked = 0
    for a, b in [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2)]:
        x = loads.index((a, b))
        B, R = [0] * 9, [0] * 90
        B[x] = 1
        expected = [F(0)] * 9
        expected[x] += F(1, 2) + 8 * (a - b) ** 2
        expected[loads.index((b, a))] -= F(1, 2)
        assert act(B, R) == (expected, [0] * 90)
        checked += 1
        for magnitude, parity, E in ((0, 1, 0), (2, 1, 32), (2, -1, 33)):
            R = [0] * 90
            for j, (_, d) in enumerate(labels):
                if abs(d) == magnitude:
                    R[9 * j + x] = parity if d < 0 else 1
            assert act([0] * 9, R) == ([0] * 9, [E * n for n in R])
            checked += 1
    return checked
