"""Exact open-chain singular-frequency brackets for the eighth-port countertest."""

from fractions import Fraction


def below(n, tick, scale):
    """Sturm variations of the tridiagonal C C^T at (tick/scale)^2."""
    unit, square = scale**2, tick**2
    previous = 1
    current = 256 * unit - square
    signs = [1]
    if current:
        signs.append(1 if current > 0 else -1)
    for _ in range(n):
        previous, current = (
            current,
            (257 * unit - square) * current - 256 * unit**2 * previous,
        )
        if current:
            signs.append(1 if current > 0 else -1)
    if current == 0:
        raise ValueError("A grid endpoint is an exact eigenfrequency")
    return sum(a != b for a, b in zip(signs, signs[1:]))


def bracket(n, index, scale):
    lo, hi = 15 * scale, 17 * scale
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if below(n, mid, scale) > index:
            hi = mid
        else:
            lo = mid
    return Fraction(lo, scale), Fraction(hi, scale)


def gap(n=49, scale=100000):
    low = [bracket(n, i, scale) for i in range(8)]
    high = [bracket(n, n - i, scale) for i in range(8)]
    return (
        sum(lower for lower, _ in high) - sum(upper for _, upper in low),
        sum(upper for _, upper in high) - sum(lower for lower, _ in low),
    )


if __name__ == "__main__":
    print(tuple(map(float, gap())))
