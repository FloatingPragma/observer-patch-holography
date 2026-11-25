"""Exact integer inertia counts certify numerical energy intervals."""

from math import floor, prod

DENOM = 100_000_000


def path_minors(diag, off):
    leading = [1, diag[0]]
    for i, v in enumerate(diag[1:]):
        leading.append(v * leading[-1] - off[i] ** 2 * leading[-2])
    return leading


def below(diag, off, corner, point):
    """Number below point / DENOM for the supplied twice-native Hamiltonian."""
    d = [DENOM * x - 2 * point for x in diag]
    f = [DENOM * x for x in off]
    m = path_minors(d, f)
    if corner:
        # The cycle differs from the open determinant by the edge transposition
        # on the endpoints and by the two oriented full cycles.
        a = DENOM * corner
        minor = path_minors(d[1:-1], f[1:-1])[-1]
        m[-1] += -a * a * minor + 2 * (-1) ** (len(d) - 1) * a * prod(f)
    if 0 in m:
        return None
    return sum((x < 0) != (y < 0) for x, y in zip(m, m[1:]))


def enclose(diag, off, corner, index, estimate):
    center = floor(estimate * DENOM)
    for delta in range(1, 12):
        a, b = center - delta, center + delta + 1
        ca, cb = below(diag, off, corner, a), below(diag, off, corner, b)
        if ca is not None and cb is not None and ca <= index < cb:
            return [a, b]
    raise AssertionError("Could not enclose the requested eigenvalue")
