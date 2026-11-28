"""Independent small-network and probability checks for unrestricted landmark records."""

from decimal import Decimal, localcontext
from fractions import Fraction as Q

import numpy as np
from sympy import Matrix

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    collective_source as source,
    landmark_entropy as entropy,
    landmark_release as ranging,
)


def test_relay_range_on_two_edges_matches_direct_grounded_inverse():
    n = 6
    adj = np.zeros((n, n), int)
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 3)):
        adj[i, j] = adj[j, i] = 1
    contacts = [(1, 4), (2, 5)]
    prices = [Q(2, 7), Q(5, 9)]
    R, den = ranging.ranges(ranging.grounded_inverse(adj), contacts, prices)
    L = Matrix(np.diag(adj.sum(0)) - adj)
    for (u, v), p in zip(contacts, prices, strict=True):
        L[u, u] += p
        L[v, v] += p
        L[u, v] -= p
        L[v, u] -= p
    G = L[:-1, :-1].inv().row_join(Matrix.zeros(n - 1, 1)).col_join(Matrix.zeros(1, n))
    assert all(
        Q(R[i, j], den) == G[i, i] + G[j, j] - 2 * G[i, j]
        for i in range(n)
        for j in range(n)
    )


def test_exponential_enclosure_includes_high_precision_reference():
    with localcontext() as ctx:
        ctx.prec = 75
        for x in (Q(-53, 7), Q(-1, 3), Q(0), Q(5, 8), Q(81, 7)):
            lo, hi = entropy.exp_bounds(x)
            ref = (Decimal(x.numerator) / Decimal(x.denominator)).exp()
            assert lo <= Q(ref) <= hi
            assert hi - lo < Q(1, 2**44) * max(1, hi)


def test_large_denominator_donor_resolves_equal_first_digit():
    class Digits:
        def __init__(self):
            self.calls = 0

        def integers(self, upper, size):
            assert upper == 2**32
            self.calls += 1
            return np.full(size, 2**31 if self.calls == 1 else 0)

    rng = Digits()
    # The leading digit ties. The true probability lies just above 1/2.
    chosen = source.fractional_donor(rng, Q(2**79 + 1, 2**80), 5)
    assert chosen.all() and rng.calls == 3
