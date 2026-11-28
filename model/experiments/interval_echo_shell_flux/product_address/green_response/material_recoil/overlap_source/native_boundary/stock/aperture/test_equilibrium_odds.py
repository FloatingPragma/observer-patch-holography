from fractions import Fraction

import numpy as np
from mpmath import iv

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    equilibrium_odds as eq,
    paired_odds_counts as counts,
)


def test_conditional_test_includes_observed_table_and_handles_fixed_margins():
    assert counts.conditional_tail([[0, 2], [2, 0]], Fraction(1)) == Fraction(1, 6)
    assert counts.conditional_tail([[1, 1], [1, 1]], Fraction(1)) == Fraction(5, 6)
    assert counts.conditional_tail([[1, 1], [1, 1]], Fraction(1, 2)) == Fraction(12, 13)
    assert counts.conditional_tail([[2, 0], [0, 0]], counts.NULL) == 1


def test_shared_lines_are_summed_after_multiplying_replica_traces():
    # An active line that attenuates strongly eventually loses thermal weight.
    # Repeating an already averaged one-copy odds factor misses this effect.
    lines = {2: {(0, 0): 1, (0, 1): 2}}
    banks = [
        {(0, 0): 1, (0, 1): 0.6, (1, 0): 0.7, (1, 1): 0.39},
        {(0, 0): 1, (0, 1): 0.9, (1, 0): 0.8, (1, 1): 0.70},
    ]
    interval = [{k: iv.mpf(v) for k, v in b.items()} for b in banks]
    floor = eq.every_inventory(interval, lines, iv.mpf(1))["odds_floor"]

    def exact(m, n):
        cells = {k: banks[0][k] ** m * banks[1][k] ** n for k in banks[0]}
        u, v, w = cells[1, 0], cells[0, 1], cells[1, 1]
        return (1 + 2 * u + 2 * v + 4 * w) / ((1 + 2 * u) * (1 + 2 * v))

    assert exact(1, 0) ** 300 < floor
    assert min(exact(i, j) for i in range(81) for j in range(61)) >= floor
    np.testing.assert_allclose(exact(300, 300), 1, rtol=0, atol=2e-15)
