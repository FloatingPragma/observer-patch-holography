"""Independent checks on finite register enumeration and physical recombination."""

from fractions import Fraction as Q

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    peer_fibers as fibers,
    peer_population as population,
)


def test_compressed_fibers_match_full_enumeration():
    for k in (1, 2):
        compressed, arms = population.histogram(k)
        direct, other = population.histogram(k, collapse=False)
        assert compressed == direct
        assert all(np.array_equal(z, other[key]) for key, z in arms.items())


def test_integer_projection_covers_an_independently_bounded_fiber():
    basis = np.array([[1, -1, 0], [0, 1, -1]])
    base = np.array([3, 1, 2])
    found = {tuple(v) for v in fibers.populations(basis, base)}
    assert found == {(i, j, 6 - i - j) for i in range(7) for j in range(7 - i)}


def test_same_peer_rule_has_both_signs_at_fixed_ordinary_readings():
    for k, sign in ((1, 0), (8, -1), (16, 1)):
        row, arms = population.histogram(k)
        original = fibers.source.records()[1]["original"]
        w = fibers.source.records()[0]
        assert np.sign(Q(row["material"]["agreed"]["blank_change"])) == sign
        _, _, _, _, codes, h, *_ = fibers.source.tidal.lab.energy.pair.setup()
        ordinary = np.vstack(
            (
                w[:, :4].T,
                h.diagonal()[codes],
                fibers.source.tidal.lab.energy.returned(w)[-1],
            )
        )
        assert np.array_equal(ordinary @ arms["agreed"], ordinary @ original)
        protected = np.isin(w[:, -1], fibers.source.PARTIAL)
        assert np.array_equal(arms["agreed"][protected], original[protected])
