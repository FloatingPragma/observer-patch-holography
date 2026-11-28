"""Independent closure, local thermal reduction and stopped-kernel checks."""

from fractions import Fraction as Q

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_cadence as cadence,
    sector_retention as retention,
    sector_spectrum as spectra,
    sector_thermal as thermal,
)


def test_native_gibbs_references_and_singletons():
    words, _, _, groups, ops = thermal.apparatus()
    h = ops[0].sum(0) / 4
    references = thermal.bath()
    raw = expm(-thermal.BETA * h)
    np.testing.assert_allclose(raw / raw.trace(), references["exchanged"], atol=1e-14)
    for g in groups:
        assert (
            abs(references["retained"][np.ix_(g, g)].trace() - len(g) / len(words))
            < 1e-14
        )
    for j in range(4):
        marginal = thermal.restrict(
            references["retained"], thermal.region_indices(j, j)
        )
        assert np.max(abs(marginal - np.diag(marginal.diagonal()))) < 1e-13


def test_four_observer_spectrum_reconstructs_one_variance():
    C, V, reads = spectra.reconstruction()
    assert list(C[:, 2]) == [0, 0, Q(8, 4225), -Q(8, 4225), 0, 0, 0, 0]
    assert all(x == 0 for x in C[:, :2])
    assert (reads[-1, :] * V).rank() == 1


def test_positive_return_resolvent_on_a_three_state_clock():
    den = retention.ReturnClock.den
    P = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]) * (den // 4)
    T = retention.ReturnClock(csr_matrix(P), [0, 2])
    result, error = T.advance(np.array([[0], [T.grid]]))
    truth = np.array([3 / 8, 5 / 8])
    assert np.all(abs(result[:, 0] / T.grid - truth) <= error)
    power, eta = cadence.mixing(T, 2)
    assert power == 256 and eta < 1e-4


def test_interval_clock_polynomial_keeps_signed_terms():
    lo, hi = cadence.evaluate(0.25, 0.5, [1, -2, 3], [1.1, -1.9, 3.1], 0.1)
    for q in np.linspace(0.25, 0.5, 17):
        assert lo <= 1 - 2 * q + 3 * q * q - 0.1 * q**3 <= hi
        assert lo <= 1.1 - 1.9 * q + 3.1 * q * q + 0.1 * q**3 <= hi
