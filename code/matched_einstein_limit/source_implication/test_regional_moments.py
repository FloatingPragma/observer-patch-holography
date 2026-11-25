"""Uniform-energy comparison with an independent finite entropy difference."""

import numpy as np
import pytest
from scipy.linalg import eigvalsh

from model.candidate_laws.interval_echo.candidate import native_law
from .regional_moments import conserved_weights_witness, embed, measure, pair_transfers


def entropy(rho):
    probabilities = eigvalsh(rho)
    assert probabilities[0] > 0
    return -probabilities @ np.log(probabilities)


@pytest.mark.parametrize(("contrast", "difference"), [(8, 0.14966922), (4, 0.07216476)])
def test_two_uniform_energy_preparations(contrast, difference):
    result = measure(contrast)
    assert result["internal_edges"] == 15
    assert result["discrepancy"] == pytest.approx(difference, abs=1e-8)
    np.testing.assert_allclose(result["total_variation"], result["slopes"], rtol=0, atol=1e-13)
    for tangent, prediction in zip(result["regional_tangents"], result["slopes"], strict=True):
        rho = result["reference"]
        step = 2**-17
        assert step < result["positive_step"]
        secant = (entropy(rho + step*tangent) - entropy(rho - step*tangent)) / (2*step)
        assert secant == pytest.approx(prediction, abs=1e-7)


def test_equal_profile_on_every_relative_seam_position():
    assert conserved_weights_witness() == -2
    ids = [native_law.SRC.index(word) for word in ((0, 0), (1, 0), (0, 1), (1, 1))]
    h = native_law.seam_term(8)[np.ix_(ids, ids)]
    for delta in pair_transfers(8):
        full = embed(delta, (0, 1), 4)
        for location, expected in (((0, 1), 1), ((0, 2), 0), ((2, 3), 0)):
            # The readout acts as identity on its spectators, independent of their state.
            observable = 4 * embed(h, location, 4)
            assert np.trace(full @ observable) == expected


def test_high_temperature_boundary_triangle_coefficient():
    beta = 1/4096
    result = measure(8, beta)
    assert result["high_temperature_coefficient"] == 25.5
    assert result["discrepancy"] / beta**2 == pytest.approx(25.5, abs=.08)
