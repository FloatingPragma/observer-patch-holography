from itertools import product
from math import atanh

import numpy as np
import pytest
from scipy.stats import binom

from .response import curvature_bound, latent_extremum, majority, majority_stationarity, native_input


def test_native_contrast_and_input_chord():
    width, z = native_input()
    assert width > 1/8
    assert z[8] == pytest.approx(-0.014370365536930434, abs=2e-14)
    assert abs(z[0]) < 1e-14


def test_score_curvature_including_extremal_bounded_outputs():
    for d in (0.125, 1.0):
        for m in (1, 2, 3, 8, 32):
            for fraction in (-0.5, 0, 0.625):
                x = d * fraction
                sharp, second = latent_extremum(m, d, x)
                bound = curvature_bound(m, d, x)
                assert second == pytest.approx(bound**2, rel=1e-12, abs=1e-11)
                assert sharp <= bound + 1e-12
    assert latent_extremum(2, 1, 0)[0] == pytest.approx(2)


def test_majority_readout_with_independent_outcomes():
    n, y = 5, -0.18
    mu, slope = majority(n, y)
    expected = sum(np.sign(sum(word)) * np.prod([(1+b*y)/2 for b in word])
                   for word in product((-1, 1), repeat=n))
    assert mu == pytest.approx(expected)
    def from_tails(u):
        return 2*binom.sf(n//2, n, (1+u)/2) - 1
    assert slope == pytest.approx((from_tails(y+1e-6)-from_tails(y-1e-6))/2e-6, rel=1e-9)
    assert majority(1, y) == pytest.approx((y, 1))


def test_two_state_stationarity_and_single_copy_floor():
    z = native_input()[1][8]
    previous = None
    for inverse in (8, 16, 32, 64):
        low, high = [majority_stationarity(inverse, z, s) for s in (.25,.75)]
        k = low["f"]*z
        assert low["entropy"] - high["entropy"] >= 0.5*k*k
        area_constant = -(low["entropy"] + high["entropy"])/2
        errors = [abs(r["entropy"] + area_constant) for r in (low, high)]
        assert min(errors) == pytest.approx((low["entropy"]-high["entropy"])/2)
        # No single-copy slope can meet the two endpoints better than this common offset.
        assert abs(k*(atanh(0.75*k)-atanh(0.25*k))) == pytest.approx(sum(errors))
        rate = max(abs(r["residual"]) for r in (low,high))
        if previous is not None:
            assert rate < 0.28*previous
        previous = rate
