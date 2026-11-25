import numpy as np
from scipy.stats import binom

from model.experiments.count_intervals import clopper_pearson


def test_large_trials_and_small_counts_against_forward_tails():
    for n, k, tau in [(2**28, 1000, 1e-7), (2**32, 999, 1e-15), (2**40, 1250000, 1e-6)]:
        lo, hi = clopper_pearson(k, n, tau)
        assert lo < k / n < hi
        for p in (binom.sf(k - 1, n, lo), binom.cdf(k, n, hi)):
            assert tau * 0.999 < p <= tau * (1 + 1e-12)


def test_endpoint_counts():
    lo, hi = clopper_pearson(np.array([0, 200]), 200, 0.001)
    assert lo[0] == 0 and hi[-1] == 1
    np.testing.assert_allclose(hi[0], -np.expm1(np.log(0.001) / 200), rtol=2e-8)
    np.testing.assert_allclose(lo[1], np.exp(np.log(0.001) / 200), rtol=2e-8)
