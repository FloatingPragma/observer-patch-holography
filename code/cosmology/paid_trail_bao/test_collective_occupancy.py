from fractions import Fraction as Q
from itertools import product

import numpy as np

from . import collective_occupancy as sky
from . import pooled_ordinary as pooled
from . import pooled_fluctuations as fluctuations


def test_count_law_dependence_envelope_and_finite_retention():
    p = np.array([0.01, 0.25, 0.75, 0.03])
    exact = np.zeros(4)
    for bit in product((0, 1), repeat=len(p)):
        prob = np.prod(np.where(bit, p, 1 - p))
        exact[min(sum(bit), 3)] += prob
    np.testing.assert_allclose(sky.count_law(p, 3), exact)
    assert sky.upper_dependent(p, 2) >= exact[2:].sum()
    assert sky.upper_dependent(np.array([0.95, 0.01, 0.01]), 2) == 0.02
    for t in (1, 8, 32):
        per = 1 - (1 - p) ** t
        np.testing.assert_allclose(
            sky.retained_pair(p, t), sky.count_law(per, 2)[-1], atol=2e-15
        )


def test_fluctuation_moments_and_independent_sphere_window():
    from math import comb

    n, w = 4, 0.125
    coeff = (0.03, -0.014, 0.211)
    active, mean, var = fluctuations.count_variance(n, w, *coeff)
    count = np.arange(n + 1)
    probability = np.array([comb(n, int(k)) * w**k * (1 - w) ** (n - k) for k in count])
    x = (48 * n * coeff[0] + coeff[2]) * (count > 0) + coeff[1] * count
    np.testing.assert_allclose(
        [active, mean, var],
        [1 - probability[0], probability @ x, probability @ (x - mean) ** 2],
    )
    # Two independently drawn top-hat locations share a random grid
    # cube with probability prod(max(1-|delta_x|/L,0)) across dimensions.
    rng = np.random.default_rng(198747)

    def ball():
        unit = rng.normal(size=(300000, 3))
        return (
            unit
            / np.linalg.norm(unit, axis=1)[:, None]
            * rng.random((len(unit), 1)) ** (1 / 3)
        )

    d = abs(ball() - ball())
    direct = np.maximum(0, 1 - d / 1.4).prod(1)
    theoretical = fluctuations.window_factor(1.4, 1, 32)
    assert abs(theoretical - direct.mean()) < 5 * direct.std() / np.sqrt(len(direct))
    assert theoretical >= 1.4**3 / (4 * np.pi / 3) * fluctuations.overlap(
        1.4 / np.sqrt(2)
    )


def test_empty_pools_and_explicit_enumeration_predict_uniform_ordinary_mass():
    p = np.array([0.25, 0.0, 0.125, 0.75])
    coeff = [float(Q(38, 1155)), float(Q(-1780, 124509)), float(Q(131792, 622545))]
    expected = 0.0
    for bits in product((0, 1), repeat=len(p)):
        m = sum(bits)
        if m:
            expected += np.prod(np.where(bits, p, 1 - p)) * (
                coeff[0] * 48 * len(p) + coeff[1] * m + coeff[2]
            )
    np.testing.assert_allclose(
        pooled.expected_per_tile(p, np.zeros(len(p), int), coeff), expected / len(p)
    )
    assert pooled.expected_per_tile(p, np.arange(len(p)), coeff)[1] == 0
    for n in (1, 5, 20):
        fill = 0.23
        w = 6 * fill**2 * (1 - fill) ** 2 * 448 / 2110
        n_eq = pooled.expected_per_tile(np.full(n, w), np.zeros(n, int), coeff)
        np.testing.assert_allclose(
            pooled.homogeneous(n, fill, coeff, 0.25), 0.25 * n_eq / (4 * fill)
        )
