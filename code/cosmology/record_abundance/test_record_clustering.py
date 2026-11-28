"""Direct enumeration of spatial cross-correlations and operator capacity."""

from itertools import product

import numpy as np

from cosmology.record_abundance import clustering


def test_all_256_patterns_on_the_cube():
    n = 8
    edges = np.asarray(
        [(i, i ^ bit) for i in range(n) for bit in (1, 2, 4) if i < (i ^ bit)]
    )
    a, b = edges.T
    incidence = (np.eye(n)[a] + np.eye(n)[b]) / 2
    configs = np.array(list(product((0, 1), repeat=n)))
    flags = configs[:, a] != configs[:, b]
    scores = flags @ incidence
    assert scores.max() == 1.5
    # A snapshot with every ordinary port filled has no record differences.
    assert np.all(scores[-1] == 0)
    for p in (0.16, 0.5, 0.84):
        weights = np.prod(np.where(configs, p, 1 - p), axis=1)
        np.testing.assert_allclose(weights @ scores, 3 * p * (1 - p), atol=1e-15)
        centered = configs - p
        cov = scores.T @ (weights[:, None] * centered)
        np.testing.assert_allclose(
            cov, (2 * p * (1 - p) * (1 - 2 * p)) * (incidence.T @ incidence), atol=5e-16
        )
        pairs = flags.T @ (weights[:, None] * flags) - (2 * p * (1 - p)) ** 2
        overlap = 4 * incidence @ incidence.T
        np.testing.assert_allclose(pairs[overlap == 0], 0, atol=2e-16)
        np.testing.assert_allclose(
            pairs[overlap == 1], p * (1 - p) * (1 - 2 * p) ** 2, atol=2e-16
        )
        np.testing.assert_allclose(
            np.max(scores / (weights @ scores)),
            clustering.cluster_caps(p)["dark_density_cap"],
            rtol=1e-13,
        )
