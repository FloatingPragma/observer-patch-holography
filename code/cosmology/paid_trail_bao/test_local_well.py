import numpy as np

from . import local_well as well


def test_density_matches_enclosure_derivative_and_interpolation():
    edges = np.array([0, 1.5, 3.0])
    weights = np.array([0.2, 0.8])
    shape = np.array([0, 2]), np.array([0.3, 0.7])
    sizes = np.linspace(0.01, 3, 151)
    r = 1.7
    spline = well.density_times_length(r, sizes, shape, (edges, weights))
    direct = np.concatenate(
        [
            well.density_times_length(r, group, shape, (edges, weights))
            for group in np.array_split(sizes, 3)
        ]
    )
    np.testing.assert_allclose(spline, direct, rtol=1e-13, atol=1e-15)

    def enclosed(x):
        zero = np.clip((x * x - edges[:-1] ** 2) / np.diff(edges**2), 0, 1)
        moving = well.profiles.annulus_fraction(
            x, 2 * sizes[:, None], edges[:-1], edges[1:]
        )
        return 0.3 * (zero @ weights) + 0.7 * (moving @ weights)

    eps = 1e-6
    estimate = (enclosed(r + eps) - enclosed(r - eps)) / (2 * eps * 4 * np.pi * r * r)
    np.testing.assert_allclose(spline / sizes, estimate, atol=1e-9)


def test_continuous_search_accepts_a_constructed_off_grid_curve():
    shape = np.array([0, 1]), np.array([0.25, 0.75])
    annuli = np.array([0, 1, 2, 3]), np.array([0.7, 0.2, 0.1])
    r = np.array([0.3, 0.5, 0.8, 1.2, 1.6, 2.2, 2.4])
    length = np.array([0.80637])
    field = np.array(
        [well.field_times_length(x, length, shape, annuli)[0] / length[0] for x in r]
    )
    obs = np.sqrt(100**2 + 2000 * field)
    rows = dict(
        radius_kpc=r,
        baryons_kms=np.ones_like(r) * 100,
        observed_kms=obs,
        error_kms=np.ones_like(r) * 1e-4,
    )
    test = well.shape_test(rows, annuli, shape)
    assert len(test["ring_indices"]) >= 3 and not test["excluded"]


def test_exact_amplitude_elimination_including_the_zero_profile():
    from fractions import Fraction as Q
    from .well_certificate import segment_feasible

    lo, hi = [Q(1), Q(2)], [Q(1), Q(2)]
    assert segment_feasible([Q(1), Q(0)], [Q(0), Q(1)], lo, hi)
    assert segment_feasible([Q(0), Q(0)], [Q(1), Q(2)], lo, hi)
    assert not segment_feasible([Q(0), Q(0)], [Q(-1), Q(-2)], lo, hi)
    assert not segment_feasible([Q(1), Q(1)], [Q(2), Q(3)], lo, hi)
