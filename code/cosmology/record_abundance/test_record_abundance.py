"""Check the sharp cut bounds exhaustively and break the homogeneous-field control."""

from itertools import combinations, product

import numpy as np

from cosmology.record_abundance import measure
from cosmology.record_abundance import stationary


def test_pair_transfer_matches_every_native_swap_word():
    initial = np.array([1, 1, 0, 0])
    matchings = [([0, 2], [1, 3]), ([0, 1], [3, 2]), ([0, 2], [1, 3])]
    mean = initial.astype(float)
    second = np.outer(mean, mean)
    for x, y in matchings:
        stationary.exchange_moments(mean, second, x, y)
    joint = np.zeros_like(second)
    for choices in product((0, 1), repeat=6):
        field = initial.copy()
        for t, (x, y) in enumerate(matchings):
            for a, b, swap in zip(x, y, choices[2 * t : 2 * t + 2], strict=True):
                if swap:
                    field[a], field[b] = field[b], field[a]
        joint += np.outer(field, field) / 64
    np.testing.assert_array_equal(second, joint)


def test_stationary_cut_variance_by_direct_enumeration():
    owner = np.asarray(measure.federation(0).external) // 12
    n, k = 20, 3
    cuts = []
    for occupied in combinations(range(n), k):
        field = np.isin(np.arange(n), occupied)
        cuts.append(np.count_nonzero(field[owner[:, 0]] != field[owner[:, 1]]))
    expected = stationary.equilibrium(n, k)
    np.testing.assert_allclose(
        [np.mean(cuts), np.var(cuts)],
        [expected["mean"], expected["variance"]],
        atol=1e-13,
        rtol=0,
    )


def test_every_complementary_half_on_the_80_cell_carrier():
    owners, region, galaxy, banks = measure.policies(1)
    size = int(owners.max()) + 1

    def cut(x):
        return (x[owners[:, 0]] != x[owners[:, 1]]).sum()

    means = []
    for half in combinations(region, len(region) // 2):
        mask = np.isin(np.arange(size), half)
        inverse = np.isin(np.arange(size), region) & ~mask
        assert not np.any(mask[owners[(owners == galaxy).any(axis=1)]])
        means.append((cut(mask) + cut(inverse)) / 2)
    assert len(means) == 70
    assert min(means) == np.mean([cut(mask) for mask in banks["whole_face"]]) == 6
    assert max(means) == np.mean([cut(mask) for mask in banks["staggered"]]) == 12


def test_moving_an_equal_injection_budget_changes_gravity():
    visits, eigen = measure.placements()
    n, s = measure.WIDTH, measure.SURVIVAL
    coord = np.arange(n) - n // 2
    mode = np.cos(2 * np.pi * coord[:, None, None] / n) * np.ones_like(visits)
    uniform = np.full_like(visits, measure.UNIFORM_PLACEMENT)
    eigenvalue = (3 + np.cos(2 * np.pi / n)) / 4
    varied = uniform * (1 + 0.5 * (1 - s) * mode / (1 - s * eigenvalue))
    driving = varied - s * measure.array.periodic_step(varied)
    np.testing.assert_allclose(
        driving, (1 - s) * uniform * (1 + 0.5 * mode), atol=1e-16, rtol=0
    )
    np.testing.assert_allclose(driving.sum(), (1 - s) * uniform.sum())
    gain = -0.014
    halo = 3 * measure.charge(visits, gain)
    a, esd, _ = measure.gravity(halo + 6 * measure.charge(uniform, gain), eigen)
    b, esd_varied, _ = measure.gravity(halo + 6 * measure.charge(varied, gain), eigen)
    assert np.max(np.abs(a - b)) > 1e-8
    assert np.max(np.abs(esd - esd_varied)) > 1e-9
