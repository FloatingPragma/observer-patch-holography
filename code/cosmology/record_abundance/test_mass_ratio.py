from itertools import product

import numpy as np
import pytest

from cosmology.record_abundance import mass_ratio as m


def check_boundary(ledger):
    k = ledger["ordinary"]
    owned = ledger["owned"]
    twice = ledger["symmetric_twice"]
    incoming = ledger["active_outside_crossings"]
    outgoing = ledger["active_inside_crossings"]
    assert owned <= 3 * k
    assert twice == 2 * owned + incoming - outgoing
    assert twice <= 6 * k + ledger["boundary_edges"]
    assert twice <= 6 * k + 3 * ledger["external_occupied_neighbors"]
    if not ledger["external_occupied_neighbors"]:
        assert twice <= 6 * k


def test_exhaustive_degree_three_boundary_and_sharpness():
    # A complete four-vertex graph is cubic. All 256 preparation/aperture pairs
    # include a lone interior mark, an exterior mark, and the full/empty aperture.
    n = 4
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    sets = [
        {i for i, on in enumerate(word) if on} for word in product((0, 1), repeat=n)
    ]
    for occ, region in product(sets, repeat=2):
        check_boundary(m.aperture(edges, occ, region))
    assert m.aperture(edges, {0}, set(range(n)))["owned"] == 3
    external = m.aperture(edges, {0}, {1, 2, 3})
    assert external["symmetric_twice"] == 3 and external["ordinary"] == 0


def test_native_board_and_unequal_volume_assignment():
    rng = np.random.default_rng(1461)
    board = m.stock.measure.federation(2)
    edges = np.asarray(board.external) // 12
    assert np.all(np.bincount(edges.ravel()) == 3)
    sizes = rng.lognormal(0, 2.5, board.cells)
    for _ in range(18):
        occ = set(np.flatnonzero(rng.random(board.cells) < 0.16))
        region = set(np.flatnonzero(rng.random(board.cells) < 0.5))
        ledger = m.aperture(edges, occ, region)
        check_boundary(ledger)
        # Integrate local counts as densities through arbitrarily unequal volumes.
        density = np.zeros(board.cells)
        for x, y in edges:
            if (x in occ) != (y in occ):
                owner = x if x in occ else y
                density[owner] += 1 / sizes[owner]
        members = list(region)
        assert np.dot(density[members], sizes[members]) == pytest.approx(
            ledger["owned"]
        )


def test_interior_exterior_and_background_thresholds():
    result = m.run()
    row = result["ngc3741"]
    known, dyn, cap = m.GAS + m.STAR, m.DYN, result["dark_to_baryons_cap"]
    assert (1 + cap) * (
        known + row["minimum_extra_interior_baryons_msun"]
    ) == pytest.approx(dyn)
    exterior = row["minimum_one_link_exterior_baryons_if_no_hidden_interior_msun"]
    assert known + cap * (known + exterior / 2) == pytest.approx(dyn)
    scale = result["required_distance_if_no_unseen_mass_mpc"] / m.DISTANCE
    assert (1 + cap) * known * scale**2 == pytest.approx(dyn * scale)
    ref = result["uniform_reference_baryons_msun_inside_upper_scale"]
    # Test signed excess form directly from a raw mass at the sharp cap.
    assert cap * known - result["external_2018_cold_to_baryons"] * ref == pytest.approx(
        cap * (known - ref) + result["maximum_background_subtraction_correction_msun"]
    )
    assert result["planck2013_clusters"]["below_at_upper_two_quoted_errors"] == [
        "2955",
        "S0296",
    ]


def test_recombination_and_dilute_surface_brightness():
    diagnostic = m.warm_gas(m.run())
    assert 7.35e38 < diagnostic["minimum_luminosity_erg_per_s"] < 7.37e38
    assert (
        7.95e-19
        < diagnostic["minimum_mean_surface_brightness_erg_s_cm2_arcsec2"]
        < 8.0e-19
    )
    assert (
        2.80e8
        < diagnostic["warm_gas_upper_msun_if_reported_light_is_complete"]
        < 2.81e8
    )
