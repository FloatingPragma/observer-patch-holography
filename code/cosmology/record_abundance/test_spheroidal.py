from itertools import combinations
import json
from math import pi, radians

import numpy as np
import pytest

from cosmology.record_abundance import spheroidal as s


def test_two_occupied_all_placements_on_native_80_board():
    board = s.stationary.measure.federation(1)
    edges = np.asarray(board.external) // 12
    adjacency = [set() for _ in range(board.cells)]
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    # With precisely two marks the second can block a directed cut whenever
    # it occupies either endpoint or any node reached in one or two steps.
    blocked = {}
    for a, b in edges:
        neighborhood = {a, b}
        for _ in range(2):
            neighborhood |= set().union(*(adjacency[i] for i in neighborhood))
        blocked[a, b] = blocked[b, a] = neighborhood
    tally = [
        sum(y not in blocked[x, z] for x, y in ((i, j), (j, i)) for z in adjacency[x])
        for i, j in combinations(range(board.cells), 2)
    ]
    assert np.mean(tally) / 2 == pytest.approx(
        s.native_census(1, 2)["expected_failed_mismatches_per_ordinary"]
    )


def test_actual_sky_coverage_against_independent_points():
    saved = json.loads(s.RECEIPT.read_text())
    for key in ("half_light", "outer"):
        radius = saved[key]["radius_pc"]
        angle = np.radians(np.arange(360))
        offsets = np.degrees(np.arcsin(radius / s.DIST)) * np.array(
            [np.cos(angle), np.sin(angle)]
        )
        lon = s.CENTER[0] + offsets[0] / np.cos(np.radians(s.CENTER[1]))
        lat = s.CENTER[1] + offsets[1]
        beams = np.array(saved[key]["neutral"]["beams_lon_lat_deg"])
        hav = (
            np.sin(np.radians(lat - beams[:, 1:2]) / 2) ** 2
            + np.cos(np.radians(lat))
            * np.cos(np.radians(beams[:, 1:2]))
            * np.sin(np.radians(lon - beams[:, 0:1]) / 2) ** 2
        )
        arc = 2 * np.degrees(np.arcsin(np.sqrt(hav)))
        response = np.exp(-4 * np.log(2) * (arc / s.FWHM) ** 2).sum(axis=0)
        assert response.min() >= saved[key]["neutral"]["minimum_beam_sum"]


def test_warm_beam_and_density_occupation_bounds():
    x = json.loads(s.RECEIPT.read_text())["half_light"]
    n = x["warm"]["warm_mass_ceiling_msun"] * s.MSUN / (1.33 * s.MP)
    volume = 4 * pi / 3 * (291 * s.PC) ** 3
    projected_area = pi * (s.DIST * s.PC * np.tan(radians(0.5))) ** 2
    rayleigh = n**2 / (volume * projected_area) / (2.5 * s.PC)
    assert 0.9 * rayleigh - 0.003 == pytest.approx(0.024)
    assert x["proper_site_density_enhancement_minimum_from_stars_at_ml_1"] > 27400
    p = x["limits"]["two_step"]["maximum_ordinary_occupation_at_ml_3"]
    p0 = s.stationary.screen.DENSITY

    def charge(p):
        return ((1 - p) ** 12 + (1 - p) ** 13) / 2

    ratio = 0.11933 / 0.02242 * charge(p) / charge(p0)
    ordinary = (
        3 * x["light_used_for_inventory_lsun"] + x["allowed_atomic_plus_warm_msun"]
    )
    assert ordinary * (1 + ratio) == pytest.approx(x["mass_proxy_msun"])
