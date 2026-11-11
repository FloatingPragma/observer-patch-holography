"""The stored receipt equals its recomputation, and the named checks pass."""
import json
import math
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from . import complexes as C
from . import exact as X
from . import fits as ft
from .certificate import certificate


@pytest.fixture(scope='module')
def value():
    """The receipt, computed once and only when a test of this module runs."""
    return certificate()


@pytest.fixture(scope='module')
def stored():
    return json.loads(Path(__file__).with_name('results.json').read_text())


def test_every_named_check_passes(value):
    failed = [name for name, ok in value['checks'].items() if not ok]
    assert not failed and value['all_checks_pass']
    assert {name.split('.')[0] for name in value['checks']} == {
        '1_one_register', '2_cartan_circle', '3_flux_densities', '4_wilson_classical', '5_second_reading',
        '6_estimator', '7_gates', '8_law', '9_sector_law', '10_doublet', '11_second_reading'}


def test_stored_receipt_matches_recomputation(value, stored):
    assert stored['checks'] == value['checks']
    assert stored['measured_points_sha256'] == value['measured_points_sha256']
    for k, e in value['scans'].items():
        assert abs(stored['scans'][k]['fit_p1']['c'] - e['fit_p1']['c']) < 1e-9


def test_one_link_law():
    tau = X.G2_U / 2
    assert abs(float(X.cartan_law(tau) / mp.exp(-X.S_U)) - 2) < 1e-12
    assert abs(float(X.register_law(2, tau) / mp.exp(-X.S_U)) - 2 * (2 * X.S_U - 1)) < 1e-9
    assert abs(X.S_U - (math.pi / 2) * X.X_U) < 1e-12


def test_weight_powers_one_and_zero_have_no_exponential():
    tau = X.G2_U / 2
    assert 0.03 < float(X.register_law(1, tau)) < 0.04
    assert 0.16 < float(X.register_law(0, tau)) < 0.18
    for c in (0.5, 1 / 3):
        assert float(X.register_law(1, tau / c)) > 0.05


def test_centre_ratio_of_the_heat_kernel():
    for g2 in (0.4, 0.5168, 1.0):
        S = 2 * math.pi ** 2 / g2
        r = math.exp(X.ln_heat_kernel(2 * math.pi, g2) - X.ln_heat_kernel(0.0, g2))
        assert abs(r / (2 * (2 * S - 1) * math.exp(-S)) - 1) < 1e-9


def test_exact_projectors(value):
    V, E, Fc = C.icosahedron()
    assert (len(V), len(E), len(Fc)) == (12, 30, 20)
    t = value['cartan_circle']
    assert t['icosahedral_surface_euclidean']['c_exact'] == '1/20'
    assert t['icosahedral_space_one_tick_gauss_law']['c_exact'] == '19/30'
    assert t['face_pair_euclidean_pinned_boundary']['c_exact'] == '1/2'
    assert t['cubic_d3_L3']['c_exact'] == str(1 - F(2, 3) * (1 - F(1, 27)))
    assert t['cubic_d4_L2']['c_exact'] == str(1 - F(1, 2) * (1 - F(1, 16)))
    assert abs(t['icosahedral_space_many_ticks']['c'] - 0.398231) < 1e-6
    assert abs(t['face_pair_space_many_ticks']['c'] - 1 / (2 * math.sqrt(2))) < 1e-6


def test_count_table(value):
    t = value['cartan_circle']
    assert abs(t['cubic_d3_infinite']['count_4_over_c'] - 12) < 1e-12
    assert abs(t['cubic_d4_infinite']['count_4_over_c'] - 8) < 1e-12
    assert abs(t['pinned_link']['count_4_over_c'] - 4) < 1e-12
    assert abs(t['icosahedral_space_many_ticks']['count_4_over_c'] - 10.0444) < 1e-3
    d = value['flux_densities']['cubic']
    assert [d[k]['count_density'] for k in ('2', '3', '4')] == [4.0, 8 / 3, 2.0]


def test_wilson_minimum(value):
    w = value['wilson_classical']
    assert abs(w['cub4_L4']['s_W'] - 8) < 1e-7 and w['cub4_L4']['flux_on_the_twisted_plaquette_over_2pi'] < 1e-4
    assert abs(w['pinned_link']['c'] - 4 / math.pi ** 2) < 1e-15


def test_measured_laws(value):
    s = value['scans']
    assert set(s) == {'cub4_L6_hk', 'cub4_L6_w', 'cub3_L8_hk', 'cub3_L8_w', 'icoT_NT8_hk', 'icoT_NT8_w'}
    for k, e in s.items():
        assert abs(e['fit_p1']['c'] - 1) > 0.3            # no lattice ensemble has the coefficient of the one-link law
        assert abs(e['fit_p1']['c'] / e['c_predicted'] - 1) < 0.05
        if e['action'] == 'heat kernel':
            assert abs(e['c_predicted'] - e['c_gaussian_twist']) < 1e-12
    assert abs(s['cub4_L6_w']['c_predicted'] - 4 / math.pi ** 2) < 1e-12
    assert abs(s['cub4_L6_w']['fit_p_free']['c'] / (4 / math.pi ** 2) - 1) < 0.03


def test_doublet_and_second_reading(value):
    ch = value['doublet_change_across_the_line']
    assert ch['heat kernel']['difference'] < -3 and ch['Wilson']['difference'] < -1
    one = value['one_link_law_lnA_at_chain_coupling']
    assert abs(one['heat kernel'] - (math.log(2 * (2 * X.S_U - 1)) - X.S_U)) < 1e-12
    sec = value['second_reading']
    assert sec['cub4_L4']['c'] > 1.2 and 0.7 < sec['icoT_NT8']['c'] < 0.75
    clock = sec['clock']
    assert abs(clock['cubic space, 3 dimensions']['2']['face_step'] - 0.944) < 2e-3
    assert abs(clock['icosahedral space']['continuous']['face_step'] - 0.540506) < 1e-5


def test_fit_function_recovers_a_known_law():
    S = np.array([20.0, 25.0, 30.0, 38.0, 44.0])
    y = math.log(3.0) + np.log(S) - 0.4 * S
    f = ft.fit_scan(S, y, np.full(5, 0.01), 1.0)
    assert abs(f['c'] - 0.4) < 1e-10 and abs(f['ln_M0'] - math.log(3.0)) < 1e-9
    g = ft.fit_scan(S, y, np.full(5, 0.01), None)
    assert abs(g['p'] - 1) < 1e-8


def test_sector_law_fit_recovers_weight_power():
    n = np.arange(1, 25)
    tau = 0.3
    p = n ** 2 * np.exp(-0.5 * (n ** 2 - 1) / 4)
    p /= p.sum()
    s = ft.sector_law(p, p * 1e-3, tau)
    assert abs(s['w'] - 2) < 1e-8 and abs(s['tau_over_t'] - 0.6) < 1e-8
