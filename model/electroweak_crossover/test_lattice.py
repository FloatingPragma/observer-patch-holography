"""Fast gates of the lattice code (numpy implementation, small lattices, fixed seeds)."""
import math

import numpy as np

from . import lattice as lt


def _binned_error(a, nb=10):
    n = len(a) // nb
    return float(np.asarray(a)[:n * nb].reshape(nb, n).mean(1).std(ddof=1) / math.sqrt(nb))


def test_relations_invert():
    for beta_G in (5.0, 8.0, 12.0, 20.0):
        for x, y in ((0.036, 0.02), (0.13, -0.04), (0.296, -0.09)):
            beta_H, beta_R = lt.lattice_couplings(x, y, beta_G)
            xx, yy = lt.continuum_parameters(beta_G, beta_H, beta_R)
            assert abs(xx - x) < 1e-12 and abs(yy - y) < 1e-9
    assert abs(lt.CONDENSATE_CONSTANT - 0.66796) < 2e-5


def test_relations_reproduce_published_critical_temperatures():
    # hep-lat/9510020, table 4: infinite-volume critical couplings and the temperatures quoted for them
    table = ((60, 5, 0.358495, 137.534), (60, 8, 0.3479735, 137.669), (60, 12, 0.3426840, 137.842),
             (60, 20, 0.3387418, 138.019), (35, 8, 0.3450806, 94.181), (70, 8, 0.3491523, 153.620))
    for m_h, beta_G, beta_H, t_c in table:
        h2 = (m_h / 80.6) ** 2
        x = -0.00550 + 0.12622 * h2
        _, y = lt.continuum_parameters(beta_G, beta_H, x * beta_H ** 2 / beta_G)
        t = math.sqrt(2.58088 * m_h ** 2 / (0.39818 + 0.15545 * h2 - 0.00190 * h2 ** 2 - y))
        assert abs(t - t_c) < 0.03


def test_gauge_invariance_of_every_observable():
    rng = np.random.default_rng(11)
    beta_H, beta_R = lt.lattice_couplings(0.296, -0.09, 8.0)
    lat = lt.Lattice(4, 8.0, beta_H, beta_R, rng, start="hot", radius=1.5)
    for _ in range(3):
        lat.sweep(1)
    before = list(lat.higgs_sums()) + [lat.plaquette(), lat.action()]
    lat.gauge_transform(lt.random_su2((4, 4, 4), rng))
    after = list(lat.higgs_sums()) + [lat.plaquette(), lat.action()]
    assert max(abs(a - b) / (1.0 + abs(a)) for a, b in zip(before, after)) < 1e-12


def test_overrelaxation_preserves_the_action():
    rng = np.random.default_rng(12)
    beta_H, beta_R = lt.lattice_couplings(0.296, -0.09, 8.0)
    lat = lt.Lattice(4, 8.0, beta_H, beta_R, rng, start="hot", radius=1.5)
    for _ in range(3):
        lat.sweep(1)
    s0 = lat.action()
    lat.overrelax_links()
    s1 = lat.action()
    assert abs(s1 - s0) < 1e-8 * abs(s0)
    # the Higgs move keeps the action where the step along the force is accepted; rejected steps are rare
    lat.accept_overrelax = [0, 0]
    lat.overrelax_higgs()
    assert lat.accept_overrelax[0] >= 0.9 * lat.accept_overrelax[1]


def test_free_field_limit():
    rng = np.random.default_rng(13)
    lat = lt.Lattice(4, 1.0, 0.30, 0.0, rng)
    lat.update_links = False
    for _ in range(50):
        lat.sweep(1)
    ser = lat.run(1500, 1)
    r2 = ser[:, 1] / lat.V
    hop = ser[:, 0] / lat.V
    exact_r2, exact_hop = lt.free_field_R2(4, 0.30)
    assert abs(r2.mean() - exact_r2) < 5 * _binned_error(r2) + 1e-3
    assert abs(hop.mean() - exact_hop) < 5 * _binned_error(hop) + 1e-3


def test_plaquette_at_zero_hopping_against_strong_coupling_series():
    rng = np.random.default_rng(14)
    lat = lt.Lattice(4, 1.0, 0.0, 0.0, rng, start="hot")
    values = []
    for i in range(700):
        lat.heatbath_links()
        if i >= 100:
            values.append(lat.plaquette())
    values = np.array(values)
    assert abs(values.mean() - lt.strong_coupling_plaquette(1.0)) < 5 * _binned_error(values) + 2e-4


def test_schwinger_dyson_identity():
    # invariance of the measure under phi -> (1 + eps) phi:
    # <R^2> - beta_H <H>/V + 2 beta_R <R^4 - R^2> = 2
    rng = np.random.default_rng(15)
    beta_H, beta_R = lt.lattice_couplings(0.296, -0.09, 8.0)
    lat = lt.Lattice(4, 8.0, beta_H, beta_R, rng)
    for _ in range(50):
        lat.sweep(2)
    ser = lat.run(600, 2)
    sd = ser[:, 1] / lat.V - beta_H * ser[:, 0] / lat.V + 2 * beta_R * (ser[:, 2] - ser[:, 1]) / lat.V
    assert abs(sd.mean() - 2.0) < 5 * _binned_error(sd) + 2e-3
    assert 0.4 < lat.accept_global[0] / lat.accept_global[1] < 0.95
