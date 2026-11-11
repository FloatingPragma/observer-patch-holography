import math

import numpy as np
import pytest

from . import centre, scale, su3
from .certificate import BETA_C_LITERATURE, certificate, result

RESULT = result()


def test_every_named_check_passes():
    value = certificate()
    failed = [name for name, ok in value['checks'].items() if not ok]
    assert not failed and value['all_checks_pass']
    assert {name.split('.')[0] for name in value['checks']} == {
        '1_action', '2_centre', '3_lattice', '4_strong_coupling', '5_kernel', '6_tables', '7_scale'}
    stored = dict(RESULT)
    fresh = dict(value)
    stored.pop('replay', None)
    fresh.pop('replay', None)
    assert stored['checks'] == fresh['checks'] and stored['scale'] == fresh['scale']
    assert stored['input_sha256'] == fresh['input_sha256']


def test_colour_factor_of_the_plaquette_action():
    rng = np.random.default_rng(3)
    for V in su3.random_su3(20, rng):
        assert abs(centre.five_trace_cost(V3=V) - 0.6 * (1 - np.trace(V).real / 3)) < 1e-14
    assert centre.thin_flux_cost()['2'] == centre.thin_flux_cost()['4'] == '9/10'
    assert centre.thin_flux_cost()['3'] == '4/5'


def test_centre_classes():
    ft = centre.flux_table()
    assert ft['2']['classes_rotated'] == ft['4']['classes_rotated'] == [1, 2, 4, 5]
    assert ft['3']['classes_rotated'] == [1, 3, 5]
    assert ft['2']['colour_phase_turns'] == '1/3' and ft['4']['colour_phase_turns'] == '2/3'
    assert ft['2']['weak_phase_turns'] == ft['4']['weak_phase_turns'] == '0'
    assert ft['3']['colour_phase_turns'] == '0' and ft['3']['weak_phase_turns'] == '1/2'


def test_centre_transformation_on_the_lattice():
    for value in centre.lattice_checks(seed=21).values():
        assert value < 1e-12
    v = centre.vortex_checks()
    assert v['flux_plaquettes'] == v['region_surface_links'] == 52
    assert v['cube_closure_violation'] == 0 and v['spatial_flux_plaquettes'] == 0


def test_gauge_invariance_and_overrelaxation():
    lat = su3.Lattice((4, 4, 4, 2), 5.5, seed=4, start='hot')
    lat.compound_sweeps(2, 1)
    p, l = lat.plaquette(), lat.polyakov()
    lat.gauge_transform(su3.random_su3(lat.volume, lat.rng).reshape(lat.dims + (3, 3)))
    assert abs(lat.plaquette() - p) < 1e-13 and abs(lat.polyakov() - l) < 1e-13
    lat.sweep(False)
    assert abs(lat.plaquette() - p) < 1e-13
    assert max(lat.unitarity_defect()) < 1e-12


def test_x0_generator():
    rng = np.random.default_rng(8)
    q = np.linspace(-1, 1, 200001)
    for a in (0.3, 2.0, 8.0):
        x = su3.draw_x0(np.full(40000, a), rng)
        w = np.sqrt(1 - q * q) * np.exp(a * q)
        assert abs(x.mean() - (w * q).sum() / w.sum()) < 5 * x.std() / math.sqrt(len(x))


def test_strong_coupling_series():
    assert [round(su3.haar_moment_re_trace(k), 12) for k in range(1, 6)] == [0.0, 0.5, 0.25, 0.75, 0.9375]
    for beta in (0.2, 0.5):
        series = beta / 18 + beta ** 2 / 216 - 5 * beta ** 4 / 93312
        assert abs(su3.one_plaquette_average(beta) - series) < 2 * beta ** 5 / 1e4


def test_numba_kernel_agrees_with_numpy_kernel():
    nb = pytest.importorskip('numba')
    from . import su3_numba as K
    dims = (4, 4, 4, 4)
    up, dn = K.neighbours(*dims)
    K.seed_rng(5)
    U = K.hot_start(256)
    K.compound_sweeps(U, up, dn, 5.7, 3, 2)
    lat = su3.Lattice(dims, 5.7)
    # site index s = i0 + 4 (i1 + 4 (i2 + 4 i3)): reshape in Fortran order of the site axes
    lat.U = np.transpose(U.reshape(4, 4, 4, 4, 4, 3, 3), (4, 3, 2, 1, 0, 5, 6)).copy()
    ps, pt = K.plaquette(U, up)
    assert abs(lat.plaquette() - 0.5 * (ps + pt)) < 1e-13
    assert abs(lat.polyakov() - K.polyakov(U, up, 4)) < 1e-13


def test_recorded_transition():
    lat = RESULT['lattice']
    assert lat['N_t'] == 4 and [s['N_s'] for s in lat['sizes']] == [8, 12, 16]
    assert abs(lat['beta_c'] - BETA_C_LITERATURE) < 1e-3 * BETA_C_LITERATURE
    chi = [s['chi_max'] for s in lat['sizes']]
    assert chi[0] < chi[1] < chi[2]
    assert RESULT['register']['Z3']['order'] == 'first' and RESULT['register']['Z2']['order'] == 'continuous'


def test_scale():
    t = scale.table()
    assert abs(t['primary']['ratio'] - 0.46302) < 1e-5
    assert abs(t['primary']['T_MeV'] - 157.426) < 1e-3
    assert abs(t['t_seconds_g_above'] - 1.2428e-5) < 1e-9 and abs(t['t_seconds_g_below'] - 2.3513e-5) < 1e-9
    assert abs(t['ticks_g_above'] / 1.805e38 - 1) < 1e-3
    assert scale.g_star(18, 50) == 61.75 and scale.g_star(5, 14) == 17.25


def test_register_ensemble_kernel():
    from . import register as rg
    assert abs(rg.dual_critical_coupling(2) - 0.761413) < 1e-6
    assert abs(rg.dual_critical_coupling(3) - 1.0843) < 1e-4
    for N, kappa in ((2, 0.15), (3, 0.2)):
        reg = rg.Register(6, N, kappa, seed=N, start='hot')
        for _ in range(20):
            reg.sweep()
        vals = []
        for _ in range(200):
            reg.sweep()
            vals.append(reg.mean_plaquette())
        err = np.std(vals, ddof=1) / math.sqrt(len(vals))
        assert abs(np.mean(vals) - rg.one_plaquette_average(N, kappa)) < 5 * err + 2e-3
    reg = rg.Register(6, 3, 1.0, seed=9, start='hot')
    for _ in range(5):
        reg.sweep()
    p = reg.plaquettes().copy()
    reg.gauge_transform(reg.rng.integers(0, 3, (6, 6, 6)))
    assert (reg.plaquettes() == p).all()


def test_series_runs_with_the_numpy_kernel():
    P, L, kernel = su3.series((4, 4, 4, 2), 5.0, 6, n_or=1, seed=1, start='hot', kernel='numpy')
    assert kernel == 'numpy' and len(P) == len(L) == 6
    assert 0.2 < P[-1] < 0.6 and abs(L[-1]) < 1.0
