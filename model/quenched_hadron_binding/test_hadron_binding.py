"""The stored receipt equals its recomputation from the recorded tables, and the named checks pass."""
import json
from pathlib import Path

import numpy as np
import pytest

from . import analysis as an
from . import lattice as lt
from .__main__ import COMPARED
from .certificate import HEADLINE, certificate, close, weighted_slope

@pytest.fixture(scope='module')
def result():
    """The receipt, computed once and only when a test of this module runs."""
    return json.loads(json.dumps(certificate()))


@pytest.fixture(scope='module')
def stored():
    return json.loads(Path(__file__).with_name('results.json').read_text())


def test_every_named_check_passes(result):
    failed = [name for name, ok in result['checks'].items() if not ok]
    assert not failed and result['all_checks_pass']
    sections = {name.split('.')[0] for name in result['checks']}
    assert sections == {'1_gamma', '2_gauge', '3_dirac', '4_hadrons', '5_potential', '5_gauge', '6_solver', '7_hadrons',
                        '8_chiral', '8_covariance', '9_volume', '10_melting', '11_recorded_gates'}


def test_stored_receipt_equals_recomputation(result, stored):
    assert stored['checks'] == result['checks']
    for key in COMPARED:
        assert not close(stored[key], result[key], path=key)


def test_headline_numbers(result):
    h = result['headline']
    assert h['lattice'] == [12, 12, 12, 24] and h['configurations'] == 128 and h['streams'] == 32
    assert abs(h['a2_sigma']['value'] - 0.1667) < 5e-5 and abs(h['r0_over_a']['value'] - 2.946) < 5e-4
    assert abs(h['kappa_c']['value'] - 0.16919) < 5e-6
    assert abs(h['mN_over_mrho_chiral']['value'] - 1.4495) < 5e-5 and abs(h['mN_over_mrho_chiral']['err'] - 0.132) < 5e-4
    assert abs(h['mN_over_mrho_phys']['value'] - 1.4544) < 5e-5 and abs(h['mN_over_mrho_phys']['err'] - 0.124) < 5e-4
    assert abs(h['mrho_over_sqrtsigma_chiral']['value'] - 1.3326) < 5e-5
    assert abs(h['mN_over_sqrtsigma_chiral']['value'] - 1.9317) < 5e-5
    assert h['nucleon_not_resolved_at'] == [0.1675]


def test_wilson_projectors_remove_backtracking_paths():
    g, g5 = lt.gammas()
    for mu in range(4):
        assert np.abs((np.eye(4) - g[mu]) @ (np.eye(4) + g[mu])).max() < 1e-14
        assert np.linalg.matrix_rank(np.eye(4) - g[mu]) == 2


def test_hopping_term_equals_its_definition():
    rng = np.random.default_rng(5)
    dims = (4, 4, 4, 6)
    uf = lt.fermion_links(lt.hot_links(dims, rng))
    psi = rng.normal(size=dims + (4, 3)) + 1j * rng.normal(size=dims + (4, 3))
    g, _ = lt.gammas()
    ref = np.zeros_like(psi)
    for mu in range(4):
        fwd = np.einsum('...ab,...sb->...sa', uf[mu], np.roll(psi, -1, axis=mu))
        ref += np.einsum('ts,...sa->...ta', np.eye(4) - g[mu], fwd)
        bwd = np.einsum('...ba,...sb->...sa', np.conj(np.roll(uf[mu], 1, axis=mu)), np.roll(psi, 1, axis=mu))
        ref += np.einsum('ts,...sa->...ta', np.eye(4) + g[mu], bwd)
    assert np.abs(lt.hop(uf, psi) - ref).max() < 1e-13
    assert np.abs(lt.hop(uf, psi, lt.hopping_blocks(uf)) - ref).max() < 1e-13


def test_potential_rises_and_string_tension_is_positive(result):
    for label, res in result['production'].items():
        p = res['potential']
        assert p['sigma_over_err'] > 5
        r = np.array(p['r'])
        v = np.array(p['V_t2'], float)
        on_axis = [k for k in range(len(r)) if abs(r[k] - round(r[k])) < 1e-9 and np.isfinite(v[k])]
        assert all(v[b] > v[a] for a, b in zip(on_axis[:-1], on_axis[1:]))


def test_mass_fit_recovers_a_synthetic_mass():
    t = np.arange(5, 12)
    big_t = 24
    y = 3.0 * (np.exp(-0.41 * t) + np.exp(-0.41 * (big_t - t)))
    m, chi2 = an.fit_mass(y, 0.01 * y, t, big_t, 'cosh')
    assert abs(m - 0.41) < 1e-9 and chi2 < 1e-12
    y = 2.0 * np.exp(-0.9 * t)
    m, chi2 = an.fit_mass(y, 0.01 * y, t, big_t, 'exp')
    assert abs(m - 0.9) < 1e-9
    c = 3.0 * (np.exp(-0.41 * np.arange(big_t + 1)) + np.exp(-0.41 * (big_t - np.arange(big_t + 1))))
    assert np.allclose(an.effective_mass(c, big_t, 'cosh')[:11], 0.41, atol=1e-9)


def test_mass_fit_is_the_minimum_of_chi2():
    rng = np.random.default_rng(7)
    t = np.arange(5, 13)
    big_t = 24
    y = 1.7 * (np.exp(-0.38 * t) + np.exp(-0.38 * (big_t - t))) * (1 + 0.02 * rng.normal(size=t.size))
    s = 0.02 * y
    m, chi2 = an.fit_mass(y, s, t, big_t, 'cosh')
    for dm in (-1e-4, 1e-4):
        f = np.exp(-(m + dm) * t) + np.exp(-(m + dm) * (big_t - t))
        a = (y * f / s ** 2).sum() / (f * f / s ** 2).sum()
        assert (((y - a * f) / s) ** 2).sum() > chi2


def test_potential_fit_recovers_synthetic_parameters():
    r = np.array([1.414, 1.732, 2.0, 2.236, 2.828, 3.0, 3.464, 4.0])
    v = 0.6 + 0.16 * r - 0.25 / r
    v0, sigma, e, chi2 = an.potential_fit(r, v, 0.01 * np.ones_like(r))
    assert abs(v0 - 0.6) < 1e-10 and abs(sigma - 0.16) < 1e-10 and abs(e - 0.25) < 1e-10


def test_normalized_mass_systematics_include_both_fit_choices():
    primary = {'mrho_over_sqrtsigma_chiral': [2.0, 0.1], 'mrho_r0_chiral': [3.0, 0.1],
               'mN_over_mrho_chiral': [1.5, 0.1], 'sqrt_sigma': [2.0, 0.1], 'r0': [3.0, 0.1]}
    alternate = {k: list(v) for k, v in primary.items()}
    alternate['mrho_over_sqrtsigma_chiral'][0] = 2.2
    alternate['mrho_r0_chiral'][0] = 3.3
    fits = {'primary': {'sigma': 4.0, 'r0': 3.0}, 'alternate': {'sigma': 1.0, 'r0': 6.0}}
    scale_only = an.chiral_systematics({'primary': primary}, fits)
    both = an.chiral_systematics({'primary': primary, 'alternate': alternate}, fits)
    assert scale_only['mrho_over_sqrtsigma_chiral'] == 2.0
    assert scale_only['mrho_r0_chiral'] == 3.0
    assert both['mrho_over_sqrtsigma_chiral'] == pytest.approx(2.4)
    assert both['mrho_r0_chiral'] == pytest.approx(3.6)
    assert both['mN_over_mrho_chiral'] == 0.0
    assert both['sqrt_sigma'] == 1.0 and both['r0'] == 3.0


def test_recorded_mass_systematics_cover_potential_fit_variations(result):
    for res in result['production'].values():
        ch = res['chiral']
        mass = ch['primary']['mrho_chiral'][0]
        central = ch['primary']['mrho_over_sqrtsigma_chiral'][0]
        for fit in res['potential']['fits'].values():
            shift = abs(mass / np.sqrt(fit['sigma']) - central)
            assert ch['systematic']['mrho_over_sqrtsigma_chiral'] >= shift - 1e-14
        assert ch['systematic_potential_only']['mrho_over_sqrtsigma_chiral'] > 0
        assert ch['systematic']['mN_over_mrho_chiral'] == ch['systematic_hadron_only']['mN_over_mrho_chiral']


def test_pion_lies_below_rho_and_resolved_nucleon_above_it(result):
    for res in result['production'].values():
        used = res['chiral']['kappas_used']
        assert len(used) >= 2
        for row in res['hadrons']:
            if row['kappa'] in used:
                assert row['mass_pion'][0] < row['mass_rho'][0]
                if row['kappa'] not in res['nucleon_not_resolved_at']:
                    assert row['mass_rho'][0] < row['mass_nucleon'][0]
    assert result['production'][HEADLINE]['nucleon_not_resolved_at'] == [0.1675]


def test_melting(result):
    lo, hi = result['melting']['below'], result['melting']['above']
    assert abs(lo['slope'][0] - 0.9896) < 5e-5 and abs(hi['slope'][0] - 0.00107) < 5e-6
    assert hi['slope'][0] < 0.01 * lo['slope'][0]
    r = np.array([2.0, 3.0, 4.0, 5.0])
    assert abs(weighted_slope(r, 0.3 + 0.7 * r, np.array([0.1, 0.2, 0.1, 0.3])) - 0.7) < 1e-12


def test_compiled_kernel_equals_numpy_reference():
    pytest.importorskip('numba')
    from ..colour_deconfinement import su3_numba as su3
    from .production import kernel as fast
    rng = np.random.default_rng(3)
    dims = (4, 4, 4, 6)
    u = lt.hot_links(dims, rng)
    psi = rng.normal(size=dims + (4, 3)) + 1j * rng.normal(size=dims + (4, 3))
    ref = lt.hop(u, psi)
    vol = int(np.prod(dims))
    flat_u = np.ascontiguousarray(np.transpose(u, (4, 3, 2, 1, 0, 5, 6)).reshape(vol, 4, 3, 3))
    flat_psi = np.ascontiguousarray(np.transpose(psi, (3, 2, 1, 0, 4, 5)).reshape(vol, 4, 3))
    up, dn = su3.neighbours(*dims)
    out = np.zeros_like(flat_psi)
    fast.hop(flat_u, up, dn, flat_psi, out, np.arange(vol, dtype=np.int64))
    back = np.transpose(out.reshape(dims[3], dims[2], dims[1], dims[0], 4, 3), (3, 2, 1, 0, 4, 5))
    assert np.abs(back - ref).max() < 1e-12


def test_shared_jackknife_covariance_cancels_a_common_scale():
    samples = np.array([[1., 2.], [2., 4.], [3., 6.], [4., 8.]])
    covariance = an.jk_cov(samples)
    contrast = np.array([-2., 1.])
    assert abs(contrast @ covariance @ contrast) < 1e-12
    assert an.jk_err(samples[:, 1]/samples[:, 0]) == 0
    assert covariance[0, 0] + covariance[1, 1] > 0


def test_joint_covariance_rejects_unpaired_or_nonfinite_replicas():
    import pytest
    for samples in ([[1., 2.]], [[1., 2.], [float('nan'), 3.]], [1., 2.]):
        with pytest.raises(ValueError):
            an.jk_cov(samples)
