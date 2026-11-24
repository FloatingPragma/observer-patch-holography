"""Independent checks on the joint phase mixture and the disk-domain calculation."""
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

sys.path.insert(0, str(Path(__file__).resolve().parent))
import disk_domain
import fit


def test_shared_phase_mean_and_dependence():
    r = np.array([0.2, 0.7, 1.2])
    weights, multipliers = fit.phase_law(r, 1.0)
    np.testing.assert_allclose(weights @ multipliers, np.ones(3))
    counts = multipliers * r
    np.testing.assert_allclose(counts[:, 2] - counts[:, 0], 1)
    np.testing.assert_allclose(weights @ (counts[:, 0] * counts[:, 1]), 0.2)


def test_galaxy_density_against_dense_gaussian_and_direct_phase_integral():
    rows = np.array([[0.25, 16., 12., 1.], [0.75, 25., 14., 2.], [1.75, 49., 23., 1.5]])
    x = [np.log(0.8), 0.03, 0.09]
    actual = fit.GalaxyLikelihood({'one': rows}, 1.0).log_scores(x)[0]
    radius, vb2, speeds, errors = rows.T
    terms = []
    for phase in (0.125, 0.375, 0.625, 0.875):
        delta = np.sqrt(vb2*radius*fit.base.KPC_M*8e-11/1e6) * np.floor(radius+phase)/radius
        mean = np.sqrt(vb2+delta)
        covariance = np.diag(errors**2+0.03*mean**2) + 0.09*np.outer(mean, mean)
        terms.append(multivariate_normal.logpdf(speeds, mean=mean, cov=covariance))
    np.testing.assert_allclose(actual, logsumexp(terms)-np.log(4.), atol=1e-12, rtol=0)


def test_independent_radial_quadrature_and_force_conversion():
    result = disk_domain.run()
    radius = result['outer']['figure8']['radius_kpc']
    theta = np.linspace(0, 220, 220001)
    surface = np.exp(-theta/0.676 + 559.231*np.log1p(theta/375.220))
    radial_mass = theta * surface
    boundary = radius / disk_domain.KPC_PER_ARCSEC
    inner = theta < boundary
    total = np.trapezoid(radial_mass, theta)
    enclosed = np.trapezoid(np.r_[radial_mass[inner], np.interp(boundary,theta,radial_mass)],
                            np.r_[theta[inner], boundary])/total
    np.testing.assert_allclose(enclosed, result['gas_mass_fraction_in_outer_cylinder'], atol=1e-8, rtol=0)
    assert enclosed < result['baryonic_mass_fraction_upper_in_outer_cylinder'] < 0.7
    for row in result['outer'].values():
        assert row['spherical_equivalent_required_excess_msun'] > 0
        np.testing.assert_allclose(row['spherical_equivalent_required_excess_msun']
                                   / row['spherical_equivalent_min_excess_msun'],
                                   row['central_excess_fraction_of_min'])


def test_exterior_selection_and_speed_with_independent_si_units():
    import json

    result = disk_domain.run()
    orbits = result['exterior_orbits']
    assert set(orbits) == {'NGC2685', 'NGC3838', 'NGC6798'}
    assert [k for k,v in orbits.items() if v['mass90_pass']] == ['NGC3838']
    galaxy = next(row for row in json.loads((fit.ROOT/'outer_orbit_transfer/inputs.json').read_text())['galaxies']
                  if row['name'] == 'NGC3838')
    gm = 1.3271244e20*(0.8*galaxy['stellar_unit_msun']+galaxy['helium_hi_msun'])
    radius = galaxy['radius_kpc']*fit.base.KPC_M
    speed = np.sqrt(gm/radius + np.sqrt(gm*result['a0_m_s2']))/1000
    np.testing.assert_allclose(speed, orbits['NGC3838']['smooth_point_mass_kms'])
    lo,hi = orbits['NGC3838']['phase_conditioned_speed_range_kms']
    assert lo < speed < hi < orbits['NGC3838']['observed_kms']
