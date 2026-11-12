"""Complete native-mass-channel time evolution and finite-spacing boost tests."""
from itertools import product
import numpy as np
from scipy.linalg import expm
from model.native_mass_channels import algebra
from . import lattice as l


def hamiltonian(momentum, spacing, masses=(2., 1., 2., 2.), speed=1.):
    if spacing <= 0 or speed <= 0:
        raise ValueError('positive source spacing and speed required')
    gamma, beta = algebra.folded_generators()
    kinetic = speed*np.sin(spacing*np.asarray(momentum))/spacing
    h = sum(float(k)*np.asarray(g, dtype=complex) for k, g in zip(kinetic, gamma))
    h += sum(float(m)*np.asarray(b, dtype=complex) for m, b in zip(masses, beta))
    return h


def run():
    gamma, beta = algebra.folded_generators()
    q = np.asarray(beta[0]*(beta[1]+2*beta[2]+2*beta[3])/3, dtype=complex)
    one = np.eye(16)
    rng = np.random.default_rng(786)
    raw = rng.normal(size=16)+1j*rng.normal(size=16)
    rows = []
    worst_time = 0.
    for a in (1/2, 1/4, 1/8, 1/16):
        for p in ((0., 0., 0.), (.4, -.3, .2), (1., 1., 1.)):
            h = hamiltonian(p, a)
            eigen = np.linalg.eigvalsh(h)
            expected = []
            for sign in (-1, 1):
                sector = (one+sign*q)/2
                mass = abs(2+sign*3)
                kinetic2 = float(np.sum((np.sin(a*np.asarray(p))/a)**2))
                energy = np.sqrt(kinetic2+mass*mass)
                assert np.linalg.norm((h@h-(kinetic2+mass*mass)*one)@sector) < 1e-12
                vector = sector@raw
                vector /= np.linalg.norm(vector)
                time = .73
                exact = np.cos(energy*time)*vector-1j*np.sin(energy*time)*(h@vector)/energy
                direct = expm(-1j*time*h)@vector
                error = float(np.linalg.norm(direct-exact)); worst_time = max(worst_time, error)
                assert error < 1e-12
                shell_defect = float(np.dot(p, p)-kinetic2)
                bound = float(a*a*np.sum(np.asarray(p)**4)/3)
                assert -1e-14 <= shell_defect <= bound+1e-14
                expected.extend([-energy]*4+[energy]*4)
                rows.append({'spacing': a, 'momentum': list(p), 'sector': sign,
                             'supplied_rest_coefficient': mass,
                             'energy': round(float(energy), 12),
                             'squared_shell_defect': round(shell_defect, 12),
                             'shell_error_upper': round(bound, 12)})
            assert np.max(np.abs(eigen-np.sort(expected))) < 1e-12
    # Compare boosted continuum momentum with the same finite-spacing source.
    velocity = .4
    lorentz = 1/np.sqrt(1-velocity*velocity)
    p0 = np.array([.5, .2, -.1]); mass = 1.
    energy0 = np.sqrt(np.dot(p0, p0)+mass*mass)
    boosted_p = p0.copy(); boosted_p[0] = lorentz*(p0[0]-velocity*energy0)
    boosted_energy = lorentz*(energy0-velocity*p0[0])
    assert abs(boosted_energy**2-np.dot(boosted_p, boosted_p)-mass*mass) < 1e-12
    boosts = []
    for a in (1/2, 1/4, 1/8, 1/16):
        measured = np.sqrt(np.sum((np.sin(a*boosted_p)/a)**2)+mass*mass)
        defect = float(abs(measured-boosted_energy))
        bound = float(a*a*np.sum(boosted_p**4)/(6*mass))
        assert defect <= bound+1e-14
        boosts.append({'spacing': a, 'boosted_momentum': [round(float(v), 12) for v in boosted_p],
                       'energy_difference': round(defect, 14), 'bound': round(bound, 14)})
    assert all(boosts[i]['energy_difference'] > 3.9*boosts[i+1]['energy_difference'] for i in range(3))
    h = hamiltonian((.4, -.3, .2), .25)
    bad = h+.1*one
    sector = (one-q)/2
    k2 = float(np.sum((np.sin(.25*np.array([.4, -.3, .2]))/.25)**2))
    assert np.linalg.norm((bad@bad-(k2+1)*one)@sector) > .1
    return {'task': 'moving_mass_shell', 'source': 'model.native_mass_channels.algebra.folded_generators',
            'full_folded_dimension': 16, 'retained_valleys': 8, 'signed_energy_rank_per_mass': 4,
            'mass_coefficients': [2, 1, 2, 2], 'rows': rows, 'boost_refinement': boosts,
            'maximum_direct_time_error': round(worst_time, 14),
            'checks': {'complete_spectrum_and_signed_branches': True,
                       'direct_exponential_matches_squared_wave_law': True,
                       'finite_spacing_shell_bound': True, 'boost_error_vanishes_with_spacing': True},
            'mutations': {'unreported_energy_shift_rejected': True},
            'scope': 'Moving free one-body modes of the declared native mass-channel family. Coefficients are supplied; no Yukawa derivation or identified moving charged defect is claimed. Finite spacing has a preferred lattice frame.'}
