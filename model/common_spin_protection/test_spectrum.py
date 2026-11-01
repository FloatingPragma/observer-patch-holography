"""Complete fixed-number diagonalization and time-domain response checks."""
import math

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from ..full_fock_field.field import FullFock
from . import spectrum as s, spin


@pytest.mark.parametrize('number', [2, 4, 5, 6, 8])
def test_highest_weight_census_matches_every_direct_number_sector_eigenvalue(number):
    # This direct matrix retains all states of the selected number; it uses
    # no spin kernel or highest-weight reduction to obtain the oracle spectrum.
    full = FullFock(12).hamiltonian(spin.GEOMETRY)
    addresses = np.array([i for i in range(4096) if i.bit_count() == number])
    direct = np.linalg.eigvalsh(full[addresses][:, addresses].toarray())-number/2
    repeated = []
    for twice_spin in range(number % 2, min(number, 12-number)+1, 2):
        chart = s.chart(number, twice_spin)
        repeated.extend(np.repeat(chart.energies, twice_spin+1))
    assert len(repeated) == math.comb(12, number)
    np.testing.assert_allclose(np.sort(repeated), direct, atol=3e-14, rtol=0)


def test_all_number_sectors_and_full_105_level_spin_one_spectrum_are_retained():
    rows = s.census()
    assert sum(row['total_sector_dimension'] for row in rows) == 4096
    for n in range(13):
        assert sum(row['total_sector_dimension'] for row in rows if row['number'] == n) == math.comb(12, n)
    triplet = s.chart(4, 2)
    assert triplet.dimension == 105
    assert len(triplet.energies) == 105
    assert triplet.largest_residual < 1e-12
    assert np.ptp(triplet.energies) > 1


@pytest.mark.parametrize('kind', ['hole', 'particle', 'pair_hole', 'pair_particle'])
def test_every_local_probe_weight_against_direct_moments_and_time_evolution(kind):
    receipt = s.probe_spectrum(kind)
    e0, ground = s.ground()
    f = FullFock(12)
    a, b = f.annihilator(0), f.annihilator(2)
    probe = {'hole': a, 'particle': a.conj().T, 'pair_hole': b@a,
             'pair_particle': a.conj().T@b.conj().T}[kind]
    state = probe@ground
    energies = np.array(receipt['hamiltonian_model_gaps'])
    weights = np.array(receipt['all_highest_weight_residues'])
    assert len(energies) == (105 if kind.startswith('pair_') else 210)
    assert abs(sum(weights)-np.vdot(state, state).real) < 1e-14
    generator = s.hubbard_generator()
    assert abs(np.dot(weights, energies)-np.vdot(state, generator@state-e0*state).real) < 1e-14
    for time in (-.37, .83):
        response = np.vdot(state, expm_multiply(-1j*time*generator, state))*np.exp(1j*time*e0)
        assert abs(response-np.dot(weights, np.exp(-1j*time*energies))) < 2e-14


def test_reference_and_particle_candidate_scope_are_explicit():
    result = s.report()
    assert result['ground_number'] == 6 and result['ground_twice_spin'] == 0
    assert result['reference_is_faithful'] is False
    assert result['absolute_physical_masses_predicted'] == 0
    assert result['interacting_spectrum_is_stationary_under_repair'] is False
    assert 0 < result['full_Hamiltonian_gap'] < .01
    for kind, probe in result['complete_local_probe_responses'].items():
        assert probe['lowest_multiplet_residue'] > 1e-4
        assert probe['physical_detector_executed'] is False
        if not kind.startswith('pair_'):
            assert probe['repair_access'] == 'outside outcome only'


@pytest.mark.parametrize('args', [(True, 2), (4, True), (4, 1), (13, 1), (4, 6)])
def test_malformed_charge_spin_domains_reject(args):
    with pytest.raises(ValueError):
        s.chart(*args)
