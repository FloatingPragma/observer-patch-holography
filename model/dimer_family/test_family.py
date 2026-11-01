import math

import numpy as np
import pytest

from model.dimer_family.family import certificate
from model.full_fock_field.field import FieldGeometry
from model.particle_response.spectrum import FieldSpectrum


def test_exact_whole_parameter_family():
    result = certificate()
    assert result['maximum_distinct_energies_per_charge'] == 2
    assert result['gap_ratio_upper_inclusive'] == 3
    assert not result['three_same_charge_distinct_mass_families_possible']


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('u,t', [(0., .2), (.001, .5), (1., .0625), (13., .4)])
def test_symbolic_levels_against_actual_all_axis_complete_field(axis, u, t):
    sides = [1, 1, 1]; sides[axis] = 2
    spectrum = FieldSpectrum(FieldGeometry(tuple(sides), 2*t, u),
                             chemical_potential=u/2, reference='ground')
    d = math.hypot(u, 4*t)
    assert abs(spectrum.energies[0]+(u+d)/2) < 2e-12
    for sector in (1, 3):
        gaps = spectrum.energies[spectrum.sectors == sector]-spectrum.energies[0]
        np.testing.assert_allclose(gaps, [d/2-t]*2+[d/2+t]*2, atol=2e-12)
        assert 1 < gaps[-1]/gaps[0] <= 3+1e-12


def test_recentring_is_a_real_hypothesis_not_hidden_freedom():
    u, t = 1., 1/16
    low = math.hypot(u, 4*t)/2-t
    spectrum = FieldSpectrum(FieldGeometry((2,1,1), 2*t, u),
                             chemical_potential=u/2+.99*low, reference='ground')
    assert spectrum.sectors[0] == 2
    gaps = spectrum.energies[spectrum.sectors == 3]-spectrum.energies[0]
    assert gaps[-1]/gaps[0] > 20
    assert len(np.unique(np.round(gaps, 10))) == 2


@pytest.mark.parametrize('u', [0., .01, 1., 8.])
def test_unique_ground_interval_and_charge_shifts_include_noninteracting_endpoint(u):
    t = .25
    d = math.hypot(u, 4*t)
    low = d/2-t
    for direction in (-1, 1):
        delta = direction*.99*low
        spectrum = FieldSpectrum(FieldGeometry((2,1,1), 2*t, u),
                                 chemical_potential=u/2+delta, reference='ground')
        assert spectrum.sectors[0] == 2
        for sector, sign in ((1, 1), (3, -1)):
            gaps = spectrum.energies[spectrum.sectors == sector]-spectrum.energies[0]
            np.testing.assert_allclose(gaps, [low+sign*delta]*2+[d/2+t+sign*delta]*2,
                                       atol=3e-12)
        outside = FieldSpectrum(FieldGeometry((2,1,1), 2*t, u),
                                chemical_potential=u/2+direction*1.01*low,
                                reference='gibbs', beta=.1)
        assert outside.sectors[0] != 2


def test_noninteracting_upper_doublet_is_dark_without_deleting_its_states():
    t = .25
    spectrum = FieldSpectrum(FieldGeometry((2,1,1), 2*t, 0.),
                             chemical_potential=0., reference='ground')
    operator, charge = spectrum.probe('single', 0)
    response = spectrum.measure(operator, charge)
    upper = np.abs(response.frequencies-3*t) < 1e-12
    lower = np.abs(response.frequencies-t) < 1e-12
    assert abs(response.forward[lower].sum()-.5) < 1e-12
    assert response.forward[upper].sum() < 1e-24
    # The full N=3 doublet remains in the spectrum despite zero probe weight.
    gaps = spectrum.energies[spectrum.sectors == 3]-spectrum.energies[0]
    assert np.count_nonzero(abs(gaps-3*t) < 1e-12) == 2
