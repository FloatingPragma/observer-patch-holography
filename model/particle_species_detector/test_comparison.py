from copy import deepcopy
from fractions import Fraction as F

import pytest

from . import comparison


def test_zero_masses_fail_positive_targets_without_a_fitted_unit():
    result = comparison.certificate()
    assert all(not row['direct_mass_preserving_identification'] for row in result['charged_leptons'])
    assert result['new_physical_species_identifications'] == 0
    assert result['absolute_energy_unit_fitted'] is False
    assert result['mass_term_or_source_coefficient_fitted'] is False
    ratios = result['observed_charged_lepton_central_mass_ratios']
    assert 206 < F(ratios[0]['ratio_of_reported_central_masses']) < 207
    assert 3477 < F(ratios[1]['ratio_of_reported_central_masses']) < 3478
    assert result['finite_energy_as_mass_shortcut']['charged_lepton_hierarchy_reproduced'] is False


def test_massless_neutrino_and_composite_photon_are_not_overexcluded():
    result = comparison.certificate()
    assert result['photon']['direct_spin_preserving_identification'] is False
    assert result['photon']['composite_or_other_sector_photon_excluded'] is False
    assert result['three_neutrino_mass_spectrum']['single_massless_neutrino_excluded'] is False
    assert result['three_neutrino_mass_spectrum']['absolute_neutrino_masses_inferred'] is False
    assert all(not row['model_zero_in_reported_range'] for row in result['three_neutrino_mass_spectrum']['comparisons'])
    assert result['whole_OPH_or_other_scalings_excluded'] is False


@pytest.mark.parametrize('key,value', [('continuum_mass','1'), ('helicity','1')])
def test_a_changed_model_requires_review_instead_of_reusing_the_exclusion(monkeypatch, key, value):
    modes = deepcopy(comparison.mode_catalog())
    modes[0][key] = value
    monkeypatch.setattr(comparison, 'mode_catalog', lambda: modes)
    with pytest.raises(ArithmeticError, match='masses or helicities changed'):
        comparison.certificate()
