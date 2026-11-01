"""Retrospective unit-independent comparisons of the constructed model modes.

The source law and detector do not consume physical_targets.json. The
comparison permits a common positive unit conversion; it never calibrates
one from a measured particle, fits a mass term, or treats finite-momentum
energy as rest mass. Conclusions concern this constructed sector only.
"""
from fractions import Fraction as F
from pathlib import Path

from tools import strict_json
from .detector import mode_catalog


TARGETS = Path(__file__).with_name('physical_targets.json')


def certificate():
    targets = strict_json.load(TARGETS)
    modes = mode_catalog()
    if len(modes) != 16 or len({row['flavor'] for row in modes}) != 8:
        raise ArithmeticError('the comparison requires the actual eight-flavor detector')
    if any(F(row['continuum_mass']) != 0 or abs(F(row['helicity'])) != F(1, 2)
           for row in modes):
        raise ArithmeticError('constructed masses or helicities changed; review the comparison')
    leptons = targets['charged_leptons']
    charged = []
    for row in leptons:
        mass, error = F(row['mass_MeV_c2']), F(row['reported_uncertainty_MeV_c2'])
        if mass <= error or error < 0:
            raise ArithmeticError('the published positive-mass target changed')
        charged.append({'name': row['name'], 'observed_mass_MeV_c2': str(mass),
            'reported_uncertainty_MeV_c2': str(error),
            'model_continuum_mass_in_any_positive_units': '0',
            'direct_mass_preserving_identification': False,
            'reason': 'a common positive unit conversion maps zero to zero, not to the observed positive mass',
            'source': row['source']})
    electron_mass = F(leptons[0]['mass_MeV_c2'])
    ratios = [{'name': row['name']+'/electron',
               'ratio_of_reported_central_masses': str(F(row['mass_MeV_c2'])/electron_mass),
               'decimal_display': float(F(row['mass_MeV_c2'])/electron_mass)}
              for row in leptons[1:]]
    magnitudes = {abs(F(row['energy_sqrt2_coefficient'])) for row in modes}
    if magnitudes != {F(1, 16)}:
        raise ArithmeticError('the selected finite-momentum energy degeneracy changed')
    neutrino = []
    for row in targets['neutrino_splittings']['entries']:
        low, high = map(F, row['reported_three_sigma_range'])
        if not 0 < low <= F(row['best_fit']) <= high:
            raise ArithmeticError('the nonzero published splitting target changed')
        neutrino.append({'name': row['name'], 'unit': row['unit'],
            'model_difference_of_squared_masses': '0',
            'reported_best_fit': row['best_fit'],
            'reported_three_sigma_range': row['reported_three_sigma_range'],
            'model_zero_in_reported_range': low <= 0 <= high})
    return {'status': 'retrospective comparison of the specified constructed sector',
        'target_edition': targets['edition'],
        'model_massive_relativistic_species_count': 0,
        'model_massless_Weyl_flavor_count': 8,
        'model_finite_detector_energy_magnitudes': 'sqrt(2)/16 in original model units',
        'model_continuum_mass_ratios': 'undefined (0/0), not a mass hierarchy',
        'charged_leptons': charged, 'observed_charged_lepton_central_mass_ratios': ratios,
        'finite_energy_as_mass_shortcut': {
            'permitted_dictionary_for_this_control': 'one common positive conversion of the selected finite-momentum energy magnitudes',
            'all_model_ratios': '1',
            'charged_lepton_hierarchy_reproduced': False,
            'finite_momentum_energy_is_not_rest_mass': True},
        'photon': {'masslessness_alone_compatible': True,
            'model_helicity_magnitude': '1/2', 'target_spin': targets['photon']['spin'],
            'direct_spin_preserving_identification': False,
            'source': targets['photon']['source'],
            'composite_or_other_sector_photon_excluded': False},
        'three_neutrino_mass_spectrum': {
            'direct_identification_of_three_of_these_massless_flavors': False,
            'comparisons': neutrino,
            'source': targets['neutrino_splittings']['source'],
            'single_massless_neutrino_excluded': False,
            'absolute_neutrino_masses_inferred': False,
            'scope': targets['neutrino_splittings']['interpretation']},
        'absolute_energy_unit_fitted': False,
        'mass_term_or_source_coefficient_fitted': False,
        'model_energies_converted_to_MeV': False,
        'new_prospective_prediction_successes': 0,
        'new_physical_species_identifications': 0,
        'whole_OPH_or_other_scalings_excluded': False,
        'remaining': 'A common physical vacuum, charge/spin and detector dictionary, mass generation and independently justified units are required for physical identification.'}
