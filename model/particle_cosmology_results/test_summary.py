import json
from pathlib import Path

import pytest

from .summary import authenticate, certificate, digest, require_false


def test_actual_source_join_and_mass_ancestry():
    result = certificate()
    assert len(result['tasks']) == 20
    assert len(result['conditional_GeV_targets']) == 17
    assert result['physical_premise_promotions'] == []
    assert result['conditional_neutrino_mass_observables_GeV']['all_masses_are_conditional_calibration'] is True
    assert result['lepton_archive_resources'][-1]['retained_quantum_archive_qubits'] == 2


def test_modified_dependency_and_forged_success_rejected(tmp_path: Path):
    directory = tmp_path / 'model' / 'control'
    directory.mkdir(parents=True)
    source = directory / 'source.py'
    source.write_text('fixed_source = 1\n')
    record = {'schema': 'oph.control.v1', 'all_checks_pass': True,
              'source_bindings': {'model/control/source.py': digest(source)}}
    path = directory / 'results.json'
    path.write_text(json.dumps(record))
    authenticate('control', tmp_path)
    source.write_text('fixed_source = 2\n')
    with pytest.raises(ValueError, match='changed input'):
        authenticate('control', tmp_path)
    record['all_checks_pass'] = False
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='failed scientific checks'):
        authenticate('control', tmp_path)
    del record['all_checks_pass']
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='failed scientific checks'):
        authenticate('control', tmp_path)


def test_physical_promotion_cannot_follow_from_receipt_success():
    require_false(False, 'control')
    for forged in (True, None, 0, 'false'):
        with pytest.raises(ValueError, match='Unsupported physical promotion'):
            require_false(forged, 'control')


def test_quantum_sector_nonclosure_and_complete_histories_cannot_be_forged(tmp_path: Path):
    source = Path(__file__).parents[1] / 'quantum_sector_measure' / 'results.json'
    original = json.loads(source.read_text())
    target = tmp_path / 'model' / 'quantum_sector_measure' / 'results.json'
    target.parent.mkdir(parents=True)
    for mutation in ('closure', 'history', 'weights'):
        result = json.loads(json.dumps(original))
        if mutation == 'closure':
            result['closure']['scalar_sector_adjoint_closure'] = True
        elif mutation == 'history':
            result['complete_histories'][0]['complete_records'].pop('000')
        else:
            result['stationarity']['family_weights_selected_by_channel'] = True
        target.write_text(json.dumps(result))
        with pytest.raises(ValueError, match='incomplete quantum sector certificate'):
            authenticate('quantum_sector_measure', tmp_path)


def test_impossible_expansion_counts_rejected_even_with_success_flag(tmp_path: Path):
    target = tmp_path / 'model' / 'expanding_quadratic_particles' / 'results.json'
    target.parent.mkdir(parents=True)
    source = target.with_name('source.py')
    source.write_text('declared_modes = 180\n')
    law = {'all181_total_number_probabilities': [1.0]+[0.0]*180,
           'expected_total_number': 0.0,
           'complete12_species_number_laws': [{'mean_number': 0.0, 'probabilities': [1.0, 0.0]} for _ in range(12)]}
    good = {'schema': 'oph.expanding_quadratic_particles.v1', 'all_checks_pass': True,
            'source_bindings': {'model/expanding_quadratic_particles/source.py': digest(source)},
            'protocols': [{'preparations': [{'complete_number_readout': law} for _ in range(2)]}
                          for _ in range(7)]}
    target.write_text(json.dumps(good))
    authenticate('expanding_quadratic_particles', tmp_path)
    for corruption in ('probabilities', 'mean', 'species', 'species_nan', 'species_cancellation', 'species_distribution'):
        bad = json.loads(json.dumps(good))
        row = bad['protocols'][0]['preparations'][0]['complete_number_readout']
        if corruption == 'probabilities':
            row['all181_total_number_probabilities'][1] = 2.0
        elif corruption == 'mean':
            row['expected_total_number'] = 29915.95
        elif corruption == 'species':
            row['complete12_species_number_laws'][0]['mean_number'] = 181.0
        elif corruption == 'species_nan':
            row['complete12_species_number_laws'][0]['mean_number'] = float('nan')
        elif corruption == 'species_cancellation':
            row['complete12_species_number_laws'][0]['mean_number'] = 181.0
            row['complete12_species_number_laws'][1]['mean_number'] = -181.0
        else:
            row['complete12_species_number_laws'][0]['probabilities'] = [1.0, 0.5]
        target.write_text(json.dumps(bad))
        with pytest.raises(ValueError, match='invalid normalized (source|species) count law'):
            authenticate('expanding_quadratic_particles', tmp_path)


def test_spectroscopy_confidence_cannot_promote_its_input_domain(tmp_path: Path):
    original = json.loads((Path(__file__).parents[1] / 'finite_particle_spectroscopy' / 'results.json').read_text())
    target = tmp_path / 'model' / 'finite_particle_spectroscopy' / 'results.json'
    target.parent.mkdir(parents=True)
    mutations = {
        'coverage_domain': 'all original pre-probe source inputs',
        'all_source_unconditional_65_over72_coverage_claimed': True,
        'vacuum_many_excitation_outcomes_supply_mass_interval': True,
        'all_input_bad_phase_effect_bound': 'F_bad <= (7/72) I',
        'complete_original_probe_bound': '',
    }
    for key, value in mutations.items():
        bad = json.loads(json.dumps(original))
        bad['uniform_confidence'][key] = value
        target.write_text(json.dumps(bad))
        with pytest.raises(ValueError, match='confidence domain or complete-probe pullback lost'):
            authenticate('finite_particle_spectroscopy', tmp_path)
