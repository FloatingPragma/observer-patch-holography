"""Pole normalization, domain boundaries and source-bound catalogue replay."""
import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'radial_catalogue', Path(__file__).with_name('derive_radial_excitation_catalogue.py'))
catalogue = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(catalogue)


def test_complete_low_energy_states_and_distinct_bulk_domain():
    data = catalogue.scientific_catalogue()
    zero, minus, plus, negative, positive = data['rows']
    assert zero['multiplicity'] == 2
    assert (minus['particle_number'], plus['particle_number']) == (26, 28)
    assert minus['multiplicity'] == plus['multiplicity'] == 1
    assert minus['energy_model_interval'] == plus['energy_model_interval']
    lo, hi = map(F, plus['energy_model_interval'])
    assert F('0.00000762310') < lo <= hi < F('0.00000762974')
    assert list(map(F, minus['signed_retarded_pole_interval'])) == [-hi, -lo]
    assert F(minus['other_charged_excitation_energy_lower']) > 61
    assert negative['multiplicity_at_each_signed_energy'] + positive['multiplicity_at_each_signed_energy'] == 16
    assert all(row['observed_species'] is None for row in data['rows'])
    assert all(row['moving_defect_rest_energy'] is None for row in data['rows'])
    assert data['scope']['complete_physical_particle_catalogue'] is False
    assert plus['physical_units']['energy_GeV_interval'] is None


def test_normalized_probe_bounds_preserve_probability_and_length_scaling():
    core = catalogue.core_bounds()
    row = catalogue.spectral_probe_bounds(core)
    center = F(row['center_probability_lower'])
    fourier = F(row['each_Fourier_probe_probability_lower'])
    assert F('0.99959') < center < 1
    assert F('0.02989') < fourier < F(1, 27)
    assert 27*fourier < 1
    assert F(row['each_Fourier_probe_probability_upper']) == 1-26*fourier
    assert F(row['RMS_position_squared_about_origin_upper_model_length_squared']) < F('0.00123')
    doubled = dict(core, spacing=2*core['spacing'])
    scaled = catalogue.spectral_probe_bounds(doubled)
    assert scaled['center_probability_lower'] == row['center_probability_lower']
    assert F(scaled['RMS_position_squared_about_origin_upper_model_length_squared']) == 4*F(row['RMS_position_squared_about_origin_upper_model_length_squared'])


@pytest.mark.parametrize('key,value', [('radius', 2), ('each_pole_matrix_rank', 3),
    ('projection_angle_bound', F(1, 6)), ('origin_weight_lower', F(1, 2))])
def test_inapplicable_visibility_premises_are_rejected(key, value):
    source = dict(catalogue.core_bounds(), **{key: value})
    with pytest.raises(ValueError):
        catalogue.spectral_probe_bounds(source)


def test_outward_decimal_display_does_not_shrink_the_interval():
    for value in (F(1, 3), F(-1, 3), F(960588587329007, 126010204829305274368)):
        assert F(catalogue.decimal_bound(value)) <= value
        assert F(catalogue.decimal_bound(value, lower=False)) >= value


def test_retained_json_markdown_and_scientific_mutations():
    value = catalogue.strict_json.load(catalogue.OUTPUT)
    catalogue.validate(value)
    assert catalogue.MARKDOWN.read_text() == catalogue.render_markdown(value)
    for mutation in ('mass', 'species', 'pin'):
        damaged = copy.deepcopy(value)
        if mutation == 'mass':
            damaged['catalogue']['rows'][1]['moving_defect_rest_energy'] = '1'
        elif mutation == 'species':
            damaged['physical_particle_identification_count'] = 1
        else:
            key = next(iter(damaged['source_bindings']))
            damaged['source_bindings'][key] = '0'*64
        with pytest.raises(ValueError):
            catalogue.validate(damaged)
