"""Reject physical promotions and omitted source sectors in the field panel."""
import copy
import pytest

from . import field_panels


@pytest.fixture(scope='module')
def data():
    return field_panels.load()


def test_source_receipts_and_domains(data):
    field_panels.validate(data)
    assert len(data['pins']) > 350


@pytest.mark.parametrize('name,change', [
    ('interacting limit', lambda d: d['field']['theorems'].__setitem__('nontrivial_interacting_relativistic_limit_proved', True)),
    ('vector identification', lambda d: d['field']['theorems'].__setitem__('Lorentz_vector_or_Proca_identification_proved', True)),
    ('hidden empirical input', lambda d: d['field']['construction']['measured_inputs'].append('electron mass')),
    ('dropped screening', lambda d: d['charge']['full_source_charged_bins']['ec']['energy_filtered_localization'].__setitem__('origin_probability_lower_inferred_by_discarding_screened_complement', True)),
    ('cycle promoted to continuum', lambda d: d['collision']['interval_certified_five_block_collision'].__setitem__('three_dimensional_continuum_claim', True)),
])
def test_invalid_interpretation_is_rejected(data, name, change):
    altered = copy.deepcopy(data)
    change(altered)
    with pytest.raises(ValueError):
        field_panels.validate(altered)
