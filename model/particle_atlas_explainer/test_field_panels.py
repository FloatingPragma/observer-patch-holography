"""Reject physical promotions and omitted source sectors in the field panel."""
import copy
import re
from fractions import Fraction
import pytest

from . import field_panels
from .canvas import Canvas


@pytest.fixture(scope='module')
def data():
    return field_panels.load()


def test_source_receipts_and_domains(data):
    field_panels.validate(data)
    assert len(data['pins']) > 350


def test_printed_energy_and_fidelity_bounds_are_outward(data):
    canvas = Canvas()
    field_panels.charged_field_panel(canvas, data, 0)
    text = ''.join(canvas.items)
    match = re.search(r'Total charged energy: \[([\d.]+), ([\d.]+)\]', text)
    lo, hi = map(Fraction, match.groups())
    actual_lo, actual_hi = map(Fraction, data['charged_field']['actual_charged_ground_energy_interval'])
    assert lo <= actual_lo <= actual_hi <= hi
    shown = Fraction(re.search(r'fidelity above ([\d.]+)\.', text).group(1))
    assert shown <= Fraction(data['charged_clock']['finite_rational_gate_record_bounds']['accepted_source_ground_fidelity_lower'])


@pytest.mark.parametrize('name,change', [
    ('interacting limit', lambda d: d['field']['theorems'].__setitem__('nontrivial_interacting_relativistic_limit_proved', True)),
    ('vector identification', lambda d: d['field']['theorems'].__setitem__('Lorentz_vector_or_Proca_identification_proved', True)),
    ('hidden empirical input', lambda d: d['field']['construction']['measured_inputs'].append('electron mass')),
    ('dropped screening', lambda d: d['charge']['full_source_charged_bins']['ec']['energy_filtered_localization'].__setitem__('origin_probability_lower_inferred_by_discarding_screened_complement', True)),
    ('cycle promoted to continuum', lambda d: d['collision']['interval_certified_five_block_collision'].__setitem__('three_dimensional_continuum_claim', True)),
    ('charged GeV promotion', lambda d: d['charged_field']['scope'].__setitem__('mass_GeV', 0.511)),
    ('missing charged complement', lambda d: d['charged_field'].__setitem__('additional_state_truncation', True)),
    ('paused clock', lambda d: d['charged_clock']['phase_apparatus'].__setitem__('common_source_evolution_retained', False)),
    ('exact projection claim', lambda d: d['charged_clock']['scope'].__setitem__('exact_ground_projection', True)),
    ('named vector boson', lambda d: d['vector']['construction'].__setitem__('Standard_Model_particle_assignment', 'Z')),
    ('silent reference conjugation', lambda d: d['conjugation']['reference_boundary'].__setitem__('same_reference_total_energy_equality_claimed', True)),
])
def test_invalid_interpretation_is_rejected(data, name, change):
    altered = copy.deepcopy(data)
    change(altered)
    with pytest.raises(ValueError):
        field_panels.validate(altered)
