"""Source, finite-reference and physical-scope mutations fail closed."""
from copy import deepcopy

import pytest
from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def retained():
    return strict_json.load(packet.OUTPUT)


def test_complete_same_source_receipt(retained):
    packet.validate(retained)
    assert retained['scope']['same_bounded_radial_one_body_and_interacting_Hamiltonian']
    assert retained['core']['retained_27_site_mass64_source']['full_Fock_dimension'] == 2**54


@pytest.mark.parametrize('key', [
    'actual_center_parity_is_exact_zero_orbital_occupation',
    'source_selects_or_prepares_radial_mass_texture',
    'Gibbs_preparation_or_repair_equilibration_executed',
    'two_to_54_Fock_matrix_executed', 'volume_uniform_interacting_gap_proved',
    'strong_mass_bound_is_fixed_mass_spatial_continuum_limit',
    'native_repair_protection_of_poles_proved', 'responsive_Einstein_geometry_instantiated',
    'completion_gate_accepted'])
def test_domain_promotions_rejected(retained, key):
    value = deepcopy(retained)
    value['scope'][key] = True
    with pytest.raises(ValueError, match='boundary'):
        packet.validate(value)


def test_radial_symmetry_cannot_be_replaced_by_plain_time_reversal(retained):
    value = deepcopy(retained)
    value['algebra']['source_symmetries']['commuting_Kramers_antiunitary'] = 'J K'
    with pytest.raises(ValueError, match='algebra'):
        packet.validate(value)


def test_exact_thermal_domain_and_json_types_are_preserved(retained):
    for key, replacement in [('stationary_under_native_repair_inferred', True),
                             ('actual_finite_interacting_Gibbs_state_faithful', 1),
                             ('modular_to_Hamiltonian_time_factor', '1')]:
        value = deepcopy(retained)
        value['core']['faithful_thermal_reference'][key] = replacement
        with pytest.raises(ValueError, match='core'):
            packet.validate(value)


def test_receipt_missing_pin_or_import_staleness_rejected(retained, monkeypatch):
    value = deepcopy(retained)
    value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(value)
    changed = dict(packet._IMPORT_BINDINGS)
    changed['model/radial_interacting_core/core.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='after.*import'):
        packet._produce()
