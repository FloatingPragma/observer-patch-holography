"""Complete payload mutation and source-custody rejection controls."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_stored_packet_matches_complete_cold_replay(fresh):
    packet.validate(strict_json.load(packet.OUTPUT), fresh)
    assert len(fresh['source_cases']) == 4


@pytest.mark.parametrize('field', ['physical_particle_species_identified',
    'physical_clock_or_GeV_calibration_derived', 'three_families_forced_by_topological_index',
    'fermion_creation_changes_texture_winding', 'unmodified_site_repair_preserves_core',
    'full_interacting_Fock_band_isolation_proved',
    'relativistic_dispersion_or_pole_mass_identification_proved', 'integrated_gate_accepted'])
def test_false_physical_promotions_rejected(fresh, field):
    wrong = deepcopy(fresh)
    wrong['scope'][field] = True
    with pytest.raises(ValueError, match='complete replay'):
        packet.validate(wrong, fresh)


@pytest.mark.parametrize('field', ['all_Nambu_rows', 'core_zero_Nambu_rows', 'all_native_port_count'])
def test_dropped_modes_or_currents_rejected(fresh, field):
    wrong = deepcopy(fresh)
    wrong['source_cases'][0][field] -= 1
    with pytest.raises(ValueError, match='complete replay'):
        packet.validate(wrong, fresh)


def test_stale_inherited_source_rejected_before_replay(monkeypatch):
    original = strict_json.load
    def corrupt(path):
        result = original(path)
        if str(path).endswith('radial_defect_gap/results.json'):
            key = next(iter(result['source_bindings']))
            result['source_bindings'][key] = '0'*64
        return result
    monkeypatch.setattr(strict_json, 'load', corrupt)
    with pytest.raises(ValueError, match='stale inherited source'):
        packet.source_bindings()
