from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_retained_receipt_reproduces_from_cold_source(fresh):
    packet.compare(strict_json.load(packet.OUTPUT), fresh)
    assert fresh['local_instrument_and_transport']['local_intertwiner']['native_raw_slots'] == 720
    assert fresh['locality_optimality']['all_port_subsets_checked'] == 4096


@pytest.mark.parametrize('kind', ['empty', 'boolean', 'physical', 'distance', 'hash', 'history'])
def test_plausible_but_false_receipts_fail_owned_replay(fresh, kind):
    value = deepcopy(fresh)
    if kind == 'empty':
        value = {}
    elif kind == 'boolean':
        value['physical_comparison_count_added'] = False
    elif kind == 'physical':
        value['remaining_obligations']['physical_particle_or_electron_identified'] = True
    elif kind == 'distance':
        value['locality_optimality']['bounds_by_core_size'][2]['sharp_sector_min_exterior_trace_distance'] = '0'
    elif kind == 'hash':
        value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    else:
        value['local_instrument_and_transport']['witness'][2]['exterior_trace_distance'] = '0'
    with pytest.raises(ValueError, match='complete source replay'):
        packet.compare(value, fresh)


def test_public_json_round_trip_and_sources_are_portable(fresh):
    packet.compare(strict_json.loads(packet.canonical(fresh)), fresh)
    assert all('\\' not in path and ':' not in path for path in fresh['source_bindings'])
    for path in ('model/local_exchange_charge/core.py',
                 'model/local_exchange_charge/test_packet.py',
                 'model/defect_field_bridge/full_domain.py',
                 'model/observer/source_scalar/compiler.py'):
        assert path in fresh['source_bindings']
