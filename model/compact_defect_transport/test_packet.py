from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_receipt_replays_all_source_states_and_gates(replay):
    packet.validate(strict_json.load(packet.OUTPUT), replay)
    controls = replay['complete_source_transport_controls']
    assert controls['all_continuing_records_retained'] == 12
    assert len(controls['nonidentity_transport_pairs']) == 3
    assert controls['native_current_count'] == 12
    assert controls['localized_annihilators'] == 180
    assert all(row['bare_quench_leakage_probability'] > 0
               for row in controls['nonidentity_transport_pairs'])


@pytest.mark.parametrize('mutation', ['missing_record', 'false_scope', 'lost_current', 'corrupt_bound', 'hash'])
def test_mutated_receipts_fail_complete_validation(replay, mutation):
    value = deepcopy(replay)
    if mutation == 'missing_record':
        value['exact_complete_repair']['all_continuing_instrument_records'].pop()
    elif mutation == 'false_scope':
        value['scope']['relativistic_dispersion_or_interacting_continuum'] = True
    elif mutation == 'lost_current':
        value['complete_source_transport_controls']['native_current_count'] = 11
    elif mutation == 'corrupt_bound':
        value['rational_transport_bounds'][0]['spectral_projection_difference_norm_upper'] = '0'
    else:
        value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError):
        packet.validate(value, replay)


def test_source_manifest_contains_implementation_and_no_research_note():
    names = packet.source_bindings()
    assert 'model/compact_defect_transport/transport.py' in names
    assert 'model/compact_defect_pair/spectral.py' in names
    assert 'model/lattice_texture_selection/shells.py' in names
    assert not any(name.startswith('docs/') for name in names)
