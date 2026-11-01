import copy
import json

import pytest

from model.thermodynamic_ground import packet


def test_fresh_exact_packet_roundtrip_and_all_source_bindings():
    value=packet.report()
    assert json.loads(packet.canonical(value))==value
    assert value['source_bindings']==packet.source_bindings()
    packet.validate(value)


@pytest.mark.parametrize('key',[
    'ground_state_computed_or_prepared',
    'complete_repair_stationarity_asserted',
    'spacing_zero_or_Lorentz_continuum_asserted',
    'state_functional_local_faithfulness_asserted',
])
def test_no_silent_upgrade_from_existence_to_execution_or_continuum(key):
    value=packet._produce()
    changed=copy.deepcopy(value)
    changed['existence_theorem'][key]=True
    with pytest.raises(ValueError,match='promoted'):packet.validate(changed)


def test_receipt_rejects_half_filling_or_lost_source_control():
    value=packet._produce()
    value['density']['particle_hole']['half_filling_concluded']=True
    with pytest.raises(ValueError,match='density'):packet.validate(value)
