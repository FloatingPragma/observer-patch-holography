from copy import deepcopy

import pytest

from . import packet


def test_receipt_replays_exactly():
    packet.verify(packet.OUTPUT)


def test_false_physical_promotion_is_rejected(tmp_path):
    value=packet.build();value['scope']['physical_alpha_s_prediction']='0.118'
    path=tmp_path/'promoted.json';path.write_text(packet.canonical(value))
    with pytest.raises(ValueError):packet.verify(path)


def test_source_defaults_and_probability_interval_are_bound(tmp_path):
    original=packet.build()
    for change in ('default','interval','hash'):
        value=deepcopy(original)
        if change=='default':
            value['native_color_and_electric_action']['committed_source_defaults']['hopping_J']='2'
        elif change=='interval':
            value['complete_two_time_readout'][0]['intervals']['electric_squared'][0]='0'
        else:
            key=next(iter(value['source_bindings']));value['source_bindings'][key]='0'*64
        path=tmp_path/(change+'.json');path.write_text(packet.canonical(value))
        with pytest.raises(ValueError):packet.verify(path)
