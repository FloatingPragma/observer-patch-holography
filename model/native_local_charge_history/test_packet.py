"""Receipt, scalar typing and explicit physical scope controls."""
import copy
import pytest
from tools import strict_json
from . import packet


def test_owned_receipt_replays():
    packet.check()


@pytest.mark.parametrize('mutation', ('physical', 'omitted_erasure', 'boolean'))
def test_changed_receipt_rejected(tmp_path, monkeypatch, mutation):
    value = copy.deepcopy(packet.produce())
    if mutation == 'physical':
        value['physical_particle_identification'] = True
    elif mutation == 'omitted_erasure':
        del value['finite_control']['erasure_words']
    else:
        value['finite_control']['input_columns'] = True
    path = tmp_path/'results.json'
    path.write_text(packet.inherited.canonical(value))
    # Keep the owned source root stable; only replace the submitted input.
    monkeypatch.setattr(packet, 'OUTPUT', path)
    with pytest.raises(ValueError, match='fresh exact replay'):
        packet.check()
