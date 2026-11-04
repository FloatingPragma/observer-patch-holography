from copy import deepcopy
import hashlib

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_cold_source_bound_receipt_and_proof_custody(replay):
    stored = strict_json.load(packet.OUTPUT)
    packet.validate(stored, replay)
    own = {str(path.relative_to(packet.ROOT)) for path in packet.OUTPUT.parent.glob('*.py')}
    inherited = strict_json.load(packet.ROOT/packet.PRIOR)['source_bindings']
    assert own | set(inherited) | {packet.PRIOR, packet.PACKAGE+'/PROOF.md', packet.PACKAGE+'/README.md'} <= set(stored['source_bindings'])
    for path, digest in stored['source_bindings'].items():
        assert hashlib.sha256((packet.ROOT/path).read_bytes()).hexdigest() == digest
    assert stored['index']['atomic']['reflection_index_antihedgehog'] == 2


def test_wrong_topological_identification_and_scope_promotions_are_rejected(replay):
    for key, value in replay['scope'].items():
        wrong = deepcopy(replay)
        wrong['scope'][key] = not value
        with pytest.raises(ValueError):
            packet.validate(wrong, replay)
    wrong = deepcopy(replay)
    wrong['index']['atomic']['reflection_index_antihedgehog'] = -2
    with pytest.raises(ValueError):
        packet.validate(wrong, replay)
    wrong = deepcopy(replay)
    wrong['localization']['examples'][0]['reference_radius_sites'] = 1
    with pytest.raises(ValueError):
        packet.validate(wrong, replay)
    wrong = deepcopy(replay)
    wrong['all_pass'] = 1
    with pytest.raises(ValueError):
        packet.validate(wrong, replay)


def test_proof_and_source_hash_corruption_is_rejected(replay):
    for name in (packet.PRIOR, packet.PACKAGE+'/PROOF.md', packet.PACKAGE+'/finite.py'):
        wrong = deepcopy(replay)
        del wrong['source_bindings'][name]
        with pytest.raises(ValueError, match='manifest'):
            packet.validate(wrong)
    wrong = deepcopy(replay)
    wrong['source_bindings'][packet.PRIOR] = '0'*64
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(wrong)


def test_parent_custody_without_reexecuting_its_spectrum(monkeypatch):
    def forbidden():
        raise AssertionError('parent proof was reexecuted')
    monkeypatch.setattr(packet.parent_packet, 'report', forbidden)
    monkeypatch.setattr(packet.parent_packet, '_produce', forbidden)
    assert packet.inherited_custody()['every_inherited_source_digest_verified'] is True
    parent = strict_json.load(packet.ROOT/packet.PRIOR)
    parent['source_bindings'][next(iter(parent['source_bindings']))] = '0'*64
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: parent)
    with pytest.raises(ValueError, match='unchanged'):
        packet.inherited_custody()
