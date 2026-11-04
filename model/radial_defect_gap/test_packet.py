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
    assert stored['bounds']['examples'][1]['origin_weight_lower'] == '7688/7691'


def test_scope_algebra_and_strong_mass_promotions_are_rejected(replay):
    for key, value in replay['scope'].items():
        wrong = deepcopy(replay)
        wrong['scope'][key] = not value
        with pytest.raises(ValueError):
            packet.validate(wrong)
    for section, key, change in (
            ('bounds', 'exact_kernel_dimension', 3),
            ('bounds', 'sharp_strong_mass_domain', 'm>0'),
            ('bounds', 'finite_volume_convergence_is_quantitative', False),
            ('algebra', 'origin_D0_charged_minus_single_contrast', 'rho_0'),
            ('algebra', 'interacting_ground_or_charge_gap_proved_in_this_package', True)):
        wrong = deepcopy(replay)
        wrong[section][key] = change
        with pytest.raises(ValueError, match='algebra|bounds'):
            packet.validate(wrong)
    wrong = deepcopy(replay)
    wrong['bounds']['examples'][1]['one_body_gap_lower'] = '64'
    with pytest.raises(ValueError, match='bounds'):
        packet.validate(wrong)
    wrong = deepcopy(replay)
    wrong['all_pass'] = 1
    with pytest.raises(ValueError):
        packet.validate(wrong)


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


def test_parent_custody_without_reexecuting_its_index(monkeypatch):
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
