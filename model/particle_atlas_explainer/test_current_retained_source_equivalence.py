"""The additive atlas service cannot weaken the retained source contract."""
import copy

import pytest

from tools import strict_json
from tools import renew_canonical_live_receipts as live
from . import current_retained_source_equivalence as current
from . import retained_source_equivalence as retained


def test_current_navigation_preserves_exact_historical_members_and_guards():
    assert current.RETAINED_NAMES == retained.RETAINED_NAMES
    assert isinstance(current.RETAINED_NAMES, frozenset)
    assert current.load_authority(current.ROOT) == retained.load_authority(current.ROOT)
    value = current.navigation()
    replay = strict_json.load(current.ROOT / retained.REPLAY)
    assert value['retained_physical_ids'] == replay['members']
    assert value['authority_sha256'] == retained.AUTHORITY_SHA256
    assert value['historical_replay_sha256'] == retained.REPLAY_SHA256
    assert value['additional_executed_receipts'] == sorted(current.NAMES)
    old = retained.load_authority()
    assert set(value['live_custody_transitions']) == set(old['allowed_live_receipt_transitions']) | current.NAMES
    assert value['physical_ids_rewritten'] is False
    assert value['physical_specification_changed'] is False
    assert value['retained_authority_bytes_changed'] is False
    assert value['live_producer_receipts_claimed_for_retained_primary_evidence'] is False


def test_retained_module_still_rejects_new_receipts_and_is_never_mutated():
    before = (retained.load_authority, retained.validate_transition, retained.navigation)
    current.service()
    assert before == (retained.load_authority, retained.validate_transition, retained.navigation)
    with pytest.raises(ValueError, match='unreviewed historical source change'):
        retained.historical_input_paths()


def test_complete_historical_input_map_has_original_bytes():
    paths = current.historical_input_paths()
    original = retained.load_authority()['retained_input_sha256']
    assert paths.keys() == original.keys()
    assert all(current.sha(path) == original[name] for name, path in paths.items())


@pytest.mark.parametrize('name', sorted(current.NAMES))
def test_each_executed_transition_rejects_science_or_unreviewed_source_scope(name):
    spec = current.authority()['additional_live_receipt_transitions'][name]
    old = strict_json.load(current.ROOT / spec['retained_path'])
    now = strict_json.load(current.ROOT / name)
    current.compare_transition(name, spec, old, now)
    wrong = copy.deepcopy(now)
    wrong['unreviewed_physical_result'] = True
    with pytest.raises(ValueError):
        current.compare_transition(name, spec, old, wrong)
    wrong = copy.deepcopy(now)
    next(live.source_maps(wrong))[1]['model/unreviewed.py'] = '0' * 64
    with pytest.raises(ValueError):
        current.compare_transition(name, spec, old, wrong)


def test_unreviewed_stale_pin_is_rejected():
    with pytest.raises(ValueError, match='unreviewed stale atlas input'):
        current.authenticate_pins({retained.AUTHORITY: '0' * 64})


def test_current_authentication_exposes_both_authorities_and_actual_owner_record():
    pins = {p: s['retained_sha256'] for p, s in
            current.authority()['additional_live_receipt_transitions'].items()}
    authenticated = current.authenticate_pins(pins)
    assert authenticated[current.AUTHORITY] == current.AUTHORITY_SHA256
    assert authenticated[retained.AUTHORITY] == retained.AUTHORITY_SHA256
    for p in pins:
        assert authenticated[p] == current.sha(current.ROOT / p)
    replay = current.authority()['actual_owner_replay']
    assert authenticated[replay] == current.sha(current.ROOT / replay)


def test_executable_and_live_receipt_tampering_are_rejected(monkeypatch):
    real_sha = current.sha
    with monkeypatch.context() as patch:
        patch.setattr(current, 'sha', lambda p: '0' * 64
                      if str(p).endswith(current.RETAINED_SOURCE) else real_sha(p))
        with pytest.raises(ValueError, match='retained atlas or actual owner evidence changed'):
            current.authority()
    name = sorted(current.NAMES)[0]
    spec = current.authority()['additional_live_receipt_transitions'][name]
    with monkeypatch.context() as patch:
        patch.setattr(current, 'sha', lambda p: '0' * 64
                      if p == current.ROOT / name else real_sha(p))
        with pytest.raises(ValueError, match='executed current atlas receipt'):
            current.validate_current_transition(name, spec)


def test_retained_receipt_values_cannot_acquire_new_physical_claims():
    name = sorted(retained.RETAINED_NAMES)[0]
    value = strict_json.load(current.ROOT / name)
    current.validate_receipt_value(name, value)
    value['new_physical_gate'] = True
    with pytest.raises(ValueError, match='retained proof science, identity or gate changed'):
        current.validate_receipt_value(name, value)
