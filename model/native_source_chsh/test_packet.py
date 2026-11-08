from copy import deepcopy
import pytest
from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.produce()


def test_retained_owner_receipt_matches_complete_fresh_replay(fresh):
    assert packet.check_payload(strict_json.load(packet.OUTPUT), fresh)


@pytest.mark.parametrize('field', (
    'source_bindings', 'exact_setting_angles', 'same_pre_setting_herald_for_all_settings',
    'abort_outputs_independent_of_settings', 'complete_unconditional_score_certified_above',
    'exact_ideal_complete_score_interval', 'one_setting_complete_halfdiamond_upper',
    'four_setting_score_error_upper', 'actual_stored_score_and_gain_checked_without_normalization',
    'herald_native_program', 'readout_native_program', 'complete_records',
    'resources_per_complete_path', 'scope',
))
def test_every_bound_record_resource_and_scope_field_is_authenticated(fresh, field):
    changed = deepcopy(fresh)
    changed[field] = None
    with pytest.raises(ArithmeticError):
        packet.check_payload(changed, fresh)


def test_payload_types_and_unknown_keys_cannot_be_silently_coerced(fresh):
    for left, right in ((True, 1), (False, 0), ('1', 1)):
        with pytest.raises(ArithmeticError):
            packet.check_payload({'x': left}, {'x': right})
    for value in (1.0, float('nan'), float('inf')):
        with pytest.raises(TypeError):
            packet.canonical({'x': value})
    changed = deepcopy(fresh)
    changed['unapproved_physical_claim'] = True
    with pytest.raises(ArithmeticError):
        packet.check_payload(changed, fresh)


def test_exact_receipt_has_no_platform_dependent_runtime_mass_pins(fresh):
    assert 'complete_CHSH_score_diagnostic' not in fresh
    assert 'complete_CHSH_score' not in fresh
    assert fresh['complete_unconditional_score_certified_above'] == '12/5'
    paths = fresh['source_bindings']
    for name in ('chsh.py', 'packet.py', 'test_chsh.py', 'test_packet.py', 'README.md', 'PROOF.md'):
        assert 'model/native_source_chsh/' + name in paths
    assert 'model/native_coherent_occupation_control/derived_control.py' in paths
    assert 'model/native_program/interpreter.py' in paths


def test_duplicate_fields_are_rejected(tmp_path):
    target = tmp_path / 'receipt.json'
    target.write_text('{"schema":"a","schema":"b"}')
    with pytest.raises(ValueError):
        strict_json.load(target)
