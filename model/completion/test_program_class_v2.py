"""The proposed class preserves source families, evidence and decision scope."""
from copy import deepcopy
import hashlib
import json
import shutil

import pytest

from tools import strict_json
from . import program_class_v2 as contract


def packet():
    return strict_json.load(contract.ROOT/contract.PATH)


def test_proposed_full_target_class_and_retrospective_exposure_are_explicit():
    result = contract.validate(packet())
    assert result['target_observation_rows'] == 71
    assert not result['prospective_claim_allowed']
    assert result['physical_obligations_discharged'] == []
    assert not result['CP00_acceptance_established_by_this_checker']
    assert not result['CP01_acceptance_established_by_this_checker']


def test_actual_continuing_observer_counts_fit_declared_finite_request_budget():
    p = packet()
    receipt = strict_json.load(contract.ROOT/p['observer_acceptance']['source_receipt'])
    actual = receipt['resources']
    ceiling = p['resource_envelope']
    assert actual['total_native_attempts'] == 31458883
    assert actual['total_native_attempts'] < int(ceiling['native_attempts_per_finite_request_max'])
    assert actual['complete_supplied_realization_pulses_upper'] < int(
        ceiling['source_pulses_per_finite_request_max'])
    assert actual['all_committed_native_record_bits_upper'] < int(
        ceiling['archive_bits_per_finite_request_max'])
    # These inequalities establish resource eligibility, not grammar membership.


@pytest.mark.parametrize('damage', [
    'retire_federation', 'change_weights', 'extra_primitive', 'free_clock',
    'false_refinement', 'false_local_sector', 'q3_executed', 'erase_exposure',
    'blind_relabel', 'drop_target', 'stale_source', 'premature_adoption',
    'resource_escape', 'false_field_success',
])
def test_proposal_mutations_are_rejected(damage):
    p = deepcopy(packet())
    if damage == 'retire_federation':
        p['carriers'] = [x for x in p['carriers'] if x['id'] != 'binary_native_federation']
    elif damage == 'change_weights':
        p['laws']['local_alternative_family']['port_weights']['max'] = 1024
    elif damage == 'extra_primitive':
        p['services']['source_word_kinds'].append('ARBITRARY_TARGET_MATRIX')
    elif damage == 'free_clock':
        p['resource_envelope']['host_limits_are_physical_clock'] = True
    elif damage == 'false_refinement':
        p['refinements']['dynamical_intertwining_assumed'] = True
    elif damage == 'false_local_sector':
        next(x for x in p['carriers'] if x['id'] == 'binary_native_federation')[
            'local_total_two_tensor_invariance_claimed'] = True
    elif damage == 'q3_executed':
        p['observer_acceptance']['q3']['historical_status'] = 'EXECUTED'
    elif damage == 'erase_exposure':
        p['exposure']['source_hashes'].pop(p['exposure']['registers'][0])
    elif damage == 'blind_relabel':
        p['targets']['prospective_arm']['prospective_claim_allowed'] = True
    elif damage == 'drop_target':
        p['hypothesis']['mandatory_obligations'].pop()
    elif damage == 'stale_source':
        p['exposure']['source_hashes'][next(iter(p['exposure']['source_hashes']))] = '0'*64
    elif damage == 'premature_adoption':
        p['adoption']['status'] = 'ADOPTED_AFTER_REVIEW'
        p['adoption']['independent_review'] = None
    elif damage == 'resource_escape':
        p['resource_envelope']['repeat_count_max'] *= 2
    else:
        p['hypothesis']['imported_field_action_counts_as_physical_success'] = True
    with pytest.raises(ValueError):
        contract.validate(p)


def test_exact_comparison_rule_retains_inconclusive_and_calibrated_cases():
    flags = dict(complete_error=True, same_quantity=True, independent_output=True)
    assert contract.compare_intervals(('1', '2'), ('0', '3'), **flags) == 'RETROSPECTIVE_INTERVAL_COMPATIBILITY'
    assert contract.compare_intervals(('1', '2'), ('2', '3'), **flags) == 'INCONCLUSIVE'
    assert contract.compare_intervals(('1', '2'), ('3', '4'), **flags) == 'CANDIDATE_DISCREPANCY'
    for premise in flags:
        missing = flags | {premise: False}
        assert contract.compare_intervals(('1', '2'), ('0', '3'), **missing) == 'INCONCLUSIVE'
    for bad in ((1.0, 2), (True, 2), ('2', '1')):
        with pytest.raises(ValueError):
            contract.compare_intervals(bad, ('0', '3'), **flags)


@pytest.mark.parametrize('damage', ['charged_menu', 'calibration', 'native_law',
                                  'review_body', 'review_map_drift', 'missing_map'])
def test_pass_review_cannot_be_reused_after_contract_change(tmp_path, damage):
    p = deepcopy(packet())
    # Bind the actual interpreter when available; this synthetic review tests
    # custody logic and asserts no scientific acceptance of its own.
    if 'interpreter_binding' not in p:
        names = ['model/native_program/interpreter.py', 'model/native_program/members.py',
                 'model/native_program/receipt.py', 'model/native_program/results.json',
                 'model/native_program/test_interpreter.py']
        p['interpreter_binding'] = {
            'source_hashes': {name: hashlib.sha256((contract.ROOT/name).read_bytes()).hexdigest()
                              for name in names},
            'member_receipt': 'model/native_program/results.json',
            'independent_check': 'model/native_program/test_interpreter.py'}
    paths = set(p['exposure']['source_hashes']) | set(p['interpreter_binding']['source_hashes'])
    paths |= set(contract.EXPOSURE_SNAPSHOTS.values())
    paths |= set(p.get('review_implementation_bindings', {}))
    map_path = 'model/completion/finite_attachment_map.json'
    paths.add(map_path)
    for name in paths:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(contract.ROOT/name, target)
    reviewed_digest = contract.reviewed_contract_digest(p)
    review_sources = {name: hashlib.sha256((tmp_path/name).read_bytes()).hexdigest()
                      for name in (map_path, 'model/native_program/results.json')}
    review = {'verdict': 'PASS', 'reviewed_contract_digest': reviewed_digest,
              'source_bindings': review_sources,
              'scope': 'synthetic test fixture; not a scientific approval'}
    name = 'review.json'
    target = tmp_path/name
    target.write_text(json.dumps(review))
    p['adoption'] = {'status': 'ADOPTED_AFTER_REVIEW', 'changes_current_registry': False,
                     'independent_review': {'verdict': 'PASS', 'path': name,
                         'sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                         'reviewed_contract_digest': reviewed_digest}}
    assert contract.validate(p, root=tmp_path)['adoption_blockers'] == []
    if damage == 'charged_menu':
        p['targets']['charged_spectrum']['quantities'] = ['ye']
    elif damage == 'calibration':
        p['calibrations']['allowed_fit_parameters'] = 'one arbitrary scale per species'
    elif damage == 'native_law':
        p['laws']['fixed'] = 'arbitrary target Hamiltonian'
    elif damage == 'review_body':
        review['reviewed_contract_digest'] = '0'*64
        target.write_text(json.dumps(review))
        p['adoption']['independent_review']['sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
    elif damage == 'review_map_drift':
        (tmp_path/map_path).write_text('{}')
    else:
        del review['source_bindings'][map_path]
        target.write_text(json.dumps(review))
        p['adoption']['independent_review']['sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match='exact class contract|stale reviewed source|actual theorem map'):
        contract.validate(p, root=tmp_path)


@pytest.mark.parametrize('name', list(contract.EXPOSURE_SNAPSHOTS))
def test_exposure_snapshot_uses_original_pin_and_rejects_damage(tmp_path, name):
    p = packet()
    frozen_name = contract.EXPOSURE_SNAPSHOTS[name]
    target = tmp_path/frozen_name
    target.parent.mkdir(parents=True, exist_ok=True)
    original = (contract.ROOT/frozen_name).read_bytes()
    target.write_bytes(original)
    # A live register is not an alternate source for the immutable exposure.
    live = tmp_path/name
    live.parent.mkdir(parents=True, exist_ok=True)
    live.write_text('{"commentary": "different live contents"}')
    expected = p['exposure']['source_hashes'][name]
    contract.bound_exposure(tmp_path, name, expected)
    target.write_bytes(original+b' ')
    with pytest.raises(ValueError, match='stale reviewed source'):
        contract.bound_exposure(tmp_path, name, expected)
    target.unlink()
    live.write_bytes(original)
    with pytest.raises(ValueError, match='missing or external bound file'):
        contract.bound_exposure(tmp_path, name, expected)


def test_live_target_replacement_is_rejected(tmp_path):
    p = packet()
    paths = set(p['exposure']['source_hashes']) | set(contract.EXPOSURE_SNAPSHOTS.values())
    for name in paths:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(contract.ROOT/name, target)
    target = tmp_path/'tracking/observation_ledger.json'
    live = json.loads(target.read_text())
    live['rows'][0]['target'] = 'a different physical prediction'
    target.write_text(json.dumps(live))
    with pytest.raises(ValueError, match='physical targets differ from frozen exposure'):
        contract.validate(p, root=tmp_path)
