"""CI partitions execution, never the canonical receipt or artifact coverage."""
from collections import Counter
import gzip
import hashlib
import json
import subprocess
import sys

import pytest
import yaml

from . import verify


def flatten(groups):
    return tuple(replay for group in groups for replay in group)


@pytest.mark.parametrize('count', range(1, len(verify.replay_groups())+1))
def test_shards_are_exhaustive_disjoint_and_keep_empirical_worker_atomic(count):
    shards = [verify.selected_replay_groups(index, count) for index in range(count)]
    assert Counter(replay for groups in shards for replay in flatten(groups)) == Counter(verify.REPLAYS)
    assert sum(group == verify.EMPIRICAL_REPLAYS for groups in shards for group in groups) == 1
    assert all(group == verify.EMPIRICAL_REPLAYS or len(group) == 1 for groups in shards for group in groups)
    if count > 1:
        assert shards[0] == (verify.EMPIRICAL_REPLAYS,)
    if count >= 3:
        assert shards[1] == ((verify.OCCUPATION_TRANSFER_REPLAY,),)
    assert all(shards)


def test_default_and_single_shard_keep_every_original_replay_in_order():
    assert flatten(verify.selected_replay_groups()) == verify.REPLAYS
    assert flatten(verify.selected_replay_groups(0, 1)) == verify.REPLAYS


def test_new_replays_are_included_without_a_parallel_allowlist():
    extra = ('model.example_new_certificate', None, 'model/example_new_certificate/results.json')
    declaration = (extra,)+verify.REPLAYS
    shards = [verify.selected_replay_groups(index, 3, replays=declaration) for index in range(3)]
    assert Counter(replay for groups in shards for replay in flatten(groups)) == Counter(declaration)
    assert extra in flatten(shards[2])


def test_two_shard_fallback_keeps_transfer_and_all_other_nonempirical_work():
    second = flatten(verify.selected_replay_groups(1, 2))
    assert verify.OCCUPATION_TRANSFER_REPLAY in second
    assert second == tuple(replay for replay in verify.REPLAYS if replay not in verify.EMPIRICAL_REPLAYS)


@pytest.mark.parametrize('count', (4, 6))
def test_shards_keep_both_expensive_workers_isolated_and_order_remaining_work(count):
    remaining = tuple(group for group in verify.replay_groups()
                      if group not in (verify.EMPIRICAL_REPLAYS, (verify.OCCUPATION_TRANSFER_REPLAY,)))
    for index in range(2, count):
        assert verify.selected_replay_groups(index, count) == remaining[index-2::count-2]


def test_every_new_research_receipt_has_a_canonical_replay():
    assert {
        ('model.band_clock_budget', None, 'model/band_clock_budget/results.json'),
        ('model.common_source.native_repair_mass', None, 'model/common_source/native_repair_mass_results.json'),
        ('model.common_source.native_gauge_matter_transport', None, 'model/common_source/native_gauge_matter_transport_results.json'),
        ('model.common_source.native_sm_current', None, 'model/common_source/native_sm_current_results.json'),
        ('model.common_source.native_matter_attachment', None, 'model/common_source/native_matter_attachment_results.json'),
        ('model.common_source.native_particle_records', None, 'model/common_source/native_particle_records_results.json'),
        ('model.common_source.native_gauge_object', None, 'model/common_source/native_gauge_object_results.json'),
        ('model.common_source.native_persistent_response', None, 'model/common_source/native_persistent_response_results.json'),
        ('model.common_source.native_virtual_exchange', None, 'model/common_source/native_virtual_exchange_results.json'),
        ('model.common_source.native_code_transport', None, 'model/common_source/native_code_transport_results.json'),
        ('model.common_source.native_charged_spectroscopy', None, 'model/common_source/native_charged_spectroscopy_results.json'),
        ('model.native_equilibration.common_variation', None, 'model/native_equilibration/variation_results.json'),
        ('model.native_interface.tomography', None, 'model/native_interface/tomography_results.json'),
        ('model.native_modular_reference.cap_local_clock', None, 'model/native_modular_reference/cap_local_clock_results.json'),
        ('model.occupation_control.transfer', None, 'model/occupation_control/transfer_results.json.gz'),
        ('model.common_source.replay', 'mass_response', 'model/common_source/mass_response/results.json'),
        ('model.common_source.replay', 'particle_identification', 'model/common_source/particle_identification/results.json'),
        ('model.common_source.replay', 'chiral_candidate', 'model/common_source/chiral_candidate/results.json'),
        ('model.common_source.hierarchy_diagnosis.audit', None, 'model/common_source/hierarchy_diagnosis/results.json'),
        ('model.common_source.full_repair_response.run', None, 'model/common_source/full_repair_response/results.json'),
        ('model.common_source.native_staggered_field', None, 'model/common_source/native_staggered_field_results.json'),
        ('model.common_source.native_field_history', None, 'model/common_source/native_field_history_results.json'),
        ('model.common_source.native_charged_spectrum', None, 'model/common_source/native_charged_spectrum_results.json'),
        ('model.native_equilibration.wedge_modular', None, 'model/native_equilibration/wedge_modular_results.json'),
        ('model.screened_recovery.experiment', None, 'model/screened_recovery/screened_recovery_results.json'),
        ('model.conserved_sector_action.sector_action', None, 'model/conserved_sector_action/sector_action_results.json'),
        ('model.clock_defect_energy.shared_event_clock.run', None, 'model/clock_defect_energy/shared_event_clock/results.json'),
        ('model.screened_transport.transport', None, 'model/screened_transport/transport_results.json'),
        ('model.exchange_sector.statistics_packet', None, 'model/exchange_sector/statistics_results.json'),
        ('model.particle_screening.mobile_packet', None, 'model/particle_screening/mobile_results.json'),
    }.issubset(verify.REPLAYS)


def test_completion_candidate_receipts_have_replay_custody_and_workflow_coverage():
    statistics_receipt = 'model/exchange_sector/statistics_results.json'
    assert ('model.exchange_sector.statistics_packet', None, statistics_receipt) in verify.REPLAYS
    statistics_sources = json.loads((verify.ROOT/statistics_receipt).read_text())['source_sha256']
    for name in ('statistics.py', 'statistics_service.py', 'statistics_topology.py',
                 'statistics_packet.py', 'statistics_contract.json'):
        assert 'model/exchange_sector/'+name in statistics_sources
    packets = {
        'particle_screening': ('screening', 'screening_results.json', 'campaign.json'),
        'active_code_energy': ('energy_bridge', 'energy_bridge_results.json', 'campaign.json'),
        'clock_defect_energy': ('clock_count_inference', 'clock_count_inference_results.json',
                                'clock_count_inference_campaign.json'),
        'local_repair_family': ('observer', 'observer_results.json', 'observer_campaign.json'),
    }
    for packet, (module, receipt_name, campaign_name) in packets.items():
        prefix = f'model/{packet}/'
        receipt_path = prefix + receipt_name
        assert (f'model.{packet}.{module}', None, receipt_path) in verify.REPLAYS
        receipt = json.loads((verify.ROOT/receipt_path).read_text())
        sources = receipt['source_sha256']
        assert prefix + module + '.py' in sources
        assert prefix + campaign_name in sources
    successors = {
        'screened_recovery': ('experiment', 'screened_recovery_results.json'),
        'conserved_sector_action': ('sector_action', 'sector_action_results.json'),
        'screened_transport': ('transport', 'transport_results.json'),
    }
    for packet, (module, receipt_name) in successors.items():
        prefix = f'model/{packet}/'
        receipt_path = prefix + receipt_name
        assert (f'model.{packet}.{module}', None, receipt_path) in verify.REPLAYS
        receipt = json.loads((verify.ROOT/receipt_path).read_text())
        sources = receipt['source_sha256']
        assert prefix + module + '.py' in sources
        assert prefix + 'campaign.json' in sources
    clock_prefix = 'model/clock_defect_energy/shared_event_clock/'
    record_receipt = 'model/local_repair_family/record_clock_results.json'
    assert ('model.local_repair_family.record_clock', None, record_receipt) in verify.REPLAYS
    record_sources = json.loads((verify.ROOT/record_receipt).read_text())['source_sha256']
    for name in ('record_clock.py', 'record_clock_campaign.json', 'inference.py',
                 'acquisition.py', 'learned_flow.py'):
        assert 'model/local_repair_family/'+name in record_sources
    clock_receipt = clock_prefix + 'results.json'
    assert ('model.clock_defect_energy.shared_event_clock.run', None, clock_receipt) in verify.REPLAYS
    clock_sources = json.loads((verify.ROOT/clock_receipt).read_text())['source_sha256']
    assert clock_prefix + 'run.py' in clock_sources
    assert clock_prefix + 'campaign.json' in clock_sources
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/completion-program.yml').read_text())
    events = workflow.get('on', workflow.get(True))
    for event in ('push', 'pull_request'):
        paths = events[event]['paths']
        assert 'model/particle_screening/**' in paths
        assert 'model/active_code_energy/**' in paths
        assert 'model/screened_recovery/**' in paths
        assert 'model/conserved_sector_action/**' in paths
        assert 'model/screened_transport/**' in paths
        assert 'model/clock_defect_energy/**' in paths
        assert 'model/local_repair_family/**' in paths
    commands = [step['run'] for step in workflow['jobs']['program']['steps'] if 'run' in step]
    pytest_commands = [command for command in commands if 'python -m pytest' in command]
    assert pytest_commands
    pytest_command = '\n'.join(pytest_commands)
    assert 'model/particle_screening' in pytest_command
    assert 'model/active_code_energy' in pytest_command
    assert 'model/screened_recovery' in pytest_command
    assert 'model/conserved_sector_action' in pytest_command
    assert 'model/screened_transport' in pytest_command
    assert 'model/clock_defect_energy' in pytest_command
    assert 'model/local_repair_family' in pytest_command
    for packet in ('moving_defect', 'record_action', 'particle_obstruction', 'interacting_field'):
        assert f'model/{packet}' in pytest_command
        for event in ('push', 'pull_request'):
            assert f'model/{packet}/**' in events[event]['paths']
    for command in ('python -m model.moving_defect.experiment --check',
                    'python -m model.record_action.action --check',
                    'python -m model.particle_obstruction.classification --check',
                    'python -m model.interacting_field --check model/interacting_field/results.json'):
        assert command in commands
    for packet in ('bound_pair', 'particle_clock', 'field_regions', 'field_stress', 'field_tastes', 'biased_reset', 'support_refinement'):
        assert f'model/{packet}' in pytest_command
        receipt_path=f'model/{packet}/results.json'
        assert (f'model.{packet}',None,receipt_path) in verify.REPLAYS
        assert f'python -m model.{packet} --check {receipt_path}' in commands
        for event in ('push','pull_request'):
            assert f'model/{packet}/**' in events[event]['paths']
        receipt=json.loads((verify.ROOT/receipt_path).read_text())
        sources=receipt['source_sha256' if packet in ('particle_clock', 'support_refinement') else 'source_bindings']
        assert sources and all((verify.ROOT/path).is_file() for path in sources)


def test_full_band_and_native_localization_receipts_are_mandatory():
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/completion-program.yml').read_text())
    commands = [step['run'] for step in workflow['jobs']['program']['steps'] if 'run' in step]
    events = workflow.get('on', workflow.get(True))
    for package, module in [('band_clock', 'model.band_clock'),
                            ('native_localization', 'model.native_localization.experiment'),
                            ('native_program', 'model.native_program.receipt')]:
        receipt = f'model/{package}/results.json'
        assert (module, None, receipt) in verify.REPLAYS
        assert f'python -m {module} --check {receipt}' in commands
        assert any('python -m pytest' in c and f'model/{package}' in c for c in commands)
        for event in ('push', 'pull_request'):
            assert f'model/{package}/**' in events[event]['paths']


@pytest.mark.parametrize('index', range(6))
def test_missing_dedicated_transfer_fails_on_every_ci_shard(index):
    declaration = tuple(replay for replay in verify.REPLAYS if replay != verify.OCCUPATION_TRANSFER_REPLAY)
    with pytest.raises(ValueError, match='occupation-transfer'):
        verify.selected_replay_groups(index, 6, replays=declaration)


@pytest.mark.parametrize('index,count', [(None,3),(0,None),(-1,3),(3,3),(0,0),
                                       (0,10000),(True,3),(0,True),(0.,3),(0,3.)])
def test_invalid_shard_arguments_fail_closed(index, count):
    with pytest.raises(ValueError):
        verify.selected_replay_groups(index, count)


def test_malformed_empirical_group_cannot_silently_run_split():
    first = verify.REPLAYS.index(verify.EMPIRICAL_REPLAYS[0])
    changed = list(verify.REPLAYS)
    changed[first], changed[first+1] = changed[first+1], changed[first]
    with pytest.raises(ValueError, match='contiguous'):
        verify.replay_groups(changed)
    with pytest.raises(ValueError):
        verify.replay_groups(tuple(x for x in verify.REPLAYS if x not in verify.EMPIRICAL_REPLAYS))
    with pytest.raises(ValueError, match='unique'):
        verify.replay_groups(verify.REPLAYS+(verify.REPLAYS[0],))


@pytest.mark.parametrize('index,count', [(None,None),(0,1),(0,2),(1,2),
                                       (0,4),(1,4),(2,4),(3,4),
                                       (0,6),(1,6),(2,6),(3,6),(4,6),(5,6)])
def test_default_and_every_shard_hashes_all_artifacts_and_executes_complete_groups(monkeypatch, index, count):
    argv = ['verify'] if count is None else ['verify', '--shard-index', str(index), '--shard-count', str(count)]
    monkeypatch.setattr(sys, 'argv', argv)
    hashed, expanded = [], []
    def hashes():
        hashed.append(True)
        return 534
    def run(argv, **kwargs):
        assert kwargs['check'] is True
        assert kwargs['cwd'] == verify.ROOT
        if argv == [sys.executable, '-m', 'model.replay_empirical']:
            expanded.extend(verify.EMPIRICAL_REPLAYS)
        else:
            position = argv.index('--check')
            expanded.append((argv[2], None if position == 3 else argv[3], argv[position+1]))
    monkeypatch.setattr(verify, 'verify_artifacts', hashes)
    monkeypatch.setattr(verify.subprocess, 'run', run)
    verify.main()
    assert hashed == [True]
    assert tuple(expanded) == flatten(verify.selected_replay_groups(index, count))


def test_invalid_cli_shard_stops_before_hashes_or_execution(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['verify', '--shard-index', '3', '--shard-count', '3'])
    monkeypatch.setattr(verify, 'verify_artifacts', lambda: pytest.fail('invalid arguments reached hashes'))
    with pytest.raises(SystemExit) as error:
        verify.main()
    assert error.value.code == 2


def test_hash_only_shard_checks_all_artifacts_without_replays(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['verify', '--hashes-only', '--shard-index', '2', '--shard-count', '3'])
    calls = []
    monkeypatch.setattr(verify, 'verify_artifacts', lambda: calls.append('all') or 534)
    monkeypatch.setattr(verify.subprocess, 'run', lambda *a, **k: pytest.fail('hash-only executed a replay'))
    verify.main()
    assert calls == ['all']


@pytest.mark.parametrize('field', ['sources', 'source_sha256', 'source_bindings',
    'input_sha256', 'source_binding_sha256', 'source_control_bindings', 'source_pins',
    'implementation_sha256', 'source_interface_hashes'])
@pytest.mark.parametrize('suffix', ['.json', '.json.gz'])
def test_hash_only_rejects_stale_inner_source_even_with_current_outer_manifest(
        monkeypatch, tmp_path, field, suffix):
    packet = tmp_path/'model/example'
    packet.mkdir(parents=True)
    source = packet/'source.py'
    source.write_text('value = 1\n', encoding='utf-8')
    retained = hashlib.sha256(source.read_bytes()).hexdigest()
    source_key = {'implementation_sha256': 'source.py',
                  'source_interface_hashes': 'example/source.py'}.get(field, 'model/example/source.py')
    receipt = packet/f'results{suffix}'
    def write(value):
        raw = json.dumps(value).encode('utf-8')
        receipt.write_bytes(gzip.compress(raw, mtime=0) if suffix.endswith('.gz') else raw)
    write({'nested': [{field: {source_key: retained}}]})
    source.write_text('value = 2\n', encoding='utf-8')
    # Both files have current outer hashes; the receipt's dependency binding
    # is stale and must fail before any scientific replay.
    manifest = {'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (source, receipt)}}
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'PACKETS', ('model/example',))
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, f'model/example/results{suffix}'),))
    monkeypatch.setattr(verify, 'selected_replay_groups', lambda *args: ())
    monkeypatch.setattr(sys, 'argv', ['verify', '--hashes-only'])
    monkeypatch.setattr(verify.subprocess, 'run', lambda *a, **k: pytest.fail('a replay executed'))
    with pytest.raises(ValueError, match=f'root.nested\\[0\\].{field}'):
        verify.main()
    write({'nested': [{field: {source_key: hashlib.sha256(source.read_bytes()).hexdigest()}}]})
    manifest['files'][receipt.name] = hashlib.sha256(receipt.read_bytes()).hexdigest()
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    verify.main()


def test_implementation_bindings_resolve_package_parent_without_root_fallback(monkeypatch, tmp_path):
    folder = tmp_path/'model/example'
    folder.mkdir(parents=True)
    source = folder.parent/'source.py'
    source.write_text('actual implementation\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = folder/'results.json'
    receipt.write_text(json.dumps({'implementation_sha256': {'../source.py': digest}}),
                       encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json'),))
    verify.verify_receipt_source_bindings()
    source.write_text('changed implementation\n', encoding='utf-8')
    # A matching file at another plausible root cannot hide the local drift.
    (tmp_path/'source.py').write_text('actual implementation\n', encoding='utf-8')
    with pytest.raises(ValueError, match='digest mismatch'):
        verify.verify_receipt_source_bindings()


@pytest.mark.parametrize('binding', [{}, [], False, {'missing.py': '0'*64},
    {'../../../outside.py': '0'*64}, {'https://example.org/source.py': '0'*64}])
@pytest.mark.parametrize('field', ['implementation_sha256', 'source_interface_hashes',
    'source_sha256', 'source_bindings', 'input_sha256', 'source_binding_sha256',
    'source_control_bindings', 'source_pins'])
def test_local_bindings_cannot_skip_malformed_or_unresolved_inputs(
        monkeypatch, tmp_path, binding, field):
    folder = tmp_path/'model/example'
    folder.mkdir(parents=True)
    receipt = folder/'results.json'
    receipt.write_text(json.dumps({field: binding}), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json'),))
    with pytest.raises(ValueError, match=field):
        verify.verify_receipt_source_bindings()


def test_endpoint_pins_are_checked_at_the_declared_nested_location(monkeypatch, tmp_path):
    source = tmp_path/'source.py'
    source.write_text('original\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = tmp_path/'results.json'
    receipt.write_text(json.dumps({'retrospective_endpoint_audit': {
        'source': {'pins': {'source.py': digest}}}, 'unrelated': {'pins': {'alias': 'owner'}}}),
        encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'results.json'),))
    verify.verify_receipt_source_bindings()
    source.write_text('changed\n', encoding='utf-8')
    with pytest.raises(ValueError, match='root.retrospective_endpoint_audit.source.pins'):
        verify.verify_receipt_source_bindings()


def test_live_source_gate_leaves_external_aliases_and_frozen_captures_to_owners(
        monkeypatch, tmp_path):
    source = tmp_path/'source.py'
    source.write_text('value = 1\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = tmp_path/'results.json'
    receipt.write_text(json.dumps({'sources': {'source.py': digest,
        'external/raw.csv': 'external custody', 'https://example.org/data': 'external',
        '../outside.py': 'not repository-relative'},
        'external_interface_hashes': {'unknown_alias.py': 'owning verifier'}}), encoding='utf-8')
    (tmp_path/'frozen_capture.json').write_text(
        json.dumps({'source_sha256': {'source.py': 'capture-time hash'}}), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'results.json'),))
    verify.verify_receipt_source_bindings()


@pytest.mark.parametrize('suffix', ['.json', '.json.gz'])
def test_receipt_source_gate_rejects_duplicate_keys_before_hash_selection(monkeypatch, tmp_path, suffix):
    receipt = tmp_path/f'results{suffix}'
    raw = b'{"prediction":{"source_bindings":{"source.py":"stale","source.py":"replacement"}}}'
    receipt.write_bytes(gzip.compress(raw, mtime=0) if suffix.endswith('.gz') else raw)
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, receipt.name),))
    with pytest.raises(ValueError, match='duplicate'):
        verify.verify_receipt_source_bindings()


def test_compressed_receipt_requires_outer_manifest_coverage(monkeypatch, tmp_path):
    packet = tmp_path/'model/example'
    packet.mkdir(parents=True)
    source = packet/'source.py'
    source.write_text('value = 1\n', encoding='utf-8')
    receipt = packet/'results.json.gz'
    receipt.write_bytes(gzip.compress(b'{}', mtime=0))
    manifest = {'files': {'source.py': hashlib.sha256(source.read_bytes()).hexdigest()}}
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'PACKETS', ('model/example',))
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json.gz'),))
    with pytest.raises(ValueError, match='Unbound artifact'):
        verify.verify_artifacts()
    manifest['files'][receipt.name] = hashlib.sha256(receipt.read_bytes()).hexdigest()
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    assert verify.verify_artifacts() == 2


def test_ci_covers_all_six_shards_both_platforms_and_complete_regression():
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/observer-model.yml').read_text())
    observer, regression = (workflow['jobs'][key] for key in ('observer', 'regression'))
    assert observer['strategy']['matrix']['os'] == ['ubuntu-latest', 'windows-latest']
    assert observer['strategy']['matrix']['shard'] == [0,1,2,3,4,5]
    assert observer['timeout-minutes'] == (
        "${{ matrix.os == 'windows-latest' && matrix.shard == 1 && 120 || 90 }}")
    assert any(step.get('run') == 'python -m model.verify --shard-index ${{ matrix.shard }} --shard-count 6'
               for step in observer['steps'])
    scalar = [step for step in observer['steps'] if 'verify_source_scalar_time_refinement.py' in step.get('run','')]
    assert len(scalar) == 1 and scalar[0]['if'] == 'matrix.shard == 2'
    assert regression['strategy']['matrix']['os'] == ['ubuntu-latest', 'windows-latest']
    assert regression['strategy']['matrix']['shard'] == [0,1,2]
    commands = [step for step in regression['steps'] if 'run' in step]
    pytest_steps = [step for step in commands if 'python -m pytest' in step['run']]
    assert len(pytest_steps) == 1 and 'if' not in pytest_steps[0]
    assert pytest_steps[0]['run'] == ('python -m pytest -q -W error -p model.pytest_shards '
        '--model-shard-index ${{ matrix.shard }} --model-shard-count 3 model')
    extra = [step for step in commands if 'model.common_source.repair_field_attachment' in step['run']
             or 'native_staggered_field' in step['run'] or 'observer_particles/run.py' in step['run']]
    assert len(extra) == 4
    assert all(step.get('if') == 'matrix.shard == 0' for step in extra)


@pytest.mark.parametrize("matches", [True, False])
def test_ci_fetches_canonical_release_and_checks_the_recorded_commit(monkeypatch, tmp_path, matches):
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/observer-model.yml').read_text())
    steps = workflow['jobs']['observer']['steps']
    fetch = next(s for s in steps if s.get('name') == "Fetch the model's recorded base release")
    assert fetch['run'] == ('git fetch --no-tags --depth=1 '
        'https://github.com/FloatingPragma/observer-patch-holography.git '
        'refs/tags/r2053:refs/tags/r2053')
    check = next(s for s in steps if s.get('name') == 'Verify the recorded base release commit')
    assert check['shell'] == 'python'
    receipt = tmp_path/'model/common_source/mass_response/results.json'
    receipt.parent.mkdir(parents=True)
    expected = 'a'*40
    receipt.write_text(json.dumps({'base_release_commit': expected}), encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    def rev_parse(command, *, text):
        assert command == ['git', 'rev-parse', 'r2053^{commit}'] and text is True
        return (expected if matches else 'b'*40) + '\n'

    monkeypatch.setattr(subprocess, 'check_output', rev_parse)
    if matches:
        exec(compile(check['run'], '<ci-base-release-check>', 'exec'), {})
    else:
        with pytest.raises(SystemExit, match='Recorded base release mismatch'):
            exec(compile(check['run'], '<ci-base-release-check>', 'exec'), {})


def test_particle_services_and_actual_data_replay_remain_mandatory():
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/completion-program.yml').read_text())
    events = workflow.get('on', workflow.get(True))
    commands = [step['run'] for step in workflow['jobs']['particle-services']['steps'] if 'run' in step]
    for package in ('reference_recovery', 'particle_response', 'contact_scattering', 'particle_comparison', 'particle_transport'):
        for event in ('push', 'pull_request'):
            assert f'model/{package}/**' in events[event]['paths']
        assert any('python -m pytest' in c and f'model/{package}' in c for c in commands)
    module, receipt = 'model.particle_comparison.compare', 'model/particle_comparison/results.json'
    assert (module, None, receipt) in verify.REPLAYS
    assert f'python -m {module} --check {receipt}' in commands
    data = json.loads((verify.ROOT/receipt).read_text())
    assert data['CP11']['eligible_physical_comparisons_added'] == 0
    assert not data['CP11']['OPH_universe_membership_confirmed']
    assert data['custody']['input_sha256']


def test_same_candidate_transport_is_replayed_without_promoting_physical_identity():
    module, receipt = 'model.particle_transport.packet', 'model/particle_transport/results.json'
    assert (module, None, receipt) in verify.REPLAYS
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/completion-program.yml').read_text())
    commands = [step['run'] for step in workflow['jobs']['particle-services']['steps'] if 'run' in step]
    assert f'python -m {module} --check {receipt}' in commands
    data = json.loads((verify.ROOT/receipt).read_text())
    assert not data['physical_particle_identification']
    assert not data['integrated_gate_accepted']
    assert not data['native_topological_protection_established']
    assert not data['failed_dispersion_benchmark_repaired']
    assert data['physical_comparison_count_added'] == 0
