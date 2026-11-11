"""Mutation controls for exact case custody and independent current validation."""
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from . import live
from .runner import (Execution, RetainedExecutionError, case_failure,
                     validate_imported_sources, validate_reports)


def reports(*names):
    return {'exitstatus': 0, 'reports': [
        {'nodeid': name, 'when': phase, 'outcome': 'passed', 'wasxfail': '', 'longrepr': ''}
        for name in names for phase in ('setup', 'call', 'teardown')]}


def test_complete_inventory_and_current_producer_coverage():
    value = live.authority()
    assert len(value['test_files']) == 81
    nodes = [n for group in value['nodeids_by_file'].values() for n in group]
    assert len(nodes) == len(set(nodes)) == 1005
    assert all(n.split('::')[0] == p for p, group in value['nodeids_by_file'].items() for n in group)
    assert value['calibration_bound_input_count'] == 1424
    assert value['occupation_bound_input_count'] == 260
    assert len(value['current_changes']) == 9
    assert all(row['path'].endswith('.json') for row in value['current_changes'])
    assert value['allowed_skips'] == []


def test_archived_full_execution_has_real_complete_phase_records_and_exact_source_scope():
    value=live.authority()
    folder=live.ROOT/'model/retained_identity_execution/verification'
    assert live.sha(folder/'authority_before_documentary_handoff.json')=='c6e188a4fc996e7218bbdc20c50e4cf0389c7c18f324295979d483553910f7ca'
    execution=json.loads((folder/'activated_1005_execution.json').read_text())
    path=folder/'activated_1005_case_reports.json'
    assert live.sha(path)==execution['case_journal_sha256']
    checked=validate_reports(execution['selected_nodeids'],json.loads(path.read_text()))
    assert len(checked)==1005 and all(case_failure(rows) is None for rows in checked.values())
    assert execution['observed_imported_python_files_unchanged']==644
    assert execution['outcomes_reused_from_receipts'] is False
    assert execution['case_properties']['current_physics_verified_by_historical_pass'] is False
    assert set(execution['selected_nodeids'])==set(value['baseline_collection_order'])


def test_completion_routing_executes_all_old_suites_and_keeps_current_validators():
    from tools import check_retained_source_records_custody as retained
    workflow=yaml.safe_load((live.ROOT/'.github/workflows/completion-program.yml').read_text())
    jobs=workflow['jobs']
    steps=jobs['program-contract']['steps']
    checkout=next(step for step in steps if step.get('uses','').startswith('actions/checkout@'))
    assert checkout['with']['fetch-depth']==0
    commands=[step['run'] for step in steps if 'run' in step]
    assert commands.count('python tools/check_retained_source_records_custody.py')==1
    assert commands.count('python -m model.documentary_release_successor --check')==1
    assert not any('python -m model.completion.source_records_' in command for command in commands)
    assert {'tools/test_retained_custody_routing.py','tools/test_retained_field_custody.py'}<=set(retained.TESTS)
    current=next(command for command in commands if 'tools/test_retained_source_records_custody.py' in command)
    assert 'tools/test_diagnostic_parent_custody.py' in current
    assert 'tools/test_numeric_evidence_medium.py' in current
    all_commands=[step['run'] for job in jobs.values() for step in job['steps'] if 'run' in step]
    for script in ('check_retained_completion_custody.py','check_retained_field_custody.py'):
        assert all_commands.count('python tools/'+script)==1


def test_actual_full_collection_preserves_all_1005_original_node_ids():
    value = live.authority()
    expected = value['baseline_collection_order']
    result = subprocess.run([sys.executable, '-B', '-W', 'error', '-m', 'pytest',
                             '--collect-only', '-q', *value['test_files']],
                            cwd=live.ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    actual = [line.strip() for line in result.stdout.splitlines()
              if line.startswith('model/') and '::' in line]
    assert actual == expected
    assert len(actual) == len(set(actual)) == 1005


@pytest.mark.parametrize('mutation', ['drop_case', 'drop_phase', 'duplicate_case', 'duplicate_phase',
                                      'extra_case', 'extra_phase', 'bad_outcome', 'bad_exit'])
def test_missing_or_fabricated_outcomes_fail_closed(mutation):
    value = reports('a', 'b')
    nodes = ['a', 'b']
    if mutation == 'drop_case': value['reports'] = value['reports'][:3]
    if mutation == 'drop_phase': value['reports'].pop()
    if mutation == 'duplicate_case': nodes.append('a')
    if mutation == 'duplicate_phase': value['reports'].append(copy.deepcopy(value['reports'][0]))
    if mutation == 'extra_case': value['reports'][0]['nodeid'] = 'unselected'
    if mutation == 'extra_phase': value['reports'][0]['when'] = 'cached'
    if mutation == 'bad_outcome': value['reports'][0]['outcome'] = 'cached_pass'
    if mutation == 'bad_exit': value['exitstatus'] = 2
    with pytest.raises(RetainedExecutionError): validate_reports(nodes, value)


@pytest.mark.parametrize('outcome,xfail', [('failed', ''), ('skipped', ''), ('passed', 'expected failure')])
def test_actual_failure_skip_and_xfail_never_become_pass(outcome, xfail):
    value = reports('a')
    value['reports'][1].update(outcome=outcome, wasxfail=xfail, longrepr='actual original diagnostic')
    value['exitstatus'] = int(outcome == 'failed')
    checked = validate_reports(['a'], value)
    assert case_failure(checked['a']) == 'actual original diagnostic'


def test_setup_failure_is_reported_without_inventing_call():
    value = reports('a')
    value['reports'] = [value['reports'][0], value['reports'][2]]
    value['reports'][0].update(outcome='failed', longrepr='original fixture failure')
    value['exitstatus'] = 1
    assert case_failure(validate_reports(['a'], value)['a']) == 'original fixture failure'


def test_selected_original_case_outcomes_preserve_order_and_count():
    nodes = ['c', 'a', 'b']
    checked = validate_reports(nodes, reports(*nodes))
    assert list(checked) == nodes
    assert all(case_failure(phases) is None for phases in checked.values())


def test_unlisted_identity_is_rejected_by_executor(tmp_path):
    with pytest.raises(RetainedExecutionError, match='unreviewed'):
        Execution(['model/unlisted/test_fake.py::test_fake'], log_directory=tmp_path)


@pytest.mark.parametrize('name', ['runner.py', 'reporter.py', 'plugin.py', 'live.py'])
def test_changed_execution_implementation_is_rejected(tmp_path, name):
    value = live.authority()
    names = set(value['implementation_sha256']) | {'model/retained_identity_execution/live.py',
                                                'model/retained_identity_execution/authority.json'}
    for path in names:
        target = tmp_path/path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(live.ROOT/path, target)
    target = tmp_path/'model/retained_identity_execution'/name
    target.write_text(target.read_text() + '\n# unreviewed implementation mutation\n')
    with pytest.raises(ValueError, match='implementation'):
        live.authority(tmp_path)


@pytest.mark.parametrize('index', range(9))
def test_exact_nine_dependency_transitions_reject_scientific_changes(index, tmp_path):
    value = live.authority()
    row = value['current_changes'][index]
    old = json.loads(subprocess.check_output(['git', 'show', value['revision']+':'+row['path']],
                                            cwd=live.ROOT, text=True))
    current = json.loads((live.ROOT/row['path']).read_text())
    for pointer in row['typed_nested_pointers']:
        target = tmp_path/pointer['target']
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(subprocess.check_output(['git', 'show', value['revision']+':'+pointer['target']],
                                                   cwd=live.ROOT))
    live.compare_transition(row, old, current, view=tmp_path)
    current['unreviewed_scientific_claim'] = {'positive_new_signal': True}
    with pytest.raises(ValueError): live.compare_transition(row, old, current, view=tmp_path)


def test_imported_executable_mutation_and_path_escape_are_rejected(tmp_path):
    old, current = tmp_path/'old', tmp_path/'current'
    for folder in (old, current):
        folder.mkdir()
        (folder/'source.py').write_text('answer = 1\n')
    value = {'imported_repository_python_sha256': {'source.py': live.sha(old/'source.py')}}
    assert validate_imported_sources(value, view=old, root=current) == 1
    (current/'source.py').write_text('answer = 2\n')
    with pytest.raises(RetainedExecutionError, match='executable'):
        validate_imported_sources(value, view=old, root=current)
    value['imported_repository_python_sha256'] = {'../escape.py': '0'*64}
    with pytest.raises(RetainedExecutionError): validate_imported_sources(value, view=old, root=current)


def test_documentary_handoff_is_exact_and_unreviewed_manifest_never_falls_back(tmp_path,monkeypatch):
    from model.completion import field_current_records_custody as custody
    from model.documentary_release_successor import model as documentary
    path=tmp_path/'paper/paper_release_manifest.json'
    path.parent.mkdir(parents=True)
    path.write_text('sealed')
    expected=live.sha(path)
    policy={'release_transition_scope':{'reviewed_release_source_sha256':{
        'paper/paper_release_manifest.json':expected}}}
    monkeypatch.setattr(custody,'load_authority',lambda root:policy)
    monkeypatch.setattr(custody,'validate',lambda root:{'original_current_checks':True})
    calls=[]
    def checked_documentary(root):
        calls.append(root)
        return {'actual_current_science_and_publication_checks':True}
    monkeypatch.setattr(documentary,'validate_current_science',checked_documentary)
    assert live.validate_scientific_custody(tmp_path)['domain']=='current_source_and_sealed_release_identity'
    assert not calls
    path.write_text('changed')
    assert live.validate_scientific_custody(tmp_path)['domain']=='current_source_with_separately_validated_documentary_successor'
    assert calls==[tmp_path.resolve()]
    def reject(root): raise ValueError('unreviewed documentary manifest')
    monkeypatch.setattr(documentary,'validate_current_science',reject)
    with pytest.raises(ValueError,match='unreviewed documentary'):
        live.validate_scientific_custody(tmp_path)


def test_passive_reporter_observes_real_pass_failure_and_skip(tmp_path):
    (tmp_path/'test_actual.py').write_text(
        'import pytest\ndef test_pass(): assert 1 == 1\n'
        'def test_fail(): assert 1 == 2\ndef test_skip(): pytest.skip("actual skip")\n')
    shutil.copyfile(Path(__file__).with_name('reporter.py'), tmp_path/'observer.py')
    output = tmp_path/'journal.json'
    bootstrap = ('import sys;sys.path.insert(0,'+repr(str(tmp_path))+');'
                 'import pytest;from observer import Reporter;'
                 'raise SystemExit(pytest.main(["-q","test_actual.py"],plugins=[Reporter('
                 +repr(str(output))+','+repr(str(tmp_path))+')]))')
    result = subprocess.run([sys.executable, '-I', '-B', '-W', 'error', '-c', bootstrap],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 1, result.stdout+result.stderr
    value = json.loads(output.read_text())
    checked = validate_reports(['test_actual.py::test_'+name for name in ('pass','fail','skip')], value)
    assert case_failure(checked['test_actual.py::test_pass']) is None
    assert 'assert 1 == 2' in case_failure(checked['test_actual.py::test_fail'])
    assert 'actual skip' in case_failure(checked['test_actual.py::test_skip'])
    assert 'test_actual.py' in value['imported_repository_python_sha256']
