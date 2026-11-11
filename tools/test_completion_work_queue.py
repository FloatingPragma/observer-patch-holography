"""Small mutation fixtures for task scope, readiness and publication rendering."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

import build_completion_work_queue as queue_builder
import build_completion_program_current as completion
import print_completion_issue as issue_printer


@pytest.fixture
def fixture(tmp_path):
    source = tmp_path / 'model/example/source.py'
    source.parent.mkdir(parents=True)
    source.write_text('def finite_kernel():\n    return 1\n')
    (source.parent / 'PROOF.md').write_text('# Finite channel\n\nComplete normalization.\n')
    (source.parent / 'test_source.py').write_text('def test_source():\n    assert True\n')
    (tmp_path / 'tracking').mkdir()
    (tmp_path / queue_builder.PREMISES).write_text(json.dumps({'rows': [{'id': 'PR-01'}]}))
    stages = []
    for gate, issue, status, dependencies in [('CP-01', 1136, 'ACCEPTED', []),
                                             ('CP-02', 1138, 'OPEN', ['CP-01']),
                                             ('CP-07', 1143, 'OPEN', ['CP-02'])]:
        stages.append(dict(id=gate, issue=issue, title='Declared source task', status=status,
                           requires=dependencies, deliverable='Construct the scoped source map.',
                           acceptance=['Bind the complete finite law.'], obstruction='Retain any counterexample.',
                           evidence=['model/example/PROOF.md'], legacy_issues=[]))
    program = dict(program_id='test-program', umbrella_issue=1134, stages=stages)
    task = dict(id='WQ-01', title='Complete the finite map', gate='CP-02', issue=1138,
                depends_on=[], status='ready', blocked_by=[], source_paths=['model/example/source.py'],
                reused_theorems=[dict(path='model/example/PROOF.md', anchor='Complete normalization.',
                                     role='Use completeness for the admitted finite channel.')],
                premise_ids=['PR-01'], exact_prompt='Construct the finite map on the stated source and retain every outcome.',
                outputs=[dict(path='model/example/result.json', description='Complete finite output law.')],
                validation=[dict(command='python -m pytest model/example/test_source.py',
                                 expected='The independently recomputed finite law agrees.')],
                closure={'pass': ['Every declared output is normalized.'],
                         'fail': ['Retain an exact nonnormalization witness.'],
                         'inconclusive': ['State the unbounded case or resource limit.']},
                scope_boundary='Finite source result; no physical gate acceptance.')
    second = deepcopy(task)
    second.update(id='WQ-02', gate='CP-07', issue=1143, status='blocked', depends_on=['WQ-01'],
                  outputs=[dict(path='model/example/scale.json', description='Declared scale interface.')])
    queue = dict(schema='oph.completion_work_queue.v1', queue_id='OPH-WORK-QUEUE-v1',
                 purpose='Construct bounded source results with explicit closure criteria.',
                 gate_acceptance_changed=False, physical_comparison_count_added=0,
                 execution_policy=['Run only tasks whose prerequisites are supplied.'],
                 premise_policy=['Keep architecture choices distinct from empirical inputs.'],
                 priority_order=['WQ-01', 'WQ-02'], first_batch=['WQ-01'], tasks=[task, second])
    (tmp_path / 'claims').mkdir()
    (tmp_path / queue_builder.QUEUE).write_text(json.dumps(queue))
    (tmp_path / queue_builder.PROGRAM).write_text(json.dumps(program))
    return tmp_path, program, queue


def test_live_shape_and_validation_never_execute_commands_or_mutate_gates(fixture):
    root, program, queue = fixture
    before = deepcopy(program)
    sentinel = root / 'command-must-not-run'
    queue['tasks'][0]['validation'][0]['command'] = f'touch {sentinel}'
    assert queue_builder.validate_queue(queue, root, program) is queue
    assert program == before
    assert not sentinel.exists()


@pytest.mark.parametrize('damage', [
    'duplicate_id', 'wrong_issue', 'unknown_gate', 'self_dependency', 'unknown_dependency',
    'cycle', 'premature_ready', 'unexplained_block', 'blocked_first_batch', 'duplicate_priority',
    'missing_source', 'escaping_source', 'missing_anchor', 'invented_premise', 'directory_output',
    'glob_output', 'escaping_output', 'missing_expected', 'missing_inconclusive', 'accept_gate',
    'physical_comparison', 'bool_comparison', 'unknown_task_promotion', 'completed_without_evidence',
])
def test_invalid_task_contracts_rejected(fixture, damage):
    root, program, original = fixture
    queue = deepcopy(original)
    first, second = queue['tasks']
    if damage == 'duplicate_id':
        second['id'] = first['id']
    elif damage == 'wrong_issue':
        first['issue'] = 1143
    elif damage == 'unknown_gate':
        first['gate'] = 'CP-99'
    elif damage == 'self_dependency':
        first['depends_on'] = ['WQ-01']
    elif damage == 'unknown_dependency':
        first['depends_on'] = ['WQ-99']
    elif damage == 'cycle':
        first.update(depends_on=['WQ-02'], status='blocked')
    elif damage == 'premature_ready':
        second['status'] = 'ready'
    elif damage == 'unexplained_block':
        first['status'] = 'blocked'
    elif damage == 'blocked_first_batch':
        queue['first_batch'] = ['WQ-02']
    elif damage == 'duplicate_priority':
        queue['priority_order'] = ['WQ-01', 'WQ-01']
    elif damage == 'missing_source':
        first['source_paths'] = ['model/example/missing.py']
    elif damage == 'escaping_source':
        first['source_paths'] = ['../outside.py']
    elif damage == 'missing_anchor':
        first['reused_theorems'][0]['anchor'] = 'theorem not_here'
    elif damage == 'invented_premise':
        first['premise_ids'] = ['PR-999']
    elif damage == 'directory_output':
        first['outputs'][0]['path'] = 'model/example'
    elif damage == 'glob_output':
        first['outputs'][0]['path'] = 'model/example/*.json'
    elif damage == 'escaping_output':
        first['outputs'][0]['path'] = '../outside.json'
    elif damage == 'missing_expected':
        first['validation'][0]['expected'] = ''
    elif damage == 'missing_inconclusive':
        del first['closure']['inconclusive']
    elif damage == 'accept_gate':
        queue['gate_acceptance_changed'] = True
    elif damage == 'physical_comparison':
        queue['physical_comparison_count_added'] = 1
    elif damage == 'bool_comparison':
        queue['physical_comparison_count_added'] = False
    elif damage == 'unknown_task_promotion':
        first['integrated_gate_accepted'] = True
    else:
        first['status'] = 'completed'
    with pytest.raises(ValueError):
        queue_builder.validate_queue(queue, root, program)


def complete_first(root, queue):
    result = root / queue['tasks'][0]['outputs'][0]['path']
    result.write_text('{"normalization": 1}\n')
    checker = root / 'model/example/test_source.py'
    queue['tasks'][0].update(status='completed', completion_outcome='positive', completion_evidence=[{
        'path': result.relative_to(root).as_posix(),
        'sha256': hashlib.sha256(result.read_bytes()).hexdigest(),
        'independent_check': 'model/example/test_source.py', 'scope': 'Only the declared finite channel.'}],
        completion_bindings={path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in (result, checker)})
    queue['tasks'][1]['status'] = 'ready'
    queue['first_batch'] = ['WQ-02']


def test_evidenced_task_completion_unlocks_dependency_without_accepting_gate(fixture):
    root, program, queue = fixture
    before = deepcopy(program)
    complete_first(root, queue)
    queue_builder.validate_queue(queue, root, program)
    assert program == before
    assert program['stages'][1]['status'] == 'OPEN'


@pytest.mark.parametrize('damage', ['stale_hash', 'missing_check', 'missing_output', 'ready_evidence',
                                     'stale_checker', 'unbound_output', 'own_checker', 'aliased_checker',
                                     'missing_outcome', 'unknown_outcome', 'missing_manifest'])
def test_completion_evidence_requires_real_files_and_current_hash(fixture, damage):
    root, program, queue = fixture
    complete_first(root, queue)
    row = queue['tasks'][0]['completion_evidence'][0]
    if damage == 'stale_hash':
        row['sha256'] = '0' * 64
    elif damage == 'missing_check':
        row['independent_check'] = 'model/example/missing_test.py'
    elif damage == 'missing_output':
        (root / row['path']).unlink()
    elif damage == 'ready_evidence':
        queue['tasks'][0]['status'] = 'ready'
    elif damage == 'stale_checker':
        (root / row['independent_check']).write_text('def test_source():\n    assert False\n')
    elif damage == 'unbound_output':
        second_output = root / 'model/example/second.json'
        second_output.write_text('{"unbound": true}\n')
        queue['tasks'][0]['outputs'].append(dict(path='model/example/second.json', description='Second required result.'))
    elif damage == 'own_checker':
        row['independent_check'] = row['path']
    elif damage == 'aliased_checker':
        alias = root / 'model/example/aliased_check.json'
        alias.symlink_to(root / row['path'])
        row['independent_check'] = alias.relative_to(root).as_posix()
    elif damage == 'missing_outcome':
        del queue['tasks'][0]['completion_outcome']
    elif damage == 'unknown_outcome':
        queue['tasks'][0]['completion_outcome'] = 'tests_passed'
    else:
        del queue['tasks'][0]['completion_bindings']
    with pytest.raises(ValueError):
        queue_builder.validate_queue(queue, root, program)


@pytest.mark.parametrize('outcome', ['negative', 'inconclusive'])
def test_nonpositive_completion_requires_explicit_dependency_outcome_admission(fixture, outcome):
    root, program, queue = fixture
    complete_first(root, queue)
    queue['tasks'][0]['completion_outcome'] = outcome
    with pytest.raises(ValueError, match='readiness disagrees'):
        queue_builder.validate_queue(queue, root, program)
    queue['tasks'][1]['status'] = 'blocked'
    queue['first_batch'] = []
    queue_builder.validate_queue(queue, root, program)
    queue['tasks'][1]['accepted_dependency_outcomes'] = {'WQ-01': ['positive', outcome]}
    queue['tasks'][1]['status'] = 'ready'
    queue['first_batch'] = ['WQ-02']
    queue_builder.validate_queue(queue, root, program)
    rendered = queue_builder.render_section(queue)
    assert f'completed: {outcome}' in rendered
    assert f'(positive/{outcome})' in rendered
    assert program['stages'][1]['status'] == 'OPEN'


@pytest.mark.parametrize('admission', [{'WQ-99': ['positive']}, {'WQ-01': []},
                                      {'WQ-01': ['tests_passed']}])
def test_unknown_or_untyped_dependency_outcome_admission_rejected(fixture, admission):
    root, program, queue = fixture
    queue['tasks'][1]['accepted_dependency_outcomes'] = admission
    with pytest.raises(ValueError):
        queue_builder.validate_queue(queue, root, program)


def test_explicit_external_blocker_and_whole_proof_are_admitted(fixture):
    root, program, queue = fixture
    queue['tasks'][0].update(status='blocked', blocked_by=['A declared external dataset is unavailable.'])
    queue['tasks'][0]['reused_theorems'][0]['anchor'] = ''
    queue['first_batch'] = []
    queue_builder.validate_queue(queue, root, program)


def test_duplicate_json_keys_and_symlink_escape_are_rejected(fixture, tmp_path):
    root, program, queue = fixture
    (root / queue_builder.QUEUE).write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError, match='duplicate JSON key'):
        queue_builder.load_queue(root, program)
    outside = tmp_path.parent / 'external-proof.md'
    outside.write_text('External proof')
    (root / 'model/example/outside.md').symlink_to(outside)
    queue['tasks'][0]['source_paths'] = ['model/example/outside.md']
    with pytest.raises(ValueError, match='escaping repository path'):
        queue_builder.validate_queue(queue, root, program)


def test_render_is_deterministic_copyable_and_uses_cross_owner_links(fixture):
    root, program, queue = fixture
    queue['tasks'][0]['exact_prompt'] += '\n```python\nassert True\n```'
    queue_builder.validate_queue(queue, root, program)
    full = queue_builder.render_section(queue)
    assert full == queue_builder.render_section(queue)
    assert '````text\n' in full
    assert full.count('### WQ-01:') == full.count('### WQ-02:') == 1
    owned = queue_builder.render_section(queue, 'CP-07')
    assert '### WQ-01:' not in owned and '### WQ-02:' in owned
    assert 'docs/COMPLETION_PROGRAM.md#work-wq-01' in owned
    assert queue['execution_policy'][0] in full and queue['execution_policy'][0] not in owned
    assert 'execution and premise policies' in owned
    assert queue['tasks'][1]['exact_prompt'] in owned
    assert queue_builder.render_section(queue, 'CP-01') == ''


def test_current_wrapper_inserts_queue_and_preserves_retained_output(fixture, monkeypatch):
    root, program, queue = fixture
    retained = '# Program\nGenerated by `tools/build_completion_program.py`\n\n## Existing foundations to extend\n\nFrozen contracts.\n'
    monkeypatch.setattr(completion.retained, 'render', lambda p, h: retained)
    result = completion.render(program, {}, root=root)
    section = queue_builder.render_section(queue_builder.load_queue(root, program))
    assert result.count(queue_builder.BEGIN) == result.count(queue_builder.END) == 1
    assert result.replace(section + '\n', '') == retained.replace(
        'Generated by `tools/build_completion_program.py`',
        'Generated by `tools/build_completion_program_current.py`')


def test_issue_printer_has_owned_tasks_by_default_and_explicit_legacy_mode(fixture):
    root, program, queue = fixture
    stage = program['stages'][1]
    body = issue_printer.issue_body(stage, program, root=root)
    old = issue_printer.issue_body(stage, program, root=root, include_work_queue=False)
    section = queue_builder.render_section(queue_builder.load_queue(root, program), 'CP-02')
    assert body.replace(section + '\n', '') == old
    assert '### WQ-01:' in body and '### WQ-02:' not in body
    assert '- [x] ' in body
    assert '## Acceptance\n' in body
