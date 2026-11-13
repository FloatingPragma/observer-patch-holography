"""Failure controls for conservative CI dependency selection."""
import json
from tools import ci_affected_model as c

H = 'a' * 64

def receipt(*paths):
    return json.dumps({'nested': {'source_sha256': {p: H for p in paths}}}).encode()

def fixture():
    return {
        'model/__init__.py': b'',
        'model/a/__init__.py': b'',
        'model/a/source.py': b'VALUE = 2\n',
        'model/a/test_source.py': b'from .source import VALUE\ndef test_value(): assert VALUE == 2\n',
        'model/a/results.json': receipt('model/a/source.py', 'model/a/test_source.py', 'docs/a.md'),
        'model/b/__init__.py': b'',
        'model/b/source.py': b'from model.a.source import VALUE\n',
        'model/b/test_source.py': b'from .source import VALUE\ndef test_value(): assert VALUE == 2\n',
        'model/b/results.json': receipt('model/b/source.py', 'model/b/test_source.py', 'model/a/results.json'),
        'model/c/__init__.py': b'',
        'model/c/source.py': b'VALUE = 3\n',
        'model/c/test_source.py': b'from .source import VALUE\ndef test_value(): assert VALUE == 3\n',
        'model/c/results.json': receipt('model/c/source.py', 'model/c/test_source.py'),
        'model/test_global.py': b'def test_global(): assert True\n',
    }

def select(files, *changed, old=None):
    return c.plan_from_files(files, changed, old if old is not None else {p: files[p] for p in changed if p in files})

def test_transitive_import_and_nested_child_receipt_edges_select_all_consumers():
    files = fixture()
    for change in ('model/a/source.py', 'model/a/results.json', 'docs/a.md'):
        plan = select(files, change)
        assert plan['mode'] == 'affected'
        assert set(plan['selected_test_modules']) == {'model/a/test_source.py', 'model/b/test_source.py', 'model/test_global.py'}
    assert c.source_paths(json.loads(files['model/b/results.json'])) == {'model/b/source.py', 'model/b/test_source.py', 'model/a/results.json'}

def test_unrelated_documents_keep_unowned_tests_and_never_invent_a_replay():
    plan = select(fixture(), 'docs/unrelated.md')
    assert plan['selected_test_modules'] == ['model/test_global.py']
    assert plan['full_source_guards_required'] and not plan['scientific_comparisons_changed']

def test_dynamic_and_unresolved_imports_are_always_selected():
    files = fixture()
    files['model/c/source.py'] = b'import importlib\nx = importlib.import_module("model.a.source")\n'
    assert 'model/c/test_source.py' in select(files, 'docs/unrelated.md')['selected_test_modules']
    files['model/c/source.py'] = b'from model.missing import value\n'
    assert 'model/c/test_source.py' in select(files, 'docs/unrelated.md')['selected_test_modules']

def test_current_and_old_import_edges_are_both_retained():
    files = fixture()
    old = {'model/b/source.py': files['model/b/source.py'], 'model/a/source.py': files['model/a/source.py']}
    files['model/b/source.py'] = b'VALUE = 2\n'
    plan = select(files, 'model/a/source.py', 'model/b/source.py', old=old)
    assert 'model/b/test_source.py' in plan['selected_test_modules']

def test_changed_input_key_set_is_full_even_if_new_map_drops_evidence():
    files = fixture(); old = {'model/b/results.json': files['model/b/results.json']}
    files['model/b/results.json'] = receipt('model/b/test_source.py')
    plan = select(files, 'model/b/results.json', old=old)
    assert plan['mode'] == 'full'
    assert any('input inventory changed' in r for r in plan['full_reasons'])

def test_new_deleted_and_renamed_modules_force_full():
    for mode in ('new', 'deleted', 'renamed'):
        files = fixture(); old = {}; changed = []
        if mode in ('deleted', 'renamed'):
            old['model/a/source.py'] = files.pop('model/a/source.py'); changed.append('model/a/source.py')
        if mode in ('new', 'renamed'):
            files['model/a/new.py'] = b'X = 3\n'; changed.append('model/a/new.py')
        assert select(files, *changed, old=old)['mode'] == 'full'

def test_config_unknown_data_and_unreadable_metadata_fail_closed():
    for path in ('pytest.ini', 'model/requirements.txt', 'model/unowned.dat'):
        assert select(fixture(), path)['mode'] == 'full'
    files = fixture(); files['model/a/results.json'] = b'broken json'
    assert select(files, 'model/a/results.json')['mode'] == 'full'

def test_source_map_rejects_escaping_inputs():
    import pytest
    with pytest.raises(ValueError, match='escaping'):
        c.source_paths({'source_bindings': {'../outside.py': H}})

def test_shared_core_and_explicit_full_keep_every_test():
    files = fixture()
    assert select(files, 'model/__init__.py')['selected_test_modules'] == select(files, 'model/__init__.py')['all_test_modules']
    plan = c.plan_from_files(files, [], force_reason='scheduled/manual/release complete coverage')
    assert plan['mode'] == 'full' and plan['shard_count'] == 16
    assert plan['selected_test_modules'] == plan['all_test_modules']

def test_unavailable_base_cannot_make_an_empty_green_selection(tmp_path):
    import subprocess
    subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
    subprocess.run(['git', '-C', str(tmp_path), '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid', 'commit', '--allow-empty', '-qm', 'fixture'], check=True)
    (tmp_path/'model').mkdir(); (tmp_path/'model/test_a.py').write_text('def test_a(): pass\n')
    plan = c.build_plan(tmp_path, base='0'*40)
    assert plan['mode'] == 'full' and plan['selected_test_modules'] == ['model/test_a.py']

def test_directory_proximity_and_downstream_receipts_do_not_establish_ownership():
    files = fixture()
    files['model/c/test_unowned.py'] = b'def test_other(): assert True\n'
    files['model/b/results.json'] = receipt('model/a/test_source.py', 'model/b/test_source.py')
    plan = select(files, 'docs/unrelated.md')
    assert 'model/c/test_unowned.py' in plan['always_test_modules']
    # A's actual package owner remains its nearest explicitly declared receipt.
    assert 'model/a/test_source.py' not in plan['selected_test_modules']


def test_importlib_aliases_function_references_and_reflection_cannot_hide_dynamic_loading():
    examples = [
        'from importlib import import_module as f\nx = f("model.a.source")\n',
        'from runpy import run_module as f\nx = f("model.a.source")\n',
        'import importlib as m\nf = m.import_module\nx = f("model.a.source")\n',
        'import importlib as m\nx = getattr(m, "import_module")("model.a.source")\n',
        'import subprocess as s\ns.run(["python", "other.py"])\n',
    ]
    for text in examples:
        files = fixture(); files['model/c/source.py'] = text.encode()
        assert 'model/c/test_source.py' in select(files, 'docs/unrelated.md')['always_test_modules']
    files = fixture(); files['model/c/source.py'] = b'import importlib.util\nx = importlib.util.resolve_name(".source", "model.c")\n'
    assert 'model/c/test_source.py' not in select(files, 'docs/unrelated.md')['always_test_modules']


def test_all_guard_map_spellings_records_and_relative_bases_are_covered():
    from tools.audit_receipt_sources import MAP_FIELDS, PACKAGE_LEVELS
    assert MAP_FIELDS <= c.SOURCE_MAPS
    assert PACKAGE_LEVELS == c.PACKAGE_LEVELS
    for field in MAP_FIELDS:
        assert c.source_paths({field: {'model/a/source.py': H}}) == {'model/a/source.py'}
    assert c.source_paths({'source_bindings': [{'path': 'model/a/source.py', 'sha256': H}]}) == {'model/a/source.py'}
    assert c.source_paths({'implementation_sha256': {'source.py': H}}, 'model/a/results.json') == {'model/a/source.py'}
    assert c.source_paths({'source_interface_hashes': {'a/source.py': H}}, 'model/a/results.json') == {'model/a/source.py'}
    value = {'retrospective_endpoint_audit': {'source': {'pins': {'model/a/source.py': H}}}}
    assert c.source_paths(value) == {'model/a/source.py'}
    value = {'source_sha256': {'oph_exact/file.py': H}}
    overrides = {('evidence/spec.json', 'root.source_sha256'): 'evidence/archive'}
    assert c.source_paths(value, 'evidence/spec.json', source_bases=overrides) == {'evidence/archive/oph_exact/file.py'}


def test_malformed_and_duplicate_record_maps_never_make_an_omission():
    import pytest
    with pytest.raises(ValueError, match='duplicate'):
        c.source_paths({'source_bindings': [{'path': 'a.py', 'sha256': H}] * 2})
    with pytest.raises(ValueError, match='malformed'):
        c.source_paths({'source_bindings': [{'path': 'a.py', 'wrong': H}]})
    files = fixture(); files['model/c/results.json'] = b'{"source_sha256": "wrong"}'
    assert 'model/c/test_source.py' in select(files, 'docs/unrelated.md')['selected_test_modules']


def test_pytest_ignored_modules_do_not_create_empty_or_invalid_shards():
    files = fixture()
    files['pytest.ini'] = b'[pytest]\naddopts = --ignore=model/c\n'
    assert 'model/c/test_source.py' not in select(files, 'docs/unrelated.md')['all_test_modules']


def test_costs_change_worker_balance_only_not_coverage():
    files = fixture(); first = select(files, 'docs/unrelated.md')
    files['model/pytest_shard_costs.json'] = json.dumps({'seconds': {'model/test_global.py': 1000}}).encode()
    second = select(files, 'docs/unrelated.md')
    assert first['selected_test_modules'] == second['selected_test_modules']
    assert 1 <= second['shard_count'] <= len(second['selected_test_modules'])

def signed_plan(tmp_path, monkeypatch):
    import hashlib
    from pathlib import Path
    plan = select(fixture(), 'docs/unrelated.md')
    plan.update(head='b'*40, base='a'*40, selector_sha256=hashlib.sha256(Path(c.__file__).read_bytes()).hexdigest())
    def write(value):
        value = dict(value); value.pop('plan_sha256', None)
        digest = hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        value['plan_sha256'] = digest
        path = tmp_path/'plan.json'; path.write_text(json.dumps(value))
        return path, digest
    monkeypatch.setattr(c, 'git', lambda *args: b'b'*40+b'\n')
    path, digest = write(plan)
    return plan, path, digest, write


def test_plan_artifact_is_bound_to_job_hash_exact_commit_base_and_selector(tmp_path, monkeypatch):
    import pytest
    plan, path, digest, write = signed_plan(tmp_path, monkeypatch)
    assert c.authenticate_plan(tmp_path, path, digest, 'a'*40)['selected_module_count'] == 1
    with pytest.raises(ValueError, match='SHA256'):
        c.authenticate_plan(tmp_path, path, 'c'*64, 'a'*40)
    with pytest.raises(ValueError, match='head/base'):
        c.authenticate_plan(tmp_path, path, digest, 'd'*40)
    bad = dict(plan); bad['selector_sha256'] = '0'*64; path, digest = write(bad)
    with pytest.raises(ValueError, match='selector source'):
        c.authenticate_plan(tmp_path, path, digest, 'a'*40)
    bad = dict(plan); bad['head'] = 'c'*40; path, digest = write(bad)
    with pytest.raises(ValueError, match='head/base'):
        c.authenticate_plan(tmp_path, path, digest, 'a'*40)


def test_authenticated_plan_rejects_duplicate_keys_paths_and_changed_matrix(tmp_path, monkeypatch):
    import pytest
    plan, path, digest, write = signed_plan(tmp_path, monkeypatch)
    path.write_text('{"x":1,"x":2}')
    with pytest.raises(ValueError, match='duplicate'):
        c.authenticate_plan(tmp_path, path, digest, 'a'*40)
    for key, value in [('matrix', [99]), ('selected_test_modules', ['model/../outside/test_a.py'])]:
        bad = dict(plan); bad[key] = value; path, digest = write(bad)
        with pytest.raises(ValueError):
            c.authenticate_plan(tmp_path, path, digest, 'a'*40)


def test_full_runner_executes_original_pytest_command_without_changed_comparison(monkeypatch):
    import sys
    plan = c.plan_from_files(fixture(), [], force_reason='full')
    plan.update(base=None, head='b'*40)
    monkeypatch.setattr(c, 'build_plan', lambda *args: dict(plan))
    monkeypatch.setattr(sys, 'argv', ['selection', '--full', '--run-shard', '3'])
    calls = []; monkeypatch.setattr(c.subprocess, 'run', lambda args, **kw: calls.append((args,kw)))
    c.main()
    assert calls[0][0] == [sys.executable, '-m', 'pytest', '-q', '-W', 'error', '--durations=30',
        '-p', 'model.pytest_shards', '--model-shard-index', '3', '--model-shard-count', '16',
        '--deselect', 'model/charged_lepton_response/test_packet.py::test_generated_results_and_runtime_custody', 'model']
    assert calls[0][1]['check'] is True


def test_replacement_preserves_full_science_commands_source_guards_and_filters():
    from pathlib import Path
    import yaml
    root = Path(c.__file__).resolve().parents[1]
    old = yaml.safe_load((root/'.github/workflows/observer-model.yml').read_text())
    new = yaml.safe_load((root/'.github/workflows/observer-affected.yml').read_text())
    runs = lambda job: [step['run'] for step in job['steps'] if 'run' in step]
    assert runs(old['jobs']['observer']) == runs(new['jobs']['observer'])
    assert old['jobs']['observer']['if'] == new['jobs']['observer']['if']
    assert old['jobs']['observer']['strategy'] == new['jobs']['observer']['strategy']
    old_reg = runs(old['jobs']['regression']); new_reg = runs(new['jobs']['regression'])
    ordinary = next(command for command in old_reg if '--model-shard-count 16' in command)
    selected = next(command for command in new_reg if '--plan ' in command)
    assert [v for v in old_reg if v != ordinary] == [v for v in new_reg if v != selected]
    assert '--expect-plan-sha256' in selected and '--run-shard' in selected
    assert 'build_plan' not in selected and '--plan "$RUNNER_TEMP/' in selected
    assert 'python -m model.verify_current --hashes-only' in runs(new['jobs']['hashes'])
    preflight = next(v for v in runs(new['jobs']['hashes']) if 'audit_receipt_sources.py' in v)
    assert '--mode check --receipt model/particle_identification/completion_results.json' in preflight
    assert set(new['jobs']['regression']['needs']) == {'scope', 'hashes', 'infrastructure'}
    assert 'needs' not in new['jobs']['scope'] and 'needs' not in new['jobs']['hashes'] and 'needs' not in new['jobs']['infrastructure']
    old_on, new_on = old.get('on', old.get(True)), new.get('on', new.get(True))
    for event in ('push', 'pull_request'):
        assert set(old_on[event]['paths']) <= set(new_on[event]['paths'])
        assert 'tools/ci_affected_model.py' in new_on[event]['paths']
    assert old_on['schedule'] == new_on['schedule']
    assert 'workflow_dispatch' in new_on and new_on['push']['tags'] == ['r*']
    cancel = new['concurrency']['cancel-in-progress']
    assert "github.event_name == 'pull_request'" in cancel and "refs/heads/" in cancel
    assert 'github.run_id' in new['concurrency']['group']

def test_precollection_partitions_match_original_sharder_without_dropping_modules():
    from model.pytest_shards import selected_modules
    modules = [f'model/p/test_{i}.py' for i in range(19)]
    costs = {modules[0]: 1000, modules[3]: 80, modules[11]: 4}
    for count in (1, 2, 4, 16):
        parts = c.partition_modules(modules, count, costs)
        assert [set(group) for group in parts] == [set(selected_modules(modules, i, count, costs)) for i in range(count)]
        assert sorted(p for group in parts for p in group) == sorted(modules)
    assert c.partition_modules([], 1, {}) == [[]]


def test_artifact_cannot_drop_or_duplicate_a_partitioned_module(tmp_path, monkeypatch):
    import pytest
    plan, path, digest, write = signed_plan(tmp_path, monkeypatch)
    for partitions in ([[]], [['model/test_global.py', 'model/test_global.py']]):
        bad = dict(plan); bad['module_partitions'] = partitions; path, digest = write(bad)
        with pytest.raises(ValueError, match='duplicates or drops'):
            c.authenticate_plan(tmp_path, path, digest, 'a'*40)


def test_affected_worker_collects_only_assigned_modules_but_keeps_original_fixture_plugin(monkeypatch):
    import sys
    plan = select(fixture(), 'model/a/source.py')
    plan.update(base=None, head='b'*40)
    monkeypatch.setattr(c, 'build_plan', lambda *args: dict(plan))
    monkeypatch.setattr(sys, 'argv', ['selection', '--run-shard', '0'])
    calls = []; monkeypatch.setattr(c.subprocess, 'run', lambda args, **kw: calls.append(args))
    c.main()
    command = calls[0]
    assert command[command.index('-p')+1] == 'model.pytest_shards'
    assert command[command.index('--model-shard-index')+1] == '0'
    assert command[command.index('--model-shard-count')+1] == '1'
    assert command[-len(plan['module_partitions'][0]):] == plan['module_partitions'][0]
    assert 'model' not in command

def test_record_serialization_is_not_confused_with_loading_code():
    safe = [
        'def serial(record): return dict(record.__dict__)\n',
        'def serial(record): return vars(record)\n',
        'def serial(record): return vars(record)["value"]\n',
        'def serial(): return sorted(globals())\n',
        'def serial(record): return locals()\n',
    ]
    for text in safe:
        assert c.imports('model/example.py', text, {'model/example.py'})[1] is False
    unsafe = [
        'globals()[name]("model.changed")\n',
        'f = globals()[name]\ng = f\ng("model.changed")\n',
        'namespace = globals()\nf = namespace["__import__"]\nf("model.changed")\n',
        'import importlib as loader\nf = loader.__dict__[name]\nf("model.changed")\n',
        'def run(f): return f.__globals__[name]("model.changed")\n',
    ]
    for text in unsafe:
        assert c.imports('model/example.py', text, {'model/example.py'})[1] is True

def test_actual_fixed_loader_boundary_covers_every_archived_library_file():
    from pathlib import Path
    root = Path(c.__file__).resolve().parents[1]
    contracts = c.reviewed_loaders(root)
    assert set(contracts) == set(c.LOADER_SHA256)
    archive = c.fixed_archive(root)
    assert len(archive) == 259
    assert archive <= contracts['model/federated_observer/source.py']
    assert archive <= contracts['model/common_source/slow_geometry.py']
    assert contracts['model/observer/p_response/experiment.py'] == {
        'code/P_derivation/paper_math.py',
        'code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json',
        'model/observer/native.py'}


def test_archive_membership_bytes_and_reverse_imports_fail_closed(tmp_path):
    import hashlib
    import pytest
    path = tmp_path/c.ARCHIVE/'module.py'; path.parent.mkdir(parents=True); path.write_text('VALUE = 1\n')
    names = [path.relative_to(tmp_path).as_posix()]
    expected_names = c.digest_json(names)
    expected_content = c.digest_json({names[0]: hashlib.sha256(path.read_bytes()).hexdigest()})
    assert c.fixed_archive(tmp_path, expected_names, expected_content) == set(names)
    path.write_text('VALUE = 2\n')
    with pytest.raises(ValueError, match='source bytes'):
        c.fixed_archive(tmp_path, expected_names, expected_content)
    path.write_text('VALUE = 1\n'); extra = path.with_name('added.py'); extra.write_text('pass\n')
    with pytest.raises(ValueError, match='membership'):
        c.fixed_archive(tmp_path, expected_names, expected_content)
    extra.unlink(); path.unlink()
    with pytest.raises(ValueError, match='membership'):
        c.fixed_archive(tmp_path, expected_names, expected_content)
    path.write_text('from model import current\n')
    changed_content = c.digest_json({names[0]: hashlib.sha256(path.read_bytes()).hexdigest()})
    with pytest.raises(ValueError, match='imports current RER'):
        c.fixed_archive(tmp_path, expected_names, changed_content)
    path.write_text('import importlib\nimportlib.import_module("tools.helper")\n')
    changed_content = c.digest_json({names[0]: hashlib.sha256(path.read_bytes()).hexdigest()})
    with pytest.raises(ValueError, match='imports current RER'):
        c.fixed_archive(tmp_path, expected_names, changed_content)


def test_fixed_loader_byte_path_and_data_pin_mutations_require_full_review(tmp_path, monkeypatch):
    from pathlib import Path
    import pytest
    root = Path(c.__file__).resolve().parents[1]
    for name in c.LOADER_SHA256:
        target=tmp_path/name; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes((root/name).read_bytes())
    geometry = tmp_path/'model/common_source/slow_geometry.py'
    original = geometry.read_bytes(); geometry.write_bytes(original.replace(c.ARCHIVE.encode(), b'changed/archive'))
    with pytest.raises(ValueError, match='loader requires renewed review'):
        c.reviewed_loaders(tmp_path)
    geometry.write_bytes(original)
    monkeypatch.setattr(c, 'fixed_archive', lambda root: frozenset())
    paper = tmp_path/'code/P_derivation/paper_math.py'; paper.parent.mkdir(parents=True); paper.write_text('changed\n')
    with pytest.raises(ValueError, match='response input pin'):
        c.reviewed_loaders(tmp_path)
    loader='model/observer/p_response/experiment.py'; text=(root/loader).read_text()
    with pytest.raises(ValueError, match='loader source changed'):
        c.imports(loader, text+'\n# mutation\n', {loader}, {loader: frozenset({'a.py'})})


def test_bounded_loader_dependencies_still_select_their_consumers():
    from pathlib import Path
    root = Path(c.__file__).resolve().parents[1]
    loader='model/observer/p_response/experiment.py'; text=(root/loader).read_text()
    dependencies={'code/P_derivation/paper_math.py','model/observer/native.py'}
    edges, unknown=c.imports(loader,text,{loader}|dependencies,{loader:frozenset(dependencies)})
    assert dependencies <= edges and not unknown
    # The same mechanism does not excuse an unreviewed arbitrary loader.
    assert c.imports('model/unknown.py','import importlib\nimportlib.import_module(name)\n',set(),{})[1]

def test_fast_route_preserves_union_and_falls_back_on_changed_scientific_contract(tmp_path, monkeypatch):
    import hashlib
    path=tmp_path/'guard.py'; path.write_text('reviewed\n')
    monkeypatch.setattr(c, 'FAST_SOURCE_HASHES', {'guard.py':hashlib.sha256(path.read_bytes()).hexdigest()})
    costs=tmp_path/'model/pytest_shard_costs.json'; costs.parent.mkdir(); costs.write_text('{"seconds":{}}')
    plan=select(fixture(),'docs/unrelated.md'); plan['selected_test_modules']=[c.FAST_TEST,'model/test_global.py']
    plan['all_test_modules']=sorted(set(plan['all_test_modules'])|{c.FAST_TEST})
    routed=c.apply_fast_route(dict(plan),tmp_path,True)
    assert routed['separate_test_modules']==[c.FAST_TEST]
    assert routed['ordinary_test_modules']==['model/test_global.py']
    assert [p for g in routed['module_partitions']for p in g]==['model/test_global.py']
    path.write_text('changed\n')
    fallback=c.apply_fast_route(dict(plan),tmp_path,True)
    assert not fallback['fast_identification'] and fallback['ordinary_test_modules']==plan['selected_test_modules']
    assert not c.apply_fast_route(dict(plan),tmp_path,False)['fast_identification']


def test_fast_route_is_not_executed_if_module_is_unselected_or_harness_changes(tmp_path,monkeypatch):
    monkeypatch.setattr(c,'FAST_SOURCE_HASHES',{})
    plan=select(fixture(),'docs/unrelated.md')
    assert not c.apply_fast_route(dict(plan),tmp_path,True)['fast_identification']
    plan['selected_test_modules']=[c.FAST_TEST]
    (tmp_path/'conftest.py').write_text('pass\n')
    assert not c.apply_fast_route(dict(plan),tmp_path,True)['fast_identification']


def test_fast_route_pins_all_scientific_implementations_and_full_events_stay_plain():
    from tools.immutable_custody_replay import TARGET_SOURCE_SHA256
    from tools.exact_walsh import SYNTHESIS_PATH,SYNTHESIS_SHA256
    assert c.FAST_SOURCE_HASHES[SYNTHESIS_PATH]==SYNTHESIS_SHA256
    for name,h in TARGET_SOURCE_SHA256.items():
        assert c.FAST_SOURCE_HASHES[name.replace('.','/')+'.py']==h
    from pathlib import Path
    import yaml
    workflow=yaml.safe_load((Path(c.__file__).resolve().parents[1]/'.github/workflows/observer-affected.yml').read_text())
    script=next(step['run']for step in workflow['jobs']['scope']['steps']if step.get('id')=='plan')
    assert 'then arguments+=(--full); else arguments+=(--fast-identification); fi' in script
    assert workflow['jobs']['exact-interface']['if']=="needs.scope.outputs.fast_identification == 'true'"
    assert set(workflow['jobs']['exact-interface']['needs'])=={'scope','hashes','infrastructure'}
