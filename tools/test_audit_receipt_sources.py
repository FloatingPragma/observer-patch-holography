"""Mutation tests for custody preflight, independent of scientific calculations."""
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from tools.audit_receipt_sources import Inventory, Owner, audit, discover_owners


def write(root, path, value):
    p = root / path
    p.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value).encode()
    p.write_bytes(gzip.compress(data) if path.endswith('.gz') else data)
    return hashlib.sha256(p.read_bytes()).hexdigest()


def hierarchy(tmp_path):
    source = 'model/child/source.py'
    (tmp_path / source).parent.mkdir(parents=True)
    (tmp_path / source).write_text('source = 1\n')
    digest = hashlib.sha256((tmp_path / source).read_bytes()).hexdigest()
    child = 'model/child/results.json'
    child_sha = write(tmp_path, child, {'source_bindings': {source: digest}, 'value': 7})
    parent = 'model/parent/results.json'
    write(tmp_path, parent, {'nested': {'source_sha256': {child: child_sha}}})
    owners = {p: Owner('model.' + n + '.packet', None, p, 'test original owner')
              for n, p in [('parent', parent), ('child', child)]}
    return source, child, parent, owners


def test_hidden_stale_child_is_rejected_before_parent_write(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    (tmp_path / source).write_text('source = 2\n')
    result = audit(tmp_path, [parent], owners, mode='write')
    assert result['verdict'] == 'FAIL'
    assert result['failures'][0]['receipt'] == child
    assert [o['receipt'] for o in result['ordered_affected_owners']] == [child, parent]
    assert not result['scientific_replay_performed']


def test_write_permits_only_named_root_and_check_rejects_it(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    data = json.loads((tmp_path / parent).read_text())
    data['nested']['source_sha256'][child] = '0' * 64
    write(tmp_path, parent, data)
    allowed = audit(tmp_path, [parent], owners, mode='write')
    assert allowed['verdict'] == 'PASS'
    assert len(allowed['permitted_selected_root_refreshes']) == 1
    assert audit(tmp_path, [parent], owners)['verdict'] == 'FAIL'


@pytest.mark.parametrize('bad', [{}, [], {'missing.py': 'a' * 64},
    {'../../escape.py': 'a' * 64}, {'/absolute.py': 'a' * 64}, {'http://outside': 'a' * 64}])
def test_malformed_local_source_maps_fail_even_in_write_mode(tmp_path, bad):
    source, child, parent, owners = hierarchy(tmp_path)
    write(tmp_path, parent, {'nested': [{'source_bindings': bad}]})
    assert audit(tmp_path, [parent], owners, mode='write')['verdict'] == 'FAIL'


def test_invalid_digest_is_not_a_permitted_root_refresh(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    write(tmp_path, parent, {'source_bindings': {source: 'not a digest'}})
    assert audit(tmp_path, [parent], owners, mode='write')['verdict'] == 'FAIL'


def test_unknown_bound_result_owner_fails_closed(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    del owners[child]
    (tmp_path / source).write_text('source = 2\n')
    result = audit(tmp_path, [parent], owners)
    assert result['verdict'] == 'FAIL'
    assert result['failures'][0]['receipt'] == child
    assert child in result['affected_dependencies_without_executable_owner']
    assert all(row['receipt'] != child for row in result['ordered_affected_owners'])


def test_historical_gzip_capture_is_reported_not_rewritten(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    historical = 'model/historical/results.json.gz'
    digest = write(tmp_path, historical, {'source_bindings': {source: '0' * 64}})
    write(tmp_path, parent, {'source_sha256': {child: hashlib.sha256((tmp_path / child).read_bytes()).hexdigest(), historical: digest}})
    result = audit(tmp_path, [parent], owners, frozen_captures={historical: digest})
    assert result['verdict'] == 'PASS'
    assert result['background_findings'][0]['classification'] == 'authenticated_fixed_capture'
    assert historical not in result['live_dependency_receipts']


def test_literal_owner_discovery_does_not_execute_producer(tmp_path):
    p = tmp_path / 'model' / 'child' / 'packet.py'
    p.parent.mkdir(parents=True)
    p.write_text('from pathlib import Path\nHERE=Path(__file__).resolve().parent\nOUTPUT=HERE/"results.json"\nraise RuntimeError("must not import")\nif __name__ == "__main__":\n OUTPUT.write_text("--write --check")\n')
    owners, ambiguous = discover_owners(tmp_path)
    assert not ambiguous
    assert owners['model/child/results.json'].module == 'model.child.packet'


def test_duplicate_owner_is_not_guessed(tmp_path):
    package = tmp_path / 'model' / 'child'
    package.mkdir(parents=True)
    for filename in ('one.py', 'two.py'):
        (package / filename).write_text('from pathlib import Path\nOUTPUT=Path(__file__).with_name("results.json")\nif __name__ == "__main__": OUTPUT.write_text("--write --check")\n')
    owners, ambiguous = discover_owners(tmp_path)
    assert not owners
    assert len(ambiguous['model/child/results.json']) == 2


def test_duplicate_json_key_is_rejected_for_live_receipt(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    (tmp_path / parent).write_text('{"source_bindings": {}, "source_bindings": {}}')
    assert audit(tmp_path, [parent], owners)['verdict'] == 'FAIL'


def test_documentary_allowance_is_exact_and_does_not_cover_other_inputs(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    write(tmp_path, child, {'source_bindings': {source: '0' * 64}})
    write(tmp_path, parent, {'source_sha256': {child: hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()}})
    exact = frozenset({(child, 'root.source_bindings', source, '0' * 64)})
    assert audit(tmp_path, [parent], owners, documentary_exemptions=exact)['verdict'] == 'PASS'
    write(tmp_path, child, {'nested': {'source_bindings': {source: '0' * 64}}})
    write(tmp_path, parent, {'source_sha256': {child: hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()}})
    assert audit(tmp_path, [parent], owners, documentary_exemptions=exact)['verdict'] == 'FAIL'


def test_selected_child_cannot_evade_parent_strictness(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    (tmp_path / source).write_text('source = 2\n')
    assert audit(tmp_path, [parent, child], owners, mode='write')['verdict'] == 'FAIL'


def test_source_mutation_during_preflight_is_rejected(tmp_path, monkeypatch):
    source, child, parent, owners = hierarchy(tmp_path)
    original = Inventory.digest
    changed = False
    def mutate(self, path):
        nonlocal changed
        result = original(self, path)
        if path == tmp_path / source and not changed:
            path.write_text('source = 2\n')
            changed = True
        return result
    monkeypatch.setattr(Inventory, 'digest', mutate)
    result = audit(tmp_path, [parent], owners)
    assert result['verdict'] == 'FAIL'
    assert source in result['files_changed_during_audit']


def test_renamed_hidden_child_keeps_its_source_validation(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    renamed = 'model/child/audit.json'
    (tmp_path / child).rename(tmp_path / renamed)
    del owners[child]
    write(tmp_path, parent, {'nested': {'source_bindings': {renamed: hashlib.sha256((tmp_path / renamed).read_bytes()).hexdigest()}}})
    (tmp_path / source).write_text('source = 2\n')
    result = audit(tmp_path, [parent], owners, mode='write')
    assert result['verdict'] == 'FAIL'
    assert renamed in result['affected_dependencies_without_executable_owner']


@pytest.mark.parametrize('roots', [[], ['../escape.json'], ['/tmp/absolute.json']])
def test_empty_or_escaping_roots_are_rejected(tmp_path, roots):
    source, child, parent, owners = hierarchy(tmp_path)
    with pytest.raises(ValueError):
        audit(tmp_path, roots, owners)


def test_fixed_capture_requires_exact_outer_bytes(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    write(tmp_path, child, {'source_bindings': {source: '0' * 64}})
    digest = hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()
    write(tmp_path, parent, {'source_bindings': {child: digest}})
    assert audit(tmp_path, [parent], owners, frozen_captures={child: digest})['verdict'] == 'PASS'
    with pytest.raises(ValueError, match='fixed capture bytes changed'):
        audit(tmp_path, [parent], owners, frozen_captures={child: '0' * 64})


def test_cli_output_cannot_overwrite_audited_source(tmp_path, monkeypatch):
    from tools import audit_receipt_sources as module
    monkeypatch.setattr(module, 'ROOT', tmp_path)
    monkeypatch.setattr(module.sys, 'argv', ['audit', '--output', str(tmp_path / 'source.json')])
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2
    assert not (tmp_path / 'source.json').exists()


def test_path_sha_source_binding_records_are_checked_strictly(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    digest = hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()
    write(tmp_path, parent, {'source_bindings': [{'path': child, 'sha256': digest}]})
    assert audit(tmp_path, [parent], owners)['verdict'] == 'PASS'
    (tmp_path / source).write_text('source = 2\n')
    assert audit(tmp_path, [parent], owners, mode='write')['verdict'] == 'FAIL'


def test_duplicate_source_binding_records_are_rejected(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    item = {'path': child, 'sha256': hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()}
    write(tmp_path, parent, {'source_bindings': [item, item]})
    assert audit(tmp_path, [parent], owners, mode='write')['verdict'] == 'FAIL'


def test_changed_archived_spec_is_not_dynamically_exempted(tmp_path):
    from tools import audit_receipt_sources as module
    path = tmp_path / module.ARCHIVED_SPEC
    path.parent.mkdir(parents=True)
    path.write_text('{}')
    with pytest.raises(ValueError, match='archived source contract changed'):
        module.archived_source_map_roots(tmp_path)


def test_declared_source_map_root_checks_inner_archive_bytes(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    vendor = tmp_path / 'archive'
    vendor.mkdir()
    (vendor / 'native.py').write_text('native = 1\n')
    digest = hashlib.sha256((vendor / 'native.py').read_bytes()).hexdigest()
    write(tmp_path, child, {'source_bindings': {'native.py': digest}})
    write(tmp_path, parent, {'source_sha256': {child: hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()}})
    bases = {(child, 'root.source_bindings'): vendor}
    assert audit(tmp_path, [parent], owners, source_map_roots=bases)['verdict'] == 'PASS'
    (vendor / 'native.py').write_text('native = 2\n')
    assert audit(tmp_path, [parent], owners, source_map_roots=bases)['verdict'] == 'FAIL'


def test_nonscalar_expected_digest_fails_without_crashing(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    write(tmp_path, parent, {'source_bindings': {source: {'unexpected': 'shape'}}})
    result = audit(tmp_path, [parent], owners, mode='write')
    assert result['verdict'] == 'FAIL'
    assert result['failures'][0]['error'] == 'invalid SHA256 in local source map'


def test_historical_directory_name_does_not_exempt_hidden_child(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    historical = 'model/historical/audit.json'
    digest = write(tmp_path, historical, {'source_bindings': {source: '0' * 64}})
    write(tmp_path, parent, {'source_sha256': {historical: digest}})
    result = audit(tmp_path, [parent], owners, mode='write')
    assert result['verdict'] == 'FAIL'
    assert historical in result['affected_dependencies_without_executable_owner']


def test_bound_malformed_unowned_child_cannot_fall_outside_closure(tmp_path):
    source, child, parent, owners = hierarchy(tmp_path)
    del owners[child]
    (tmp_path / child).write_text('{invalid JSON')
    write(tmp_path, parent, {'source_sha256': {child: hashlib.sha256((tmp_path / child).read_bytes()).hexdigest()}})
    result = audit(tmp_path, [parent], owners)
    assert result['verdict'] == 'FAIL'
    assert child in result['live_dependency_receipts']
    assert any(row['receipt'] == child for row in result['failures'])


@pytest.mark.parametrize('guard', ["__name__ != '__main__'", "__name__ == '__main__' == False"])
def test_nonentrypoint_guard_does_not_claim_executable_owner(tmp_path, guard):
    package = tmp_path / 'model' / 'child'
    package.mkdir(parents=True)
    (package / 'packet.py').write_text('from pathlib import Path\nOUTPUT=Path(__file__).with_name("results.json")\nif ' + guard + ': OUTPUT.write_text("--write --check")\n')
    owners, ambiguous = discover_owners(tmp_path)
    assert not owners and not ambiguous
