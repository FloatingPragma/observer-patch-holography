"""Separate current scientific validation from a sealed documentary manifest."""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
AUTHORITY = 'model/documentary_release_successor/authority.json'
AUTHORITY_SHA256 = 'c9615eac98fba6e93ee7b43fef13ffb7a05cb23e81e48395b3720229649cee12'
MANIFEST = 'paper/paper_release_manifest.json'
RETAINED = 'model/documentary_release_successor/retained_manifest.json'
ALLOWED_ROWS = (('book',),
    ('papers', 'from_observer_consensus_to_standard_physics'),
    ('papers', 'deriving_the_particle_zoo_from_observer_consistency'),
    ('papers', 'executable_observer_construction'))
ALLOWED_LEAVES = frozenset(row + (key,) for row in ALLOWED_ROWS for key in ('sha256', 'size_bytes'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def differences(old, current, path=()):
    if type(old) is not type(current):
        return [(path, old, current)]
    if isinstance(old, dict):
        require(old.keys() == current.keys(), 'unlisted documentary key inventory change')
        return [change for key in old for change in differences(old[key], current[key], path+(key,))]
    if isinstance(old, list):
        require(len(old) == len(current), 'unlisted documentary list scope change')
        return [change for i, (a, b) in enumerate(zip(old, current)) for change in differences(a, b, path+(i,))]
    return [] if old == current else [(path, old, current)]


def compare_manifests(old, current, reviewed_changes):
    changes = differences(old, current)
    require(all(path in ALLOWED_LEAVES for path, _, _ in changes), 'unlisted documentary manifest change')
    expected = [(tuple(row['json_path']), row['before'], row['after']) for row in reviewed_changes]
    require(len(expected) == len(set(path for path, _, _ in expected)), 'duplicate reviewed documentary leaf')
    require(sorted(changes) == sorted(expected), 'documentary changes differ from exact reviewed values')
    for path, _, value in changes:
        if path[-1] == 'size_bytes':
            require(type(value) is int and value > 0, 'invalid documentary artifact size')
        else:
            require(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None,
                    'invalid documentary artifact digest')
    require(old['release_id'] == current['release_id'] == 'r2064', 'documentary release identity changed')
    return changes


def template_sha(path):
    source, count = re.subn(r"AUTHORITY_SHA256 = '[0-9a-f]{64}'", "AUTHORITY_SHA256 = '"+'0'*64+"'", Path(path).read_text())
    require(count == 1, 'documentary authority slot changed')
    return hashlib.sha256(source.encode()).hexdigest()


def authority(root=ROOT):
    root = Path(root).resolve()
    require(sha(root/AUTHORITY) == AUTHORITY_SHA256, 'documentary authority changed')
    value = strict_json.load(root/AUTHORITY)
    tests = value['original_test_sha256']
    expected_tests = {'model/completion/test_field_current_records_custody.py',
                      'model/completion/test_field_current_records_navigation.py'}
    nodes = value['original_test_nodeids']
    require(set(tests) == expected_tests and len(nodes) == len(set(nodes)) == 35
            and {node.split('::')[0] for node in nodes} == expected_tests,
            'documentary current-custody test scope changed')
    for name, expected in tests.items():
        require(sha(root/name) == expected, 'original current-custody test changed: '+name)
    require(template_sha(root/'model/documentary_release_successor/model.py') == value['validator_template_sha256'],
            'documentary validator implementation changed')
    for name, expected in value['input_sha256'].items():
        require(sha(root/name) == expected, 'documentary source changed: '+name)
    require(value['retained_manifest_sha256'] == 'a68eefc1d92a87bd30ee2fbc11d2a38d41bb9062b69e83877f91accbf02f971c',
            'retained r2064 documentary identity changed')
    require(sha(root/RETAINED) == value['retained_manifest_sha256'], 'retained documentary manifest changed')
    require(sha(root/MANIFEST) == value['current_manifest_sha256'], 'current documentary manifest is unreviewed')
    compare_manifests(strict_json.load(root/RETAINED), strict_json.load(root/MANIFEST), value['reviewed_changes'])
    return value


def validate_actual_release(root=ROOT):
    root = Path(root).resolve()
    value = authority(root)
    result = subprocess.run([sys.executable, '-B', '-W', 'error',
                             str(root/'tools/validate_paper_release_manifest.py')],
                            cwd=root, capture_output=True, text=True)
    require(result.returncode == 0, 'actual current publication artifacts failed strict validation:\n'+result.stdout+result.stderr)
    return {'current_manifest_sha256': value['current_manifest_sha256'],
            'actual_strict_manifest_and_artifact_validation': True,
            'strict_validator_output': result.stdout.strip(),
            'release_id_and_date_unchanged': True,
            'scientific_identity_or_physical_status_changed': False}


def current_paths(root):
    root = Path(root).resolve()
    result = subprocess.run(['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
                            cwd=root, check=True, capture_output=True)
    paths = sorted(set(result.stdout.decode().split('\0'))-{''})
    require(paths and MANIFEST in paths, 'current source inventory has no publication manifest')
    for name in paths:
        path = (root/name).resolve()
        require(path.is_relative_to(root) and path.is_file() and '.git' not in Path(name).parts,
                'missing or escaping current source path: '+name)
    return paths


def assert_view_equal(root, view, snapshot, retained_sha):
    root, view = Path(root).resolve(), Path(view).resolve()
    actual_paths = {path.relative_to(view).as_posix() for path in view.rglob('*') if path.is_file() or path.is_symlink()}
    require(actual_paths == set(snapshot), 'private current-source file inventory changed')
    for name, expected in snapshot.items():
        require(sha(root/name) == expected, 'current source changed during documentary validation: '+name)
        actual_expected = retained_sha if name == MANIFEST else expected
        require(sha(view/name) == actual_expected, 'private current-source view changed: '+name)


@contextmanager
def scientific_view(root=ROOT):
    """All current files, with exactly one historical documentary substitution."""
    root = Path(root).resolve()
    approved = authority(root)
    paths = current_paths(root)
    snapshot = {name: sha(root/name) for name in paths}
    with tempfile.TemporaryDirectory(prefix='oph-current-source-retained-release-') as folder:
        view = Path(folder).resolve()/'source'
        view.mkdir()
        for name in paths:
            destination = view/name
            destination.parent.mkdir(parents=True, exist_ok=True)
            # Never alias the actual source: even deliberate mutation controls
            # in an unchanged test must be unable to write through to it.
            shutil.copy2((root/name).resolve(), destination)
        # Replace the one documentary input only after copying all current files.
        (view/MANIFEST).unlink()
        shutil.copyfile(root/RETAINED, view/MANIFEST)
        assert_view_equal(root, view, snapshot, approved['retained_manifest_sha256'])
        try:
            yield view, {'current_file_count': len(snapshot),
                         'current_source_snapshot_sha256': hashlib.sha256(json.dumps(snapshot,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                         'all_current_scientific_files_preserved': True,
                         'only_substituted_path': MANIFEST,
                         'retained_manifest_sha256': approved['retained_manifest_sha256'],
                         'current_manifest_sha256': approved['current_manifest_sha256'],
                         'publication_validity_of_retained_manifest_view_claimed': False}
        finally:
            require(current_paths(root) == paths, 'current source inventory changed during documentary validation')
            assert_view_equal(root, view, snapshot, approved['retained_manifest_sha256'])


def validate_current_science(root=ROOT):
    root = Path(root).resolve()
    release = validate_actual_release(root)
    with scientific_view(root) as (view, report):
        bootstrap = ('import json,sys;sys.path.insert(0,'+repr(str(view))+');'
                     'from model.completion.field_current_records_custody import validate;'
                     'print(json.dumps(validate(),sort_keys=True))')
        result = subprocess.run([sys.executable, '-I', '-B', '-W', 'error', '-c', bootstrap],
                                cwd=view, capture_output=True, text=True)
        require(result.returncode == 0, 'unchanged current scientific custody validation failed:\n'+result.stdout+result.stderr)
        report['actual_unchanged_current_custody_validation'] = json.loads(result.stdout)
    require(validate_actual_release(root) == release, 'actual publication changed during scientific validation')
    return {'current_publication_validation': release, 'current_scientific_validation': report}
