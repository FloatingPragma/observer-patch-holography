"""Historical inputs cannot masquerade as current tests or changed history."""
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from tools import check_retained_completion_custody as retained


def test_current_collection_excludes_only_named_history_and_keeps_current_tests(tmp_path):
    shutil.copyfile(retained.ROOT / 'pytest.ini', tmp_path / 'pytest.ini')
    sources = {
        'model/completion/historical/snapshot/test_archive.py': 'raise RuntimeError("archived bytes are not current modules")\n',
        'model/completion/test_custody_renewal.py': 'raise RuntimeError("dated suite must use its own checkout")\n',
        'model/completion/test_field_custody.py': 'raise RuntimeError("field suite needs its c7 view")\n',
        'model/completion/test_field_navigation.py': 'raise RuntimeError("field navigation needs its c7 view")\n',
        'model/completion/test_source_records_custody.py': 'def test_current_source_record_custody(): pass\n',
        'model/another_package/test_history.py': 'def test_unrelated_history_named_module_is_retained(): pass\n',
    }
    for name, text in sources.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    result = subprocess.run([sys.executable, '-m', 'pytest', '--collect-only', '-q', 'model'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'test_current_source_record_custody' in result.stdout
    assert 'test_unrelated_history_named_module_is_retained' in result.stdout
    assert '2 tests collected' in result.stdout


@pytest.mark.parametrize('name', list(retained.ORIGINAL_SHA256))
def test_historical_authentication_rejects_modified_checkers_tests_and_config(tmp_path, name):
    for source in retained.ORIGINAL_SHA256:
        path = tmp_path / source
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(subprocess.check_output(['git', 'show', retained.CHECKPOINT + ':' + source], cwd=retained.ROOT))
    retained.authenticate(tmp_path)
    path = tmp_path / name
    path.write_bytes(path.read_bytes() + b'\nchanged\n')
    with pytest.raises(ValueError, match='retained checker, test, authority or configuration changed'):
        retained.authenticate(tmp_path)


def test_retained_checkpoint_is_full_commit_and_original_suite_is_exactly_named():
    assert retained.CHECKPOINT == '0194fc685be3f7ba5f732681713c9cfd04e5290f'
    assert retained.TEST == 'model/completion/test_custody_renewal.py'
    assert retained.TEST_COUNT == 53
