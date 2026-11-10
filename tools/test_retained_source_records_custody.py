"""Exact historical routing preserves old suites and live successor coverage."""
import shutil
import subprocess
import sys
import pytest
from tools import check_retained_source_records_custody as retained


@pytest.mark.parametrize('name', list(retained.ORIGINAL_SHA256))
def test_retained_7fc_inputs_reject_every_modified_file(tmp_path, name):
    for path in retained.ORIGINAL_SHA256:
        target = tmp_path/path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(subprocess.check_output(['git', 'show', retained.CHECKPOINT+':'+path], cwd=retained.ROOT))
    retained.authenticate(tmp_path)
    path = tmp_path/name
    path.write_bytes(path.read_bytes()+b'\nchanged\n')
    with pytest.raises(ValueError, match='retained checker, test, authority or configuration changed'):
        retained.authenticate(tmp_path)


def test_current_collection_routes_only_the_explicit_retained_files(tmp_path):
    shutil.copyfile(retained.ROOT/'pytest.ini', tmp_path/'pytest.ini')
    for name in retained.TESTS:
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('raise RuntimeError("requires sealed 7fc source")\n')
    current = ['model/completion/test_field_current_records_navigation.py',
               'model/unrelated/test_history.py', 'tools/test_field_current_documentary_owners.py']
    for i, name in enumerate(current):
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('def test_live_'+str(i)+'(): pass\n')
    result = subprocess.run([sys.executable, '-m', 'pytest', '-q', '--collect-only', 'model', 'tools'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert '3 tests collected' in result.stdout
    assert all('test_live_'+str(i) in result.stdout for i in range(3))


def test_retained_suite_is_exact_and_distinct_from_integration_baseline():
    assert retained.CHECKPOINT == '7fc48ae176316cb6d2b7ef6decb2060dab9f7db8'
    assert retained.TEST_COUNT == 73
    assert len(retained.TESTS) == 8
    config = (retained.ROOT/'pytest.ini').read_text()
    assert all('--ignore='+name in config for name in retained.TESTS)
    assert '--ignore=model/completion/test_field_current_records_navigation.py' not in config
    assert '--ignore=model/completion/test_field_current_records_custody.py' not in config
