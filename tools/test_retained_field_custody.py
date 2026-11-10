"""The field checkpoint uses its own sources, configuration and test inventory."""
import subprocess
from pathlib import Path
import pytest
from tools import check_retained_field_custody as retained


@pytest.mark.parametrize('name',list(retained.ORIGINAL_SHA256))
def test_c7_retained_inputs_cannot_be_changed(tmp_path,name):
    for relative in retained.ORIGINAL_SHA256:
        path=tmp_path/relative;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(subprocess.check_output(['git','show',retained.CHECKPOINT+':'+relative],cwd=retained.ROOT))
    retained.authenticate(tmp_path)
    path=tmp_path/name;path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='retained checker, test, authority or configuration changed'):
        retained.authenticate(tmp_path)


def test_c7_suite_is_bounded_and_distinct_from_current_suite():
    assert retained.CHECKPOINT=='c7b09261853f5a3c5a67e21294c3d9511ed59fc7'
    assert retained.TEST_COUNT==31
    assert retained.TESTS==('model/completion/test_field_custody.py','model/completion/test_field_navigation.py')
    config=(retained.ROOT/'pytest.ini').read_text()
    for test in retained.TESTS:assert '--ignore='+test in config
    assert '--ignore=model/completion/test_source_records_custody.py' not in config
