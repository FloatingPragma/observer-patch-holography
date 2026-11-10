"""Artifact slot renewal retains every scientific inventory value."""
import copy,hashlib
import pytest
from tools import renew_source_record_manifests as current


def fixture(tmp_path):
    path=tmp_path/'model/pkg/result.json';path.parent.mkdir(parents=True);path.write_text('{"mass":"7/3","physical":false}\n')
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    old={'schema':'owned','files':{'result.json':'0'*64,'other.py':'b'*64},'scope':'finite only','gates':[False]}
    approved={'only_live_artifact_entries':[{'manifest':'model/pkg/manifest.json','artifact_key':'result.json',
        'target':'model/pkg/result.json','retained_pointer_sha256':'0'*64,'executed_current_sha256':sha}]}
    return old,approved,path


def test_exact_slot_only_and_retained_input_unchanged(tmp_path):
    old,approved,path=fixture(tmp_path);before=copy.deepcopy(old)
    value=current.project('model/pkg/manifest.json',old,approved,root=tmp_path)
    assert old==before
    value['files']['result.json']='0'*64
    assert value==old


@pytest.mark.parametrize('mutation',['retained_pointer','target_bytes','target_path'])
def test_unexecuted_artifact_or_path_is_rejected(tmp_path,mutation):
    old,approved,path=fixture(tmp_path)
    if mutation=='retained_pointer':old['files']['result.json']='f'*64
    elif mutation=='target_bytes':path.write_text('{"mass":"8/3"}')
    else:approved['only_live_artifact_entries'][0]['target']='other.json'
    with pytest.raises(ValueError):current.project('model/pkg/manifest.json',old,approved,root=tmp_path)


def test_current_exact_seven_slots_bind_actual_unchanged_owner_replays():
    result=current.execute()
    assert result['exact_artifact_slots']==7 and not result['scientific_statuses_changed']
