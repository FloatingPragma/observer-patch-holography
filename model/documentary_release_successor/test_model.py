import copy
from contextlib import contextmanager
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

from . import model


def fixture():
    old = json.loads((model.ROOT/model.RETAINED).read_text())
    new = copy.deepcopy(old)
    new['book']['sha256'] = 'f'*64
    new['book']['size_bytes'] += 1
    reviewed = [{'json_path': list(path), 'before': a, 'after': b}
                for path, a, b in model.differences(old, new)]
    return old, new, reviewed


def test_exact_documentary_change_keeps_release_membership_paths_and_date():
    old, new, reviewed = fixture()
    before = copy.deepcopy(old)
    assert len(model.compare_manifests(old, new, reviewed)) == 2
    assert old == before


@pytest.mark.parametrize('damage', ['release', 'date', 'extra_key', 'removed_paper', 'unlisted_paper',
                                    'pdf_path', 'extra_field', 'false_hash', 'false_size', 'false_review'])
def test_unlisted_manifest_change_or_fabricated_review_is_rejected(damage):
    old, new, reviewed = fixture()
    if damage == 'release': new['release_id'] = 'r2065'
    if damage == 'date': new['released_at'] = 'Unreviewed date'
    if damage == 'extra_key': new['physical_gate_accepted'] = True
    if damage == 'removed_paper': new['papers'].pop(next(iter(new['papers'])))
    if damage == 'unlisted_paper': new['papers']['observers_are_all_you_need']['sha256'] = 'a'*64
    if damage == 'pdf_path': new['book']['pdf_path'] = 'wrong.pdf'
    if damage == 'extra_field': new['book']['physical_source_verified'] = True
    if damage == 'false_hash': new['book']['sha256'] = 'fabricated'
    if damage == 'false_size': new['book']['size_bytes'] = -1
    if damage == 'false_review': reviewed[0]['after'] = 'e'*64
    with pytest.raises(ValueError): model.compare_manifests(old, new, reviewed)


def test_only_four_artifact_rows_and_two_leaf_types_are_authorized():
    assert len(model.ALLOWED_LEAVES) == 8
    assert {path[-1] for path in model.ALLOWED_LEAVES} == {'sha256', 'size_bytes'}
    assert {path[:-1] for path in model.ALLOWED_LEAVES} == set(model.ALLOWED_ROWS)


def test_view_mutation_or_live_mutation_is_rejected(tmp_path):
    root, view = tmp_path/'root', tmp_path/'view'
    for folder in (root, view):
        (folder/'paper').mkdir(parents=True)
        (folder/'source.py').write_text('answer = 1\n')
        (folder/model.MANIFEST).write_text('current')
    snapshot = {p: model.sha(root/p) for p in ('source.py',model.MANIFEST)}
    (view/model.MANIFEST).write_text('retained')
    retained = model.sha(view/model.MANIFEST)
    model.assert_view_equal(root, view, snapshot, retained)
    (view/'source.py').write_text('answer = 2\n')
    with pytest.raises(ValueError, match='private current-source'):
        model.assert_view_equal(root, view, snapshot, retained)
    (view/'source.py').write_text('answer = 1\n')
    (root/'source.py').write_text('answer = 3\n')
    with pytest.raises(ValueError, match='current source changed'):
        model.assert_view_equal(root, view, snapshot, retained)


def test_manifest_is_unlinked_before_replacement_and_current_science_stays_exact(tmp_path, monkeypatch):
    root = tmp_path/'root'
    for name, value in ((model.MANIFEST,'current'), (model.RETAINED,'retained'), ('source.py','answer = 1\n')):
        target = root/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(value)
    paths = sorted([model.MANIFEST,model.RETAINED,'source.py'])
    approved = {'retained_manifest_sha256':model.sha(root/model.RETAINED),
                'current_manifest_sha256':model.sha(root/model.MANIFEST)}
    monkeypatch.setattr(model,'authority',lambda path:approved)
    monkeypatch.setattr(model,'current_paths',lambda path:paths)
    with model.scientific_view(root) as (view, report):
        assert (view/model.MANIFEST).read_text() == 'retained'
        assert (root/model.MANIFEST).read_text() == 'current'
        assert (view/'source.py').read_bytes() == (root/'source.py').read_bytes()
        assert report['only_substituted_path'] == model.MANIFEST
        assert report['publication_validity_of_retained_manifest_view_claimed'] is False
    assert (root/model.MANIFEST).read_text() == 'current'
    assert not view.exists()


def test_private_scientific_mutation_cannot_write_through_to_root(tmp_path,monkeypatch):
    root = tmp_path/'root'
    for name,value in ((model.MANIFEST,'current'),(model.RETAINED,'retained'),('source.py','answer = 1\n')):
        path=root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(value)
    paths=sorted([model.MANIFEST,model.RETAINED,'source.py'])
    approved={'retained_manifest_sha256':model.sha(root/model.RETAINED),
              'current_manifest_sha256':model.sha(root/model.MANIFEST)}
    monkeypatch.setattr(model,'authority',lambda path:approved)
    monkeypatch.setattr(model,'current_paths',lambda path:paths)
    with pytest.raises(ValueError,match='private current-source view changed'):
        with model.scientific_view(root) as (view,_):
            (view/'source.py').write_text('answer = 2\n')
            assert (root/'source.py').read_text() == 'answer = 1\n'
    assert (root/'source.py').read_text() == 'answer = 1\n'


def test_extra_non_python_data_file_in_private_view_is_rejected(tmp_path):
    root,view=tmp_path/'root',tmp_path/'view'
    for directory in (root,view):
        (directory/'paper').mkdir(parents=True)
        (directory/model.MANIFEST).write_text('manifest')
        (directory/'input.dat').write_text('original datum')
    snapshot={name:model.sha(root/name) for name in (model.MANIFEST,'input.dat')}
    (view/'unreviewed.dat').write_text('injected value')
    with pytest.raises(ValueError,match='file inventory changed'):
        model.assert_view_equal(root,view,snapshot,model.sha(root/model.MANIFEST))


def test_publication_change_between_initial_check_and_view_is_rejected(tmp_path,monkeypatch):
    observations=iter([{'manifest':'first'},{'manifest':'changed'}])
    monkeypatch.setattr(model,'validate_actual_release',lambda root:next(observations))
    @contextmanager
    def protected_view(root):
        yield tmp_path,{}
    monkeypatch.setattr(model,'scientific_view',protected_view)
    monkeypatch.setattr(model.subprocess,'run',lambda *args,**kwargs:SimpleNamespace(returncode=0,stdout='{}',stderr=''))
    with pytest.raises(ValueError,match='actual publication changed'):
        model.validate_current_science(tmp_path)


@pytest.mark.parametrize('leaf,bad', [('size_bytes',True),('size_bytes',-1),('sha256','fabricated')])
def test_even_reviewed_invalid_artifact_metadata_is_rejected(leaf,bad):
    old, new, _ = fixture()
    new['book'][leaf] = bad
    reviewed = [{'json_path':list(path),'before':a,'after':b} for path,a,b in model.differences(old,new)]
    with pytest.raises(ValueError): model.compare_manifests(old,new,reviewed)
