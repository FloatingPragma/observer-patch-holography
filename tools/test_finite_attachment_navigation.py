"""Navigation maintenance cannot renew scientific acceptance or other sources."""
from copy import deepcopy
import pytest
from tools import strict_json
from tools import renew_finite_attachment_navigation as owner


def documents():return {p:strict_json.load(owner.ROOT/p)for p in owner.FILES}


def test_actual_map_and_class_cold_navigation():
    expected=owner.build()
    for path,text in expected.items():assert (owner.ROOT/path).read_text()==text
    owner.validate_class(strict_json.load(owner.ROOT/owner.CLASS),owner.ROOT)


@pytest.mark.parametrize('target',owner.FILES)
def test_every_scientific_document_rejects_unreviewed_change(target):
    changed=documents();changed[target]['schema']='HOSTILE_CHANGED_SCIENTIFIC_CONTRACT'
    with pytest.raises(ValueError,match='scientific contents changed'):owner.unchanged_science(changed)


def test_only_declared_navigation_and_pointer_fields_are_renewable():
    changed=documents()
    for key in owner.FIELDS:changed[owner.MAP]['existing_navigation_coverage'][key]='stale'
    changed[owner.REVIEW]['source_bindings'][owner.MAP]='0'*64
    changed[owner.CLASS]['adoption']['independent_review']['sha256']='1'*64
    owner.pointer_update(changed[owner.PROGRAM],{p:'2'*64 for p in (owner.MAP,owner.REVIEW,owner.CLASS)})
    owner.unchanged_science(changed)
    source=next(p for p in changed[owner.REVIEW]['source_bindings']if p!=owner.MAP)
    changed[owner.REVIEW]['source_bindings'][source]='3'*64
    with pytest.raises(ValueError,match='scientific contents changed'):owner.unchanged_science(changed)


def test_gate_statuses_and_finite_classifications_are_not_refreshable():
    changed=documents();changed[owner.MAP]['canonical_premises_discharged']=['PR-23']
    with pytest.raises(ValueError,match='scientific contents changed'):owner.unchanged_science(changed)
    changed=documents();changed[owner.PROGRAM]['stages'][0]['status']='UNREVIEWED'
    with pytest.raises(ValueError,match='scientific contents changed'):owner.unchanged_science(changed)


def _presentation_root(tmp_path):
    for name in (owner.CANDIDATE,owner.PRESENTATION_PARENT):
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes((owner.ROOT/name).read_bytes())
    return tmp_path


def test_exact_inverse_and_operative_ast_control():
    owner.candidate_presentation_custody(documents())
    raw=(owner.ROOT/owner.CANDIDATE).read_bytes()
    import hashlib
    assert hashlib.sha256(raw.replace(owner.CANDIDATE_PRESENTATION_ADDITION,b'',1)).hexdigest()==owner.CANDIDATE_OLD_SHA256


@pytest.mark.parametrize('target',('candidate','parent_scope','parent_physics'))
def test_changed_source_or_presentation_scope_rejects(tmp_path,target):
    root=_presentation_root(tmp_path)
    if target=='candidate':
        path=root/owner.CANDIDATE;path.write_bytes(path.read_bytes()+b'\nOPERATIVE_CHANGE = 1\n')
    else:
        path=root/owner.PRESENTATION_PARENT;raw=path.read_bytes()
        if target=='parent_scope':
            raw=raw.replace(b"('extra/yang_mills_gap_clay_problem.tex',)",b"('model/particle_yang_mills_gap_contract/PROOF.md',)")
        else:raw=raw.replace(b"ratio='1/2'",b"ratio='1/3'",1)
        path.write_bytes(raw)
    with pytest.raises(ValueError,match='presentation.*changed'):
        owner.unchanged_science(documents(),root)


@pytest.mark.parametrize('mutation',('stale_hash','wrong_row','wrong_scope','physical_status'))
def test_only_exact_current_row_hash_is_renewable(mutation):
    changed=documents();stage=changed[owner.PROGRAM]['stages'][6];row=stage['progress_evidence'][38]
    if mutation=='stale_hash':row['sha256']=owner.CANDIDATE_OLD_SHA256
    elif mutation=='wrong_row':stage['progress_evidence'][37],stage['progress_evidence'][38]=row,stage['progress_evidence'][37]
    elif mutation=='wrong_scope':row['scope']='All physical particle masses are proved.'
    else:stage['status']='ACCEPTED'
    with pytest.raises(ValueError):owner.unchanged_science(changed)
