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
