"""Reference damage and scope promotion must fail the conditional interface."""
from copy import deepcopy
import json

import pytest

from . import contract


def test_conditional_theorem_does_not_require_source_selection_or_gate_acceptance():
    value = contract.theorem_contract()
    checked = contract.validate_contract(json.loads(json.dumps(value)))
    assert checked['theorem_rows'] == 9
    assert checked['conditional_interfaces'] == 10
    assert checked['canonical_discharges_added'] == 0
    assert value['conditional_derivation_requires_source_selection'] is False
    assert value['conditional_derivation_requires_all_physical_gates_accepted'] is False
    assert value['verification']['lean_modules_recompiled_by_contract'] is False
    assert value['verification']['python_scientific_proofs_reexecuted_by_contract'] is False


@pytest.mark.parametrize('damage', [
    'missing_path', 'absolute_path', 'escaping_path', 'wrong_path_type',
    'missing_anchor', 'empty_anchors', 'unknown_premise', 'duplicate_premise',
    'omitted_index', 'unreviewed_claim', 'source_selection', 'boolean_as_integer',
    'omitted_matching', 'assume_standard_model', 'einstein_promotion'])
def test_reference_and_scope_mutations_fail_closed(damage):
    value = deepcopy(contract.theorem_contract())
    first = value['theorem_chain'][0]
    reference = first['references'][0]
    if damage == 'missing_path':
        reference['path'] = 'Lean/missing_conditional_theorem.lean'
    elif damage == 'absolute_path':
        reference['path'] = str(contract.ROOT / reference['path'])
    elif damage == 'escaping_path':
        reference['path'] = '../outside/theorem.lean'
    elif damage == 'wrong_path_type':
        reference['path'] = 17
    elif damage == 'missing_anchor':
        reference['anchors'][0] = 'theorem unproved_physical_standard_model'
    elif damage == 'empty_anchors':
        reference['anchors'] = []
    elif damage == 'unknown_premise':
        first['canonical_premises'].append('PR-999')
    elif damage == 'duplicate_premise':
        first['canonical_premises'].append(first['canonical_premises'][0])
    elif damage == 'omitted_index':
        value['canonical_premise_register']['ids'].pop()
    elif damage == 'unreviewed_claim':
        first['conclusion'] = 'The finite source proves every physical Standard Model particle.'
    elif damage == 'source_selection':
        value['source_selection_derived'] = True
    elif damage == 'boolean_as_integer':
        value['conditional_derivation_requires_source_selection'] = 0
    elif damage == 'omitted_matching':
        value['conditional_interfaces'] = [r for r in value['conditional_interfaces']
                                           if r['id'] != 'same_source_continuum_matching']
    elif damage == 'assume_standard_model':
        row = next(r for r in value['conditional_interfaces']
                   if r['id'] == 'same_source_continuum_matching')
        row['requirements'] = ['Assume that the source limit is the Standard Model.']
    else:
        value['einstein_regional_premises_instantiated'] = True
    with pytest.raises(ValueError):
        contract.validate_contract(value)


def test_contract_returns_detached_data():
    first = contract.theorem_contract()
    first['theorem_chain'][0]['references'][0]['anchors'].clear()
    assert contract.theorem_contract()['theorem_chain'][0]['references'][0]['anchors']


def test_register_ownership_and_duplicate_canonical_ids_are_checked(tmp_path):
    value = contract.theorem_contract()
    target = tmp_path / contract.REGISTER
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps({'rows': [{'id': 'PR-09'}, {'id': 'PR-09'}]}))
    with pytest.raises(ValueError, match='unique strings'):
        contract.validate_contract(value, tmp_path)


def test_symlink_cannot_import_a_theorem_outside_the_repository(tmp_path):
    outside = tmp_path / 'outside.lean'
    outside.write_text('theorem fabricated')
    root = tmp_path / 'repository'
    root.mkdir()
    (root / 'theorem.lean').symlink_to(outside)
    with pytest.raises(ValueError, match='external referenced path'):
        contract._local_file(root.resolve(), 'theorem.lean')
