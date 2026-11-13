import json
import pytest
from .check import verify,ROOT,HERE,traced_uses


def test_three_exact_roles_authenticate_original_scientific_and_provenance_bytes():
    result=verify()
    assert len(result['edges'])==3
    assert result['raw_legacy_documentary_policy_findings_retained']==3
    assert result['catch_all_docs_exception'] is False


def test_operative_provider_or_original_data_mutation_fails():
    registry=json.loads((HERE/'roles.json').read_text())
    row=registry['edges'][0]
    with pytest.raises(AssertionError,match='provider changed'):
        verify(override={row['provider']:(ROOT/row['provider']).read_bytes()+b'\n# mutation\n'})
    owner=json.loads((ROOT/row['owner']).read_text())
    operative=next(p for p in owner[row['legacy_map']] if p.endswith('.json') and not p.startswith('docs/'))
    with pytest.raises(AssertionError,match='custody changed'):
        verify(override={operative:b'{}'})


def test_legacy_prose_is_authenticated_and_no_wildcard_is_admitted():
    registry=json.loads((HERE/'roles.json').read_text());row=registry['edges'][0]
    with pytest.raises(AssertionError,match='presentation bytes changed'):
        verify(override={row['presentation_path']:b'changed presentation'})
    registry['edges'].append(dict(row,presentation_path='docs/another.md'))
    with pytest.raises(AssertionError):verify(registry)


def test_hidden_physical_use_of_note_is_rejected_by_read_trace():
    registry=json.loads((HERE/'roles.json').read_text());row=registry['edges'][-1]
    # Select the chiral owner explicitly rather than relying on registry order.
    row=next(x for x in registry['edges'] if x['provider'].endswith('chiral_candidate/run.py'))
    body=(ROOT/row['provider']).read_text()+'\ndef hidden_mass():\n    return float(NOTE.read_text())\n'
    with pytest.raises(AssertionError):traced_uses(body,row['presentation_path'],row['provider'])
