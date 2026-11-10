"""Reject unreviewed source, state, reference, clock or status replacements."""
import copy
import json
import pytest
from . import source_records_navigation as current


def fixture():
    approved=current.authority()
    baseline={p:json.loads((current.ROOT/approved['baseline_files'][p]['retained_path']).read_text())
              for p in current.JSON_FILES}
    return baseline,approved


def test_current_projection_and_original_scope_are_exact():
    baseline,approved=fixture();before=copy.deepcopy(baseline)
    result=current.project(baseline,approved)
    assert baseline==before
    current.normalize(result,approved)
    assert current.validate()['claim_count']==392
    assert len(approved['package_receipt_sha256'])==10
    assert len(approved['new_assumption_definitions'])==2
    assert not approved['physical_status_changes']
    assert not approved['release_id_transition_included']
    assert all(not p.endswith('.raw') for p in approved['new_evidence_sha256'])
    for index in (8, 18):
        sources=result[current.TRACKING]['model_supply']['groups'][index]['sources']
        assert all(isinstance(row, dict) and isinstance(row['label'], str)
                   and isinstance(row['path'], str) for row in sources)


@pytest.mark.parametrize('mutation',['status','prior_evidence','extra_evidence','assumption',
                                    'binding','scope','species','selector','history','raw_artifact'])
def test_unlisted_semantics_are_rejected(mutation):
    baseline,approved=fixture();value=current.project(baseline,approved)
    claim=value[current.REGISTRY]['claims'][326]
    group=value[current.TRACKING]['model_supply']['groups'][18]
    if mutation=='status':group['status']='closed'
    elif mutation=='prior_evidence':claim['evidence'].pop(0)
    elif mutation=='extra_evidence':claim['evidence'].append('unreviewed.py')
    elif mutation=='assumption':claim['assumptions'].pop(0)
    elif mutation=='binding':value[current.TRACKING]['model_supply']['contract_coverage']['bindings'][1]['assumptions'].pop()
    elif mutation=='scope':group['scope']='physical'
    elif mutation=='species':claim['statement']+=' The Standard Model mass is established.'
    elif mutation=='selector':claim['statement']=claim['statement'].replace('separate choices','the same eigenvalue')
    elif mutation=='history':claim['statement']=claim['statement'].replace('conditional-state','unconditional-state')
    elif mutation=='raw_artifact':claim['evidence'].append('model/particle_construction_scope/retained/CM/result.json.raw')
    with pytest.raises(ValueError):current.normalize(value,approved)


def test_retained_authority_bytes_are_authenticated(tmp_path):
    baseline,approved=fixture()
    for path in [current.AUTHORITY]+[approved[k]['path'] for k in ('predecessor_custody_authority','predecessor_navigation_authority')]:
        dest=tmp_path/path;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((current.ROOT/path).read_bytes())
    path=tmp_path/approved['predecessor_navigation_authority']['path'];path.write_bytes(path.read_bytes()+b' ')
    with pytest.raises(ValueError,match='retained field authority changed'):current.authority(tmp_path)
