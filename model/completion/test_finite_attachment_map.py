"""The model map must retain complete ownership, proof domains and controls."""
from copy import deepcopy
from pathlib import Path
import pytest
from tools import strict_json
from .finite_attachment_map import PATH,ROOT,validate,predecessor_review


def packet():
    return strict_json.load(ROOT/PATH)


def test_complete_canonical_ownership_and_realized_finite_domains():
    result=validate(packet())
    assert result['canonical_premises']==84
    assert result['finite_obligations_realized']==34
    assert result['canonical_premises_discharged']==0
    assert result['gate_acceptance_owner']=='claims/completion_program.json'


@pytest.mark.parametrize('damage',[
    'missing_owner','wrong_issue','missing_premise','stale_theorem','missing_anchor',
    'missing_clock','missing_control','canonical_promotion','unproved_closure',
    'missing_attachment','wrong_classification','changed_existing_coverage'])
def test_review_mutations_fail_closed(damage):
    p=deepcopy(packet())
    if damage=='missing_owner': p['canonical_ownership'][0]['owners']=[]
    elif damage=='wrong_issue': p['canonical_ownership'][0]['issues']=[9999]
    elif damage=='missing_premise': p['canonical_ownership'].pop()
    elif damage=='stale_theorem':
        key=next(iter(p['source_bindings']));p['source_bindings'][key]='0'*64
    elif damage=='missing_anchor': next(iter(p['reviewed_sources'].values()))['anchors']=['not an actual theorem']
    elif damage=='missing_clock': del p['attachments'][0]['interface']['clock']
    elif damage=='missing_control': p['attachments'][0]['countercontrols']=[]
    elif damage=='canonical_promotion': p['canonical_premises_discharged']=['PR-23']
    elif damage=='unproved_closure': p['issue_closure_review']['1137']['close_issue']=True
    elif damage=='missing_attachment': p['attachments'].pop()
    elif damage=='wrong_classification': next(a for a in p['attachments'] if a['id']=='FA-12')['classification']='realized_on_claimed_domain'
    else: p['existing_navigation_coverage']['contract_coverage_sha256']='0'*64
    with pytest.raises(ValueError):validate(p)


def test_algebraic_cp_interface_does_not_claim_a_new_lean_inhabitant_or_spacetime():
    p=packet()
    regional=next(a for a in p['attachments'] if a['id']=='FA-02')
    assert regional['new_lean_inhabitant'] is False
    assert 'nonempty' in regional['interface']['algebra']
    assert 'Kraus' in regional['proof_kind']
    assert 'positivity' in p['reviewed_sources']['cp_net']['boundary']
    assert all(a.get('physical_identification') is False for a in p['attachments'])


def test_no_implicit_deletion_or_theorem_error_from_failed_transfers():
    p=packet()
    assert p['cleanup']['deleted_sources']==[]
    assert p['cleanup']['frozen_receipts_rewritten']==[]
    assert p['cleanup']['false_general_theorems_found']==[]
    corrected=next(a for a in p['attachments'] if a['id']=='FA-09')
    assert corrected['classification']=='needs_correction'
    assert 'identification' in corrected['title']
    assert p['issue_closure_review']['1137']['required_predecessors']==['CP-00','CP-01']


def test_open_reference_gate_can_have_all_predecessors_accepted():
    p = packet()
    program = strict_json.load(ROOT/'claims/completion_program.json')
    for issue, row in p['issue_closure_review'].items():
        expected = predecessor_review(program, issue)
        assert all(row[key] == value for key, value in expected.items())
    row = p['issue_closure_review']['1138']
    assert row['required_predecessors'] == ['CP-01', 'CP-13']
    assert row['unaccepted_predecessors'] == []
    assert row['unfulfilled_acceptance'] and row['close_issue'] is False
    validate(p)


def test_new_finite_interfaces_keep_reference_scattering_and_physical_scope_separate():
    p = packet()
    rows = {a['id']: a for a in p['attachments']}
    assert 'reference_recovery' in rows['FA-18']['reviewed_sources']
    assert 'pure half-filled ground' in rows['FA-21']['interface']['state']
    assert 'supplied logical access' in rows['FA-21']['interface']['operation']
    assert 'eta>0' in rows['FA-22']['interface']['clock']
    assert 'on-shell' in rows['FA-22']['interface']['refinement']
    assert p['physical_comparison_count_added'] == 0
    assert p['issue_closure_review']['1147']['close_issue'] is False
    assert 'particle_comparison' in p['reviewed_sources']


@pytest.mark.parametrize('field,value', [
    ('unaccepted_predecessors', ['CP-01', 'CP-13']),
    ('required_predecessors', []),
])
def test_stale_predecessor_status_and_removed_dependencies_are_rejected(field,value):
    p = deepcopy(packet())
    p['issue_closure_review']['1138'][field] = value
    with pytest.raises(ValueError, match='prerequisite status or ownership'):
        validate(p)
