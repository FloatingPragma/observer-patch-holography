"""Exact source derivation and adversarial source/response-boundary controls."""
from __future__ import annotations
import copy,json,shutil,sys
from pathlib import Path
import pytest
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import source_w5_response_constraints as producer
import verify_source_w5_response_constraints as verifier

@pytest.fixture(scope='module')
def packet():
    value=producer.build()
    assert value==json.loads(producer.DEFAULT_OUTPUT.read_text())
    return value

def test_exact_producer_and_independent_verifier(packet):
    assert verifier.verify(packet)

@pytest.mark.parametrize('path,value',[
 (('carrier','dimension'),6),
 (('carrier','D3_line_complement_dimension'),5),
 (('carrier','normalized_native_projector_weights',0),'1/6'),
 (('carrier','normalized_native_projector_weights',1),'4/15'),
 (('carrier','D3_multiplicities','sign'),1),
 (('carrier','character_norm'),'2'),
 (('carrier','D3_line_projector',0,0),'1/5'),
 (('carrier','quadrupole_Gram',0,0),'1'),
 (('carrier','invariant_Hermitian_functional_dimensions','A5'),5),
 (('operators','native_W5_laplacian_eigenvalues'),{'0':1,'3':4}),
 (('operators','D3_transposition_EndW5_spectrum'),{'0':1,'3':4,'6':1}),
 (('operators','channel_Gram',0,0),'5'),
 (('operators','channel_Gram_inverse',0,1),'0'),
 (('operators','central_rate_eigenvalue_coefficients',1,1),0),
 (('operators','faithful_group_heat_channel_identity'),False),
 (('proposal','symmetrized_covariance',0,0),'2/3'),
 (('proposal','raw_moment_coefficients','fourth_cumulant_S2_squared'),'0'),
 (('proposal','symmetrized_unit_direction_examples','one_high_five_low','sixth_cumulant'),'0'),
 (('proposal','native_accepted_step_witnesses',0,'accepted_seams'),30),
 (('proposal','native_accepted_step_witnesses',0,'covariance',0,0),'1/3'),
 (('proposal','native_accepted_step_witnesses',1,'mean',0),'0'),
 (('color','A_rank'),4),
 (('color','A_nonzero_eigenvalue'),'1'),
 (('color','normalized_trace'),'4/5'),
 (('color','normalized_HS_square'),'4/15'),
 (('color','scalar_cumulants',1),'16/45'),
 (('color','twirled_rank_r_HS_square_coefficient'),'1'),
 (('color','conditional_Gaussian_coefficients','d432_rank1_twirl_before_response'),'1/432'),
 (('exact_checks','A5_irreducibility'),1),
 (('status',),'SOURCE_SELECTED_QUARK_MASSES'),
])
def test_independent_rejection_of_mathematical_corruption(packet,path,value):
    changed=copy.deepcopy(packet);cursor=changed
    for key in path[:-1]:cursor=cursor[key]
    cursor[path[-1]]=value
    with pytest.raises(ValueError):verifier.verify(changed)

@pytest.mark.parametrize('guard',list(producer.read_spec()['guards']))
def test_rejects_each_physical_promotion(packet,guard):
    changed=copy.deepcopy(packet);changed['specification']['guards'][guard]=True
    with pytest.raises(ValueError,match='scope|specification'):verifier.verify(changed)

def test_false_guard_cannot_be_replaced_by_numeric_zero(packet):
    changed=copy.deepcopy(packet);changed['specification']['guards']['quark_mass_prediction']=0
    with pytest.raises(ValueError):verifier.verify(changed)

@pytest.mark.parametrize('field',[
 'D3_channel_interpretation','D3_scope','trace_law','coordinate','proof_grade'])
def test_response_assumptions_cannot_silently_change(packet,field):
    changed=copy.deepcopy(packet);changed['specification'][field]='source-selected without further assumptions'
    with pytest.raises(ValueError):verifier.verify(changed)

def test_rejects_source_pin_and_hidden_receipt_input(packet):
    changed=copy.deepcopy(packet);changed['source_pins']['code/particles/data/quark_targets.json']='0'*64
    with pytest.raises(ValueError):verifier.verify(changed)

@pytest.mark.parametrize('injection',[
 '\nimport conditional_quark_mass_replay\n',
 '\nfrom compare_conditional_quark_masses import build\n',
 '\nreference = "code/particles/data/pdg_2025_conditional_quark_comparison.json"\n'])
def test_target_import_and_read_surface_rejected(injection):
    with pytest.raises(ValueError):verifier.audit_source((producer.REPO/producer.PRODUCER_REL).read_text()+injection)

def test_changed_spec_and_added_fitted_parameter_are_rejected(tmp_path):
    d=producer.read_spec();d['quark_scale_fit']='1.01';p=tmp_path/'spec.json';p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='reviewed law'):producer.read_spec(p)

def test_tampered_geometry_cannot_be_repaired_by_rehashing(packet,tmp_path):
    for name in verifier.SOURCE_FILES+[verifier.SPEC_REL,verifier.PRODUCER_REL]:
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(verifier.REPO/name,target)
    p=tmp_path/verifier.SOURCE_FILES[0];d=json.loads(p.read_text());d['carrier']['edges'].pop();p.write_text(json.dumps(d))
    changed=copy.deepcopy(packet);changed['source_pins'][verifier.SOURCE_FILES[0]]=verifier.sha(p)
    with pytest.raises(ValueError,match='incidence'):verifier.verify(changed,repo=tmp_path)
