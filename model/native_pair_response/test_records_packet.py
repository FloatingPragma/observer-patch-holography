from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import pytest
from . import response as r,records,packet


@pytest.mark.parametrize('q',(2,80))
def test_actual_same_q_body_and_global_tree_readback(q):
    value=records.snapshot_program(0,q,2)
    initial=sum(1<<(i*(240//q)) for i in range(q))
    result=records.trace_snapshots(value,initial,((0,1),(0,1)),(0,1)+(0,)*239)
    assert result['decoded']['body_masks']==(initial,(initial&~1)|2,initial)
    assert records.decode_snapshots(value,result['records'])==result['decoded']
    assert result['weight']==F(1,1260)**2*F(1,2)**239
    assert len(result['archive'])==240 and len(result['memory'])==1
    wrong=deepcopy(value);wrong['parameters']['repair_probability']='1/2';wrong['declared_resources']=records.source._native_cost(wrong)
    with pytest.raises(ValueError):records.decode_snapshots(wrong,result['records'])
    altered=(replace(result['records'][0],coin=False),)+result['records'][1:]
    with pytest.raises(ValueError):records.decode_snapshots(value,altered)
    with pytest.raises(ValueError):records.decode_snapshots(value,result['records'][:-1])


def test_native_admission_and_host_stops_are_distinct_from_kernel_queries():
    f,probes=r.federation_menu(0);value=records.program(0,80,3,probes[0],probes[2])
    assert value['declared_resources']['attempts_upper']==7
    initial=sum(1<<i for i in range(80));rho={(initial,0,initial,0):F(1)}
    stopped=records.source.execute(value,rho,max_branches=1)
    assert stopped.status=='RESOURCE_STOP' and stopped.attempted_steps==0
    assert stopped.branches[0].density==rho
    with pytest.raises(ValueError):records.program(0,80,records.source.MAX_ATTEMPTS,probes[0],probes[2])
    with pytest.raises(ValueError):records.snapshot_program(0,80,records.source.MAX_ATTEMPTS)
    value=records.snapshot_program(2,1280,1)
    assert value['declared_resources']['AST_nodes']<=records.source.MAX_NODES
    assert value['declared_resources']['attempts_upper']==3840


@pytest.mark.parametrize('bad',(True,1.0,-1))
def test_body_count_rejects_ambiguous_types(bad):
    with pytest.raises(ValueError):records.snapshot_program(0,2,bad)


@pytest.fixture(scope='module')
def candidate():return packet.produce()


def test_owner_receipt_replays_and_keeps_scientific_scope(candidate):
    assert packet.verify()
    assert candidate['scope']['same_working_q_reference_and_native_attempt_clock']
    assert not candidate['scope']['physical_clock_Lorentzian_QFT_or_particle_masses_derived']
    assert candidate['independent_controls']['path4_one_attempt']['difference']=='1/12'
    assert candidate['native_record_access']['complete_pair_slots_reused']==625


@pytest.mark.parametrize('kind',('boolean','float','collision','hash','scope'))
def test_strict_receipt_science_and_type_mutations_fail_closed(candidate,tmp_path,monkeypatch,kind):
    monkeypatch.setattr(packet,'produce',lambda:candidate);value=deepcopy(candidate)
    if kind=='boolean':value['native_record_access']['expanded_control_attempts']=False
    elif kind=='float':value['native_record_access']['complete_pair_slots_reused']=625.0
    elif kind=='collision':value['independent_controls']['path4_one_attempt']['difference']='0'
    elif kind=='hash':del value['source_bindings'][next(iter(value['source_bindings']))]
    else:value['scope']['collision_survivors_postselected']=True
    path=tmp_path/'bad.json';path.write_text(packet.canonical(value))
    with pytest.raises(ValueError):packet.verify(path)


def test_duplicate_json_keys_fail_before_owner_replay(candidate,tmp_path):
    path=tmp_path/'duplicate.json';path.write_text('{"schema":"wrong",'+packet.canonical(candidate)[1:])
    with pytest.raises(ValueError):packet.verify(path)
