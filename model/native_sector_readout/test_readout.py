from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import pytest

from . import readout as r,packet


@pytest.fixture(scope='module')
def local():return r.program()


@pytest.fixture(scope='module')
def trace(local):return r.trace_word(local,3,tuple(i%2 for i in range(11)))


def test_all_small_number_sectors_and_unknown_total_complete_instrument():
    value=packet.small_full_instrument()
    assert value['formal_words']==625 and value['unknown_total_rank_two_words']==4


def test_unknown_total_empty_full_reference_coherence_is_preserved():
    family=r.native.Family((2,3,5),((0,1),(1,2)),(1,1))
    slots=tuple(next(i for i,(rec,w) in enumerate(r.alphabet(family,e))
                     if w and rec.coin==0 and rec.label==0) for e in r.tree(family))
    rho={(0,0,0,0):F(1,2),(7,1,7,1):F(1,2),(0,0,7,1):F(1,2),(7,1,0,0):F(1,2)}
    weight,pairs,_=r.branch(family,slots)
    assert pairs==((0,0),(7,7))
    assert r.apply_branch(family,slots,rho)=={k:weight*v for k,v in rho.items()}
    assert r.branch(family,slots,0)[1]==((0,0),)
    assert r.branch(family,slots,3)[1]==((7,7),)


def test_known_total_reads_all_local_columns_and_coin_patterns(local):
    family=r.cone.source_family(local)
    for initial in (3,5,65,(1<<10)|(1<<11)):
        for coins in ((0,)*11,(1,)*11,tuple(i%2 for i in range(11))):
            h=r.trace_word(local,initial,coins)
            assert r.decode(family,h['records'],2)==((initial,h['terminal_mask']),)
            assert r.decode(family,h['records'])==((initial,h['terminal_mask']),)
            assert r.authenticate_trace(local,h)


def test_sparse_reference_map_matches_literal_native_kernel(local,trace):
    rho={(3,0,3,0):F(1,3),(3,0,3,1):F(-1,7),(3,1,3,0):F(-1,7),(5,1,5,1):F(2,3)}
    out=r.continue_raw(local,trace['slots'],rho)
    assert out['density']==r.apply_branch(r.cone.source_family(local),trace['slots'],rho,2)
    assert len(out['archive'])==10 and len(out['memory'])==1


def test_detector_authentication_retains_actual_fifo_capacity(local):
    value=r.source.program(commands=local['commands'],capacity=3)
    h=r.trace_word(value,3,(0,)*11)
    assert len(h['archive'])==8 and len(h['memory'])==3
    assert r.authenticate_trace(value,h)


@pytest.mark.parametrize('total',[0,1,17,120,239,240])
def test_full_connected_source_readout_without_global_basis_enumeration(total):
    value=r.program(level=0,total=total);mask=(1<<total)-1
    h=r.trace_word(value,mask,tuple(i%2 for i in range(239)))
    assert r.authenticate_trace(value,h)
    assert h['weight']==F(1,2**239)


@pytest.mark.parametrize('field',['program_digest','weight','initial_mask','terminal_mask','archive','memory','records'])
def test_authenticated_detector_trace_rejects_mutations(local,trace,field):
    altered=deepcopy(trace)
    if field=='program_digest':altered[field]='0'*64
    elif field=='weight':altered[field]=F(1)
    elif field in ('initial_mask','terminal_mask'):altered[field]=True
    elif field in ('archive','memory'):altered[field]=()
    else:altered[field]=(replace(altered[field][0],coin=bool(altered[field][0].coin)),)+altered[field][1:]
    with pytest.raises(ValueError):r.authenticate_trace(local,altered)


def test_zero_addresses_and_wrong_total_are_retained(local,trace):
    family=r.cone.source_family(local);slots=(0,)+trace['slots'][1:]
    assert r.branch(family,slots,2)[1]==()
    assert r.continue_raw(local,slots,{(3,0,3,0):F(1)})['density']=={}
    assert r.decode(family,trace['records'],1)==()
    with pytest.raises(ValueError):r.decode(family,trace['records'],True)
    with pytest.raises(ValueError):r.branch(family,(True,)+slots[1:],2)
    with pytest.raises(ValueError):r.branch(family,slots,True)


def test_terminal_readback_recovers_full_history_and_rejects_forged_currents(local,trace):
    family=r.cone.source_family(local)
    states=r.reverse_currents(family,trace['terminal_mask'],trace['records'])
    assert states[0]==trace['initial_mask'] and len(states)==12
    broken=(replace(trace['records'][0],current=999),)+trace['records'][1:]
    with pytest.raises(ValueError):r.reverse_currents(family,trace['terminal_mask'],broken)


def test_tree_record_types_and_source_digest_cannot_be_relabelled(local,trace):
    family=r.cone.source_family(local)
    broken=(replace(trace['records'][0],seam=list(trace['records'][0].seam)),)+trace['records'][1:]
    with pytest.raises(ValueError):r.decode(family,broken,2)
    other=r.source.program(commands=local['commands'],capacity=3)
    with pytest.raises(ValueError):r.authenticate_trace(other,trace)
    broken=(replace(trace['records'][0],raw=True),)+trace['records'][1:]
    with pytest.raises(ValueError):r.reverse_currents(family,trace['terminal_mask'],broken)


def test_actual_mixed_body_readback_needs_no_initial_census():
    value=r.source.program(commands=[
        {'op':'repair','selector':{'kind':'fixed','seam':[0,1]},'law':'E','repeat':1},
        {'op':'repair','selector':{'kind':'fixed','seam':[1,7]},'law':'C','repeat':1},
        {'op':'hold','repeat':1}])
    rho=r.source_density(value,{(3,0,3,0):F(1)});records=[];memory=archive=()
    for command in value['commands']:
        out=next(out for out in r.source._outcomes(value,command,rho,memory,archive)
                 if r.source.trace(out.density)>0)
        rho,memory,archive=out.density,out.memory,out.archive;records.append(out.record)
    terminal=next(iter(r.cone.source_masks(value,rho)))[0]
    states=r.reverse_currents(r.cone.source_family(value),terminal,tuple(records))
    assert states[0]==3 and states[-1]==terminal and len(states)==4


def test_no_motion_lower_bound_covers_actual_adaptive_controller():
    value=r.source.program(commands=[{'op':'repair','selector':{'kind':'last_current','positive':[0,5],'otherwise':[0,1]},'law':'E','repeat':3}])
    out=r.source.execute(value,{(0,0,0,0):F(1)})
    stay=sum(block.density.get((0,0,0,0),F()) for block in out.branches)
    bound=r.no_reset_bound(value)
    assert bound['minimum_coin']=='8/33' and stay>=F(bound['no_motion_probability_lower'])>0
    assert bound['exact_reset_on_nontrivial_sector_possible'] is False


def test_actual_resource_domains_not_silently_enlarged(local):
    stopped=r.source.execute(local,{(0,0,0,0):F(1)})
    assert stopped.status=='RESOURCE_STOP' and stopped.attempted_steps==3
    for level,n in ((0,239),(1,959),(2,3839)):
        value=r.program(level=level,total=0)
        assert len(value['commands'])==n<=4096
    with pytest.raises(ValueError):r.program(level=3,total=0)
    with pytest.raises(ValueError):r.program(level=0,total=True)
    with pytest.raises(ValueError):r.program(total=3)
    altered=deepcopy(local);altered['commands'][0]['law']='C'
    with pytest.raises(ValueError):r.specification(altered)
    for value in (True,1.0):
        altered=deepcopy(local);altered['commands'][0]['repeat']=value
        with pytest.raises(ValueError):r.specification(altered)
    with pytest.raises(ValueError):r.program(total=2.0)


def test_receipt_replay():assert packet.verify()


def test_receipt_bool_int_substitution_rejected(tmp_path,monkeypatch):
    original=packet.strict_json.load(packet.OUTPUT);bad=deepcopy(original)
    bad['scope']['all_binary_number_sectors_in_general_tree_proof']=1
    p=tmp_path/'receipt.json';p.write_text(packet.canonical(bad))
    monkeypatch.setattr(packet,'produce',lambda:original)
    with pytest.raises(ValueError):packet.verify(p)
