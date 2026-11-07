from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import json
import sympy as sp
import pytest
from . import geometry as g,propagation as p,records as r,packet


def test_exact_base_radius_and_crossing_energy():
    q=sp.sqrt(5);R2=(5+q)/2
    assert sp.simplify(4/R2-(2-2*q/5))==0
    assert sp.simplify(30*(6+2*q)/(9*R2)-(20+4*q)/3)==0
    assert g.edge_energy(0)==(F(20,3),F(4,3))
    x,den=g.coordinates(0)
    high=x[12*6];low=x[12*8]
    # These opposite original face means justify the all-level diameter proof.
    dz=(0,0,0,0,high[4]-low[4],high[5]-low[5])
    assert g.sign(g.add(g.normalized_square(dz,den),(F(-16,9),F())))>0


@pytest.mark.parametrize('ell',[0,1,2,3])
def test_literal_graph_all_halos_and_exact_energy(ell):
    graph=g.graph(ell);assert g.audit(graph)==g.census(ell)
    for cell in range(graph.cells):
        root,window=g.halo(graph,cell)
        assert len(root)==12 and len(window)==15
    assert g.edge_energy(ell)==g.energy_formula(ell)
    assert g.at_most(g.edge_energy(ell),F(12))


def test_routing_collision_and_mismatched_neighbours_rejected():
    graph=g.graph(0);edges=list(graph.seams)
    first=graph.inter[0];other=next(e for e in graph.inter if e[0]//12==first[0]//12 and e!=first)
    edges[edges.index(other)]=(first[0],other[1]);edges=tuple(sorted(edges))
    with pytest.raises(ValueError):g.audit(replace(graph,seams=edges))
    bad=list(graph.neighbours);bad[0]=bad[0][1:]
    with pytest.raises(ValueError,match='adjacency'):g.audit(replace(graph,neighbours=tuple(bad)))


def test_subdivision_counts_do_not_allocate_huge_source():
    c=g.census(200)
    assert c['ports']==240*4**200 and c['seams']==630*4**200
    assert g.at_most(g.energy_formula(200),F(12))
    with pytest.raises(ValueError,match='zero through three'):g.graph(200)


@pytest.mark.parametrize('ell',[0,1,2])
@pytest.mark.parametrize('time',[0,1,2,4])
def test_exact_native_moment_obeys_uniform_bound(ell,time):
    value=p.exact_moment(ell,time)
    assert g.sign(value)>=0
    assert g.at_most(value,p.displacement_bound(ell,time)['mean_squared_displacement_upper'])
    if time==1:
        census=g.census(ell)
        assert value==tuple(x/F(census['ports']*census['seams']) for x in g.edge_energy(ell))


def test_exact_moment_zero_and_base_value():
    assert p.exact_moment(0,0)==(0,0)
    assert p.exact_moment(0,1)==(F(1,22680),F(1,113400))


def test_ballistic_and_incident_bounds_have_distinct_clock_scalings():
    for ell in (0,1,2,20,200):
        result=p.ballistic_bound(ell)
        assert result['mean_squared_displacement_upper']==F(1,10*2**ell)
        assert p.incident_bound(ell,2,30)['mean_squared_displacement_upper']==F(1,75*4**ell)
        assert p.incident_bound(ell,0,30)['mean_squared_displacement_upper']==0
        assert p.clock_counts(ell,630*4**ell,5)['receiver_incident_mean']==5


def test_exact_stopped_kernel_and_invalid_mean_time_shortcut():
    # A tiny endpoint-law algebra control, not an alternative source geometry.
    I=sp.eye(3);swap01=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
    swap12=sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
    U=(I+swap01)/2;V=(I+swap12)/2;a=sp.Rational(1,2)
    C=a*(I-(1-a)*U).inv();K=V*C;P=(U+V)/2
    f=sp.Matrix([0,1,4]);epsilon=lambda Q:sp.Rational(2,3)*(f.T*(I-Q)*f)[0]
    assert K*sp.ones(3,1)==sp.ones(3,1)
    assert sp.ones(1,3)*K==sp.ones(1,3)
    assert epsilon(K)>epsilon(P)/a
    assert epsilon(K)<=2*epsilon(P)/a
    for k in (0,1,2,3,7):assert epsilon(K**k)<=2*k*k*epsilon(P)/a


def test_timeout_retains_complete_mass():
    a=F(30,630)
    assert p.timeout(0,0,30,5)==0
    assert p.timeout(0,1,30,5)==(1-a)**5
    assert p.timeout(0,4,30,3)==1
    assert p.timeout(0,2,30,3)==(1-a)**3+3*a*(1-a)**2


@pytest.fixture(scope='module')
def complete_trace():
    value=r.program(0,1);family=r.cone.source_family(value)
    coins=(1,)+(0,)*len(r.detector.tree(family))
    return value,r.trace(value,0,((0,1),),coins)


def test_body_readback_uses_records_and_preserves_full_terminal_state(complete_trace):
    value,result=complete_trace
    assert r.authenticate(value,result)
    decoded=r.decode(value,result['records'])
    assert decoded['body_masks']==(1,1)
    assert len(result['archive'])==len(result['records'])-1 and len(result['memory'])==1
    query=r.query(value,result['slots'],{(1,0,1,0):F(1)})
    assert r.same(query['density'],result['density'])
    assert r.source.trace(query['density'])==F(1,1260)*F(1,2)**239


@pytest.mark.parametrize('kind',['decoded_bool','weight_float','extra','record_coin_bool','record_seam_bool','density_float','archive_list','missing_record'])
def test_trace_authentication_rejects_type_aliases_and_incomplete_outputs(complete_trace,kind):
    value,original=complete_trace;result=deepcopy(original)
    if kind=='decoded_bool':result['decoded']['body_masks']=(True,1)
    elif kind=='weight_float':result['weight']=float(result['weight'])
    elif kind=='extra':result['extra']=False
    elif kind=='record_coin_bool':result['records']=(replace(result['records'][0],coin=True),)+result['records'][1:]
    elif kind=='record_seam_bool':result['records']=(replace(result['records'][0],seam=(False,1)),)+result['records'][1:]
    elif kind=='density_float':result['density']={key:float(v) for key,v in result['density'].items()}
    elif kind=='archive_list':result['archive']=list(result['archive'])
    elif kind=='missing_record':result['records']=result['records'][:-1]
    with pytest.raises(ValueError):r.authenticate(value,result)


def test_body_program_probability_and_trace_edge_types_are_bound(complete_trace):
    value,result=complete_trace;altered=deepcopy(value)
    altered['parameters']['repair_probability']='1/2'
    altered['declared_resources']=r.source._native_cost(altered)
    with pytest.raises(ValueError):r.authenticate(altered,result)
    family=r.cone.source_family(value)
    with pytest.raises(ValueError):r.trace(value,0,((False,1),),(1,)+(0,)*len(r.detector.tree(family)))


def test_zero_raw_word_is_not_normalized_into_success(complete_trace):
    value,result=complete_trace;slots=(0,)+result['slots'][1:]
    answer=r.query(value,slots,{(1,0,1,0):F(1)})
    assert not answer['density'] and len(answer['records'])==len(slots)
    assert answer['records'][0]==r.source.Record('hold')


@pytest.mark.parametrize('ell',[0,1,2,3])
def test_receiver_first_tick_positive_far_path(ell):
    w=r.arrival_witness(ell)
    assert w['graph_distance']>2**ell and w['receiver_ticks']==1
    assert w['specified_complete_body_word_probability']==F(1,2*len(g.graph(ell).seams))**len(w['edges'])
    assert all(w['receiver'] not in edge for edge in w['edges'][:-1])


def test_complete_prefix_and_crossing_record_attachment():
    value=r.certificate()
    assert value['actual_body_plus_detector']['expanded_status']=='RESOURCE_STOP'
    assert value['actual_body_plus_detector']['expanded_completed_attempts']==1
    assert value['actual_body_plus_detector']['actual_expanded_prefix_matches_one_unit_transition']
    assert value['actual_crossing_body_plus_detector']['body_path_decoded_from_complete_records']


@pytest.mark.parametrize('call',[
    lambda:g.census(True),lambda:g.census(-1),lambda:p.displacement_bound(0,True),
    lambda:p.displacement_bound(0,1,F()),lambda:p.ballistic_bound(0,-F(1)),
    lambda:p.incident_bound(0,1,0),lambda:p.incident_bound(0,1,631),
    lambda:p.timeout(0,1,5,257),lambda:p.exact_moment(0,9),lambda:r.program(3,1),
    lambda:r.program(0,10**100),
])
def test_domains_and_actual_host_geometry_limits(call):
    with pytest.raises(ValueError):call()


def test_retained_receipt_matches_replay():assert packet.verify()


def test_receipt_typed_alias_and_omission_fail_closed(tmp_path,monkeypatch):
    rebuilt=packet.produce();monkeypatch.setattr(packet,'produce',lambda:rebuilt)
    path=tmp_path/'receipt.json'
    bad=deepcopy(rebuilt);bad['scope']['host_limits_changed']=0
    path.write_text(json.dumps(bad))
    with pytest.raises(ValueError):packet.verify(path)
    bad=deepcopy(rebuilt);bad['geometry']['formal_level200_census']['level']=200.0
    path.write_text(json.dumps(bad))
    with pytest.raises(ValueError):packet.verify(path)
    bad=deepcopy(rebuilt);del bad['source_bindings'][next(iter(bad['source_bindings']))]
    path.write_text(json.dumps(bad))
    with pytest.raises(ValueError):packet.verify(path)
