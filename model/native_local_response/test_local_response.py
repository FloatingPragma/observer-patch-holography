from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations,product
import json
import pytest
from . import cone as c,bounds as b,packet


def test_actual_source_endpoint_attachment():
    report=packet.source_attachment()
    assert report['weighted_H66_fixed_raw_maps_checked']==750
    assert report['fair_federation_complete_raw_maps']==25


@pytest.mark.parametrize('federation',[False,True])
def test_complete_native_history_cone_comparison(federation):
    report=packet.actual_comparison(federation)
    assert report['actual_expanded_status']=='COMPLETE'
    assert report['source_raw_histories']==15625
    assert report['observed_raw_histories']==625
    assert report['observed_zero_blocks']==593
    assert report['exact_initial_cone']==[0,1,5]


@pytest.mark.parametrize('law',['C','E','mixture'])
def test_all_fixed_two_matrix_units_with_reference_on_disjoint_fixture(law):
    # An algebra control, separately attached to the actual source above.
    family=c.native.Family((25,8,3,7,1,9),tuple((i,i+1) for i in range(5)),(1,)*5)
    schedule=((0,1),(4,5),(0,1));final=(0,1,2,3);observed=(2,)
    masks=tuple(sum(1<<v for v in pair) for pair in combinations(range(6),2))
    for i,j in product(masks,repeat=2):
        rho={(i,0,j,1):F(-1,7)}
        full=c.instrument(family,rho,schedule,final,observed,law=law)
        reduced=c.instrument(family,rho,schedule,final,observed,reduce=True,law=law)
        assert full==reduced
    assert c.backward(schedule,final,observed)['retained_commands']==(0,2)


def test_nonzero_quantum_coherence_survives_exact_cone_readout():
    family=c.native.local_family();schedule=((0,1),(2,3),(0,1))
    final=(0,1,5,7);rho={(1<<5,0,1<<7,1):F(1,2)}
    # These two local masks can occur with the same remote spectator in H66.
    rho={(i|(1<<4),a,j|(1<<4),b):v for (i,a,j,b),v in rho.items()}
    actual=c.instrument(family,rho,schedule,final,(2,),reduce=True,law='E')
    assert sum((sum(block.values(),F()) for block in actual.values()),F())==F(1,2)
    assert any(any(i!=j and a!=b for i,a,j,b in block) for block in actual.values())


def test_observed_remote_endpoints_are_mandatory_roots():
    schedule=((0,1),(2,3))
    unobserved=c.backward(schedule,(0,1),())
    observed=c.backward(schedule,(0,1),(1,))
    assert unobserved['initial_support']==(0,1)
    assert observed['initial_support']==(0,1,2,3)
    assert observed['observation_roots']==(0,1,2,3)


def test_backward_order_and_propagation_are_not_an_undirected_union():
    assert c.backward(((0,1),(1,2)),(0,),())['initial_support']==(0,1)
    assert c.backward(((1,2),(0,1)),(0,),())['initial_support']==(0,1,2)


def test_same_reduced_cone_reference_inputs_give_identical_output():
    family=c.native.local_family();schedule=((0,1),(2,3),(0,1));final=(0,1)
    first={(5,0,5,0):F(1,2),(6,1,6,1):F(1,2)}
    second={(9,0,9,0):F(1,2),(10,1,10,1):F(1,2)}
    assert c.partial_trace(first,final,12)==c.partial_trace(second,final,12)
    assert c.instrument(family,first,schedule,final,(2,))==c.instrument(family,second,schedule,final,(2,))


def test_unreferenced_equal_marginals_do_not_bound_correlated_reference():
    family=c.native.local_family()
    first={(5,0,5,0):F(1,2),(6,1,6,1):F(1,2)}
    second={(5,1,5,1):F(1,2),(6,0,6,0):F(1,2)}
    a=c.instrument(family,first,(),(0,),())[()]
    z=c.instrument(family,second,(),(0,),())[()]
    assert sum(a.values())==sum(z.values())==1
    assert a!=z
    assert sum(abs(a.get(key,0)-z.get(key,0)) for key in a.keys()|z.keys())/2==1


def test_actual_controller_and_memory_halt_boundaries():
    result=packet.boundaries()
    assert result['actual_adaptive_receiver_probabilities']==['1/4','0']
    assert result['four_step_expanded_status']=='RESOURCE_STOP'
    assert result['four_step_expanded_attempts']==3
    with pytest.raises(ValueError,match='record-adaptive'):
        c.audit_program(packet.dataflow.feedback_program())
    halted=packet.fixed_program(((0,1),(2,3)),overflow='halt')
    with pytest.raises(ValueError,match='FIFO'):c.audit_program(halted)
    halted['memory']['capacity']=2;halted['declared_resources']=c.source._native_cost(halted)
    assert c.audit_program(halted)['attempts']==2


def test_unknown_record_and_zero_address_cannot_be_dropped(monkeypatch):
    original=c.fixed_instrument
    monkeypatch.setattr(c,'fixed_instrument',lambda *args,**kwargs:original(*args,**kwargs)[1:])
    with pytest.raises(ArithmeticError,match='address'):packet.source_attachment()


def test_history_must_bind_exact_source_probability_and_fifo_partition():
    value=packet.fixed_program(((0,1),(2,3)))
    history=c.source.execute(value,{(0,0,0,0):F(1)})
    altered=deepcopy(value);altered['parameters']['repair_probability']='1/2'
    altered['declared_resources']=c.source._native_cost(altered)
    with pytest.raises(ValueError,match='exact admitted source program'):
        c.observe_history(altered,history,(0,1),(0,))
    with pytest.raises(ValueError,match='every admitted command'):
        c.observe_history(value,replace(history,attempted_steps=1),(0,1),(0,))
    first=history.branches[0]
    corrupt=replace(first,archive=(),memory=first.archive+first.memory)
    with pytest.raises(ValueError,match='partition'):
        c.observe_history(value,replace(history,branches=(corrupt,)+history.branches[1:]),(0,1),(0,))


def test_admitted_weighted_federation_is_not_invented():
    value=packet.fixed_program(((0,1),),federation=True)
    value['parameters']['port_weights']=[25,25]+[8]*10
    value['declared_resources']=c.source._native_cost(value)
    with pytest.raises(ValueError,match='bias lift'):c.source_family(value)


@pytest.mark.parametrize('final,observed',[
    ((0,1),(-1,)),((0,1),(0,0)),((0,1),(2,)),((0,1),(True,)),
    ((-1,),(0,)),((0,0),(0,)),((12,),(0,)),((True,),(0,)),
])
def test_history_observation_indices_are_canonical_before_address_lookup(final,observed):
    value=packet.fixed_program(((0,1),(2,3)))
    history=c.source.execute(value,{(0,0,0,0):F(1)})
    with pytest.raises(ValueError,match='canonical source vertex subset'):
        c.observe_history(value,history,final,observed)


def test_all_literal_cell_halos_and_census():
    count=0
    for level in (0,1,2):
        graph=c.native.geometry.federation(level)
        for cell in range(graph.cells):
            row=b.halo_geometry(level,cell)
            assert len(row['window'])==15
            assert row['ordered_two_edge_patterns']==15
            assert row['seams']==630*4**level
            count+=1
    assert count==420


def test_exact_subset_cone_against_every_small_selector_word():
    family=c.native.Family((25,8,3,7),((0,1),(1,2),(2,3)),(1,2,3))
    root=(0,);window=(0,1)
    for attempts in range(5):
        probability=F()
        for tape in product(family.seams,repeat=attempts):
            cone=c.backward(tape,root,())['initial_support']
            if not set(cone)<=set(window):
                weight=F(1)
                for edge in tape:weight*=family.selection(edge)
                probability+=weight
        exact=b.escape(family,root,window,attempts)['escape_probability']
        assert probability==exact
        assert exact<=b.walk_bound(family,root,window,attempts)['error_upper']
        assert exact<=b.radius_bound(family,root,1,attempts)


def test_four_state_lumping_from_all_eight_actual_halo_subsets():
    graph=c.native.geometry.federation(0);halo=b.halo_geometry(0,0)
    root=set(halo['root']);window=set(halo['window']);M=len(graph.seams)
    matrix=b.halo_transition(M)
    for bits in range(8):
        selected={v for i,v in enumerate(halo['halo']) if bits>>i&1}
        present=root|selected;n=len(selected);counts=[0]*5
        for edge in graph.seams:
            target=present|set(edge) if present.intersection(edge) else present
            counts[4 if not target<=window else len(target-root)]+=1
        assert tuple(F(x,M) for x in counts[:4])==tuple(matrix[i][n] for i in range(4))
        assert counts[4]==5*n


@pytest.mark.parametrize('attempts',[0,1,2,3,6,12,30])
def test_global_exact_halo_and_ordered_pair_bound(attempts):
    family=c.native.federation_family(0,(1,)*12);halo=b.halo_geometry(0,0)
    subset=b.escape(family,halo['root'],halo['window'],attempts)
    four=b.halo_global_error(630,attempts)
    assert subset['escape_probability']==four['cone_escape_probability']
    assert four['cone_escape_probability']<=four['ordered_pair_upper']
    assert subset['maximum_subset_states']<=8


def test_halo_first_escape_exact_coefficient_and_memory_control():
    assert b.halo_global_error(630,0)['cone_escape_probability']==0
    assert b.halo_global_error(630,1)['cone_escape_probability']==0
    assert b.halo_global_error(630,2)['cone_escape_probability']==F(15,630**2)
    assert b.halo_internal_error(1)['cone_escape_probability']==F(677,50061)
    # Omitting the initial failure gap would give R**0 and the false zero.
    assert b.halo_internal_error(0)['cone_escape_probability']!=b.halo_internal_error(1)['cone_escape_probability']


@pytest.mark.parametrize('M',[630,2520,10080])
def test_geometric_stopping_resolvent_against_full_global_attempt_equations(M):
    R=b.internal_block();a=F(30,M)
    for i,j in product(range(4),repeat=2):
        # Independently form Q_failure from its physical event counts.
        lhs=R[i][j]*(1-(1-F(33+4*i,M)))
        if i:lhs-=F(4-i,M)*R[i-1][j]
        assert lhs==a*int(i==j)
    assert 1-sum(row[0] for row in R)==F(677,50061)


@pytest.mark.parametrize('k',[0,1,2,4,8,16,64])
def test_local_clock_ordered_pair_envelope(k):
    result=b.halo_internal_error(k)
    assert 0<=result['cone_escape_probability']<=min(F(1),F(k*(k+1),120))
    assert result['global_elapsed_attempt_record_retained']


@pytest.mark.parametrize('q',[F(0),F(1,4),F(1,2),F(3,4),F(99,100)])
def test_abel_response_resolvent_and_series_tail(q):
    result=b.abel_response_bound(q);R=b.internal_block();v=result['resolvent_vector']
    assert tuple(v[i]-q*b.matrix_apply(R,v)[i] for i in range(4))==(F(1),F(),F(),F())
    partial=2*(1-q)*sum(q**k*b.halo_internal_error(k)['cone_escape_probability'] for k in range(32))
    assert partial<=result['bounded_score_response_error_upper']<=partial+2*q**32
    assert result['bounded_score_response_error_upper']<=result['coarse_error_upper']


def test_generic_local_clock_timeout_is_explicit():
    family=c.native.Family((1,1,1),((0,1),(1,2)),(1,1))
    value=b.fixed_timeout_bound(family,(0,1),(0,1,2),2,3)
    assert value['timeout_probability']==F(1,2)
    assert value['cone_error_at_timeout']==0
    assert value['error_upper']==F(1,2)


def test_evaluator_and_probability_domain_guards():
    family=c.native.local_family()
    with pytest.raises(ValueError):b.escape(family,(0,),(0,),257)
    with pytest.raises(ValueError):b.escape(family,(0,),(0,1),3,max_states=1)
    with pytest.raises(ValueError):b.walk_bound(family,(0,),(0,1),65)
    with pytest.raises(ValueError):b.abel_response_bound(F(1))
    with pytest.raises(ValueError):b.abel_response_bound(.5)
    with pytest.raises(ValueError):b.halo_internal_error(257)
    with pytest.raises(ValueError):c.backward(((0,1),),(0,),(1,))


def test_source_bindings_and_hostile_json(tmp_path):
    bindings=packet.source_bindings()
    assert 'model/native_program/interpreter.py' in bindings
    assert 'model/observer/source_discrimination/circuit.py' in bindings
    assert 'model/observer/source_scalar/compiler.py' in bindings
    assert not any(path.startswith('docs/') for path in bindings)
    path=tmp_path/'duplicate.json';path.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError):packet.strict_json.load(path)


def test_retained_receipt_and_semantic_mutation(tmp_path):
    packet.verify()
    value=packet.strict_json.load(packet.OUTPUT)
    value['scope']['physical_spatial_dimension_or_Lorentz_clock_derived']=True
    path=tmp_path/'bad.json';path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='receipt differs'):packet.verify(path)


@pytest.mark.parametrize('mutation',['bool_as_int','int_as_float'])
def test_receipt_exact_json_types_cannot_alias(mutation,tmp_path):
    value=packet.strict_json.load(packet.OUTPUT)
    if mutation=='bool_as_int':
        value['scope']['physical_spatial_dimension_or_Lorentz_clock_derived']=0
    else:
        value['source_attachment']['fair_federation_endpoint_patterns']=4.0
    path=tmp_path/'typed.json';path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='receipt differs'):packet.verify(path)
