"""Controls for phase, collisions, complete outcomes and monitoring."""
from fractions import Fraction as F
import numpy as np
import pytest
from . import source as s

def test_actual_collision_sign_and_rank_two_contact():
    j=s.contact_identity();assert j['all_matrix_entries_exact']==4356 and j['contact_rank']==30
    g=s.carrier();e=g.seams.index((0,1));i=g.index[(0,1)]
    assert s.event_swap_column(g,i,2*e+1,3)==(i,2*e+1,(1,0))
    assert not j['free_exterior_identification']

def test_wrong_magnetic_square_term_is_not_a_probability():
    g=s.carrier();b=[s.ZERO]*g.states;b[g.index[(0,1)]]=(1,0);a=[s.ZERO]*g.states
    true=s.terminal_numerators(g,a,b,1);assert sum(true)==g.events and min(true)>=0
    wrong=s.apply_T_num(g,[(s.norm2(x),0)for x in b],1)
    assert any(x[1]!=0 for x in wrong)

def test_complete_event_inverse_and_conjugate_sectors():
    g=s.carrier()
    for k in range(6):
        for i in range(g.states):
            for e in range(g.events):
                j,ee,p=s.event_swap_column(g,i,e,k);ii,eee,q=s.event_swap_column(g,j,ee,k)
                assert(ii,eee)==(i,e)and s.mul(p,q)==(1,0)
    assert s.packet(g,1,0,(0,1,3,6))['records']==s.packet(g,5,0,(0,1,3,6))['records']
    b=[(int(i==0),0)for i in range(g.states)]
    assert s.apply_T_num(g,b,2,False)==s.apply_T_num(g,b,0)

def test_absorbing_protocol_is_complete_and_differs_from_terminal_tail():
    g=s.carrier();m=s.first_exit(g,1,horizon=8);assert all(r['complete_total_probability_exact']=='1'for r in m['records'])
    a=[(int(i==0),0)for i in range(g.states)];b=[s.ZERO]*g.states
    for _ in range(8):a,b=s.cyclic_step(g,a,b,1)
    nums=s.terminal_numerators(g,a,b,1);dist=s.distances(g)
    terminal_tail=F(sum(v for v,d in zip(nums,dist,strict=True)if d>1),g.events**17)
    assert m['records'][-1]['cumulative_exit_interval_exact']!=s.interval(terminal_tail)
    assert m['infinite_time_tail_assumed_zero']is False

def test_matching_code_failure_is_retained():
    j=s.matching_retention(1,(0,1));assert j['records'][0]['code_population_interval_exact']==['1','1']
    assert F(j['records'][1]['source_code_failure_interval_exact'][0])>0
    g=s.carrier();psi=np.zeros((g.states,g.events),complex);psi[g.index[(0,1)],0]=1
    assert np.sum(abs(s.matching_code_projection(g,psi))**2)==0

def test_enlargement_and_full_complement_not_removed():
    g=s.carrier(0);assert(g.ports,g.states,len(g.edges))==(240,28680,149940)
    for k in range(6):
        j=s.spectrum_census(g,k);assert j['fixed_plus_one_dimension']+j['fixed_minus_one_dimension']+j['dispersive_dimension']==36136800
    assert len(s.validate_loops(g))==2
    with pytest.raises(ValueError):s.carrier(1)

def test_compact_idle_and_edge_modes_have_exact_whole_carrier_eigenvalues():
    g=s.carrier()
    for k in range(6):
        v={(0,0):(1,0),(0,2):(-1,0)}
        assert s._compact_eigencheck(g,v,k,-1)==2
        i,j,c=g.seam_edges[0][0];p=s.phase(k,c)
        v={(i,1):(1,0),(j,1):p,(i,0):(-1,0),(j,0):s.scale(-1,p)}
        assert s._compact_eigencheck(g,v,k,-1)==4
    assert g.states*g.events-len(g.edges)-g.states==sum(g.events-d-1 for d in g.degree)+len(g.edges)
