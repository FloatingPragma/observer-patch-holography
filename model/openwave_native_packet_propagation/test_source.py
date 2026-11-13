from fractions import Fraction as F
from math import comb
import pytest
from model.openwave_native_process import source as native
from model.openwave_native_response.source import apply_population
from .source import (SCALE,power_interval,binomial_intervals,local_numerator_step,
                     distances,summarize_packet)

def test_outward_binomial_weights_include_complete_tail():
    for t in (0,1,7,19):
        weights,tail=binomial_intervals(t,1260,order=8)
        for k,(lo,hi) in enumerate(weights):
            value=comb(t,k)*F(6,1260)**k*F(1254,1260)**(t-k)
            assert F(lo,SCALE)<=value<=F(hi,SCALE)
        exact_tail=1-sum(comb(t,k)*F(6,1260)**k*F(1254,1260)**(t-k) for k in range(len(weights)))
        assert 0<=exact_tail<=F(tail,SCALE)
    with pytest.raises(ValueError):binomial_intervals(2,6)

def test_uniformization_equals_original_complete_population_columns():
    g=native.graph(0);n=g.ports
    original=[F(i==0) for i in range(n)]
    for _ in range(4):original=apply_population(g,original)
    v=[int(i==0) for i in range(n)];den=1;result=[F()]*n
    for k in range(5):
        weight=comb(4,k)*F(6,1260)**k*F(1254,1260)**(4-k)
        result=[a+weight*F(b,den) for a,b in zip(result,v,strict=True)]
        v=local_numerator_step(g,v);den*=6
    assert result==original and sum(result)==1

def test_first_arrival_is_a_retained_flag_and_coherence_decays():
    row=summarize_packet(0,2)
    for t in row['times']:
        first=list(map(F,t['first_current_record_arrival_distance_four_interval_exact']))
        terminal=list(map(F,t['terminal_distance_at_least_four_interval_exact']))
        assert first[1]>=terminal[0]
        assert F(t['omitted_binomial_tail_upper_exact'])<F(1,10**20)
    assert row['times'][-1]['intervals_decimal_for_display']['first_arrival'][0]>0.5
    assert row['times'][-1]['intervals_decimal_for_display']['coherent_fidelity'][1]<0.04

def test_power_interval_and_killed_block_do_not_renormalize():
    lo,hi=power_interval(7,9,23)
    assert F(lo,SCALE)<=F(7,9)**23<=F(hi,SCALE)
    g=native.graph(0);d=distances(g)
    boundary=next(i for i in range(g.ports) if d[i]==3 and any(d[j]==4 for j in g.neighbours[i]))
    vector=[int(i==boundary) for i in range(g.ports)]
    full=local_numerator_step(g,vector)
    nohit=local_numerator_step(g,vector,[x<4 for x in d])
    assert sum(full)==6 and sum(nohit)<6
    assert sum(nohit)!=sum(full)  # normalizing a successful branch is rejected.
