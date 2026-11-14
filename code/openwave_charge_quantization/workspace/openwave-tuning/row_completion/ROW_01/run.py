"""Bind native topology to the complete observer-only row-1 apparatus.

Uses standard-library exact arithmetic and writes only this row directory.
Existing full probability laws are consumed, not represented as fresh runs.
"""
from collections import Counter,defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import topology

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PROVIDER=ROOT/'reverse-engineering-reality/code/openwave_row1_integrity'
AUDIT=ROOT/'openwave-tuning/next/CHAIN3/row1_integrity_audit'
sys.path.insert(0,str(PROVIDER));sys.path.insert(0,str(AUDIT))
import integrity as it
from observer_only import admit_record_only

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def save(name,value):(HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def independent_event(n,edge,c):
    a,b=edge;S=n[a]+n[b];na=S//2+c*(S%2);j=n[a]-na
    out=list(n);out[a]=na;out[b]=S-na
    return tuple(out),j
def states(total):
    return [tuple(sum(v==x for x in group) for v in range(12))
            for group in combinations_with_replacement(range(12),total)]
def rejected(fn):
    try:fn()
    except (ValueError,AssertionError,IndexError):return True
    return False
def pins():
    output={}
    for directory,name,wanted in (
       (PROVIDER,'MANIFEST.json','1f8c91d4820fa2fe5f439f0a8e344645b4c29b11e3541ac32b7c52a71f1d73c4'),
       (ROOT/'reverse-engineering-reality/code/openwave_row1_native','manifest.json','d5e692b0597dbbeafabb88e6f7d417da8f7d7d1fd39a334c08fa73ffe40e8d90')):
        m=directory/name;assert sha(m)==wanted
        output[str(m.relative_to(ROOT))]=wanted
        data=json.loads(m.read_text())
        for relative,h in data['files_sha256'].items():
            p=directory/relative;assert sha(p)==h
            output[str(p.relative_to(ROOT))]=h
        if 'result_file_sha256' in data:
            p=directory/'RESULTS.json';assert sha(p)==data['result_file_sha256'];output[str(p.relative_to(ROOT))]=sha(p)
    for p in [AUDIT/'observer_only.py',AUDIT/'NATIVE_POINTER_RECORDS.json',AUDIT/'EXACT_TREE_COUNTS.json',
              HERE/'run.py',HERE/'topology.py',HERE/'CONTRACT.json',*sorted((HERE/'upstream').iterdir())]:
        if p.is_file():output[str(p.relative_to(ROOT))]=sha(p)
    return output

def main():
    before=pins();geometry=topology.construct(it.EDGES);save('GEOMETRY.json',geometry)
    closure=0
    for total in (1,2):
        for n in states(total):
            for edge in it.EDGES:
                for c in (0,1):
                    out,j=independent_event(n,edge,c);actual,(_,actualj)=it.native(n,edge,c)
                    assert (out,j)==(actual,actualj) and min(out)>=0 and sum(out)==total
                    if total==1:assert sorted(out)==[0]*11+[1]
                    closure+=1
    # All admissible preparations/reconstructed configurations, whether reachable
    # or not, obey the same finite-counter injection and regional range.
    pairs=0
    for total in (1,2):
        for G in states(total):
            for n in states(total):
                q=[g-v for g,v in zip(G,n)]
                assert sum(q)==0 and sum(v for v in q if v>0)<=total and sum(v for v in q if v<0)>=-total
                for w in range(-total,total+1):assert it.lift(w%6,total)==w
                pairs+=1
    print('Exact geometry and all native closure branches passed',closure,flush=True)
    examples=[];controls={};seen=set();max_phase_error=0.;latch_changes=0
    starts=[states(1)[i] for i in (0,3,8)]+[tuple(int(v==0)+int(v==b) for v in range(12)) for b in (0,1,5)]
    order=list(range(7,30))+list(range(7))
    fixture=json.loads((AUDIT/'NATIVE_POINTER_RECORDS.json').read_text())
    fixture_by=defaultdict(list)
    for item in fixture['records']:fixture_by[item['command']].append(item)
    for index,G in enumerate(starts):
        for pattern in (0,1):
            prefix=[(0,pattern,'C'),(7,1-pattern,'E')]
            coins=[pattern if e%2==0 else 1-pattern for e in range(30)]
            received,audit=it.source_execute(G,prefix,coins,order)
            # Independent source replay, including the disturbing closing sweep.
            n=G;current=[0]*30
            for records in (received['source_records'],received['sweep_records']):
                for record in records:
                    e=record['seam'];out,j=independent_event(n,it.EDGES[e],record['coin'])
                    assert j==record['current']
                    if record['kind']=='E':assert record['contrast']==n[it.EDGES[e][0]]-n[it.EDGES[e][1]]
                    n=out;current[e]+=j
            assert list(n)==audit['source_L1']
            q=topology.divergence(it.EDGES,current);assert q==tuple(g-v for g,v in zip(G,n))
            for region in ((0,),tuple(range(1,12)),(0,1,2),tuple(range(12))):
                w=topology.discrete_curve_degree(geometry,current,region)
                assert w==sum(q[v] for v in region)
                for deformed in (False,True):
                    err=abs(topology.phase_degree(geometry,current,region,deformed)-w)
                    max_phase_error=max(max_phase_error,err);assert err<1e-12
                seen.add(w)
                latches=it.observe_source(G,received,sum(G));latch_changes+=sum(latches['L0']['counters'][v] for v in region)%6 != w%6
                # One successful replay and each available failure stage. They
                # are conditional fixtures; their fine weights are never called
                # the complete law. Full laws are consumed separately below.
                selected={}
                for item in fixture_by[w]:
                    reason=item['native_executor_result']['reason']
                    if reason not in selected:selected[reason]=item
                for reason,item in selected.items():
                    pointer=deepcopy(item['native_executor_result'])
                    x=item['initial_column_for_audit_only']
                    actual=it.meter.execute(x,(w,),[r['coin'] for r in pointer['records']],2)
                    assert actual==pointer
                    claim=it.assemble_claim(G,received,region,pointer,2)
                    del claim['meter']['final_loads']
                    verdict=admit_record_only(claim,G,region,2,producer=it)
                    assert verdict['admitted'],verdict
                    decoded=topology.raw_pointer_degree(pointer,it.meter.CYCLE)
                    if pointer['success']:
                        assert decoded==w and verdict['readout']==-w
                    else:assert decoded is None and verdict['readout'] is None
                    row={'G':list(G),'source_pattern':pattern,'region':list(region),'actual_source_charge':-w,
                         'topological_degree':w,'claim':claim,'admission':verdict,
                         'raw_pointer_record_degree':decoded,'source_audit':audit}
                    examples.append(row)
                    if pointer['success'] and w and 'wrong_target' not in controls:
                        bad=deepcopy(claim);bad['decoded_relative_charge']=0
                        controls['wrong_target']=not admit_record_only(bad,G,region,2,producer=it)['admitted']
                        bad=deepcopy(claim);bad['meter']['records']=[]
                        controls['empty_pointer']=not admit_record_only(bad,G,region,2,producer=it)['admitted']
                        bad=deepcopy(claim);bad['sweep_records'][0]['contrast']+=1
                        controls['false_closing_contrast']=not admit_record_only(bad,G,region,2,producer=it)['admitted']
                        bad=deepcopy(claim);bad['meter']['final_loads']=['malicious']
                        controls['hidden_snapshot_ignored']=admit_record_only(bad,G,region,2,producer=it)==verdict
                        controls['double_winding_encoding_detected']=topology.discrete_curve_degree(geometry,[2*j for j in current],region)!=decoded
                        controls['constant_zero_meter_detected']=decoded!=0
    assert {-2,-1,0,1,2}<=seen,seen
    print('Signed source/topology/pointer examples passed',len(examples),flush=True)
    # A genuine half-current control still conserves its real-valued charge,
    # but its dual edge is not a circle map with common endpoint values.
    half=[0]*30;half[0]=F(1,2)
    G=(1,)+(0,)*11;n=[F(x) for x in G];a,b=it.EDGES[0];n[a]-=half[0];n[b]+=half[0]
    assert topology.divergence(it.EDGES,half)==tuple(g-v for g,v in zip(G,n))
    controls['half_current_gluing']=rejected(lambda:topology.discrete_curve_degree(geometry,half,(a,)))
    controls['phase_alias_outside_declared_cap']=rejected(lambda:topology.phase_degree(geometry,[64]+[0]*29,(a,)))
    controls['missing_spatial_edge']=rejected(lambda:topology.construct(it.EDGES[:-1]))
    controls['mod6_unrestricted_alias_exposed']=6%6==0 and topology.divergence(it.EDGES,[6]+[0]*29)!=tuple([0]*12)
    # Actual double-load split: no sign-only interpretation survives.
    double=tuple(2*int(v==a) for v in range(12));split,j=independent_event(double,(a,b),0)
    assert j==1 and split[a]==split[b]==1
    controls['double_load_can_split']=it.native(double,(a,b),0)[0]==split
    controls['sign_only_readout_detected']=any(abs(r['topological_degree'])==2 and abs(r['raw_pointer_record_degree'] or 0)==2 for r in examples)
    assert all(controls.values()),controls
    frozen=json.loads((PROVIDER/'RESULTS.json').read_text())
    complete=frozen['complete_bounded_laws']
    assert complete['instrument_count']==768 and complete['total_mass_error']=='0' and complete['failures_retained_without_normalization']
    for row in complete['rows']:
        assert sum(F(x) for x in row['outcome_mass'].values())==1
    assert pins()==before
    save('EXAMPLES.json',examples)
    result={'schema':'oph.row01.native_topological_charge.candidate.v1',
      'candidate_id':'literal-native-relative-spatial-charge-AU5-record-only-v1',
      'criterion':'integer topological charge; a split e → e/2 + e/2 impossible',
      'native_branches':closure,'all_admissible_G_n_pairs':pairs,'signed_degrees_observed':sorted(seen),
      'complete_join_examples':len(examples),'different_pre_and_post_sweep_regional_classes':latch_changes,
      'maximum_deformed_phase_readout_error':max_phase_error,'controls':controls,
      'complete_laws_consumed':{'instruments':768,'source_coin_words_per_setting':2**32,
          'all_twelve_unknown_pointer_columns':True,'failure_outcomes_retained':True,
          'result_sha256':sha(PROVIDER/'RESULTS.json'),'fresh_full_law_rerun':False},
      'checks':{'exact_integer_spatial_isomorphism':True,'native_closure_and_nonfragmentation':True,
          'signed_primitive_degrees_and_pair_zero':True,'actual_continuous_deformed_field_readout':True,
          'actual_native_pointer_dictionary':True,'observer_only_admission':True,
          'complete_probability_law_and_failures':True,'all_sources_unchanged':True},
      'source_sha256':before,'full_criterion_candidate':True,'full_criterion_met':False,
      'pending':'Frozen observer-only opposing audit custody and two independent full-criterion reviews.',
      'excluded':'Physical units, species mass/spin/force, coherent source transfer and unrestricted mod6 charges. These are not added row1 gates.',
      'official_rows_changed':[]}
    save('RESULT.json',result)
    print(json.dumps({k:result[k] for k in ('native_branches','all_admissible_G_n_pairs','signed_degrees_observed','complete_join_examples','maximum_deformed_phase_readout_error')}),flush=True)

if __name__=='__main__':main()
