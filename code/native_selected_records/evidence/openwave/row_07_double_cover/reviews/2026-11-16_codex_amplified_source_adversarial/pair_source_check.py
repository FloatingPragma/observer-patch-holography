"""Independent exact classical/joint-source audit; no producer imports.

The inspected source execution is a disclosed unregistered pilot followed by
registered replay. This review does not relabel it prospective discovery.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import random
import time

HERE=Path(__file__).resolve().parent
PACKET=HERE.parents[2]/'row_02_electron_rest_energy/work/2026-11-16_codex_amplified_common_pair_v1'
EDGES=((0,1),(0,1),(0,2),(0,2),(1,2))
STREAMS=((72365001,72366001,72365002,72366002),(72365011,72366011,72365012,72366012))
PATTERNS=(([2,-1,2,-1]+[0]*8,[2,2,-1]+[0]*9),([1,-2,2,2,-1]+[0]*7,[2,-2,1,-1,2]+[0]*7))


def need(ok,message):
    if not ok:raise AssertionError(message)


def endpoint(q,edge,c):
    a,b=edge;d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2
    out=q.copy();out[a]-=j;out[b]+=j
    return out,[a,b,c,j]


def census(records,n):
    relative=[None]*n;relative[0]=0;changes=[0]*n;index=0
    while index<len(records):
        a,b,c,j=records[index];first=(c,j);correction=changes[a]-changes[b]
        need(relative[a]is not None and relative[b]is None,'Census tree ordering')
        changes[a]-=j;changes[b]+=j;index+=1
        while True:
            need(index<len(records),'Truncated current census')
            aa,bb,cc,jj=records[index];index+=1
            need((aa,bb)==(a,b),'Census address changed')
            changes[a]-=jj;changes[b]+=jj
            if cc!=c:break
        relative[b]=relative[a]-(2*first[1]+jj-correction)
    need(all(x is not None for x in relative),'Incomplete census')
    return [x-relative[12]for x in relative]


def calculate(path):
    raw=path.read_bytes();data=json.loads(gzip.decompress(raw));g=data['graph']
    need(g['N']==223 and g['core_N']==24 and g['target']==[0,24,25] and g['blank']==[29,30],'Wrong common graph interface')
    allowed={tuple(sorted(e))for e in g['edges']};target=g['target'];total_calls=0;columns=0
    domain=Counter();joint_columns=0;coin_stream_checks=0
    for row in data['rows']:
        initials=[];signs=[];tapes=[];tables=[];backgrounds=[]
        for side in (0,1):
            observer=row['observer_records'][side];audit=row['source_audits'][side]
            initial=census(observer['admission']['probe'],24);sign=1 if sum(initial)>0 else-1
            need(initial==observer['admission']['initial_core'],'Current-derived initial state differs')
            q=row['initial_states_audit_only'][side].copy()
            need(q[:24]==initial and len(q)==223,'Initial full-state audit differs from currents')
            before_calibration=None;actual_calibration=[];actual_blank=[];coins=[[],[]]
            for stage,event in audit['events']:
                a,b,c,j=event;need(tuple(sorted((a,b)))in allowed,'Non-graph endpoint event')
                need(type(c)is int and c in(0,1),'Invalid coin')
                iscal=stage.startswith('calibration_');coins[int(iscal)].append(c)
                if iscal and before_calibration is None:before_calibration=q.copy()
                after,expected=endpoint(q,[a,b],c);need(event==expected,'Continuing native current differs')
                if stage in('calibration_route','calibration_probe'):
                    before=[q[v]for v in target];now=[after[v]for v in target]
                    need(before.count(sign)==now.count(sign)==1 and before.count(0)==now.count(0)==2,'Calibration lost isolated signed unit')
                    actual_calibration.append((stage,event,before.index(sign),now.index(sign)))
                    spectator=after.copy();spectator[target[now.index(sign)]]-=sign;spectator[target[1]]+=sign
                    need(spectator==before_calibration,'Calibration changed spectator background')
                elif stage=='calibration_blank':
                    need(after==q and j==0,'Blank changed continuing source');actual_blank.append(event)
                q=after;total_calls+=1
            need(q==audit['final']==observer['admission']['observed_final']==before_calibration,'Calibration did not return the same source')
            need([q[v]for v in target]==[0,sign,0] and observer['outcome']=='ready','Source not ready')
            need(sum(audit['literal_endpoint_calls'].values())==len(audit['events']),'Wrong continuing call count')
            filed=observer['calibration']['events']
            need(len(filed)==len(actual_calibration),'Missing actual calibration event')
            observed={}
            for record,(stage,event,before,after)in zip(filed,actual_calibration):
                need((record['stage'],record['record'],record['before_position'],record['after_position'])==(stage,event,before,after),'Calibration event reader differs')
                if stage=='calibration_probe':
                    edge=tuple(target.index(v)for v in event[:2]);key=(before,edge,event[2]);value=(event[3],after)
                    need(key not in observed or observed[key]==value,'Inconsistent observed column');observed[key]=value
            expected_keys={(p,e,c)for p in range(3)for e in((0,1),(1,2),(0,2))for c in(0,1)}
            need(set(observed)==expected_keys,'Missing actual conditional column')
            filed_columns=observer['calibration']['columns']
            table={(r['input_position'],tuple(r['edge']),r['coin']):(r['current'],r['decoded_output'])for r in filed_columns}
            need(len(filed_columns)==18 and table==observed,'Filed table differs from actual probes')
            need(observer['blank']['records']==actual_blank and {r[2]for r in actual_blank}=={0,1},'Blank coupon incomplete')
            initials.append(initial[:12]);signs.append(sign);tapes.append(coins);tables.append(table);backgrounds.append(q)
            columns+=18
        pattern=[i for i,pair in enumerate(PATTERNS)if all(initials[s]==[signs[s]*x for x in pair[s]]for s in(0,1))]
        stream=[]
        for i,seeds in enumerate(STREAMS):
            matches=True
            for side in(0,1):
                for phase in(0,1):
                    rng=random.Random(seeds[2*side+phase])
                    matches&=tapes[side][phase]==[rng.randrange(2)for _ in tapes[side][phase]]
                    coin_stream_checks+=1
            if matches:stream.append(i)
        need(len(pattern)==len(stream)==1,'Input/current-tape provenance ambiguous')
        domain[pattern[0],*signs,stream[0]]+=1
        catalog=row['joint_columns'];keys=[(tuple(r['input']),tuple(r['order']))for r in catalog]
        expected_keys={((x,y),order)for x in range(33)for y in range(33)for order in((0,1),(1,0))}
        need(len(keys)==2178 and set(keys)==expected_keys,'Missing/duplicate joint source column')
        results={}
        for entry in catalog:
            inputs=entry['input'];q=backgrounds[0]+backgrounds[1]
            for side,x in enumerate(inputs):
                off=side*223;q[off+target[1]]-=signs[side];q[off+target[x//11]]+=signs[side]
            outputs=list(inputs);events=[]
            for side in entry['order']:
                off=side*223;p,h=divmod(outputs[side],11)
                if h==10:physical=[off+v for v in g['blank']];c=0
                else:
                    e,c=divmod(h,2);physical=[off+target[v]for v in EDGES[e]]
                q,event=endpoint(q,physical,c);events.append(event);joint_columns+=1
                if h==10:need(event[3]==0,'Inactive column not identity')
                else:
                    measured=tables[side][p,EDGES[e],c]
                    need(event[3]==measured[0],'Joint current differs from actual source table')
                    outputs[side]=11*measured[1]+2*e+(c^bool(event[3]))
            need(outputs==entry['output'] and events==entry['events'],'Joint primitive column differs')
            for side,x in enumerate(outputs):
                off=side*223;q[off+target[x//11]]-=signs[side];q[off+target[1]]+=signs[side]
            need(q==backgrounds[0]+backgrounds[1],'Joint spectator changed')
            key=tuple(inputs);value=tuple(outputs)
            need(key not in results or results[key]==value,'Source address order changed operator');results[key]=value
        need(len(set(results.values()))==1089,'Joint map not a complete permutation')
    expected=Counter({(p,a,b,s):1 for p in(0,1)for a in(-1,1)for b in(-1,1)for s in(0,1)})
    need(domain==expected,'Actual source/current-tape grid incomplete')
    need(total_calls==4020==data['literal_continuing_source_calls'],'Continuing call census differs')
    need(joint_columns==69696==data['literal_global_column_evaluations'],'Column evaluation count differs')
    f=F(1219,2**64)
    need(F(data['classical_timeout_union_bound_each'])==f and F(data['classical_joint_success_lower'])==(1-f)**2,'Declared independent-copy bound differs')
    need(data['idle_connecting_edge']==[30,252] and data['idle_connecting_edge_events']==0,'Bridge interface differs')
    return {'passed':True,'same_group':True,'external_G6':False,'no_producer_imports':True,
            'input_sha256':hashlib.sha256(raw).hexdigest(),'actual_pairs':16,'continuing_sources':32,
            'continuing_endpoint_events':total_calls,'observed_calibration_columns':columns,
            'global_endpoint_column_evaluations':joint_columns,'complete_joint_columns_each_order':1089,
            'current_tape_stream_checks':coin_stream_checks,'exact_case_grid':True,
            'joint_dimension':1089,'graph_vertices':446,
            'source_execution_provenance':'Disclosed unregistered source pilot; registered replay, not prospective discovery.',
            'scope':'Classical continuation and complete primitive tensor map. Fresh quantum preparations/costs and pair readout remain separate.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record-to',type=Path)
    args=parser.parse_args()
    if args.record_to:need(not args.record_to.exists(),'First-record destination exists')
    started=time.perf_counter();result=calculate(PACKET/'SOURCE_JOIN.json.gz')
    if args.record_to:
        with args.record_to.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:need(result==json.loads((HERE/'PAIR_SOURCE_REVIEW.json').read_text()),'Saved pair-source review differs')
    print(json.dumps({'passed':True,'events':result['continuing_endpoint_events'],
                     'columns':result['global_endpoint_column_evaluations'],'seconds':time.perf_counter()-started}))
