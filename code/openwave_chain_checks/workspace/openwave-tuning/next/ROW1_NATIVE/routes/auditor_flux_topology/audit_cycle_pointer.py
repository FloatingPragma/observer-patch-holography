"""Independent source-rule and full-probability review of the native cycle meter."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict, deque
import ast
import hashlib
import importlib.util
import json
import random

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
METER=HERE.parent/'auditor_bridge/native_cycle_meter.py'
RECEIPT=METER.with_name('NATIVE_CYCLE_METER.json')


def load_target():
    spec=importlib.util.spec_from_file_location('reviewed_native_cycle_meter',METER)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def independent_endpoint(position,edge,coin):
    """Floor/ceiling repair, with no producer or native helper reused."""
    a,b=edge
    if position not in edge:return position,0
    before_a=int(position==a)
    after_a=coin
    return (a if after_a else b),before_a-after_a


def shortest(edges,start,end):
    queue=deque([start]);previous={start:None}
    while queue:
        x=queue.popleft()
        if x==end:break
        for a,b in edges:
            if x not in (a,b):continue
            y=b if x==a else a
            if y not in previous:previous[y]=x;queue.append(y)
    path=[];x=end
    while x!=start:path.append(x);x=previous[x]
    return path[::-1]


def move_law(position,destination,weight,retries):
    pending=weight;arrived=Q()
    for _ in range(retries):
        following=Q()
        for coin in (0,1):
            output,current=independent_endpoint(position,tuple(sorted((position,destination))),coin)
            if current:
                assert output==destination;arrived+=pending/2
            else:
                assert output==position;following+=pending/2
        pending=following
    return arrived,pending


def complete_law(edges,initial,q,retries,cycle):
    unresolved={initial:Q(1)};found=defaultdict(Q)
    for edge in edges:
        for _ in range(retries):
            following=defaultdict(Q)
            for position,weight in unresolved.items():
                for coin in (0,1):
                    output,current=independent_endpoint(position,edge,coin)
                    (found if current else following)[output]+=weight/2
            unresolved=following
    discovery=sum(unresolved.values(),Q());routing=Q();measurement=Q();reset=Q();success=Q()
    for position,weight in found.items():
        current=position
        for destination in shortest(edges,position,0):
            weight,failed=move_law(current,destination,weight,retries)
            routing+=failed;current=destination
        reset+=weight
        cursor=0
        for _ in range(6*abs(q)):
            next_cursor=(cursor+(1 if q>0 else -1))%6
            weight,failed=move_law(cycle[cursor],cycle[next_cursor],weight,retries)
            measurement+=failed;cursor=next_cursor
        success+=weight
    assert discovery+routing+measurement+success==1
    return {'discovery_timeout':discovery,'routing_timeout':routing,
            'measurement_timeout':measurement,'success':success,'reset_success':reset}


def replay_actual(row,initial,q,edges,cycle,retries):
    position=initial;known=None;measure_steps=0;wraps=0;cursor=0
    for record in row['records']:
        assert set(record)=={'stage','seam','coin','current','requested_cycle_direction'}
        edge=tuple(record['seam']);assert edge in edges
        output,current=independent_endpoint(position,edge,record['coin'])
        assert current==record['current']
        stage=record['stage']
        if stage=='discover' and current:known=output
        if stage=='reset_route':
            assert known is not None
            desired=shortest(edges,known,0)[0]
            assert edge==tuple(sorted((known,desired)))
            if current:known=output
        if stage=='measure':
            direction=1 if q>0 else -1
            assert record['requested_cycle_direction']==direction
            nxt=(cursor+direction)%6
            assert edge==tuple(sorted((cycle[cursor],cycle[nxt])))
            if current:
                assert output==cycle[nxt]
                if cursor==5 and nxt==0:wraps+=1
                if cursor==0 and nxt==5:wraps-=1
                cursor=nxt;measure_steps+=direction
        position=output
    assert row['final_loads']==[int(v==position) for v in range(12)]
    assert Q(row['fine_history_probability'])==Q(1,2**len(row['records']))
    if row['success']:
        assert position==0 and measure_steps==6*q and wraps==q
        assert row['winding']==q and row['retained_signed_cycle_steps']==6*q
    else:assert row['winding'] is None
    assert len(row['records'])<=39*retries


def run():
    paths=[METER,RECEIPT,ROOT/'reverse-engineering-reality/model/observer/native.py',
           ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py']
    before={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    target=load_target();reported=json.loads(RECEIPT.read_text())
    edges=tuple(target.EDGES);cycle=(0,1,8,6,2,10);retries=reported['retries_per_commanded_edge']
    assert tuple(target.CYCLE)==cycle and retries==16
    assert max(len(shortest(edges,v,0)) for v in range(12))==3
    laws=[]
    for initial in range(12):
        for q in (-1,0,1):
            law=complete_law(edges,initial,q,retries,cycle)
            assert law['reset_success']==Q(reported['exact_reset_success_by_unknown_initial_port'][initial])
            assert law['discovery_timeout']==Q(1,2**(5*retries))
            laws.append({'initial_port':initial,'winding':q,**{k:str(v) for k,v in law.items()}})
    worst=max(1-Q(row['success']) for row in laws)
    assert worst==Q(reported['compressed_primary_meter_worst_timeout'])
    rng=random.Random(112031)
    replays=0;reasons=defaultdict(int)
    for initial in range(12):
        for q in (-1,0,1):
            patterns=[[0]*624,[1]*624,[i%2 for i in range(624)],[(i//3)%2 for i in range(624)]]
            patterns += [[rng.randrange(2) for _ in range(624)] for _ in range(12)]
            for coins in patterns:
                row=target.execute(initial,(q,),iter(coins),retries)
                replay_actual(row,initial,q,edges,cycle,retries)
                replays+=1;reasons[row['reason']]+=1
    # Every moving local native C branch is rankone; every hold branch is
    # diagonal. Hence first conclusive discovery destroys no hidden pointer
    # location information needed by the controller and later reset maps it
    # to0. This proves reset on reference matrix units as well as basis inputs.
    rankone=0
    for edge in edges:
        for coin in (0,1):
            for current in (-1,0,1):
                columns=[]
                for position in range(12):
                    output,j=independent_endpoint(position,edge,coin)
                    if j==current:columns.append((position,output))
                if current:
                    assert len(columns)<=1;rankone+=bool(columns)
                else:assert all(a==b for a,b in columns)
    assert rankone==60
    assert before=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return {'schema':'oph.row1.independent_cycle_pointer_review.v1','checks':{
        'all36_exact_unconditional_output_laws_agree':True,'all_timeout_channels_sum_with_success_to_one':True,
        'unknown_pointer_discovery_effects_and_reset_reference_scope':True,
        'independent_native_rule_replays_every_observed_record':True,
        'all_successes_count_actual_oriented_pointer_wraps':True,
        'no_timeout_or_refusal_decodes_as_zero_charge':True,'attempt_cap624_checked':True,
        'no_hidden_E_contrast_record':True,'source_pins_unchanged':True},
        'complete_output_laws':laws,'actual_executor_replays':replays,'observed_reasons':dict(reasons),
        'compressed_worst_timeout':str(worst),'compressed_worst_timeout_float':float(worst),
        'maximum_attempts':624,'rankone_moving_C_branches':rankone,
        'scope':'Primary one-load pointer, bounded q=0,+1,-1 commands. Independent floor/ceiling reference and exact absorbing probabilities; actual target execute() called only for comparison. Joint source-to-counter executor is reviewed separately.',
        'source_sha256':before,'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'official_rows_changed':[]}


if __name__=='__main__':
    result=run();(HERE/'CYCLE_POINTER_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'replays':result['actual_executor_replays'],'timeout':result['compressed_worst_timeout_float']}))
