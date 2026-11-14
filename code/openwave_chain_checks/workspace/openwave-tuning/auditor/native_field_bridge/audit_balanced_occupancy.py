"""Opposing review of equal-efficiency native occupancy, using independent repair."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import random

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=ROOT/'openwave-tuning/auditor/row1_public_integration/balanced_occupancy.py'
RECEIPT=TARGET.with_name('BALANCED_OCCUPANCY.json')
PRIOR=ROOT/'openwave-tuning/next/ROW1_NATIVE/routes/auditor_flux_topology/audit_cycle_pointer.py'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def favorable(review,edges,cycle,initial,degree,retries):
    pos=initial;coins=[];found=False
    for edge in edges:
        for _ in range(retries):
            coin=int(pos==edge[1]) if pos in edge else 0
            pos,j=review.independent_endpoint(pos,edge,coin);coins.append(coin)
            if j:found=True;break
        if found:break
    for dest in review.shortest(edges,pos,0):
        coin=int(pos>dest);pos,_=review.independent_endpoint(pos,tuple(sorted((pos,dest))),coin);coins.append(coin)
    path=(1,0,1,0,1,0) if degree==0 else tuple(cycle[(degree*k)%6] for k in range(1,7))
    for dest in path:coins.append(int(pos>dest));pos=dest
    return coins


def replay(review,row,initial,degree,edges,cycle,retries):
    pos=initial;measure_moves=0;signed=0;wrap=0
    expected=(1,0,1,0,1,0) if degree==0 else tuple(cycle[(degree*k)%6] for k in range(1,7))
    for event in row['records']:
        edge=tuple(event['seam']);assert edge in edges
        out,j=review.independent_endpoint(pos,edge,event['coin']);assert j==event['current']
        if event['stage']=='balanced_measure':
            target=expected[measure_moves]
            assert edge==tuple(sorted((pos,target)))
            if j:
                assert out==target
                a,b=cycle.index(pos),cycle.index(out)
                step=1 if b==(a+1)%6 else -1
                assert step==event['requested_cycle_direction']
                signed+=step;measure_moves+=1
                if (a,b)==(5,0):wrap+=1
                if (a,b)==(0,5):wrap-=1
        pos=out
    assert row['final_loads']==[int(v==pos) for v in range(12)]
    assert F(row['fine_history_probability'])==F(1,2**len(row['records']))
    assert len(row['records'])<=39*retries
    if row['success']:
        assert pos==0 and measure_moves==6 and signed==6*degree and wrap==degree
        assert row['winding']==wrap
    else:assert row['winding'] is None


def run():
    monitored=[TARGET,RECEIPT,PRIOR]
    before={str(p.relative_to(ROOT)):sha(p) for p in monitored}
    assert sha(TARGET)=='a23385a0bdea821ec47b41c353ddbcbc45088ffb99ae6b78d259e782105b9732'
    target=load('opposing_balanced',TARGET);review=load('prior_independent_endpoint',PRIOR)
    edges=tuple(target.meter.EDGES);cycle=(0,1,8,6,2,10)
    rng=random.Random(611347);laws=[];runs=0;outcomes=Counter()
    for retries in (1,2,16):
        for initial in range(12):
            # Six commanded hops; the independently enumerated reset and
            # geometric hop law do not depend on the path's orientation.
            law=review.complete_law(edges,initial,1,retries,cycle)
            reported=target.complete_law(initial,retries)
            assert reported=={'success':law['success'],'timeout':1-law['success']}
            for degree,cls in ((-1,5),(0,0),(1,1)):
                tapes=[favorable(review,edges,cycle,initial,degree,retries),
                       [0]*(39*retries),[1]*(39*retries)]
                tapes += [[rng.randrange(2) for _ in range(39*retries)] for _ in range(4)]
                for coins in tapes:
                    row=target.native_meter(initial,cls,iter(coins),retries)
                    replay(review,row,initial,degree,edges,cycle,retries)
                    runs+=1;outcomes[row['reason']]+=1
            laws.append({'pointer_initial':initial,'retries':retries,
                         'all_three_class_success':str(law['success']),
                         'stage_law':{k:str(v) for k,v in law.items() if k!='reset_success'}})
    # Source hold/move plus each pointer column, with fresh fair pointer
    # randomness: both occupancy effects have exactly the same efficiency.
    joined=0;G=[1]+[0]*11
    for initial in range(12):
        for source_coin in (0,1):
            degree=1-source_coin
            coins=favorable(review,edges,cycle,initial,degree,2)
            row=target.execute_occupancy(G,((0,1),),(source_coin,),1,initial,iter(coins),2)
            assert row['success'] and row['occupancy']==source_coin
            replay(review,row['pointer'],initial,degree,edges,cycle,2)
            assert F(row['combined_history_probability'])==F(1,2**len(row['combined_native_tape']))
            assert row['clock']['source_during_pointer']=='idle under the declared addressed scheduler'
            joined+=1
    # Correlating the apparatus input with source output is a concrete
    # counterexample to an unqualified conclusive-normalization claim.
    eta0=review.complete_law(edges,0,1,1,cycle)['success']
    eta1=review.complete_law(edges,1,1,1,cycle)['success']
    biased=eta0/(eta0+eta1)
    assert eta0!=eta1 and biased!=F(1,2)
    # Declared independent arbitrary mixture has common averaged efficiency.
    for numerator in range(8):
        mix=F(numerator,7);eta=mix*eta0+(1-mix)*eta1
        assert (eta/2)/(eta/2+eta/2)==F(1,2)
    after={str(p.relative_to(ROOT)):sha(p) for p in monitored};assert before==after
    result={'schema':'oph.native.balanced_occupancy.opposing_review.v1',
      'checks':{'independent_endpoint_all_meter_replays':True,'equal_success_all_three_classes':True,
                'complete_failure_partitions_sum_one':True,'zero_path_actual_degree_zero':True,
                'actual_occupancy_output_from_pointer':True,'source_and_pointer_clock_retained':True,
                'correlated_pointer_counterexample_retained':True,'source_bytes_unchanged':True},
      'actual_pointer_replays':runs,'joined_occupancy_replays':joined,'exact_class_laws':len(laws)*3,
      'outcomes':dict(outcomes),'laws':laws,
      'independence_scope':'Success balancing is conditional on fixed pointer input and fresh conditionally fair native pointer coins. For an averaged pointer input, its preparation must be independent of the source history; marginal fairness or an unspecified correlated apparatus does not suffice.',
      'correlated_pointer_control':{'prior_source_survival':'1/2','pointer0_if_survives_pointer1_if_moves':True,
                                    'eta0':str(eta0),'eta1':str(eta1),'conclusive_source_survival':str(biased)},
      'source_pins':before,'auditor_sha256':sha(Path(__file__)),'official_rows_changed':[]}
    (HERE/'BALANCED_OCCUPANCY_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'runs':runs,'joined':joined,'correlated_control':result['correlated_pointer_control']}))


if __name__=='__main__':run()
