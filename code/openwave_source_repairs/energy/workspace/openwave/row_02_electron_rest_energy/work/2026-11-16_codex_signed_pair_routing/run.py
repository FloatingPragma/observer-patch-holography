"""Continue the actual signed source using only current-controlled native moves."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,random
from native_event import native
HERE=Path(__file__).resolve().parent
def encoded(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
EDGES=[tuple(e) for e in json.loads((HERE/'graph.json').read_text())['seams']]
CELLS=[[e for e in EDGES if 12*c<=min(e) and max(e)<12*c+12]for c in(0,1)]
CALLS=0
def path(start,target,cell):
    todo=[start];paths={start:[start]}
    for a in todo:
        if a==target:return paths[a]
        for u,v in CELLS[cell]:
            b=v if a==u else u if a==v else None
            if b is not None and b not in paths:paths[b]=paths[a]+[b];todo.append(b)
    raise AssertionError('disconnected')
def execute(initial,law,seed,forced=False):
    global CALLS
    q=tuple(initial);observed=list(initial);rng=random.Random(seed);moves=[]
    pos=observed.index(1);neg=observed.index(-1)
    route=[(a,b,1)for a,b in zip(path(pos,0,0),path(pos,0,0)[1:])]
    route +=[(a,b,-1)for a,b in zip(path(neg,14,1),path(neg,14,1)[1:])]
    route +=[(0,1,1),(1,0,1)]
    assert len(route)<=8
    for source,target,sign in route:
        assert observed[source]==sign and observed[target]==0
        edge=tuple(sorted((source,target)));assert edge in EDGES
        tape=[];arrived=False
        for _ in range(32):
            c=(1 if sign==1 else 0)if forced else rng.randrange(2)
            if law=='native':q,(_,j)=native(q,edge,c);CALLS+=1
            else:j=0
            observed[edge[0]]-=j;observed[edge[1]]+=j;tape.append([*edge,c,j])
            assert observed==list(q)
            if j:
                assert observed[source]==0 and observed[target]==sign
                arrived=True;break
        moves.append({'source':source,'target':target,'sign':sign,'records':tape,'arrived':arrived})
        if not arrived:break
    good=len(moves)==len(route)and all(m['arrived']for m in moves)
    wanted=[0]*240;wanted[0]=1;wanted[14]=-1
    if good:assert observed==wanted
    return {'law':law,'moves':moves,'planned_moves':route,'success':good,
            'terminal':'routed'if good else'timeout','source_audit_only':list(q),
            'observed_source':observed,'pair_code':14 if good else None,
            'positive_frame':[0,1]if good else None}
def calculate():
    global CALLS
    CALLS=0;inputs=json.loads((HERE/'HANDOFFS.json').read_text());rows=[]
    for h in inputs['handoffs']:
        assert h['actual_source']==h['observed_source']
        for law in('native','hold'):
            r=execute(h['observed_source'],law,72351000+h['index']);r['index']=h['index'];r['input_pair_code']=h['pair_code'];rows.append(r)
    # Deliberately choose32 holding outcomes of the actual native law.
    blank=[0]*240;blank[0]=1;blank[14]=-1
    forced=execute(blank,'native',0,True);assert forced['terminal']=='timeout'
    maximum=max(len(path(s,target,c))-1 for c,target in[(0,0),(1,14)]for s in range(12*c,12*c+12))
    assert maximum<=3
    good=[r for r in rows if r['law']=='native'];controls=[r for r in rows if r['law']=='hold']
    gates={'actual_continuations12':len(good)==12 and all(r['success']for r in good),
        'same_reader_hold_rejected12':len(controls)==12 and not any(r['success']for r in controls),
        'all_final_sources_and_frames_checked':all(r['pair_code']==14 and r['positive_frame']==[0,1]for r in good),
        'all_histories_retained':len(rows)==24,'forced_native_timeout_retained':not forced['success']}
    evidence={'trials':rows,'forced_native_timeout':forced}
    result={'candidate':'A-old240-native-current-routed-pair-v1','gates':gates,'all_gates_pass':all(gates.values()),
        'literal_native_calls':CALLS,'hold_body_calls':sum(len(m['records'])for r in controls for m in r['moves']),
        'max_target_distance':maximum,'max_requested_moves':8,'per_trial_timeout_upper':str(1-(1-F(1,2**32))**8),
        'actual_output_pair_code':14,'positive_oriented_frame':[0,1],'reset_used':False,
        'coherent_loader_or_energy_executed_here':False,'predecessor_manifest_sha256':inputs['manifest_sha256']}
    return result,evidence
def main():
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');a=p.parse_args()
    for n,h in json.loads((HERE/'PREREG.json').read_text())['sha256'].items():assert sha(HERE/n)==h
    result,evidence=calculate();raw=gzip.compress(encoded(evidence),mtime=0)
    result['evidence_sha256']=hashlib.sha256(raw).hexdigest();result['prereg_sha256']=sha(HERE/'PREREG.json')
    if a.record:
        for n,b in [('RECEIPT.json',encoded(result)),('RECORDS.json.gz',raw)]:
            with(HERE/n).open('xb')as f:f.write(b)
    else:
        assert(HERE/'RECORDS.json.gz').read_bytes()==raw
        assert json.loads((HERE/'RECEIPT.json').read_text())==result
    print(json.dumps(result,indent=2));assert result['all_gates_pass']
if __name__=='__main__':main()
