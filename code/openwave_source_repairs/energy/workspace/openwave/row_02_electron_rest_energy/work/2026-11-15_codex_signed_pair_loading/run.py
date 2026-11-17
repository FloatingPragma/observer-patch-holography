"""One actual signed source: native selection, census and code-action join."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import argparse,gzip,hashlib,json,random
from pathlib import Path
from native_event import native

HERE=Path(__file__).resolve().parent
GRAPH=json.loads((HERE/'graph.json').read_text());EDGES=[tuple(e) for e in GRAPH['seams']]
CELLS=[[e for e in EDGES if 12*c<=min(e) and max(e)<12*c+12] for c in (0,1)]
CALLS=Counter()

def call(x,e,c,kind):
    CALLS[kind]+=1
    return native(x,e,c)

def tree(cell):
    root=12*cell;seen={root};queue=[root];out=[]
    for a in queue:
        for u,v in CELLS[cell]:
            b=v if a==u else u if a==v else None
            if b is not None and b not in seen:seen.add(b);queue.append(b);out.append((a,b))
    assert len(out)==11
    return out

TREES=[tree(c) for c in (0,1)]

def selector(aux,rng,cell):
    tape=[]
    for _ in range(16):
        label=0
        for _ in range(5):
            c=rng.randrange(2);aux,(_,j)=call(aux,(0,1),c,'auxiliary')
            tape.append([c,j,aux.index(1)]);label=2*label+abs(j)
        if label<30:return aux,CELLS[cell][label],tape
    return aux,None,tape

def stage(x,aux,rng,cell,law):
    blocks=[]
    for _ in range(1024):
        aux,e,tape=selector(aux,rng,cell)
        if e is None:blocks.append({'reference':tape,'source':None});break
        c=rng.randrange(2)
        if law=='native':x,(_,j)=call(x,e,c,'selection')
        else:
            a,b=e;j=c*(x[a]-x[b]);y=list(x);y[a]-=j;y[b]+=j;x=tuple(y)
        blocks.append({'reference':tape,'source':[*e,c,j]})
    return x,aux,blocks

def census(x,cell,seed):
    rng=random.Random(seed);tape=[]
    for a,b in TREES[cell]:
        first=None
        for _ in range(32):
            c=rng.randrange(2);x,(_,j)=call(x,(a,b),c,'census');tape.append([a,b,c,j])
            if first is None:first=c
            elif c!=first:break
        else:return x,tape
    return x,tape

def read(tape,cell):
    base=12*cell;delta=[0]*12;v={0:0};cursor=0
    for aa,bb in TREES[cell]:
        a,b=aa-base,bb-base;correction=delta[a]-delta[b];first=None
        for _ in range(32):
            if cursor==len(tape):return None
            x,y,c,j=tape[cursor];cursor+=1;assert(x,y)==(aa,bb)
            delta[a]-=j;delta[b]+=j
            if first is None:first,j0=c,j
            elif c!=first:v[b]=v[a]-(2*j0+j-correction);break
        else:return None
    assert cursor==len(tape)
    total=1 if cell==0 else -1;offset=F(total-sum(v.values()),12)
    pre=[F(v[i])+offset for i in range(12)];assert all(z.denominator==1 for z in pre)
    pre=list(map(int,pre));post=[a+b for a,b in zip(pre,delta)]
    return pre,post

def pair_code(x):
    positive=[i for i,v in enumerate(x) if v==1];negative=[i for i,v in enumerate(x) if v==-1]
    assert len(positive)==len(negative)==1 and sum(abs(v) for v in x)==2
    return 256*positive[0]+negative[0]

def histories():
    results=[]
    for i in range(12):
        initial=[0]*240;initial[i]=8;initial[(i+1)%12]=-7
        initial[12+(i+5)%12]=-17;initial[12+(i+6)%12]=16
        for law in ('native','swap_hold'):
            x=tuple(initial);aux=(1,0);rng=random.Random(72341000+i);stages=[];observed=[0]*240
            for cell in (0,1):
                x,aux,blocks=stage(x,aux,rng,cell,law)
                x,tape=census(x,cell,72342000+2*i+cell);decoded=read(tape,cell)
                stage_record={'cell':cell,'blocks':blocks,'census':tape,'decoded':decoded}
                if decoded is not None:
                    stage_record['pre_census_span']=max(decoded[0])-min(decoded[0])
                    observed[12*cell:12*cell+12]=decoded[1]
                stages.append(stage_record)
            good=all(s['decoded'] is not None and s['pre_census_span']<=1 and len(s['blocks'])==1024
                     and all(b['source'] is not None for b in s['blocks']) for s in stages)
            if good:assert list(x)==observed and 0<=pair_code(observed)//256<12 and 12<=pair_code(observed)%256<24
            results.append({'index':i,'law':law,'initial_audit_only':initial,'stages':stages,
                'source_audit_only':x,'observed_source':observed,'selected':good,
                'pair_code':pair_code(observed) if good else None})
    return results

def conjugation():
    rows=[]
    for a,b,c in product(range(-4,5),range(-4,5),(0,1)):
        y,(_,j)=call((a,b),(0,1),1-c,'conjugation')
        z,(_,k)=call((-a,-b),(0,1),c,'conjugation')
        assert z==tuple(-v for v in y) and k==-j
        rows.append([a,b,c,list(y),j,list(z),k])
    return rows

def action_join():
    table={}
    for a,b,c in product((-1,0,1),(-1,0,1),(0,1)):
        out,(_,j)=call((a,b),(0,1),c,'local_action_table');table[(a,b,c)]=(*out,j)
    rows=[]
    for p in range(12):
        for n in range(12,24):
            initial=[0]*240;initial[p]=1;initial[n]=-1;initial=tuple(initial)
            for e,(a,b) in enumerate(EDGES):
                for c in (0,1):
                    out,(_,j)=call(initial,(a,b),c,'loaded_action_columns')
                    qa=int(p==a)-int(n==a);qb=int(p==b)-int(n==b)
                    left,right,current=table[(qa,qb,c)]
                    sparse={p:1,n:-1};sparse[a]=left;sparse[b]=right;sparse={i:v for i,v in sparse.items() if v}
                    encoded=0 if not sparse else 256*next(i for i,v in sparse.items() if v==1)+next(i for i,v in sparse.items() if v==-1)
                    literal=0 if not any(out) else pair_code(out)
                    assert encoded==literal and j==current
                    rows.append([256*p+n,e,c,j,encoded])
    vacuum=[]
    for e,edge in enumerate(EDGES):
        for c in (0,1):
            out,(_,j)=call((0,)*240,edge,c,'vacuum');assert not any(out) and j==0
            vacuum.append([e,c,j,0])
    labels=[]
    for code in range(65536):
        p,n=divmod(code,256)
        flag='outside_range' if max(p,n)>=240 else 'collision' if p==n else 'accepted' if p<12 and 12<=n<24 else 'outside_cells'
        labels.append(flag)
    counts=dict(Counter(labels));assert counts=={'collision':240,'outside_cells':57216,'accepted':144,'outside_range':7936}
    assert labels[0]=='collision'
    return {'local_table':[[a,b,c,*v] for(a,b,c),v in sorted(table.items())],
            'loaded_branches':rows,'vacuum_branches':vacuum,'product_word_flags':labels,'flag_counts':counts}

def encoded(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()

def calculate():
    CALLS.clear();conj=conjugation();trials=histories();join=action_join()
    native_trials=[r for r in trials if r['law']=='native'];controls=[r for r in trials if r['law']=='swap_hold']
    gates={'conjugation_exact':len(conj)==162,'native_actual_pairs12':all(r['selected'] for r in native_trials),
           'control_selected0':not any(r['selected'] for r in controls),'loaded_action_mismatches0':len(join['loaded_branches'])==181440,
           'all_product_words_tagged':len(join['product_word_flags'])==65536,'vacuum_distinguished_from_collision':len(join['vacuum_branches'])==1260 and join['product_word_flags'][0]=='collision'}
    evidence={'conjugation':conj,'trials':trials,'source_join':join}
    receipt={'candidate':'A-old240-native-signed-pair-loading-v1','gates':gates,'all_gates_pass':all(gates.values()),
        'literal_native_calls_by_kind':dict(CALLS),'literal_native_calls':sum(CALLS.values()),
        'actual_pair_codes':[r['pair_code'] for r in native_trials],'product_chart_flags':join['flag_counts'],
        'complete_quantum_sector_AND_implemented_here':False,'source_data_loading':'identity concatenation on admitted disjoint-cell chart',
        'quantum_loader_consumer':'energy packet native sector/AND instrument with terminalcomplements',
        'full_row_solved':False,'physical_energy_or_spin_claimed':False}
    return receipt,evidence

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args()
    reg=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in reg['sha256'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    receipt,evidence=calculate();raw=gzip.compress(encoded(evidence),mtime=0)
    receipt['evidence_sha256']=hashlib.sha256(raw).hexdigest();receipt['prereg_sha256']=hashlib.sha256((HERE/'PREREG.json').read_bytes()).hexdigest()
    if args.record:
        for name,data in [('RECEIPT.json',encoded(receipt)),('RECORDS.json.gz',raw)]:
            with(HERE/name).open('xb')as f:f.write(data)
    else:
        assert(HERE/'RECORDS.json.gz').read_bytes()==raw
        assert json.loads((HERE/'RECEIPT.json').read_text())==receipt
    print(json.dumps(receipt,indent=2));assert receipt['all_gates_pass']

if __name__=='__main__':main()
