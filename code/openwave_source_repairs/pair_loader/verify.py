"""Separate floor/ceiling replay. Imports no producer or native endpoint code."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import gzip,hashlib,json,random,time

HERE=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def step(x,a,b,c):
    y=list(x);s=x[a]+x[b];y[a]=s//2+c*(s%2);y[b]=s-y[a]
    return y,x[a]-y[a]
def sparse_step(p,n,a,b,c):
    v={p:1,n:-1};left=v.get(a,0);right=v.get(b,0);s=left+right
    v[a]=s//2+c*(s%2);v[b]=s-v[a];j=left-v[a]
    v={i:q for i,q in v.items() if q}
    if not v:return 0,j
    assert sorted(v.values())==[-1,1]
    return 256*next(i for i,q in v.items() if q==1)+next(i for i,q in v.items() if q==-1),j

def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG.json').read_text())
    for name,h in reg['sha256'].items():assert digest(HERE/name)==h
    evidence=json.loads(gzip.decompress((HERE/'RECORDS.json.gz').read_bytes()))
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    assert digest(HERE/'RECORDS.json.gz')==receipt['evidence_sha256']
    assert digest(HERE/'PREREG.json')==receipt['prereg_sha256']
    edges=json.loads((HERE/'graph.json').read_text())['seams']
    cells=[[e for e in edges if min(e)>=12*k and max(e)<12*k+12] for k in (0,1)]
    counts=Counter();accepted=[];forced_reader_checks=0
    for a,b,c,y,j,z,k in evidence['conjugation']:
        yp,jp=step([a,b],0,1,1-c);zp,kp=step([-a,-b],0,1,c)
        assert (yp,jp,zp,kp)==(y,j,z,k) and zp==[-v for v in yp] and kp==-jp
        counts['conjugation']+=2
    assert len(evidence['conjugation'])==162
    for trial in evidence['trials']:
        i=trial['index'];law=trial['law'];rng=random.Random(72341000+i)
        x=list(trial['initial_audit_only']);aux=[1,0];observer=[0]*240;good=True
        expected=[0]*240;expected[i]=8;expected[(i+1)%12]=-7
        expected[12+(i+5)%12]=-17;expected[12+(i+6)%12]=16;assert x==expected
        for stage in trial['stages']:
            cell=stage['cell'];base=12*cell
            for block in stage['blocks']:
                ref=block['reference'];assert len(ref)%5==0 and 5<=len(ref)<=80
                labels=[]
                for start in range(0,len(ref),5):
                    value=0
                    for c,j,occupied in ref[start:start+5]:
                        assert c==rng.randrange(2);aux,k=step(aux,0,1,c)
                        assert j==k and aux.index(1)==occupied
                        value=2*value+abs(k);counts['auxiliary']+=1
                    labels.append(value)
                assert all(v>=30 for v in labels[:-1])
                body=block['source']
                if body is None:assert len(labels)==16 and labels[-1]>=30;good=False;continue
                a,b,c,j=body;assert labels[-1]<30 and [a,b]==cells[cell][labels[-1]]
                assert c==rng.randrange(2)
                if law=='native':x,k=step(x,a,b,c);counts['selection']+=1
                else:
                    k=c*(x[a]-x[b]);x[a]-=k;x[b]+=k
                assert j==k
            before=x[:];delta=[0]*240;constraints={};reader_rng=random.Random(72342000+2*i+cell)
            for a,b,c,j in stage['census']:
                assert c==reader_rng.randrange(2)
                # Independently invert every current into possible initial contrasts.
                offset=delta[a]-delta[b]
                options={2*j-offset,2*j+(2*c-1)-offset}
                key=(a,b);constraints[key]=options if key not in constraints else constraints[key]&options
                assert constraints[key]
                x,k=step(x,a,b,c);assert k==j;counts['census']+=1
                delta[a]-=j;delta[b]+=j
            assert all(len(v)==1 for v in constraints.values()) and len(constraints)==11
            potential={base:0}
            while len(potential)<12:
                old=len(potential)
                for (a,b),s in constraints.items():
                    d=next(iter(s))
                    if a in potential:potential[b]=potential[a]-d
                    if b in potential:potential[a]=potential[b]+d
                assert len(potential)>old
            total=1 if cell==0 else -1
            shift=Fraction(total-sum(potential.values()),12);assert shift.denominator==1
            pre=[potential[k]+int(shift) for k in range(base,base+12)]
            post=[pre[k-base]+delta[k] for k in range(base,base+12)]
            assert [pre,post]==stage['decoded'] and pre==before[base:base+12] and post==x[base:base+12]
            span=max(pre)-min(pre);assert span==stage['pre_census_span']
            good &= span<=1 and len(stage['blocks'])==1024
            observer[base:base+12]=post;forced_reader_checks+=1
        assert x==trial['source_audit_only']==trial['observed_source']==observer
        assert good==trial['selected']
        if good:
            pos=[i for i,v in enumerate(x) if v==1];neg=[i for i,v in enumerate(x) if v==-1]
            assert len(pos)==len(neg)==1 and sum(abs(v) for v in x)==2
            code=256*pos[0]+neg[0];assert code==trial['pair_code'];accepted.append(code)
        else:assert trial['pair_code'] is None
        assert good==(law=='native')
    join=evidence['source_join']
    for a,b,c,oa,ob,j in join['local_table']:
        y,k=step([a,b],0,1,c);assert [oa,ob,j]==y+[k];counts['local_action_table']+=1
    seen=set();contact=outside=0
    for code,e,c,j,out in join['loaded_branches']:
        p,n=divmod(code,256);assert p<12 and 12<=n<24
        assert 0<=e<630 and c in (0,1);key=(code,e,c);assert key not in seen;seen.add(key)
        computed,k=sparse_step(p,n,*edges[e],c);assert (out,j)==(computed,k)
        contact+=out==0
        if out:pp,nn=divmod(out,256);outside+=not(pp<12 and 12<=nn<24)
        counts['loaded_action_columns']+=1
    assert len(seen)==144*1260
    assert len(join['vacuum_branches'])==1260
    for e,c,j,out in join['vacuum_branches']:
        y,k=step([0]*240,*edges[e],c);assert not any(y) and k==j==out==0
        counts['vacuum']+=1
    flags=[]
    for code in range(65536):
        p,n=divmod(code,256)
        if p>=240 or n>=240:flag='outside_range'
        elif p==n:flag='collision'
        elif p<12 and 12<=n<24:flag='accepted'
        else:flag='outside_cells'
        flags.append(flag)
    assert flags==join['product_word_flags'] and dict(Counter(flags))==join['flag_counts']==receipt['product_chart_flags']
    assert flags[0]=='collision' and contact>0 and outside>0
    assert dict(counts)==receipt['literal_native_calls_by_kind'] and sum(counts.values())==receipt['literal_native_calls']
    assert accepted==receipt['actual_pair_codes'] and len(accepted)==12 and forced_reader_checks==48
    assert receipt['all_gates_pass'] and all(receipt['gates'].values())
    result={'independent_floor_ceiling_replay':True,'imports_producer_code':False,
        'complete_loaded_branches':len(seen),'loaded_branches_to_vacuum':contact,
        'loaded_branches_outside_initial_cells':outside,'actual_pair_preparations':len(accepted),
        'current_only_stage_readers':forced_reader_checks,'checked_literal_call_records':sum(counts.values()),
        'seconds':time.monotonic()-started}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
