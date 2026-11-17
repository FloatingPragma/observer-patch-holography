#!/usr/bin/env python3
"""Exact native branch census and transfer-derived finite record moments."""
import importlib.util
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(100000)
HERE = Path(__file__).resolve().parent

def inv(a):
    n=len(a)
    q=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k=next(k for k in range(j,n) if q[k][j])
        q[j],q[k]=q[k],q[j]
        d=q[j][j]; q[j]=[x/d for x in q[j]]
        for i in range(n):
            if i!=j:
                d=q[i][j]; q[i]=[x-d*y for x,y in zip(q[i],q[j])]
    return [row[n:] for row in q]

def moments(T, sizes):
    """All event words summed exactly; no literal long-tape execution claim."""
    m=len(T); labels=sum(T[0]); g=[1,-1]+[0]*(m-2)
    z=[[1]*m]+[[0]*m for _ in range(4)]
    den=m; out={}
    for n in range(1,max(sizes)+1):
        nxt=[[0]*m for _ in range(5)]
        for i,row in enumerate(T):
            for j,c in enumerate(row):
                if c:
                    for k in range(5):
                        nxt[k][j]+=c*sum(math.comb(k,l)*z[l][i]*g[j]**(k-l) for l in range(k+1))
        z=nxt; den*=labels
        if n in sizes:
            raw=[F(sum(v),den) for v in z]
            y2=raw[2]/(10*n); y4=raw[4]/(100*n*n)
            v2=F(1,60); v4=F(3,3600)-F(1,3600*n)
            f2=y2+v2; f4=y4+6*y2*v2+v4
            out[n]={'mean':raw[1], 'variance':f2, 'fourth':f4,
                    'fourth_cumulant':f4-3*f2*f2}
    return out

def run():
    source=HERE/'inputs/native.py'
    assert hashlib.sha256(source.read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('occupation_native_source',source)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    faces=json.loads((HERE/'inputs/faces.json').read_text())
    edges=sorted({tuple(sorted((f[i],f[(i+1)%3]))) for f in faces for i in range(3)})
    assert len(edges)==30
    T=[[0]*12 for _ in range(12)]; records=[]; L=[[0]*12 for _ in range(12)]
    for a,b in edges:
        L[a][a]+=1; L[b][b]+=1; L[a][b]-=1; L[b][a]-=1
    for x in range(12):
        initial=tuple(int(i==x) for i in range(12))
        for e,(a,b) in enumerate(edges):
            for c in (0,1):
                final,(coin,current)=mod.native(initial,(a,b),c)
                read=list(initial);read[a]-=current;read[b]+=current
                assert tuple(read)==final and sum(final)==1 and max(final)==1
                y=read.index(1);T[x][y]+=1
                records.append([x,e,coin,current,y])
    coins=[]
    for address in range(12):
        for c in (0,1):
            final,record=mod.native((0,0),(0,1),c)
            assert final==(0,0) and record==(c,0)
            coins.append([address,*record])
    kernel_ok=all(T[i][j]==60*int(i==j)-L[i][j] for i in range(12) for j in range(12))
    # (L+J/M)^-1-J/M is the exact pseudoinverse, independently computed.
    B=inv([[F(L[i][j])+F(1,12) for j in range(12)] for i in range(12)])
    G=[[B[i][j]-F(1,12) for j in range(12)] for i in range(12)]
    target=G[0][0]+G[1][1]-2*G[0][1]
    ground=inv([row[:-1] for row in L[:-1]])
    root_ok=all(ground[i][j]==G[i][j]-G[i][11]-G[11][j]+G[11][11] for i in range(11) for j in range(11))
    sizes=(1,2,8,64,512,2048)
    native=moments(T,sizes)
    hold=moments([[60*int(i==j) for j in range(12)] for i in range(12)],sizes)
    checks={'native_branch_kernel':kernel_ok,'rooted_dirichlet_identity':root_ok,
            'source_effective_resistance':target==F(11,30),
            'mean_zero':all(v['mean']==0 for v in native.values()),
            'finite_covariance_bound':all(0<=target-v['variance']<=F(92,5*n) for n,v in native.items()),
            'variance_error_decreases':all(target-native[b]['variance']<target-native[a]['variance'] for a,b in zip(sizes,sizes[1:])),
            'fourth_gaussian_approach':abs(native[2048]['fourth']-3*target**2)<3*target**2/8,
            'hold_changes_same_variance':hold[512]['variance']>10*target,
            'doubled_record_detected':4*native[512]['variance']>3*target}
    def serial(v):
        if isinstance(v,F):return {'exact':str(v),'decimal':float(v)}
        if isinstance(v,dict):return {str(k):serial(x) for k,x in v.items()}
        if isinstance(v,list):return [serial(x) for x in v]
        return v
    return serial({'candidate':'A-native-Gaussian-v1','native_calls':744,
                   'literal_scope':'720 one-unit branches and 24 auxiliary coin branches only; block moments are exact transfer compression',
                   'records':records,'auxiliary_records':coins,'native_count_matrix':T,
                   'green_target':target,'native_moments':native,'hold_moments':hold,
                   'checks':checks,'full_physical_rows_passed':False})

if __name__=='__main__':
    result=run()
    (HERE/'RECEIPT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'native_calls':result['native_calls']},sort_keys=True))
    sys.exit(0 if all(result['checks'].values()) else 1)
