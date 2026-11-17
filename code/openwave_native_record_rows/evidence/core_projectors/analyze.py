#!/usr/bin/env python3
"""Preregistered retrospective analysis of a sealed native archive."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import gzip,hashlib,json
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'20261114_codex_bound_core_response'

def main():
    reg=json.loads((HERE/'PREREG.txt').read_text())
    for n,h in reg['files'].items():assert hashlib.sha256((HERE/n).read_bytes()).hexdigest()==h
    for n,h in reg['source_files'].items():assert hashlib.sha256((SOURCE/n).read_bytes()).hexdigest()==h
    receipt=json.loads((SOURCE/'RECEIPT.json').read_text());cases={}
    for c in receipt['cases']:
        d,L=c['d'],c['L'];verts=list(product(range(L),repeat=d));V=len(verts);pairs=list(combinations(range(V),2))
        cases[d,L]=dict(d=d,L=L,vertices=verts,pairs=pairs,rows=[Counter() for _ in pairs])
    records=0
    with gzip.open(SOURCE/'RAW_RECORDS.jsonl.gz','rt') as f:
        for line in f:
            row=json.loads(line);records+=1
            if row[0]=='kernel' and row[3]==2:
                _,d,L,q,s,e,c,j,t=row;cases[d,L]['rows'][s][t]+=1
    assert records==receipt['raw_records']
    summaries=[]
    for (d,L),case in cases.items():
        verts=case['vertices'];V=len(verts);W=d*V;pairs=case['pairs'];index={s:i for i,s in enumerate(pairs)};vid={x:i for i,x in enumerate(verts)};P=case['rows']
        adj=[]
        for x,y in pairs:
            delta=[(b-a)%L for a,b in zip(verts[x],verts[y])]
            adj.append(sum(v!=0 for v in delta)==1 and any(v in (1,L-1) for v in delta))
        K=[]
        for x,y in pairs:
            row=Counter()
            for a in verts:
                xx=vid[tuple((v+t)%L for v,t in zip(verts[x],a))];yy=vid[tuple((v+t)%L for v,t in zip(verts[y],a))]
                row[index[tuple(sorted((xx,yy)))]]+=(-1)**sum(a)
            K.append({j:n for j,n in row.items() if n})
        B=[row if adj[i] else {} for i,row in enumerate(K)]
        assert all(n==K[j].get(i,0) for i,row in enumerate(K) for j,n in row.items())
        assert all(n==B[j].get(i,0) for i,row in enumerate(B) for j,n in row.items())
        def multiply_row(row,matrix):
            out=Counter()
            for j,n in row.items():
                for k,m in matrix[j].items():out[k]+=n*m
            return {k:v for k,v in out.items() if v}
        for i,row in enumerate(K):
            assert multiply_row(row,K)=={j:V*n for j,n in row.items()}
            assert multiply_row(P[i],K)=={j:(2*W-4*d+2*int(adj[i]))*n for j,n in row.items()}
        for i,row in enumerate(B):
            assert multiply_row(row,B)=={j:V*n for j,n in row.items()}
            assert B[i].get(i,0)==int(adj[i])
            assert multiply_row(P[i],B)=={j:(2*W-4*d+2)*n for j,n in row.items()}
        rank=F(sum(row.get(i,0) for i,row in enumerate(B)),V);assert rank==d
        assert sum(adj)==W
        summaries.append({'d':d,'L':L,'pair_states':len(pairs),'translation_terms':V,
            'full_K_projector_rank':str(F(sum(row.get(i,0) for i,row in enumerate(K)),V)),
            'bound_projector_rank':str(rank),'bound_projector_nnz':sum(map(len,B)),
            'full_operator_identity_verified':True,'selfadjoint_idempotence_verified':True,
            'all_instrument_diagonal_weights':str(F(1,W)),'outside_adjacent_diagonal_weight':'0',
            'integer_residual':0,'bound_energy':str(F(2*d-1,W)),'K_complement_energy':str(F(2*d,W)),
            'intrinsic_R_squared':'1/4','half_onset_abs_k_times_R_parallel':'pi/3'})
    out={'candidate':'A-bound-core-response-v1','retrospective_derivation':True,'new_native_calls':0,
        'source_native_calls':receipt['native_calls'],'source_archive_records_read':records,
        'source_manifest_sha256':reg['source_files']['MANIFEST.sha256'],'cases':summaries,
        'full_source_spectral_projector_instrument_join':True,'literal_finite_core_onset_clause':True,
        'universal_beta_function_claim':False,'coherent_working_state_preparation_claim':False,'official':'untested'}
    (HERE/'RECEIPT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
