#!/usr/bin/env python3
"""Independent complete reader; imports neither native nor producer."""
import gzip,hashlib,json,math
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
P=Path(__file__).resolve().parent

def main():
    receipt=json.loads((P/'RECEIPT.json').read_text());seen=0;calls=0;summary=[]
    with gzip.open(P/'RAW_RECORDS.jsonl.gz','rt') as stream:
        def read(expected):
            nonlocal seen
            actual=json.loads(next(stream));assert actual==expected,(actual,expected);seen+=1
        def current(s,a,b,c):
            # For a binary pair, the resulting first endpoint is c on mixed inputs.
            return s[a]-c if s[a]!=s[b] else 0
        for case in receipt['cases']:
            d,L=case['d'],case['L'];verts=list(product(range(L),repeat=d));ids={v:i for i,v in enumerate(verts)}
            V=len(verts);edges=[]
            for a,v in enumerate(verts):
                for axis in range(d):
                    w=list(v);w[axis]=(w[axis]+1)%L;edges.append((a,ids[tuple(w)],axis))
            W=len(edges);assert case['W']==W
            assert case['vertices']==[list(v) for v in verts] and case['bonds']==[list(e) for e in edges]
            kernels={};populations={};coeffs=[]
            for q in (1,2):
                pop=list(combinations(range(V),q));index={s:i for i,s in enumerate(pop)};rows=[];populations[q]=pop
                for si,s in enumerate(pop):
                    start=[int(v in s) for v in range(V)];row=Counter()
                    for ei,(a,b,axis) in enumerate(edges):
                        for c in (0,1):
                            j=current(start,a,b,c);after=start.copy();after[a]-=j;after[b]+=j
                            ti=index[tuple(v for v,n in enumerate(after) if n)]
                            read(['kernel',d,L,q,si,ei,c,j,ti]);calls+=1;row[ti]+=1
                    rows.append(row)
                assert all(sum(r.values())==2*W for r in rows)
                assert all(n==rows[j][i] for i,row in enumerate(rows) for j,n in row.items())
                kernels[q]=rows
            bondlabels={tuple(sorted((a,b))):(a,i) for a,b,i in edges}
            for axis in range(d):
                g=[(-1)**sum(verts[bondlabels[s][0]]) if s in bondlabels and bondlabels[s][1]==axis else 0 for s in populations[2]]
                assert g==case['pair_mode_coefficients'][axis] and sum(v*v for v in g)==V
                assert all(sum(g[j] for j,s in enumerate(populations[2]) if x in s)==0 for x in range(V))
                assert all(2*W*g[i]-sum(c*g[j] for j,c in row.items())==2*(2*d-1)*g[i] for i,row in enumerate(kernels[2]))
                coeffs.append(g)
            # Check the whole one-source matrix independently, sufficient for exact free-sector cancellation.
            for a,row in enumerate(kernels[1]):
                expected=Counter({a:2*W-2*d})
                for x,y,_ in edges:
                    if x==a:expected[y]+=1
                    if y==a:expected[x]+=1
                assert row==expected
            assert F(case['pair_J_energy'])==F(2*d-1,W)
            assert F(case['free_entire_K_threshold'])==F(2*d,W)
            assert F(case['transfer_lower_bound'])==min(F(row[i],W)-1 for i,row in enumerate(kernels[2]))>0
            core_sums=[];point_sums=[]
            for bi,(x,y,axis) in enumerate(edges):
                zz=list(verts[y]);zz[axis]=(zz[axis]+1)%L;z=ids[tuple(zz)]
                for label,sx,sy in [('core',(x,y),(y,z)),('point',(x,),(y,))]:
                    X=[int(v in sx) for v in range(V)];Y=[int(v in sy) for v in range(V)]
                    vec=[];div=[0]*V
                    for ei,(a,b,_) in enumerate(edges):
                        total=0
                        for c in (0,1):
                            jx=current(X,a,b,c);jy=current(Y,a,b,c);calls+=2
                            read(['response',d,L,label,bi,ei,c,jx,jy]);total+=jx-jy
                        assert total==(X[a]-Y[a])-(X[b]-Y[b])
                        vec.append(total);div[a]+=total;div[b]-=total
                    # This record-derived divergence is the point-measured graph Laplacian of the independently prepared source.
                    (core_sums if label=='core' else point_sums).append(vec)
            assert core_sums==case['response_pair_edge_coin_sums']
            assert point_sums==case['response_point_edge_coin_sums']
            for kind in ('core','point'):
                for bi in range(W):
                    for ei in range(W):
                        for c in (0,1):read(['hold',d,L,kind,bi,ei,c,0,0])
            assert F(case['intrinsic_radius_squared'])==F(1,4)
            assert F(case['directional_radius_squared'])==F(1,4*d)
            summary.append({'d':d,'L':L,'all_pair_rows_checked':len(populations[2]),'all_native_response_inputs':2*W,
                            'bound_mode_residual':0,'pair_response_coefficient_residual':0,'point_response_coefficient_residual':0})
        assert stream.read()==''
    assert seen==receipt['raw_records'] and calls==receipt['native_calls']
    # Exact radius and cyclotomic values at executed L=12 momenta. Irrational values get rational enclosures.
    def sqrt_interval(n,bits=100):
        den=1<<bits;t=math.isqrt(n*den*den);return F(t,den),F(t+1,den)
    s2=sqrt_interval(2);s3=sqrt_interval(3);s6=sqrt_interval(6)
    intervals=[((s6[0]+s2[0])/4,(s6[1]+s2[1])/4),(s3[0]/2,s3[1]/2),
               (s2[0]/2,s2[1]/2),(F(1,2),F(1,2)),((s6[0]-s2[1])/4,(s6[1]-s2[0])/4),(F(0),F(0))]
    assert intervals[2][0]>F(1,2) and intervals[4][1]<F(1,2)
    out={'pass':True,'raw_records':seen,'literal_native_calls_in_producer':calls,
         'receipt_sha256':hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest(),'cases':summary,
         'L12_effective_coupling_intervals':[[str(a),str(b)] for a,b in intervals],
         'full_native_and_current_replay':True,'imports_producer_or_native':False,
         'scope':'Exact current response and spectral compact core. No coherent source-state preparation, continuum physical core or universal running claim.'}
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'pass':True,'records':seen,'producer_native_calls':calls}))

if __name__=='__main__':main()
