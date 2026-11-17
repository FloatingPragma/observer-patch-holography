#!/usr/bin/env python3
"""Retrospective exact integer source-graph audit; no producer import."""
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
source=OW/'shared/compact_logarithmic_flow/source_columns.npz'
data=np.load(source);edges=[tuple(map(int,e)) for e in data['seams']]
body=sorted((a,b) for a,b in edges if b<12)
assert len(body)==30
for cell in range(20):
    actual=sorted((a-12*cell,b-12*cell) for a,b in edges if a//12==cell==b//12)
    assert actual==body
B=np.zeros((12,12),dtype=np.int64)
for a,b in body:B[a,a]+=1;B[b,b]+=1;B[a,b]-=1;B[b,a]-=1
seen={0}
while True:
    nxt=seen|{j for i in seen for j in range(12) if B[i,j]==-1}
    if nxt==seen:break
    seen=nxt
assert len(seen)==12 and np.all(np.diag(B)==5)
B2=B@B;B3=B2@B;B4=B3@B
assert np.array_equal(B4-16*B3+80*B2-120*B,np.zeros((12,12),dtype=np.int64))
assert np.trace(B)==60 and np.trace(B2)==360 and np.trace(B3)==2280
# Connectedness gives multiplicity one at zero; rationality and the degree
# trace then force multiplicities3,5,3 at5-sqrt5,6,5+sqrt5.
external=[(a,b) for a,b in edges if a//12!=b//12]
degrees=np.zeros(240,dtype=int)
for a,b in external:degrees[a]+=1;degrees[b]+=1
assert len(external)==30 and max(degrees)==1
result={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'body_count':20,'body_vertices':12,'body_edges':30,'external_matching_edges':30,
 'exact_body_minimal_polynomial':'x(x-6)(x^2-10x+20)',
 'body_spectral_multiplicities':{'0':1,'5-sqrt(5)':3,'6':5,'5+sqrt(5)':3},
 'exterior_laplacian_bounds':'0<=L_external<=2I',
 'minmax_first80_upper':'7-sqrt(5)<6','minmax_remaining160_lower':'6',
 'exact_modes_strictly_below6':80,'not_a_new_experiment':True}
target=HERE/'LOWER_MODES_RESULT.json'
with target.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps(result))
