"""Prospective alternative: native two-coin AND gate, r=1/4."""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from inputs.native import native

HERE=Path(__file__).resolve().parent
graph=json.loads((HERE/'inputs/graph.json').read_text())
edges=[tuple(e) for e in graph['seams'] if max(e)<12]


def run():
    records=[]
    kernels=[]
    works=[]
    for edgeid,(a,b) in enumerate(edges):
        for h,t in ((a,b),(b,a)):
            for sign in (-1,1):
                K=[[F(0),F(0)],[F(0),F(0)]]
                for row,p in enumerate((h,t)):
                    for u,v,c in product((0,1),repeat=3):
                        for g in (u,v):
                            out,record=native((0,0),(0,1),g)
                            records.append(['gate',h,t,sign,p,[u,v,c],g,record[1],None])
                            assert out==(0,0)
                        active=p==t or (u==0 and v==0)
                        if active:
                            q=tuple(sign if i==p else 0 for i in range(12))
                            out,(coin,j)=native(q,(a,b),c)
                            # Reconstruct from public initial pointer/current.
                            observed=list(q);observed[a]-=j;observed[b]+=j
                            y,=[i for i,x in enumerate(observed) if x]
                            assert tuple(observed)==out
                            records.append(['active',h,t,sign,p,[u,v,c],c,j,y])
                        else:
                            _,(coin,j)=native((0,0),(0,1),c)
                            y=p
                            records.append(['idle',h,t,sign,p,[u,v,c],c,j,y])
                        K[row][int(y==t)]+=F(1,8)
                pi=[K[1][0]/(K[0][1]+K[1][0]),K[0][1]/(K[0][1]+K[1][0])]
                work=defaultdict(F)
                for p in (h,t):
                    for coin in (0,1):
                        q=tuple(sign if i==p else 0 for i in range(12))
                        out,(_,j)=native(q,(a,b),coin)
                        y,=[i for i,x in enumerate(out) if x]
                        records.append(['work',h,t,sign,p,[],coin,j,y])
                        work[int(y==t)-int(p==t)]+=pi[int(p==t)]/2
                kernels.append([[str(x) for x in row] for row in K])
                works.append({str(k):str(v) for k,v in sorted(work.items())})
    receipt={'candidate':'A-native-discovered-bond-quarter-v1','parent':'A-native-discovered-bond-v1',
             'scheduler_rate':'1/4','native_calls':len(records),'contexts':len(kernels),
             'kernels':kernels,'stationary_reference':['4/5','1/5'],'modular_gap':'log(4)',
             'work_distributions':works,'physical_rest_energy_established':False,
             'gates':{'kernels':all(k==[['7/8','1/8'],['1/2','1/2']] for k in kernels),
                      'work':all(w=={'-1':'1/10','0':'1/2','1':'2/5'} for w in works)}}
    return receipt,records


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    receipt,records=run()
    for name,obj in [('TUNABILITY_RECEIPT.json',receipt),('TUNABILITY_RECORDS.json',records)]:
        text=json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n'
        if a.write:
            with (HERE/name).open('x') as f:f.write(text)
        else:assert (HERE/name).read_text()==text
    print(json.dumps({'native_calls':len(records),'contexts':receipt['contexts'],'gates':receipt['gates'],'physical_rest_energy_established':False},sort_keys=True))
    raise SystemExit(0 if all(receipt['gates'].values()) else 1)
