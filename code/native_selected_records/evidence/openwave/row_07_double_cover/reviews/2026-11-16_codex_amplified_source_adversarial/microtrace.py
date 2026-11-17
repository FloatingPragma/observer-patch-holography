"""Independent exact audit of the declared V,V,V-adjoint current eraser.

This is a finite coherent-source access assumption, not three forward events.
No producer code is imported. Source spectator backgrounds are actual filed
admission outputs; their basis-column audit is not a new preparation run.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
OPENWAVE=HERE.parents[2]
ADMISSION=OPENWAVE/'row_01_charge_quantization/work/2026-11-16_codex_record_only_admission_v2'
EDGES=((0,1),(0,2),(1,2))


def require(ok,message):
    if not ok:raise AssertionError(message)


def endpoint(q,edge,coin,law):
    a,b=edge;d=q[a]-q[b]
    j=(d-(2*coin-1)*(abs(d)%2))//2 if law=='native' else d
    out=list(q);out[a]-=j;out[b]+=j
    return out,j


def dilation(edge,sign,law):
    # Label = (position*2+coin)*3+current_mod3. Only J=0 fixes
    # the physical native columns; the other columns are a declared extension.
    forward={}
    for p in range(3):
        for c in (0,1):
            q=[0]*3;q[p]=sign
            out,j=endpoint(q,edge,c,law)
            y=out.index(sign)
            forward[(p*2+c)*3]=(y*2+c)*3+j%3
    require(len(set(forward.values()))==6,'Native isometry lost orthogonality')
    unused_inputs=[i for i in range(18) if i not in forward]
    unused_outputs=[i for i in range(18) if i not in forward.values()]
    forward.update(zip(unused_inputs,unused_outputs))
    require(sorted(forward.values())==list(range(18)),'Incomplete unitary extension')
    inverse={y:x for x,y in forward.items()}
    return forward,inverse


def unpack(label):
    pc,j=divmod(label,3);p,c=divmod(pc,2)
    return p,c,j


def circuit(p,c,j1,j2,forward,inverse,trace=False):
    states=[(p,c,j1,j2)]
    p,c,j1=unpack(forward[(p*2+c)*3+j1]);states.append((p,c,j1,j2))
    c^=int(j1!=0);states.append((p,c,j1,j2))
    p,c,j2=unpack(forward[(p*2+c)*3+j2]);states.append((p,c,j1,j2))
    j1=(j1+j2)%3;states.append((p,c,j1,j2))
    p,c,j2=unpack(inverse[(p*2+c)*3+j2]);states.append((p,c,j1,j2))
    return states if trace else states[-1]


def calculate():
    evidence=ADMISSION/'records/EVIDENCE.json.gz'
    data=json.loads(gzip.decompress(evidence.read_bytes()))
    catalog={};all_states=[(p,c,j1,j2) for p in range(3) for c in (0,1)
                          for j1 in range(3) for j2 in range(3)]
    for edge in EDGES:
        for sign in (-1,1):
            for law in ('native','swap'):
                f,inv=dilation(edge,sign,law)
                out=[circuit(*s,f,inv) for s in all_states]
                require(set(out)==set(all_states),'Whole eraser circuit is not unitary')
                traces=[]
                for p in range(3):
                    for c in (0,1):
                        q=[0]*3;q[p]=sign;y,j=endpoint(q,edge,c,law)
                        expected=(y.index(sign),c^int(j!=0),0,0)
                        tr=circuit(p,c,0,0,f,inv,True)
                        require(tr[-1]==expected,'Current eraser did not retain transported unit')
                        # At inverse entry, the SECOND source tuple is exactly
                        # its own earlier output, so V-adjoint stays on its image.
                        require(tr[4][0:2]+(tr[4][3],)==tr[3][0:2]+(tr[3][3],),
                                'Current addition modified inverse source image')
                        traces.append({'input':[p,c],'states':tr})
                catalog[str((edge,sign,law))]={'V':[f[i] for i in range(18)],
                     'full54permutation':out,'code_traces':traces}
    digest=hashlib.sha256();columns=0;native_calls=0
    for row in data['rows']:
        if row['control'] or row['controller']['outcome']!='signed_target_unit':continue
        background=row['audit']['final'];n=12*row['B'];ports=(0,n,n+1)
        sign=row['controller']['sign']
        require([background[p] for p in ports]==[0,sign,0],'Not an actual admitted source')
        for edge in EDGES:
            physical=tuple(ports[a] for a in edge)
            for law in ('native','swap'):
                f,inv=dilation(edge,sign,law)
                for p in range(3):
                    for c in (0,1):
                        initial=list(background)
                        for port in ports:initial[port]=0
                        initial[ports[p]]=sign
                        first,j=endpoint(initial,physical,c,law)
                        cc=c^int(j!=0)
                        second,j2=endpoint(first,physical,cc,law);native_calls+=2
                        require(second==initial and j+j2==0,'Actual inverse image fails')
                        final=circuit(p,c,0,0,f,inv)
                        reconstructed=list(background)
                        for port in ports:reconstructed[port]=0
                        reconstructed[ports[final[0]]]=sign
                        require(reconstructed==first and final[1:]==(cc,0,0),
                                'Global source/current workspace failed')
                        require(all(first[i]==background[i] for i in range(len(first)) if i not in ports),
                                'Spectator reservoir or markers changed')
                        digest.update(json.dumps([row['B'],row['M'],initial,physical,c,law,
                                      first,j,second,j2,final],separators=(',',':')).encode())
                        columns+=1
    require(columns==864*3*2*3*2,'Missing global signed/background columns')
    return {'passed':True,'same_group':True,'external_G6':False,
            'access_assumption':'Recorded native isometry plus coherent inverse on its image and explicit unitary extension',
            'per_R_calls':{'forward':2,'inverse':1},'physical_forward_only_claim':False,
            'admission_evidence_sha256':hashlib.sha256(evidence.read_bytes()).hexdigest(),
            'complete_extensions':catalog,'complete54state_permutations':12,
            'global_basis_columns':columns,'literal_forward_endpoint_evaluations':native_calls,
            'inverse_table_applications':columns,'global_trace_sha256':digest.hexdigest(),
            'new_physical_preparations':0,'spectator_backgrounds':864}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record-to',type=Path)
    args=parser.parse_args()
    if args.record_to:require(not args.record_to.exists(),'First-record destination exists')
    start=time.perf_counter();out=calculate()
    if args.record_to:
        with args.record_to.open('x') as stream:json.dump(out,stream,indent=2);stream.write('\n')
    else:require(json.dumps(out,sort_keys=True)==json.dumps(json.loads((HERE/'MICROTRACE.json').read_text()),sort_keys=True),'Saved microtrace audit differs')
    print(json.dumps({'passed':True,'columns':out['global_basis_columns'],'seconds':time.perf_counter()-start}))
