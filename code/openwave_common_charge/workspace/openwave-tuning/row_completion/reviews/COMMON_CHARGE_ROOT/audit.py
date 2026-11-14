"""Independent coupled-charge geometry, Gauss and complete native-meter audit.

No producer Python module is imported. The six-state response generators,
declared graph data and frozen pulse tapes are data inputs, not verdicts.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import sympy as s
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PACK = ROOT/'openwave-tuning/next/CHAIN5/common_charge_topology'
EXPECTED = 'b45e64648b2379cdc0629a8789f74ca9a798a211e58e078d5c9f573e1513eb4f'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def custody():
    m = read(PACK/'MANIFEST.json')
    assert sha(PACK/'MANIFEST.json') == EXPECTED
    pins = {str((PACK/'MANIFEST.json').relative_to(ROOT)):EXPECTED, **m['files'], **m['source_pins']}
    for name,digest in pins.items(): assert sha(ROOT/name) == digest, name
    return pins

def geometry(pins):
    base = ROOT/'openwave-tuning/next/CHAIN3/coulomb_response_connection/inputs'
    wp = base/'openwave-tuning/auditor/native_walk_bridge/COMMON_WALK_RECEIPT.json'
    gp = base/'reverse-engineering-reality/model/native_spatial_refinement/results.json'
    for p in (wp,gp): pins[str(p.relative_to(ROOT))] = sha(p)
    internal = [tuple(e) for matching in read(wp)['descriptor']['internal_matchings'] for e in matching]
    links = read(gp)['intercell_oriented_port_pairs']
    coords = list(product(range(2),repeat=3)); cid = {c:i for i,c in enumerate(coords)}
    edges = set()
    for cell,c in enumerate(coords):
        edges.update(tuple(sorted((12*cell+a,12*cell+b))) for a,b in internal)
        for axis,pairs in enumerate(links):
            dest = list(c); dest[axis] ^= 1; other = cid[tuple(dest)]
            edges.update(tuple(sorted((12*cell+a,12*other+b))) for a,b in pairs)
    edges = sorted(edges); assert len(edges)==288
    B = s.zeros(96,288)
    for e,(a,b) in enumerate(edges): B[a,e],B[b,e]=1,-1
    assert all(sum(abs(B[v,e]) for e in range(288))==6 for v in range(96))
    inv = (B*B.T)[1:,1:].inv(method='DM')
    q=s.zeros(96,8)
    for c in range(8): q[0,c]+=1; q[12*c,c]-=1
    potential=s.zeros(1,8).col_join(inv*q[1:,:]); P=B.T*potential
    assert B*P==q
    em=read(PACK/'EMBEDDING.json'); D=em['coordinate_denominator']
    vertices=em['vertex_integer_coordinates']; start=em['oriented_edge_start_lifts'];end=em['oriented_edge_end_lifts']
    assert len(set(map(tuple,vertices)))==96
    for i,(a,b) in enumerate(edges):
        assert start[i]==vertices[a]
        assert all((end[i][j]-vertices[b][j])%(2*D)==0 for j in range(3))
    # Independent exact segment test, via determinant and a 2x2 minor.
    sub=lambda a,b:tuple(int(x)-int(y) for x,y in zip(a,b))
    cross=lambda a,b:(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    tested=0
    for shift in product(range(-2,3),repeat=3):
        delta=tuple(2*D*x for x in shift)
        for i,(a,b) in enumerate(zip(start,end)):
            low=np.minimum(a,b); high=np.maximum(a,b)
            # This filter only discards disjoint bounding boxes.
            cs=np.asarray(start,dtype=np.int64)+delta; ds=np.asarray(end,dtype=np.int64)+delta
            js=np.flatnonzero(np.all(np.minimum(cs,ds)<=high,axis=1)&np.all(np.maximum(cs,ds)>=low,axis=1))
            for j in js:
                if shift==(0,0,0) and j<=i: continue
                c,d=cs[j],ds[j];u=sub(b,a);v=sub(d,c);w=sub(c,a);n=cross(u,v);tested+=1
                if dot(n,w): continue
                if any(n):
                    axis=next(k for k in range(3) if n[k]); t=F(cross(w,v)[axis],n[axis]);r=F(cross(w,u)[axis],n[axis])
                    if 0<=t<=1 and 0<=r<=1:
                        assert t in (0,1) and r in (0,1)
                        assert tuple(a if t==0 else b)==tuple(c if r==0 else d)
                elif not any(cross(w,u)):
                    k=next(k for k in range(3) if u[k]);left=max(min(a[k],b[k]),min(c[k],d[k]));right=min(max(a[k],b[k]),max(c[k],d[k]))
                    assert left>=right
    topo=read(PACK/'TOPOLOGY.json');sphere=topo['sphere_complex']; faces=sphere['oriented_faces']
    occurrences=defaultdict(list)
    for f,tri in enumerate(faces):
        for a,b in zip(tri,tri[1:]+tri[:1]): occurrences[tuple(sorted((a,b)))].append((f,1 if a<b else -1))
    adjacent=defaultdict(set)
    for pairs in occurrences.values():
        assert len(pairs)==2 and pairs[0][1]+pairs[1][1]==0
        a,b=pairs[0][0],pairs[1][0];adjacent[a].add(b);adjacent[b].add(a)
    seen={0}
    while True:
        expanded=seen.union(*(adjacent[x] for x in seen))
        if expanded==seen:break
        seen=expanded
    assert len(seen)==12 and 8-len(occurrences)+12==2
    # Connected dual graph forces all integral 2-cycle coefficients equal;
    # all-one coefficients therefore give the primitive integral generator.
    routes=read(PACK/'DYNAMICS.json')['routes']
    for row in routes:
        J=s.zeros(288,1)
        for e,value in row['integer_current_nonzero']: J[e]=value
        old,new=row['cell'],row['destination'];path=row['native_path'];ind=s.zeros(288,1)
        lookup={e:i for i,e in enumerate(edges)}
        for a,b in zip(path,path[1:]):ind[lookup[tuple(sorted((a,b)))]]+=1 if a<b else -1
        assert ind==J and B*J==q[:,new]-q[:,old]
        longitudinal=P[:,new]-P[:,old]; transverse=J-longitudinal
        assert B*transverse==s.zeros(96,1) and longitudinal.dot(transverse)==0
    for region in topo['regions'].values():
        indicator=s.zeros(96,1)
        for v in region['vertices']:indicator[v]=1
        flux=indicator.T*q;cut=B.T*indicator
        assert cut.T*P==flux
        assert list(flux)==region['flux_eigenvalues_by_cell']
    # Fresh overlapping/disconnected regions, plus arbitrary integer cycles.
    for seed in range(1,58):
        indicator=s.Matrix([int((seed*(v+7)+v*v)%17<6) for v in range(96)])
        flux=indicator.T*q
        assert (B.T*indicator).T*P==flux and all(x in (-1,0,1) for x in flux)
    return {'independent_periodic_segment_tests':tested,'primitive_sphere':True,'exact_routes':len(routes),'fresh_regions':57,'exact_source_edges':len(edges)}

def native_audit(pins):
    gp=ROOT/'openwave-tuning/next/CHAIN3/native_field_connection/inputs/reverse-engineering-reality/code/source_selection_model/response.json'
    assert str(gp.relative_to(ROOT)) in pins
    generators=[]
    for blocks in read(gp)['generators']:
        mat=np.zeros((6,6),complex)
        for k,block in enumerate(blocks):
            for i,row in enumerate(block):
                for j,x in enumerate(row):mat[3*k+i,3*k+j]=float(F(x[0]))+math.sqrt(5)*float(F(x[1]))+1j*(float(F(x[2]))+math.sqrt(5)*float(F(x[3])))
        assert np.linalg.norm(mat+mat.conj().T)<1e-13;generators.append(mat)
    @lru_cache(None)
    def response(controls):return expm(sum((float(F(c))*g for c,g in zip(controls,generators)),np.zeros((6,6),complex)))
    def execute(tape,bits,stored):
        dim=1<<bits;M=np.eye(dim,dtype=complex);failure=np.zeros((dim,dim),complex);bound=F()
        for num,event in enumerate(tape):
            targets=event['targets'];other=[q for q in range(bits) if q not in targets];local=1<<len(targets)
            ix=np.array([[sum(((a>>(len(targets)-1-j))&1)<<q for j,q in enumerate(targets))+sum(((b>>j)&1)<<q for j,q in enumerate(other)) for b in range(1<<len(other))]for a in range(local)])
            old=M[ix].copy();U=response(tuple(event['controls_exact']))
            good=(U[:local,:local]@old.reshape(local,-1)).reshape(old.shape)
            bad=(U[local:,:local]@old.reshape(local,-1)).reshape(-1,dim);effect=bad.conj().T@bad
            failure+=effect;M[ix]=good;bound+=F(event['full_response_operator_error_exact'])
            assert np.max(abs(np.diag(effect).real-stored['first_export_records'][num]['probability_on_every_input_basis_column']))<2e-13
        expected=np.array(stored['failure_effect_real'])+1j*np.array(stored['failure_effect_imag'])
        assert np.linalg.norm(failure-expected)<2e-12
        gram=float(np.linalg.norm(M.conj().T@M+failure-np.eye(dim),2));assert gram<3e-10
        return M,failure,bound,gram
    tapes=read(PACK/'NATIVE_PROGRAMS.json');failures=read(PACK/'NATIVE_FIRST_FAILURES.json');results={};pulses=0;columns=0
    for name,tape in tapes.items():
        bits=3 if name=='native_READY_pointer' else (6 if '+' in name else 5)
        results[name]=execute(tape,bits,failures[name]);pulses+=len(tape);columns+=1<<bits
    pointer,pfailure,pbound,_=results['native_READY_pointer'];selected=pointer[:,[0,4]];effects={}
    # Decode the four actually reached native total-two configurations.
    loads=[(0,0,1,1)+(0,)*8,(0,1,1)+(0,)*9,(1,0,1)+(0,)*9,(0,2)+(0,)*10]
    for grammar in ('C','E'):
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)));meta={}
        for code,state in enumerate(loads):
            K=selected[[code,code+4]]
            for coin in (0,1):
                d=state[0]-state[1];current=(d-(2*coin-1)*(abs(d)%2))//2
                out=list(state);out[0]-=current;out[1]+=current
                key=(coin,current if grammar=='C' else d)
                groups[key][tuple(out)]+=K/math.sqrt(2);meta[key]={'coin':coin,'label':key[1],'current':current,'bit':int(2*current+2*coin-1<0)}
        effects[grammar]=[(meta[k],sum((K.conj().T@K for K in groups[k].values()),np.zeros((2,2),complex)))for k in sorted(groups)]
    records=read(PACK/'NATIVE_RECORDS.json');topo=read(PACK/'TOPOLOGY.json')['regions'];finecount=0;maxerror=0.
    for law,rows in effects.items():
        for (meta,E),archive in zip(rows,records['fine_pointer_effects'][law]):
            assert meta==archive['record']
            assert np.linalg.norm(E-(np.array(archive['effect_real'])+1j*np.array(archive['effect_imag'])))<2e-13
    for row in records['experiments']:
        names=row['regions'];n=len(names);clock=3+n;label='+'.join(names)+'_'+row['clock_setting'];M,failed,wordbound,_=results[label]
        phase=row['macros']/64;cv=np.array([np.exp(-.5j*phase),np.exp(.5j*phase)])/math.sqrt(2)
        # The ideal clock is used here, independent of its producer runner.
        initial=np.zeros((1<<(4+n),8),complex)
        for c in range(8):initial[c,c]=cv[0]/math.sqrt(8);initial[c+(1<<clock),c]=cv[1]/math.sqrt(8)
        out=M@initial;rho=sum((out[c::8]@out[c::8].conj().T for c in range(8)),np.zeros((1<<(n+1),1<<(n+1)),complex))
        for law,archive in row['laws'].items():
            actual=[];coarse=np.zeros(1<<(n+1))
            for terms in product(*(effects[x] for x in law)):
                E=np.ones((1,1),complex);code=0
                for j,(meta,local) in enumerate(terms):E=np.kron(local,E);code|=meta['bit']<<j
                probability=float(np.trace(E@rho).real);actual.append(probability);coarse[code]+=probability
            err=max(abs(x-y['probability'])for x,y in zip(actual,archive['fine_records']));maxerror=max(maxerror,err);assert err<5e-12;finecount+=len(actual)
            assert np.max(abs(coarse-archive['all_joint_binary_probabilities']))<5e-12
            if n==2:assert sum(coarse[i]for i in range(len(coarse))if i&3==3)<1e-12
            combined=sum(actual)
            for f in range(n+1):
                E=np.ones((1,1),complex)
                for j in range(n+1):
                    local=sum((a for _,a in effects[law[j]]),np.zeros((2,2),complex)) if j<f else (pfailure[np.ix_([0,4],[0,4])] if j==f else np.eye(2))
                    E=np.kron(local,E)
                combined+=float(np.trace(E@rho).real)
            assert abs(combined-np.trace(rho).real)<3e-11
        reset=row['extra_blank_reset_contract'];timeout=F(reset['complete_timeout_union_bound_exact'])
        assert reset['registers']==3*n+2 and timeout==reset['registers']*F(reset['base_uniform_failure_exact'])**reset['blocks_per_register']
        assert F(row['additional_meter_bound_exact'])==wordbound+(n+1)*pbound+timeout
        assert F(row['complete_common_source_record_TV_bound_exact'])==min(F(1),F(row['all_source_prefix_bound_exact'])+F(row['additional_meter_bound_exact']))
    assert finecount==12844
    integer=results['integer_2pi_loop'][0][:,1];half=results['half_charge_2pi_mutant'][0][:,1]
    assert sum(abs(integer[i])**2 for i in range(32) if i&8)<1e-12
    assert sum(abs(half[i])**2 for i in range(32) if i&8)>1-1e-12
    return {'programs':len(tapes),'pulses':pulses,'input_columns':columns,'fine_C_E_probabilities':finecount,'maximum_probability_residual':maxerror,'maximum_Gram_residual':max(x[3]for x in results.values()),'complete_failure_and_reset_budgets_checked':True}

def main():
    pins=custody();geom=geometry(pins);print('geometry',geom,flush=True)
    native=native_audit(pins);print('native',native,flush=True)
    assert custody().items()<=pins.items()
    result={'schema':'oph.common_charge.root.audit.v1','checks_pass':True,'geometry':geom,'native':native,'source_pins':pins,'code_sha256':sha(Path(__file__)),'full_criterion_met':False,'scope':'Scientific checks only; full interpretation is in the separate review. Complete huge source prefix not executed.'}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
