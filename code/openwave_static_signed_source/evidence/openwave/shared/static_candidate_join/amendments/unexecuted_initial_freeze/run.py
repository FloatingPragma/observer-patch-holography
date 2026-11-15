#!/usr/bin/env python3
"""Frozen new signed-source observations. --record writes only this packet."""
from __future__ import annotations
import argparse,gzip,hashlib,importlib.util,json,random,sys,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form
import source as C
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
LOCAL=HERE.parent/'localized_rest_join';COMPACT=HERE.parent/'compact_static_adjudication'
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
SEEDS=(714131,714137)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n'
def module(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def pointers(arrays):
    from model.observer.source_scalar import compiler
    from model.native_interface import pointer
    from model.native_clock import configuration
    local=pointer.local_matrix();target=pointer.target();cert=pointer.certificate()
    assert np.linalg.norm(local-target,'fro')<1e-10
    effects={};errors=[]
    for law in ('C','E'):
        count=pointer.native_effects(law);sup=np.zeros((2,2,2,2,2),complex)
        for b in (0,1):
            for z in range(4):
                K=local[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
                sup[b]+=count[b][z]/2*np.einsum('ac,bd->abcd',K,K.conj())
        wanted=np.zeros_like(sup);wanted[0,0,0,0,0]=wanted[1,1,1,1,1]=1
        err=float(np.max(abs(sup-wanted)));errors.append(err);assert err<1e-10
        effects[law]=np.array([np.einsum('aacd->dc',S) for S in sup])
        arrays['pointer_superoperators_'+law]=sup
    prep={};position={law:[] for law in ('C','E')};maxerr=0.
    for p in range(12):
        word=[compiler.pulse('X',np.pi/2,[b]) for b in range(8) if (p>>b)&1]
        U=compiler.execute(word,8);v=U[:,0]
        wanted=np.zeros(256,complex);wanted[p]=(1j)**p.bit_count()
        maxerr=max(maxerr,float(np.linalg.norm(v-wanted)));assert maxerr<1e-10
        assert np.linalg.norm(U.conj().T@U-np.eye(256),'fro')<1e-10
        prep[p]=v;arrays[f'basis_preparation_{p}']=U
        for law in ('C','E'):
            probs=[]
            for z in range(256):
                E=np.array([[1]],complex)
                for bit in reversed(range(8)):E=np.kron(E,effects[law][(z>>bit)&1])
                probs.append(float(np.vdot(v,E@v).real))
            assert abs(sum(probs)-1)<1e-9 and abs(probs[p]-1)<1e-9
            position[law].append(probs)
    for law in ('C','E'):arrays['basis_position_law_'+law]=np.array(position[law])
    # Full numerical work/reset/read composition, not a conditional-success rate.
    t=configuration.budget(full=True,mode='reset',target=pointer.ready_state())['failure_upper']
    e=F(1,10**11);blocks=64;gain=(1+e)**(2*blocks)
    bound=(gain-1)/2+gain*32*t
    assert bound<F(1,10**6)
    def read(label):
        # Ideal complete source-word instrument; actual words are certified above.
        # Observer receives only native current records and decodes these bits.
        bits=[];records=[]
        for bit in range(16):
            private=(label>>bit)&1
            column=pointer.READY_CODE+4*private
            out=int(np.argmax(target[:,column]));assert target[out,column]==1 and out//4==private
            aux=pointer.compiled.circuit.decode(out%4);alternatives=[]
            for coin in (0,1):
                _,(_,j)=C.S.native.native(aux,pointer.SEAM,coin)
                reported=int(2*j+2*coin-1<0)
                alternatives.append([coin,j,reported])
            assert alternatives[0][2]==alternatives[1][2]
            bits.append(alternatives[0][2]);records.append([bit,alternatives])
        return sum(b<<i for i,b in enumerate(bits)),records
    return {'basis_words':12,'complete_columns_per_word':256,'pointer_matrix_unit_max_error':max(errors),
        'basis_preparation_norm_error':maxerr,'auxiliary_uniform_timeout':str(t),
        'complete_prepare_read_error_upper':str(bound),'conservative_error_upper':'1/1000000',
        'all_position_labels':256,'raw_pointer_outcome_alternatives':'both native coins per bit, not sequential events',
        'body_semantics':'exact admitted source on ideal service inputs; complete actual-word comparison includes failure'},read

def incidence():
    B=Matrix.zeros(12,30);D=Matrix.zeros(30,20);index={e:i for i,e in enumerate(C.CELL)}
    for i,(a,b) in enumerate(C.CELL):B[a,i]=1;B[b,i]=-1
    for k,face in enumerate(FACES):
        for a,b in zip(face,face[1:]+face[:1]):D[index[tuple(sorted((a,b)))],k]+=1 if a<b else -1
    assert B*D==Matrix.zeros(12,20)
    smith=smith_normal_form(D,domain=ZZ);diag=[abs(int(smith[i,i])) for i in range(20)]
    assert diag==[1]*19+[0] and B.rank()==11
    L=B*B.T;J=Matrix.ones(12,12)/12;plus=(L+J).inv()-J
    # Spanning trees certify the primitive boundary image on both actual graphs.
    trees={}
    for name,n,edges in [('cell',12,C.CELL),('federation',240,C.SEAMS)]:
        seen={0};chosen=[]
        while len(seen)<n:
            e=next(e for e in edges if (e[0] in seen)!=(e[1] in seen));chosen.append(e);seen.update(e)
        trees[name]=[list(e) for e in chosen]
    return B,D,plus,{'cell_boundary_rank':11,'cell_face_smith':diag,'cell_kernel_equals_face_image':True,
        'federation_primitive_image':'C1/ker(B) = zero-total Z^240; no global sphere-face identification','spanning_trees':trees}

def source_branches(plus,B):
    labels=[256*p+r for p in range(12) for r in range(12) if p!=r];rows=[];lookup=set(labels)|{0};partial={};effect=Counter();bell=0
    expected={}
    for z in labels:
        p,r=divmod(z,256);q=Matrix([int(i==p)-int(i==r) for i in range(12)])
        expected[z]=[F(v) for v in B.T*plus*q]
    expected[0]=[F(0)]*30
    for z in labels+[0]:
        rhs=[F(0)]*30
        for e,seam in enumerate(C.CELL):
            for coin in (0,1):
                out,j,flag=C.endpoint(z,seam,coin);assert flag=='body' and out in lookup
                key=e,coin,j;bucket=partial.setdefault(key,set());assert out not in bucket;bucket.add(out)
                effect[z]+=1;rows.append([z,e,coin,j,out])
                for k in range(30):rhs[k]+=expected[out][k]/60
                rhs[e]+=F(j,60)
        assert rhs==expected[z];bell+=30
    assert set(effect.values())=={60}
    return rows,{'inputs_including_vacuum':133,'branches':len(rows),'exact_bellman_entries':bell,
        'all_fine_branches_partial_injections':True,'complete_effects_identity':True},labels

def incidence_update(q,seam,current):
    out=q.copy();a,b=seam
    out[a]=out.get(a,0)-current;out[b]=out.get(b,0)+current
    return {k:v for k,v in out.items() if v}

def trajectory(z,seed,mode,cap,read):
    rng=random.Random(seed);observed,pointer=read(z);assert observed==z
    hidden=z;q=C.represented(observed);q0=q.copy();flux=[0]*30;records=[]
    for _ in range(cap):
        if not q:break
        e=rng.randrange(30);coin=rng.randrange(2);a,b=C.CELL[e]
        if mode=='native':out,j,flag=C.endpoint(hidden,(a,b),coin);assert flag=='body'
        else:
            private=C.represented(hidden);d=private.get(a,0)-private.get(b,0)
            j=d if mode=='swap' and coin==0 else 0
            private[a]=private.get(a,0)-j;private[b]=private.get(b,0)+j;out=C.basis(private)
        hidden=out;records.append([e,coin,j]);q=incidence_update(q,(a,b),j);flux[e]+=j
    assert C.represented(hidden)==q
    boundary=[0]*12
    for (a,b),j in zip(C.CELL,flux):boundary[a]+=j;boundary[b]-=j
    assert boundary==[q0.get(i,0)-q.get(i,0) for i in range(12)]
    return {'initial_pointer':pointer,'source_readout':observed,'seed':seed,'mode':mode,
        'cap':cap,'records':records,'conclusive':not q,'terminal_source':hidden,'flux':flux}

def moments(rows):
    M=module('static_moment_arithmetic',HERE.parent/'paired_source_moments/run.py')
    points,anti,rotation,cross=M.geometry(C.CELL);stats=[];bad=0
    for p in range(12):
        r=anti[p];z=256*p+r;vectors=[];fixed=[]
        for source,e,coin,j,out in rows:
            if source!=z:continue
            a,b=C.CELL[e];v=[M.mul(M.q(j),M.sub(y,x)) for x,y in zip(points[a],points[b])];vectors.append(v)
            if C.CELL[e]==(0,5):fixed.append(v)
        stat=M.moments(vectors);audit=M.analyze(stat,points[p]);control=M.analyze(M.moments(fixed),points[p])
        assert audit['mean_nonzero'] and audit['covariance_nonzero']
        assert M.iszero(audit['mean_axis_residual']) and M.iszero(audit['uniaxial_residual'])
        assert audit['parallel_variance']==M.div(M.sub(M.q(3),M.q(0,1)),M.q(9))
        assert audit['transverse_variance']==M.div(M.add(M.q(0,1),M.q(2)),M.q(15))
        expected=[M.mul(M.div(M.sub(M.q(0,1),M.q(3)),M.q(15)),v) for v in points[p]]
        assert stat['mean']==expected
        bad+=not (M.iszero(control['mean_axis_residual']) and M.iszero(control['uniaxial_residual']))
        stats.append({'p':p,'r':r,'moments':stat,'audit':audit,'fixed_seam_control':control})
    a=stats[0]['moments'];b=stats[rotation[0]]['moments']
    assert b['mean']==M.rotate_vector(a['mean']) and b['covariance']==M.rotate_tensor(a['covariance'])
    assert b['mean']==[M.neg(x) for x in a['mean']] and b['covariance']==a['covariance']
    assert M.rotate_vector(b['mean'])==a['mean'] and a['covariance']!=stats[1]['moments']['covariance']
    assert bad>0
    plant=M.analyze(M.moments([[M.q(1),M.q(0),M.q(0)],[M.q(-1),M.q(0),M.q(0)]]),points[0])
    assert not M.iszero(plant['uniaxial_residual'])
    return M.encode(stats),{'source_inputs':12,'branch_alternatives':720,'fixed_seam_failed_inputs':bad,
        'proper_half_turn_tensor_return':True,'polar_frame_reversal':True,'squared_frame_return':True,
        'covariance_nonconstant':True,'off_axis_plant_detected':True,'continuous_rotation_status':'fitted representation; finite proper rotations executed'},anti

def old_pair_map():
    raw=C.S.native.native;different=0;rows=[];counts=Counter()
    for p in range(12):
      for r in range(12):
       if p==r:continue
       for s in range(12):
        if s in (p,r):continue
        x=tuple(int(i in (p,s)) for i in range(12));y=tuple(int(i in (r,s)) for i in range(12));q={p:1,r:-1}
        for e,(a,b) in enumerate(C.CELL):
         for c in (0,1):
          xo,(_,jx)=raw(x,(a,b),c);yo,(_,jy)=raw(y,(a,b),c);k=jx-jy;d=q.get(a,0)-q.get(b,0)
          effective=(1+d-2*k)//2 if abs(d)%2 else c
          assert effective in (0,1)
          out,j,flag=C.endpoint(256*p+r,(a,b),effective);observed={i:u-v for i,(u,v) in enumerate(zip(xo,yo)) if u!=v}
          assert C.represented(out)==observed and j==k and flag=='body'
          direct,jdirect,_=C.endpoint(256*p+r,(a,b),c);different+=int(direct!=out or jdirect!=k)
          counts[p,r,e,effective]+=1;rows.append([p,r,s,e,c,jx,jy,effective,out])
    assert set(counts.values())=={10} and different>0
    # Same effective flag; old original coin, common-background and individual
    # currents distinguish the components. This kills their reduced coherence.
    witness=[]
    for p,r,c in ((0,2,0),(3,4,1)):
        s=1;x=tuple(int(i in (p,s)) for i in range(12));y=tuple(int(i in (r,s)) for i in range(12))
        xo,(_,jx)=raw(x,(0,1),c);yo,(_,jy)=raw(y,(0,1),c);d=int(p==0)-int(r==0)-int(p==1)+int(r==1);k=jx-jy
        eff=(1+d-2*k)//2 if abs(d)%2 else c;assert eff==1 and k==0
        assert tuple(a-b for a,b in zip(xo,yo))==tuple(int(i==p)-int(i==r) for i in range(12))
        witness.append({'q':[p,r],'old_coin':c,'old_currents':[jx,jy],'effective_coin':eff,'relative_current':k,'common_output':[i for i in range(12) if xo[i]==yo[i]==1]})
    assert witness[0]['old_coin']!=witness[1]['old_coin'] and witness[0]['old_currents']!=witness[1]['old_currents']
    return rows,{'paired_inputs':1320,'paired_branches':len(rows),'original_coin_disagreements':different,
        'effective_coin_multiplicity_per_relative_context':10,'coherent_counterexample':witness,
        'coarse_branch_probability':'1/2','new_signed_unnormalized_offdiagonal':'1/4','old_fine_record_reduced_offdiagonal':'0','coherent_instrument_equivalence_claimed':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args();start=time.perf_counter()
    if args.record and (HERE/'RECEIPT.json').exists():raise RuntimeError('refuse to overwrite evidence')
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for path,h in prereg['files'].items():assert digest(ROOT/path)==h,path
    calls=0;native=C.S.native.native
    def counted(*a):
        nonlocal calls
        calls+=1;return native(*a)
    C.S.native.native=counted
    assert len(C.SEAMS)==630 and len(C.CELL)==30
    old_seams=json.loads((ROOT/'openwave/row_01_charge_quantization/work/2026-11-14_codex_relative_source/source/seams.json').read_text())
    assert [list(e) for e in C.CELL]==old_seams
    arrays={};pointer,read=pointers(arrays);B,D,plus,topology=incidence();branches,branch_result,labels=source_branches(plus,B)
    mom,moment_result,anti=moments(branches)
    trajectories=[trajectory(z,seed+10000*(z//256)+100*(z%256),'native',2048,read) for seed in SEEDS for z in labels]
    controls=[trajectory(256*p+anti[p],714139+100*p+c,mode,256,read) for mode in ('swap','hold') for p in range(12) for c in (0,1)]
    assert not any(t['conclusive'] for t in controls)
    paired,paired_result=old_pair_map()
    q0={0:F(1),2:F(-1)};q=incidence_update(q0,(0,1),F(1,2))
    integer_mutant=incidence_update(q0,(0,1),int(F(1,2)))
    detected=integer_mutant!=q
    assert q[0]==q[1]==F(1,2) and detected
    plant={'initial':[0,2],'seam':[0,1],'current':'1/2','regions':{'0':str(q[0]),'1':str(q[1])},'native':False,'constant_integer_decoder_rejected':detected,'mutant_regions':{'0':str(integer_mutant.get(0,0)),'1':str(integer_mutant.get(1,0))}}
    compact=json.loads((COMPACT/'RECEIPT.json').read_text());assert compact['all_pass']
    raw={'branches':branches,'trajectories':trajectories,'controls':controls,'moments':mom,'paired_ancestry':paired,'fractional_plant':plant}
    compressed=gzip.compress(canon(raw).encode(),mtime=0)
    result={'schema':'oph.common_static_source.v1','candidate':'A-static-signed-source-v1','all_pass':True,'native_calls_total':calls,
        'addressed_contexts':{'cell_compare':60,'federation_energy':1260},'source':branch_result,'pointer':pointer,'topology':topology,'moment':moment_result,'paired_ancestry':paired_result,
        'charge_trials':len(trajectories),'native_conclusive':sum(t['conclusive'] for t in trajectories),'native_timeouts':sum(not t['conclusive'] for t in trajectories),
        'native_body_records':sum(len(t['records']) for t in trajectories),'control_trials':len(controls),'control_body_records':sum(len(t['records']) for t in controls),
        'controls_never_synchronize':True,'fractional_plant':plant,'compact_provider_sha256':digest(COMPACT/'MANIFEST.sha256'),
        'same_source_compact_energy_and_regional_charge':'provider branch maps are exactly the federation_energy context; compact ray sharp(1,-1)',
        'single_energy_anchor_keV':511,'pair_energy_keV':1024.4508487023415,'physical_dictionary':'explicit proposed static signed-source interpretation; not a platform decision',
        'official_rows_solved':[],'G6_confirmed':False,'records_uncompressed_sha256':hashlib.sha256(canon(raw).encode()).hexdigest()}
    if args.record:
        (HERE/'RECORDS.json.gz').write_bytes(compressed);np.savez_compressed(HERE/'ARRAYS.npz',**arrays);(HERE/'RECEIPT.json').write_text(canon(result))
    else:
        assert result==json.loads((HERE/'RECEIPT.json').read_text())
        assert gzip.decompress((HERE/'RECORDS.json.gz').read_bytes())==canon(raw).encode()
        saved=np.load(HERE/'ARRAYS.npz');assert set(saved.files)==set(arrays)
        for k,v in arrays.items():assert np.array_equal(saved[k],v),k
    print(canon({'all_pass':True,'elapsed_seconds':time.perf_counter()-start,'native_calls':calls,'charge_trials':len(trajectories),'moment_branches':720,'old_coin_disagreements':paired_result['original_coin_disagreements']}))
    return result
if __name__=='__main__':main()
