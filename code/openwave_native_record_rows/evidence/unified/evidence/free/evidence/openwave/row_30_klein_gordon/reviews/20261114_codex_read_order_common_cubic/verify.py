#!/usr/bin/env python3
"""Independent source/common-tape reader, no import of producer or native."""
import gzip,hashlib,itertools,json,math,random
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent.parent/'work/2026-11-14_codex_common_cubic'

def main():
    raw=gzip.decompress((PACKET/'RECORDS.json.gz').read_bytes())
    r=json.loads(raw);receipt=json.loads((PACKET/'RECEIPT.json').read_text())
    assert hashlib.sha256(raw).hexdigest()==receipt['raw_uncompressed_sha256']
    reg=json.loads((PACKET/'PREREG.json').read_text())
    for name,h in reg['files'].items():assert hashlib.sha256((PACKET/name).read_bytes()).hexdigest()==h
    N=3;nv=972;W=6318;n=32768;cells=list(itertools.product(range(N),repeat=4));ci={z:i for i,z in enumerate(cells)}
    faces=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
           (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in faces for j in range(3)})
    edges=[]
    for z in cells:
        c=ci[z];edges += [(12*c+a,12*c+b,-1) for a,b in internal]
        for axis in range(4):
            zz=list(z);zz[axis]=(zz[axis]+1)%N
            edges += [(12*c+p,12*ci[tuple(zz)]+p,axis) for p in range(12)]
    assert len(edges)==W and r['graph']=={'cells':[list(z) for z in cells],'edges':[list(e) for e in edges]}
    branch_iter=iter(r['branches']);lap=np.zeros((nv,nv),dtype=np.int64)
    for eid,(a,b,axis) in enumerate(edges):
        lap[a,a]+=1;lap[b,b]+=1;lap[a,b]-=1;lap[b,a]-=1
        for slot in range(3):
            for coin in (0,1):
                j=(int(slot==0)-coin) if slot<2 else 0
                terminal=[int(i==slot) for i in range(3)];terminal[0]-=j;terminal[1]+=j
                assert next(branch_iter)==[eid,slot,coin,j,terminal.index(1)]
    assert next(branch_iter,None) is None
    lint=np.zeros((12,12),dtype=np.int64)
    for a,b in internal:lint[a,a]+=1;lint[b,b]+=1;lint[a,b]-=1;lint[b,a]-=1
    lmacro=np.zeros((81,81),dtype=np.int64)
    for z in cells:
        x=ci[z]
        for axis in range(4):
            zz=list(z);zz[axis]=(zz[axis]+1)%N;y=ci[tuple(zz)]
            lmacro[x,x]+=1;lmacro[y,y]+=1;lmacro[x,y]-=1;lmacro[y,x]-=1
    assert np.array_equal(lap,np.kron(lmacro,np.eye(12,dtype=np.int64))+np.kron(np.eye(81,dtype=np.int64),lint))
    assert np.all(np.diag(lap)==13)
    eigenvalues=np.linalg.eigvalsh(lint.astype(float));targets=np.array([0]+[5-math.sqrt(5)]*3+[6]*5+[5+math.sqrt(5)]*3)
    assert max(abs(eigenvalues-targets))<1e-12
    rng=random.Random(2026111432);counts=[];ticknums=[];censored=[]
    for channel in range(4):
        reader=rng.randrange(nv);assert reader==r['channel_starts'][channel]
        origin=cells[reader//12][3];occ=[0]*nv;ticks=[];count=0
        assert len(r['tapes'][channel])==n
        for step,record in enumerate(r['tapes'][channel],1):
            eid=rng.randrange(W);coin=rng.randrange(2);a,b,_=edges[eid]
            j=int(reader==a)-coin if reader in (a,b) else 0
            assert record==[eid,coin,j]
            if j:reader=b if j==1 else a
            occ[reader]+=1
            if cells[reader//12][3]==(origin+count+1)%3:count+=1;ticks.append([step,count,reader])
        assert occ==r['occupation_counts'][channel] and ticks==r['same_tape_clock_records'][channel]
        counts.append(occ);ticknums.append(count);censored.append(n-(ticks[-1][0] if ticks else 0))
    assert ticknums==[4,5,3,0]
    y=(math.sqrt(nv/(4*W*n))*(np.array(counts)-n/nv)).reshape(4,81,12)
    actual_lams=[x[0] for x in receipt['native_internal_bands']]
    bands=[]
    for lam in actual_lams:
        proj=np.eye(12)
        for other in actual_lams:
            if abs(lam-other)>1e-8:proj=proj@(lint-other*np.eye(12))/(lam-other)
        bands.append(proj)
    max_band=max(float(np.max(abs(y@p-np.array(saved)))) for p,saved in zip(bands,r['internal_band_readouts']))
    assert max_band<1e-14
    af=(y@np.ones(12)/math.sqrt(12)).reshape(4,3,3,3,3);dirs=(3,0,1,2);fields=[]
    for mu in range(4):
        for nu in range(mu+1,4):fields.append((np.roll(af[nu],-1,dirs[mu])-af[nu]-np.roll(af[mu],-1,dirs[nu])+af[mu]).ravel())
    max_field=float(np.max(abs(np.array(fields)-r['field_strength_readouts'])));assert max_field<1e-14
    covariance_residual=0.;clock_residual=0.;negative_electric=[]
    for case in receipt['population']:
        lam=case['internal_lambda'];qs=sum(4*math.sin(math.pi*j/3)**2 for j in case['spatial']);A=lam+qs
        vals=[A,A+3,A+3];population_n=receipt['population_block_length'];c=[]
        for val in vals:
            z=1-val/(2*W);c.append(1/val-1/(4*W)-z*(1-z**population_n)/(2*W*population_n*(1-z)**2))
        cov=[sum(c)/3,(c[0]-.5*c[1]-.5*c[2])/3]
        covariance_residual=max(covariance_residual,max(abs(x-y) for x,y in zip(cov,case['finite_block_covariance'])))
        # The two non-target planes have the same boundary value, so (A+1)h=1.
        clock_residual=max(clock_residual,abs((A+1)*case['first_hit_transform']-1))
        assert abs(case['finite_cycle_clock_energy']-math.log1p(A))<1e-12
        assert abs(case['changed_iid_separated_covariance'])<1e-12
        assert max(abs(x) for x in case['changed_iid_nonconstant_hit_transform'])<1e-12
        if 'reflected_electric_eigenvalues_with_contact' in case:
            negative_electric.append(case['reflected_electric_eigenvalues_with_contact'][0])
            assert case['reflected_electric_eigenvalues_with_contact'][0]<0
            assert case['electric_magnetic_relative_error_including_contact']<.001
    assert covariance_residual<1e-12 and clock_residual<1e-12
    tail=(1-1/(2*W))**n
    assert receipt['exact_expected_clock_attempts']==2*W and abs(tail-receipt['first_clock_tail_at_tape_cap'])<1e-15
    assert receipt['actual_native_calls']==len(r['branches'])+64*6+4*n==169364
    out={'pass':True,'review':'same-group; no G6','source_branch_records_checked':len(r['branches']),
         'actual_tape_records_checked':4*n,'all_rng_choices_reproduced':True,'integer_kronecker_identity':True,
         'ticks':ticknums,'final_open_interval_lengths':censored,'maximum_band_replay_error':max_band,
         'maximum_curl_replay_error':max_field,'maximum_population_covariance_error':covariance_residual,
         'maximum_first_hit_equation_error':clock_residual,'negative_longitudinal_electric_contact_values':negative_electric,
         'exact_mean_clock_attempts':2*W,'first_tick_tail_at_cap':tail,
         'actual_long_population_tape_executed':False,'producer_full_vector_calls_reviewed_as_code_not_reexecuted':384,
         'packet_receipt_sha256':hashlib.sha256((PACKET/'RECEIPT.json').read_bytes()).hexdigest(),
         'packet_records_sha256':hashlib.sha256((PACKET/'RECORDS.json.gz').read_bytes()).hexdigest()}
    (HERE/'VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
