#!/usr/bin/env python3
"""One cubic common-v1 native member: field, internal bands and clock records."""
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
N=3
STEPS=32768
POPULATION_N=100_000_000
SEED=2026111432
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
       (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
       (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))


def encoded(x):
    return (json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def main():
    path=HERE/'inputs/native.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()== \
        'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('common_cubic_native',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod
    spec.loader.exec_module(mod);native=mod.native
    cells=list(itertools.product(range(N),repeat=4));ci={z:i for i,z in enumerate(cells)}
    nc=len(cells);nv=12*nc
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})
    edges=[]
    for z in cells:
        c=ci[z];edges.extend((12*c+a,12*c+b,-1) for a,b in internal)
        for axis in range(4):
            zz=list(z);zz[axis]=(zz[axis]+1)%N
            edges.extend((12*c+p,12*ci[tuple(zz)]+p,axis) for p in range(12))
    w=len(edges);assert (nc,nv,w)==(81,972,6318)
    moves=np.zeros((nv,nv),dtype=np.int64);branches=[];full_checks=0
    for eid,(a,b,axis) in enumerate(edges):
        spectator=next(i for i in range(3) if i not in (a,b))
        for slot in range(3):
            loads=[0,0,0];loads[slot]=1
            for coin in (0,1):
                out,rec=native(loads,(0,1),coin)
                built=list(loads);built[0]-=rec[1];built[1]+=rec[1]
                assert tuple(built)==out and rec[0]==coin
                target=out.index(1);branches.append([eid,slot,coin,rec[1],target])
                if slot<2 and target!=slot:
                    moves[(a,b)[slot],(a,b)[target]]+=1
                if eid<64:
                    mapping=(a,b,spectator);x=[0]*nv;x[mapping[slot]]=1
                    y,r=native(x,(a,b),coin)
                    assert r==rec and y.index(1)==mapping[target]
                    full_checks+=1
    lap=np.diag(moves.sum(axis=1))-moves
    assert np.array_equal(lap,lap.T) and np.all(np.diag(lap)==13)
    lint=lap[:12,:12].copy()-8*np.eye(12,dtype=np.int64)
    lmacro=np.zeros((nc,nc),dtype=np.int64)
    for c in range(nc):
        for d in range(nc):
            if c!=d:lmacro[c,d]=lap[12*c,12*d]
    lmacro-=np.diag(lmacro.sum(axis=1))
    assert np.array_equal(lap,np.kron(lmacro,np.eye(12,dtype=np.int64))+
                          np.kron(np.eye(nc,dtype=np.int64),lint))
    eigenvalues,_=np.linalg.eigh(lint.astype(float))
    bands=[]
    for value in eigenvalues:
        if not bands or abs(value-bands[-1][0])>1e-8:bands.append([float(value),1])
        else:bands[-1][1]+=1
    targets=(0.,5-math.sqrt(5),6.,5+math.sqrt(5))
    assert [d for _,d in bands]==[1,3,5,3]
    assert max(abs(x-y) for (x,_),y in zip(bands,targets))<1e-10
    projectors=[]
    for lam,_ in bands:
        proj=np.eye(12)
        for other,_ in bands:
            if abs(lam-other)>1e-8:proj=proj@(lint-other*np.eye(12))/(lam-other)
        assert np.max(abs(proj@proj-proj))<1e-10
        projectors.append(proj)
    assert np.max(abs(sum(projectors)-np.eye(12)))<1e-10

    rng=random.Random(SEED);tapes=[];starts=[];counts=[];clock_records=[]
    for channel in range(4):
        source=reader=rng.randrange(nv);starts.append(source)
        origin=cells[reader//12][3];ticks=0;clock=[];tape=[];occ=[0]*nv
        for step in range(1,STEPS+1):
            eid=rng.randrange(w);coin=rng.randrange(2);a,b,_=edges[eid]
            out,rec=native((int(source==a),int(source==b)),(0,1),coin)
            j=rec[1]
            if j:
                assert source==reader==(a if j==1 else b)
                source=reader=b if j==1 else a
            assert out==(int(source==a),int(source==b))
            tape.append([eid,coin,j]);occ[reader]+=1
            if cells[reader//12][3]==(origin+ticks+1)%N:
                ticks+=1;clock.append([step,ticks,reader])
        tapes.append(tape);counts.append(occ);clock_records.append(clock)
    assert sum(map(len,clock_records))>=1
    for start,tape,occ,clock in zip(starts,tapes,counts,clock_records):
        reader=start;recount=[0]*nv;rebuilt_clock=[];origin=cells[reader//12][3];ticks=0
        for step,(eid,coin,j) in enumerate(tape,1):
            a,b,_=edges[eid]
            if j:
                assert reader==(a if j==1 else b);reader=b if j==1 else a
            recount[reader]+=1
            if cells[reader//12][3]==(origin+ticks+1)%N:
                ticks+=1;rebuilt_clock.append([step,ticks,reader])
        assert recount==occ and rebuilt_clock==clock
    y=math.sqrt(nv/(4*w*STEPS))*(np.asarray(counts)-STEPS/nv)
    y=y.reshape(4,nc,12)
    band_readouts=[(y@proj).tolist() for proj in projectors]
    af=(y@np.ones(12)/math.sqrt(12)).reshape((4,N,N,N,N))
    directions=(3,0,1,2);field=[]
    for mu in range(4):
        for nu in range(mu+1,4):
            field.append((np.roll(af[nu],-1,directions[mu])-af[nu]-
                          np.roll(af[mu],-1,directions[nu])+af[mu]).ravel().tolist())

    population=[];max_eigen_residual=0.;max_cov_error=0.;max_clock_error=0.
    for (lam,multiplicity),proj in zip(bands,projectors):
        u=proj[:,0];u=u/np.linalg.norm(u)
        for spatial in ((0,0,0),(1,0,0),(1,1,0),(1,1,1)):
            qs=float(sum(4*math.sin(math.pi*j/N)**2 for j in spatial));aop=lam+qs
            if aop<1e-10:continue
            symbol=[];finite=[]
            for temporal in range(N):
                momentum=2*np.pi/N*np.asarray(spatial+(temporal,))
                v=np.kron(np.asarray([np.exp(1j*np.dot(momentum,z))/math.sqrt(nc)
                                     for z in cells]),u)
                val=float(np.vdot(v,lap@v).real)
                residual=float(np.linalg.norm(lap@v-val*v));max_eigen_residual=max(max_eigen_residual,residual)
                assert residual<1e-10
                assert abs(val-(aop+4*math.sin(math.pi*temporal/N)**2))<1e-10
                symbol.append(val)
                z=1-val/(2*w)
                correction=z*(-math.expm1(POPULATION_N*math.log(z)))/(2*w*POPULATION_N*(1-z)**2)
                finite.append(1/val-1/(4*w)-correction)
            cov=[float(sum(math.cos(2*math.pi*k*t/N)*finite[k] for k in range(N))/N) for t in (0,1)]
            green=[sum(math.cos(2*math.pi*k*t/N)/symbol[k] for k in range(N))/N for t in (0,1)]
            cov_error=max(abs((cov[t]+(1/(4*w) if t==0 else 0))/green[t]-1) for t in (0,1))
            max_cov_error=max(max_cov_error,cov_error);assert cov_error<0.001
            # Compressed exact first-hit block from the same Fourier transfer.
            q=np.array([[1-(aop+2)/(2*w),1/(2*w)],
                        [1/(2*w),1-(aop+2)/(2*w)]])
            b=np.ones(2)/(2*w)
            hit=float(np.linalg.solve(np.eye(2)-q,b)[0]);target=1/(1+aop)
            clock_error=abs(hit-target);max_clock_error=max(max_clock_error,clock_error)
            assert clock_error<1e-10 and abs(green[1]/green[0]-hit)<1e-10
            control_cov=sum(math.cos(2*math.pi*k/N)/(4*w) for k in range(N))/N
            control_phi=sum(np.exp(2j*np.pi/N*np.dot(spatial,z))
                            for z in itertools.product(range(N),repeat=3))/N**3
            control_hit=control_phi*float(np.mean(u))/float(u[0])
            assert abs(control_cov)<1e-12 and abs(control_hit)<1e-12
            assert cov[1]>1e-5 and hit>1e-3
            scalar=dict(internal_lambda=lam,multiplicity=multiplicity,spatial=list(spatial),
                transfer_eigenvalues=[1-x/(2*w) for x in symbol],
                finite_block_covariance=cov,green_covariance=green,contact=1/(4*w),
                covariance_relative_error=cov_error,first_hit_transform=hit,
                finite_cycle_target=target,infinite_slice_energy=2*math.asinh(math.sqrt(aop)/2),
                finite_cycle_clock_energy=-math.log(hit),
                changed_iid_separated_covariance=float(control_cov),
                changed_iid_nonconstant_hit_transform=[float(control_hit.real),float(control_hit.imag)])
            if lam<1e-10 and qs>0:
                k=2*np.sin(np.asarray(spatial)*np.pi/N);pt=np.eye(3)-np.outer(k,k)/qs
                scalar['magnetic_covariance_eigenvalues']=np.linalg.eigvalsh(qs*cov[1]*pt).tolist()
                ee=(cov[0]-cov[1])*np.eye(3)-np.outer(k,k)*cov[1]
                scalar['reflected_electric_eigenvalues_with_contact']=np.linalg.eigvalsh(ee).tolist()
                error=max(abs(np.asarray(scalar['reflected_electric_eigenvalues_with_contact'])-
                              np.asarray(scalar['magnetic_covariance_eigenvalues'])))/scalar['magnetic_covariance_eigenvalues'][-1]
                scalar['electric_magnetic_relative_error_including_contact']=float(error)
                assert error<0.001
                assert abs(scalar['magnetic_covariance_eigenvalues'][0])<1e-10
                assert scalar['magnetic_covariance_eigenvalues'][1]>0
            population.append(scalar)
    raw=encoded(dict(graph=dict(cells=cells,edges=edges),branches=branches,
        channel_starts=starts,tapes=tapes,occupation_counts=counts,
        internal_band_readouts=band_readouts,field_strength_readouts=field,
        same_tape_clock_records=clock_records))
    (HERE/'RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    tail=(1-1/(2*w))**STEPS
    receipt=dict(schema='oph.common_native_field.cubic.v1',candidate='MFR-common-native-field-v1',
        class_='N+',member=dict(N=N,a=1,protocol='occupation-only',vertices=nv,seams=w),
        branch_calls=len(branches),full_vector_checks=full_checks,
        actual_tape_calls=4*STEPS,actual_native_calls=len(branches)+full_checks+4*STEPS,
        native_internal_bands=bands,same_tape_completed_clock_ticks=[len(x) for x in clock_records],
        final_open_clock_intervals_retained=4,first_clock_tail_at_tape_cap=tail,
        exact_expected_clock_attempts=2*w,record_only_field_and_clock_replay=True,
        population_block_length=POPULATION_N,population_block_is_an_actual_long_tape=False,
        max_source_eigenvector_residual=max_eigen_residual,max_covariance_relative_error=max_cov_error,
        max_first_hit_population_error=max_clock_error,population=population,
        unchanged_clock_and_field_decoder_on_iid_control='nonconstant clock and separated field both zero',
        all_registered_checks='pass',official='untested',
        raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest())
    (HERE/'RECEIPT.json').write_bytes(encoded(receipt))
    print(json.dumps({k:v for k,v in receipt.items() if k!='population'},indent=2,sort_keys=True))


if __name__=='__main__':main()
