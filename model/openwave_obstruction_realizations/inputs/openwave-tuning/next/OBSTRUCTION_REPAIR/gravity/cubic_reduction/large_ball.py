"""Complete G1 matter response in large cubic balls, without a full ball matrix."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools,math,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger
from group import orbit_bases

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,default=96)
    ap.add_argument('--R',type=float,default=1.5)
    ap.add_argument('--radii',type=int,nargs='+',default=[6,9])
    ap.add_argument('--widths',type=float,nargs='+',default=[27,36])
    ap.add_argument('--output',default=str(HERE/'large_R1p5'))
    args=ap.parse_args();outdir=Path(args.output);outdir.mkdir(parents=True,exist_ok=True)
    L=args.L;N=L**3;started=time.process_time()
    grid=np.stack(np.meshgrid(*([np.arange(L)]*3),indexing='ij'),-1).reshape(-1,3)
    k=2*np.pi*(grid+.5)/L;t=np.sin(k/2)**2
    d=np.sin(k)*(1+2*t/3);lap=np.sum(4*t*(1+t/3),axis=1)
    mass=args.R**3*lap**2;E=np.sqrt(np.sum(d*d,axis=1)+mass*mass)
    coeff=np.column_stack([d,mass])/E[:,None]
    gammas=np.array([*ALPHA,BETA]);I4=np.eye(4)
    kernels=[np.fft.ifftn(coeff[:,a].reshape(L,L,L)) for a in range(4)]
    maxr=max(args.radii)
    cells=np.array([x for x in itertools.product(range(-maxr,maxr+1),repeat=3) if sum(z*z for z in x)<=maxr**2])
    orbits=orbit_bases(cells);at=0
    for o in orbits:o['columns']=np.arange(at,at+o['B1'].shape[1]);at+=o['B1'].shape[1]
    size=at;C=np.zeros((size,size),complex)
    metadata=dict(L=L,R=args.R,radii=args.radii,widths=args.widths,full_largest_ball_modes=4*len(cells),
                  reduced_G1_multiplicity=size,G1_dimension=2,orbits=len(orbits),
                  parity_copies_retained=True,largest_dense_matrix_bytes=int(C.nbytes),
                  source_difference_order=4,geometry_response=False)
    (outdir/'METADATA.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('metadata',metadata,flush=True)
    for i,oi in enumerate(orbits):
        for j in range(i,len(orbits)):
            oj=orbits[j];delta=oi['cells'][:,None,:]-oj['cells'][None,:,:];ix=delta%L
            phase=np.exp(1j*np.pi*np.sum(delta,axis=-1)/L)
            Ct=np.zeros((len(oi['cells']),len(oj['cells']),4,4),complex)
            if i==j:Ct+=np.eye(len(oi['cells']))[:,:,None,None]*I4/2
            for a in range(4):Ct-=.5*kernels[a][ix[:,:,0],ix[:,:,1],ix[:,:,2]][:,:,None,None]*phase[:,:,None,None]*gammas[a]
            full=Ct.transpose(0,2,1,3).reshape(4*len(oi['cells']),4*len(oj['cells']))
            block=dagger(oi['B1'])@full@oj['B1']
            C[np.ix_(oi['columns'],oj['columns'])]=block
            if j!=i:C[np.ix_(oj['columns'],oi['columns'])]=dagger(block)
        if i%20==0:print('covariance row',i,'of',len(orbits),'cpu',time.process_time()-started,flush=True)
    np.save(outdir/'C_G1.npy',C)
    print('covariance complete cpu',time.process_time()-started,flush=True)
    ix=cells%L;phase=np.exp(1j*np.pi*np.sum(cells,axis=1)/L)
    kw=(k+np.pi)%(2*np.pi)-np.pi;prepared=[]
    for width in args.widths:
        f=np.exp(-width**2*np.sum(kw*kw,axis=1)/4);f2=f*f
        M0=np.sum(f2)*I4;MG=np.einsum('k,ka,aij->ij',f2,coeff,gammas)
        E0=np.sum(f2*E)*I4;EG=np.einsum('k,ka,aij->ij',f2*E,coeff,gammas)
        modes={};energy=0.;grams={}
        for name,sign in [('particle',1),('hole',-1)]:
            M=(M0+sign*MG)/2;ev,u=np.linalg.eigh(M);inv=(u/np.sqrt(ev))@dagger(u)
            grams[name]=float(np.max(abs(inv@M@inv-I4)))
            energy+=float(np.trace((inv@inv)@((E0+sign*EG)/2)).real)
            local=[];localH=[]
            for s in range(4):
                col=inv[:,s]
                wave=f[:,None]*(col[None,:]+sign*np.einsum('ka,aij,j->ki',coeff,gammas,col))/2
                field=np.fft.ifftn(wave.reshape(L,L,L,4),axes=(0,1,2))*math.sqrt(N)
                local.append(field[ix[:,0],ix[:,1],ix[:,2]]*phase[:,None]);del field
                field=np.fft.ifftn((sign*E[:,None]*wave).reshape(L,L,L,4),axes=(0,1,2))*math.sqrt(N)
                localH.append(field[ix[:,0],ix[:,1],ix[:,2]]*phase[:,None]);del field,wave
            X=np.stack(local,axis=-1).reshape(4*len(cells),4)
            XH=np.stack(localH,axis=-1).reshape(4*len(cells),4)
            reduced=np.zeros((size,4),complex)
            for o in orbits:reduced[o['columns']]=dagger(o['B1'])@X[o['rows']]
            modes[name]=X;modes[name+'_H']=XH;modes[name+'_reduced']=reduced
        localE=np.real(np.sum(modes['particle'].conj()*modes['particle_H']-modes['hole'].conj()*modes['hole_H'],axis=1)).reshape(len(cells),4).sum(1)
        prepared.append(dict(width=width,energy=energy,gram_errors=grams,modes=modes,local_energy=localE))
        print('prepared',width,energy,'cpu',time.process_time()-started,flush=True)
    records=[]
    for radius in args.radii:
        selected_orbits=[o for o in orbits if np.sum(o['cells'][0]**2)<=radius**2]
        columns=np.concatenate([o['columns'] for o in selected_orbits])
        keep=np.sum(cells*cells,axis=1)<=radius**2
        c=C[np.ix_(columns,columns)];ev,u=np.linalg.eigh(c);ec=np.clip(ev,1e-13,1-1e-13)
        S0=-np.sum(ec*np.log(ec)+(1-ec)*np.log(1-ec))
        rows=[]
        for p in prepared:
            Xp=p['modes']['particle_reduced'][columns];Xh=p['modes']['hole_reduced'][columns]
            delta=Xp@dagger(Xp)-Xh@dagger(Xh)
            diag=np.sum(abs(dagger(u)@Xp)**2,axis=1)-np.sum(abs(dagger(u)@Xh)**2,axis=1)
            first=2*float(np.log((1-ec)/ec)@diag)
            ee=np.linalg.eigvalsh(c+delta);eec=np.clip(ee,1e-13,1-1e-13)
            Sf=-np.sum(eec*np.log(eec)+(1-eec)*np.log(1-eec));change=2*float(Sf-S0)
            weight=(radius**2-np.sum(cells[keep]**2,axis=1))/(2*radius)
            bulk=2*np.pi*np.sum(weight*p['local_energy'][keep])
            floors=[]
            for floor in (1e-11,1e-13,1e-15):
                c0=np.clip(ev,floor,1-floor);c1=np.clip(ee,floor,1-floor)
                s0=-np.sum(c0*np.log(c0)+(1-c0)*np.log(1-c0));s1=-np.sum(c1*np.log(c1)+(1-c1)*np.log(1-c1))
                floors.append(dict(floor=floor,entropy_change=2*float(s1-s0),first_derivative=2*float(np.log((1-c0)/c0)@diag)))
            row=dict(width=p['width'],actual_total_energy=p['energy'],entropy_change=change,first_derivative=first,
                     relative_first_law_remainder=(first-change)/first,continuum_ball_modular_energy=float(bulk),
                     modular_ratio=first/bulk,finite_entropy_ratio=change/bulk,floor_checks=floors)
            rows.append(row);print('ball',radius,row,flush=True)
        rec=dict(radius=radius,ell_over_R=radius/args.R,cells=int(keep.sum()),G1_multiplicity=len(columns),records=rows)
        records.append(rec);(outdir/f'BALL_r{radius}.json').write_text(json.dumps(rec,indent=2)+'\n')
    out=dict(status='LARGE_BALL_EXACT_G1_MATTER_COMPONENT',metadata=metadata,balls=records,
             source_preparations=[{k:v for k,v in p.items() if k not in ('modes','local_energy')} for p in prepared],
             all_unchanged_sectors_cancel_by_validated_symmetry=True,geometry_response_reduced=False,
             global_periodic_equilibrium=False,lapse_response=False,cpu_seconds=time.process_time()-started)
    (outdir/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
