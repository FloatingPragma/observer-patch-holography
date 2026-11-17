"""Exact pure four-particle/four-hole Slater preparation and regional entropy."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools,math,argparse
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--L',type=int,default=60)
    parser.add_argument('--R',type=float,default=1.5)
    parser.add_argument('--order',type=int,choices=[2,4],default=4)
    parser.add_argument('--widths',type=float,nargs='+',default=[9,12,15])
    parser.add_argument('--r2',type=float,nargs='+',default=[16,20.25])
    args=parser.parse_args();L=args.L;N=L**3
    n=np.arange(L);grid=np.stack(np.meshgrid(n,n,n,indexing='ij'),-1).reshape(-1,3)
    k=2*np.pi*(grid+.5)/L;t=np.sin(k/2)**2
    if args.order==4:d=np.sin(k)*(1+2*t/3);lap=np.sum(4*t*(1+t/3),axis=1)
    else:d=np.sin(k);lap=np.sum(4*t,axis=1)
    mass=args.R**3*lap**2;energy=np.sqrt(np.sum(d*d,axis=1)+mass*mass)
    coeff=np.column_stack([d,mass])/energy[:,None]
    gammas=np.array([*ALPHA,BETA]);identity=np.eye(4)
    spectral_identity_error=float(np.max(abs(np.sum(coeff*coeff,axis=1)-1)))
    kernels=[np.fft.ifftn(coeff[:,a].reshape(L,L,L)) for a in range(4)]
    maxr=math.sqrt(max(args.r2));r=int(maxr)
    cells=np.array([x for x in itertools.product(range(-r,r+1),repeat=3) if sum(z*z for z in x)<=max(args.r2)])
    ix=cells%L;phase=np.exp(1j*np.pi*np.sum(cells,axis=1)/L)
    cell_index={tuple(x):i for i,x in enumerate(cells)}
    rotations=[]
    for name,Rmat,Sgen in [('z90',np.array([[0,-1,0],[1,0,0],[0,0,1]]),-1j*ALPHA[0]@ALPHA[1]),
                           ('x90',np.array([[1,0,0],[0,0,-1],[0,1,0]]),-1j*ALPHA[1]@ALPHA[2])]:
        permutations=np.array([cell_index[tuple(Rmat@x)] for x in cells])
        rotations.append((name,permutations,(identity-1j*Sgen)/np.sqrt(2)))
    kw=(k+np.pi)%(2*np.pi)-np.pi
    prepared=[]
    for width in args.widths:
        f=np.exp(-width**2*np.sum(kw*kw,axis=1)/4);f2=f*f
        M0=np.sum(f2)*identity
        MG=np.einsum('k,ka,aij->ij',f2,coeff,gammas)
        E0=np.sum(f2*energy)*identity
        EG=np.einsum('k,ka,aij->ij',f2*energy,coeff,gammas)
        modes={};grams={};ex_energy=0.
        for name,sign in [('particle',1),('hole',-1)]:
            M=(M0+sign*MG)/2
            ev,u=np.linalg.eigh(M);invsqrt=(u/np.sqrt(ev))@dagger(u)
            gram=invsqrt@M@invsqrt
            grams[name]=dict(error=float(np.max(abs(gram-identity))),min_mode_norm=float(ev.min()))
            ex_energy+=float(np.trace((invsqrt@invsqrt)@((E0+sign*EG)/2)).real)
            regional=[];regional_H=[]
            for s in range(4):
                col=invsqrt[:,s]
                wave=f[:,None]*(col[None,:]+sign*np.einsum('ka,aij,j->ki',coeff,gammas,col))/2
                field=np.fft.ifftn(wave.reshape(L,L,L,4),axes=(0,1,2))*math.sqrt(N)
                regional.append(field[ix[:,0],ix[:,1],ix[:,2]]*phase[:,None])
                del field
                fieldH=np.fft.ifftn((sign*energy[:,None]*wave).reshape(L,L,L,4),axes=(0,1,2))*math.sqrt(N)
                regional_H.append(fieldH[ix[:,0],ix[:,1],ix[:,2]]*phase[:,None])
                del fieldH,wave
            modes[name]=np.stack(regional,axis=-1).reshape(4*len(cells),4)
            modes[name+'_H']=np.stack(regional_H,axis=-1).reshape(4*len(cells),4)
        localE=np.real(np.sum(modes['particle'].conj()*modes['particle_H']-modes['hole'].conj()*modes['hole_H'],axis=1)).reshape(len(cells),4).sum(1)
        rotation_errors={}
        for name,permutation,spin in rotations:
            for sector in ('particle','hole'):
                X=modes[sector].reshape(len(cells),4,4)
                expected=spin@X@dagger(spin)
                rotation_errors[name+'_'+sector]=float(np.max(abs(X[permutation]-expected)))
        prepared.append(dict(width=width,energy=ex_energy,gram_checks=grams,
                             cubic_rotation_mode_errors=rotation_errors,modes=modes,local_energy=localE))
        print('prepared',width,ex_energy,grams,flush=True)
    results=[]
    for radius2 in args.r2:
        keep=np.sum(cells*cells,axis=1)<=radius2;bc=cells[keep];nb=len(bc)
        select=np.flatnonzero(np.repeat(keep,4))
        delta=bc[:,None,:]-bc[None,:,:];idx=delta%L
        ph=np.exp(1j*np.pi*np.sum(delta,axis=-1)/L)
        Ct=np.eye(nb)[:,:,None,None]*identity/2
        Ct=Ct.astype(complex)
        for a in range(4):Ct-=.5*kernels[a][idx[:,:,0],idx[:,:,1],idx[:,:,2]][:,:,None,None]*ph[:,:,None,None]*gammas[a]
        C=Ct.transpose(0,2,1,3).reshape(4*nb,4*nb);del Ct
        e,u=np.linalg.eigh(C);ec=np.clip(e,1e-13,1-1e-13)
        S0=float(-np.sum(ec*np.log(ec)+(1-ec)*np.log(1-ec)))
        records=[]
        for prep in prepared:
            Xp=prep['modes']['particle'][select];Xh=prep['modes']['hole'][select]
            change=Xp@dagger(Xp)-Xh@dagger(Xh)
            modular_diagonal=(np.sum(abs(dagger(u)@Xp)**2,axis=1)-np.sum(abs(dagger(u)@Xh)**2,axis=1))
            first=float(np.sum(np.log((1-ec)/ec)*modular_diagonal))
            ev=np.linalg.eigvalsh(C+change);evc=np.clip(ev,1e-13,1-1e-13)
            S=float(-np.sum(evc*np.log(evc)+(1-evc)*np.log(1-evc)))
            ell=math.sqrt(radius2);weight=(radius2-np.sum(bc*bc,axis=1))/(2*ell)
            geometric=2*np.pi*np.sum(weight*prep['local_energy'][keep])
            floor_checks=[]
            for floor in (1e-11,1e-13,1e-15):
                c0=np.clip(e,floor,1-floor);c1=np.clip(ev,floor,1-floor)
                s0=-np.sum(c0*np.log(c0)+(1-c0)*np.log(1-c0))
                s1=-np.sum(c1*np.log(c1)+(1-c1)*np.log(1-c1))
                floor_checks.append(dict(floor=floor,entropy_change=float(s1-s0),
                                         first_derivative=float(np.sum(np.log((1-c0)/c0)*modular_diagonal))))
            rec=dict(width=prep['width'],actual_total_energy=prep['energy'],entropy_change=S-S0,
                     exact_first_derivative=first,relative_entropy_remainder=first-(S-S0),
                     relative_first_law_error=(first-S+S0)/first,
                     continuum_ball_modular_energy=float(geometric),
                     modular_ratio=float(first/geometric),
                     min_covariance_eigenvalue=float(ev.min()),max_covariance_eigenvalue=float(ev.max()),
                     local_covariance_change_norm=float(np.linalg.norm(change,ord='fro')),
                     floor_checks=floor_checks)
            records.append(rec);print('ball',radius2,rec,flush=True)
        results.append(dict(radius_squared=radius2,cells=nb,vacuum_entropy=S0,records=records))
    serial_preps=[{key:value for key,value in p.items() if key not in ('modes','local_energy')} for p in prepared]
    out=dict(status='EXACT_PURE_NEUTRAL_PREPARATION_COMPONENT_TEST',L=L,R=args.R,source_difference_order=args.order,
             particles=4,holes=4,exact_charge_every_state=True,globally_gaussian=True,
             covariance_projector_proof='P−XhXh†+XpXp†; normalized orthogonal filled/empty spectral subspaces',
             spectral_projector_identity_error=spectral_identity_error,preparations=serial_preps,balls=results,
             geometry_prepared=False,stationarity_tested=False,physical_G_used=False,
             source_amplitude_is_arbitrary_fixed_profile_epsilon=False)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
