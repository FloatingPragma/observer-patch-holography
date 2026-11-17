"""Actual curved-source Gaussian entropy response, own Bloch implementation."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
from pathlib import Path
import sys,json,itertools,math
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import ALPHA,BETA,dagger

HERE=Path(__file__).resolve().parent


def run(L,radius,R=1,exact_radial=False):
    ng=np.arange(L)
    grid=np.stack(np.meshgrid(ng,ng,ng,indexing='ij'),-1).reshape(-1,3)
    k=2*np.pi*(grid+.5)/L
    kin=sum(np.sin(k[:,j])[:,None,None]*ALPHA[j] for j in range(3))
    lap=2*np.sum(1-np.cos(k),axis=1)
    H=kin+R**3*lap[:,None,None]**2*BETA
    E,U=np.linalg.eigh(H);occ=E<0
    P=(U*occ[:,None,:])@dagger(U)
    ker=np.fft.ifftn(P.reshape(L,L,L,4,4),axes=(0,1,2))
    cells=np.array([c for c in itertools.product(range(-radius,radius+1),repeat=3) if sum(x*x for x in c)<=radius**2])
    N=len(cells)
    delta=cells[:,None,:]-cells[None,:,:]
    assert abs(delta).max()<L
    ix=delta%L
    phase=np.exp(1j*np.pi*np.sum(delta,axis=-1)/L)
    C=(ker[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*phase[:,:,None,None]).transpose(0,2,1,3).reshape(4*N,4*N)
    ev,vec=np.linalg.eigh(C)
    floor=1e-13
    ec=np.clip(ev,floor,1-floor)
    entropy=float(-np.sum(ec*np.log(ec)+(1-ec)*np.log(1-ec)))
    K=(vec*np.log((1-ec)/ec))@vec.conj().T
    K4=K.reshape(N,4,N,4)

    def kappa(shift):
        q=2*np.pi*np.array(shift)/L
        ph=np.exp(-1j*np.pi*np.sum(delta,axis=-1)/L)*np.exp(1j*cells@q)[None,:]
        values=(K4.transpose(0,2,1,3)*ph[:,:,None,None]).reshape(N*N,4,4)
        acc=np.zeros((L,L,L,4,4),complex)
        di=ix.reshape(N*N,3)
        np.add.at(acc,(di[:,0],di[:,1],di[:,2]),values)
        return np.fft.fftn(acc,axes=(0,1,2)).reshape(L**3,4,4)

    def response(shift,full_scale=False):
        idx=(grid+np.array(shift))%L
        ids=idx[:,0]*L*L+idx[:,1]*L+idx[:,2]
        if full_scale:
            V=.5*(H+H[ids])
        else:
            V=-.5*(kin+kin[ids])-R**3*(lap+lap[ids])[:,None,None]**2*BETA
        Ub=U[ids];eb=E[ids];fb=occ[ids].astype(float)
        df=fb[:,:,None]-occ.astype(float)[:,None,:]
        de=eb[:,:,None]-E[:,None,:]
        F=np.zeros_like(de);msk=df!=0;F[msk]=df[msk]/de[msk]
        Vb=dagger(Ub)@V@U
        dP=Ub@(F*Vb)@dagger(U)
        return complex(np.einsum('kab,kba->',kappa(shift),dP)/len(k))

    z0=response((0,0,0))
    zm=[response(tuple(int(j==a) for j in range(3))) for a in range(3)]
    q=2*np.pi/L
    cosmean=np.mean(np.cos(q*cells),axis=0)
    profile=np.sum(cosmean[None,:]-np.cos(q*cells),axis=1)
    dS=float((cosmean.sum()*z0-sum(zm)).real)
    ell=(3*N/(4*np.pi))**(1/3)
    eta=.06779595762226844 # independent parent R1 planar diagnostic; finite quadrature, not fit to this ball
    z=q*ell
    fvol=3*(math.sin(z)-z*math.cos(z))/z**3
    fsurf=math.sin(z)/z
    geom=8*np.pi*ell**2*3*(fvol-fsurf)
    geom_discrete_volume=8*np.pi*ell**2*(cosmean.sum()-3*fsurf)
    # Complete pair, lowest positive and negative energy, from each nearest
    # zero momentum point. kappa is the untranslated plane-wave regional sum.
    kap0=kappa((0,0,0))
    tau=np.real(np.diagonal(dagger(U)@kap0@U,axis1=1,axis2=2))
    near=np.all((grid==0)|(grid==L-1),axis=1)
    bulk=(tau[:,2]-tau[:,1])/(E[:,2]-E[:,1])
    matter=float(np.mean(bulk[near]))
    expected=2*np.pi*4*np.pi*ell**4/15
    out=dict(L=L,R=R,radius=radius,cells=N,effective_radius=ell,entropy=entropy,
             unresolved_eigenvalues=int(np.sum((ev<floor)|(ev>1-floor))),eigenvalue_floor=floor,
             uniform_coframe_response=[z0.real,z0.imag],axis_responses=[[x.real,x.imag] for x in zm],
             proper_volume_first_derivative=float(3*np.sum(profile)),
             curved_entropy_response=dS,planar_eta=eta,
             smooth_sphere_area_response=geom,discrete_volume_sphere_area_response=geom_discrete_volume,
             smooth_ratio=dS/(eta*geom),discrete_volume_ratio=dS/(eta*geom_discrete_volume),
             same_ball_bulk_coefficient=matter,continuum_bulk_coefficient=expected,
             bulk_ratio=matter/expected,uniform_full_H_scaling_control=abs(response((0,0,0),True)),
             full_Einstein_criterion=False)
    if exact_radial:
        # Added after cosine diagnostics, to eliminate their local q^4
        # angular distortion. Geometry is exactly radial near the region;
        # the signed-coordinate quadratic is periodically continued far away.
        xs=np.where(np.arange(L)<=L//2,np.arange(L),np.arange(L)-L)
        coeff=np.fft.fft(xs.astype(float)**2)/L
        mean_r2=float(np.mean(np.sum(cells**2,axis=1)))
        radial=(3*coeff[0].real-mean_r2)*z0.real
        rows=[]
        for n in range(1,L//2+1):
            zn=zm[0] if n==1 else response((n,0,0))
            multiplier=1 if n==L//2 else 2
            radial+=3*multiplier*(coeff[n]*zn).real
            rows.append(dict(n=n,coefficient=[float(coeff[n].real),float(coeff[n].imag)],response=[zn.real,zn.imag]))
        out['exact_radial']={
            'introduced_after_cosine_diagnostics':True,
            'profile':'periodic signed-coordinate r^2 minus regional mean r^2',
            'axis_symmetry_defect':max(abs(x-zm[0]) for x in zm),
            'regional_mean_r2':mean_r2,
            'proper_volume_first_derivative':float(3*np.sum(np.sum(cells**2,axis=1)-mean_r2)),
            'entropy_response':float(radial),
            'smooth_fixed_volume_area_prediction':float(eta*16*np.pi*ell**4/5),
            'discrete_volume_area_prediction':float(eta*8*np.pi*ell**2*(ell**2-mean_r2)),
            'smooth_ratio':float(radial/(eta*16*np.pi*ell**4/5)),
            'discrete_volume_ratio':float(radial/(eta*8*np.pi*ell**2*(ell**2-mean_r2))),
            'modes':rows,
        }
    suffix='_radial' if exact_radial else ''
    (HERE/f'RESULT_L{L}_r{radius}_R{R}{suffix}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
    return out


if __name__=='__main__':run(int(sys.argv[1]),int(sys.argv[2]),exact_radial='--radial' in sys.argv)
