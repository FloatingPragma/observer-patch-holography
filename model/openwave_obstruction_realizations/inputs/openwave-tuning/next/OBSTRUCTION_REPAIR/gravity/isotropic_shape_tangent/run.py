"""Source kernels and independent digital-ball shape estimates."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
from pathlib import Path
import sys,json,itertools,math,time
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import ALPHA,BETA,dagger

HERE=Path(__file__).resolve().parent
L=40;R=1;WIDTH=6.
SHIFT=[(0,0,0),(1,0,0),(2,0,0),(0,1,0),(0,0,1),(1,1,0)]


class Source:
    def __init__(self):
        n=np.arange(L)
        self.grid=np.stack(np.meshgrid(n,n,n,indexing='ij'),-1).reshape(-1,3)
        self.k=2*np.pi*(self.grid+.5)/L
        self.kin=sum(np.sin(self.k[:,j])[:,None,None]*ALPHA[j] for j in range(3))
        self.lap=2*np.sum(1-np.cos(self.k),axis=1)
        self.H=self.kin+R**3*self.lap[:,None,None]**2*BETA
        self.E,self.U=np.linalg.eigh(self.H);self.occ=self.E<0
        self.P=(self.U*self.occ[:,None,:])@dagger(self.U)
        self.kernel=np.fft.ifftn(self.P.reshape(L,L,L,4,4),axes=(0,1,2))
        kw=(self.k+np.pi)%(2*np.pi)-np.pi
        envelope=np.exp(-WIDTH**2*np.sum(kw**2,axis=1)/4)
        self.p=(np.eye(4)-self.P)*envelope[:,None,None]
        self.h=self.P*envelope[:,None,None]
        self.p/=np.sqrt(np.sum(abs(self.p)**2,axis=(0,1)))[None,None,:]
        self.h/=np.sqrt(np.sum(abs(self.h)**2,axis=(0,1)))[None,None,:]
        self.packet_energy=float(np.real(np.einsum('kas,kab,kbs->',self.p.conj(),self.H,self.p)-np.einsum('kas,kab,kbs->',self.h.conj(),self.H,self.h))/4)

    def shifted(self,shift):
        idx=(self.grid+np.array(shift))%L
        return idx[:,0]*L*L+idx[:,1]*L+idx[:,2]

    def ball(self,r2,primary=False):
        r=int(math.sqrt(r2))
        cells=np.array([c for c in itertools.product(range(-r,r+1),repeat=3) if sum(x*x for x in c)<=r2])
        N=len(cells);delta=cells[:,None,:]-cells[None,:,:];idx=delta%L
        phase=np.exp(1j*np.pi*np.sum(delta,axis=-1)/L)
        C=(self.kernel[idx[:,:,0],idx[:,:,1],idx[:,:,2]]*phase[:,:,None,None]).transpose(0,2,1,3).reshape(4*N,4*N)
        if primary:ev,vec=np.linalg.eigh(C)
        else:ev=np.linalg.eigvalsh(C)
        ec=np.clip(ev,1e-13,1-1e-13)
        S=float(-np.sum(ec*np.log(ec)+(1-ec)*np.log(1-ec)))
        out=dict(r2=r2,cells=N,effective_radius=(3*N/(4*np.pi))**(1/3),entropy=S,
                 clipped_eigenvalues=int(np.sum((ev<1e-13)|(ev>1-1e-13))))
        if not primary:return out
        K=(vec*np.log((1-ec)/ec))@vec.conj().T
        K4=K.reshape(N,4,N,4)
        del C,vec,K
        values0=K4.transpose(0,2,1,3)
        di=idx.reshape(N*N,3)
        out['kernels']=[]
        for shift in SHIFT:
            q=2*np.pi*np.array(shift)/L;ids=self.shifted(shift)
            ph=np.exp(-1j*np.pi*np.sum(delta,axis=-1)/L)*np.exp(1j*cells@q)[None,:]
            values=(values0*ph[:,:,None,None]).reshape(N*N,4,4)
            acc=np.zeros((L,L,L,4,4),complex)
            np.add.at(acc,(di[:,0],di[:,1],di[:,2]),values)
            kap=np.fft.fftn(acc,axes=(0,1,2)).reshape(L**3,4,4)
            del acc,values
            V=-.5*(self.kin+self.kin[ids])-R**3*(self.lap+self.lap[ids])[:,None,None]**2*BETA
            eb=self.E[ids];Ub=self.U[ids]
            df=self.occ[ids].astype(float)[:,:,None]-self.occ.astype(float)[:,None,:]
            de=eb[:,:,None]-self.E[:,None,:]
            divided=np.zeros_like(de);mask=df!=0;divided[mask]=df[mask]/de[mask]
            dP=Ub@(divided*(dagger(Ub)@V@self.U))@dagger(self.U)
            Z=np.einsum('kab,kba->',kap,dP)/L**3
            packet=(np.einsum('kas,kbs->kab',self.p[ids],self.p.conj())-np.einsum('kas,kbs->kab',self.h[ids],self.h.conj()))/4
            J=np.einsum('kab,kba->',kap,packet)
            T=np.einsum('kab,kba->',.5*(self.H+self.H[ids]),packet)
            Wvol=np.sum(np.exp(1j*cells@q))
            out['kernels'].append(dict(shift=list(shift),q2=float(q@q),
                 coframe_Z=[float(Z.real),float(Z.imag)],matter_J=[float(J.real),float(J.imag)],
                 energy_T=[float(T.real),float(T.imag)],volume_transform=[float(Wvol.real),float(Wvol.imag)]))
            print('kernel',r2,shift,float(Z.real),float(J.real),flush=True)
        return out


def main():
    source=Source()
    packet=dict(L=L,R=R,envelope_width=WIDTH,branches=4,each_branch_exact_particle_hole=True,
                particle_norms=np.sum(abs(source.p)**2,axis=(0,1)).tolist(),
                hole_norms=np.sum(abs(source.h)**2,axis=(0,1)).tolist(),energy=source.packet_energy,
                spinor_seed_average=True,cubic_invariant_envelope=True,lapse=1)
    (HERE/'PACKET.json').write_text(json.dumps(packet,indent=2)+'\n')
    results={}
    for r2 in (6,8,9,10,12,14,16,17,20):
        rec=source.ball(r2,r2 in (9,16))
        results[r2]=rec
        (HERE/f'BALL_r2_{r2}.json').write_text(json.dumps(rec,indent=2)+'\n')
        print('ball',r2,rec['cells'],rec['entropy'],flush=True)
    eta=json.loads((HERE.parent/'curved_ball/PLANAR_NORMALIZATION.json').read_text())['entropy_per_coordinate_face']
    assembled=[]
    for r2,windows in [(9,[(8,10),(6,12)]),(16,[(14,17),(12,20)])]:
        rec=results[r2];ell=rec['effective_radius']
        slopes={f'secant_{lo}_{hi}':(results[hi]['entropy']-results[lo]['entropy'])/(results[hi]['cells']-results[lo]['cells']) for lo,hi in windows}
        slopes['continuum_dilation_proxy_only']=rec['kernels'][0]['coframe_Z'][0]/(3*rec['cells'])
        rows=[]
        for mode in rec['kernels']:
            q=math.sqrt(mode['q2']);Z=complex(*mode['coframe_Z']);J=complex(*mode['matter_J']);T=complex(*mode['energy_T']);WV=complex(*mode['volume_transform'])
            WM=4*np.pi*quad(lambda x:x*x*(ell*ell-x*x)/(2*ell)*np.sinc(q*x/np.pi),0,ell,epsabs=1e-11)[0]
            methods={}
            for name,derivative in slopes.items():
                K=Z-3*derivative*WV
                sig=-J/K if abs(K)>1e-14 else None
                methods[name]=dict(dS_dV=derivative,fixed_volume_tangent_K=[K.real,K.imag],
                         inferred_sigma=None if sig is None else [sig.real,sig.imag],
                         area_kernel_ratio=None if q==0 else K.real/(-2*eta*q*q*WM),
                         sigma_over_energy=None if sig is None or abs(T)<1e-14 else (sig/T).real)
            rows.append(dict(shift=mode['shift'],q2=mode['q2'],J_real=J.real,T_real=T.real,
                 continuum_modular_weight=WM,bulk_ratio=J.real/(2*np.pi*WM*T.real) if abs(T.real)>1e-14 else None,
                 methods=methods,postsolve_tower_sigma_over_energy=None if q==0 else np.pi/(eta*q*q)))
        assembled.append(dict(r2=r2,cells=rec['cells'],shape_slopes=slopes,modes=rows))
    out=dict(status='SOURCE_KERNELS_WITH_DIGITAL_SHAPE_DISCRETIZATION_UNCERTAINTY',packet=packet,
             independent_shape_estimates=True,shape_spread_is_certified_bound=False,
             continuum_dilation_proxy_is_independent_shape_estimate=False,
             prepared_common_global_geometry=False,stationarity_proved=False,physical_G_used=False,
             eta_for_postsolve_comparison=eta,balls=assembled)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
