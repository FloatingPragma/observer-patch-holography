"""Connection-complete regulator, preserving the exact same flat source.

No Einstein coefficient or entropy slope is used to define the operators.
The optional CAR pullback is measured separately from the metric response.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools,math,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
spec=importlib.util.spec_from_file_location('_isotropic_shape_source',HERE.parent/'isotropic_shape_tangent/run.py')
prior=importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
from audit import BETA,dagger

STENCILS={'central2':[.5],'central4':[2/3,-1/12],
          'central6':[3/4,-3/20,1/60],'central8':[4/5,-1/5,4/105,-1/280]}

def compute(source,r2):
    L=prior.L
    r=int(math.sqrt(r2))
    cells=np.array([c for c in itertools.product(range(-r,r+1),repeat=3) if sum(x*x for x in c)<=r2])
    N=len(cells);delta=cells[:,None,:]-cells[None,:,:]
    def covariance(d):
        ix=d%L
        phase=np.exp(1j*np.pi*np.sum(d,axis=-1)/L)
        return source.kernel[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*phase[:,:,None,None]
    C=covariance(delta).transpose(0,2,1,3).reshape(4*N,4*N)
    ev,vec=np.linalg.eigh(C); ec=np.clip(ev,1e-13,1-1e-13)
    K=(vec*np.log((1-ec)/ec))@vec.conj().T
    K4=K.reshape(N,4,N,4)
    del C,vec,K
    # Aggregate by separation before applying the translation-invariant shape
    # kernel. This is exactly the same trace, with O(r^3) rather than O(N^2)
    # covariance blocks per stencil and much smaller temporary allocations.
    diameter=2*r;side=4*r+1
    separations=np.stack(np.meshgrid(*([np.arange(-diameter,diameter+1)]*3),indexing='ij'),-1)
    aggregate=np.zeros((side,side,side,4,4),complex)
    di_small=(delta+diameter).reshape(-1,3)
    np.add.at(aggregate,(di_small[:,0],di_small[:,1],di_small[:,2]),
              K4.transpose(0,2,1,3).reshape(-1,4,4))
    def separated_covariance(d):
        ix=d%L
        phase=np.exp(1j*np.pi*np.sum(d,axis=-1)/L)
        return source.kernel[ix[...,0],ix[...,1],ix[...,2]]*phase[...,None,None]
    pullbacks={}
    for name,weights in STENCILS.items():
        dC=np.zeros_like(aggregate)
        reverse=-separations
        for j in range(3):
            for m,w in enumerate(weights,1):
                shift=np.zeros(3,dtype=int);shift[j]=m
                dC-=w*((reverse[...,j]+m)[...,None,None]*separated_covariance(reverse+shift)
                       -(reverse[...,j]-m)[...,None,None]*separated_covariance(reverse-shift))
        dS=np.einsum('xyzab,xyzba->',aggregate,dC)
        pullbacks[name]=dict(entropy_inward_unit_dilation=[float(dS.real),float(dS.imag)],
                             max_mode_hop=len(weights),antihermitian_generator=True)
    idx=delta%L;di=idx.reshape(N*N,3)
    values0=K4.transpose(0,2,1,3)
    out=dict(r2=r2,cells=N,entropy=float(-np.sum(ec*np.log(ec)+(1-ec)*np.log(1-ec))),
             clipped_eigenvalues=int(np.sum((ev<1e-13)|(ev>1-1e-13))),pullbacks=pullbacks,modes=[])
    for shift in prior.SHIFT:
        q=2*np.pi*np.array(shift)/L;ids=source.shifted(shift)
        ph=np.exp(-1j*np.pi*np.sum(delta,axis=-1)/L)*np.exp(1j*cells@q)[None,:]
        values=(values0*ph[:,:,None,None]).reshape(N*N,4,4)
        acc=np.zeros((L,L,L,4,4),complex)
        np.add.at(acc,(di[:,0],di[:,1],di[:,2]),values)
        kap=np.fft.fftn(acc,axes=(0,1,2)).reshape(L**3,4,4)
        del acc,values
        A=source.kin;Ap=A[ids];la=source.lap;lap=la[ids]
        oldQ=-(la+lap)[:,None,None]*np.eye(4)
        connection=.5*(A@A+Ap@Ap)-Ap@A
        newQ=oldQ+connection
        eb=source.E[ids];Ub=source.U[ids]
        df=source.occ[ids].astype(float)[:,:,None]-source.occ.astype(float)[:,None,:]
        de=eb[:,:,None]-source.E[:,None,:]
        divided=np.zeros_like(de);mask=df!=0;divided[mask]=df[mask]/de[mask]
        responses={}
        for name,dQ in [('old_vLv',oldQ),('covariant_Dsquare',newQ)]:
            V=-.5*(A+Ap)+prior.R**3*(la+lap)[:,None,None]*(BETA@dQ)
            dP=Ub@(divided*(dagger(Ub)@V@source.U))@dagger(source.U)
            Z=np.einsum('kab,kba->',kap,dP)/L**3
            responses[name]=[float(Z.real),float(Z.imag)]
        Wv=np.sum(np.exp(1j*cells@q))
        item=dict(shift=list(shift),q2=float(q@q),responses=responses,
                  volume_transform=[float(Wv.real),float(Wv.imag)],moving_algebra={})
        for name,pullback in pullbacks.items():
            shape=complex(*pullback['entropy_inward_unit_dilation'])*Wv/N
            item['moving_algebra'][name]={key:[float((complex(*z)+shape).real),float((complex(*z)+shape).imag)] for key,z in responses.items()}
        out['modes'].append(item)
        print('kernel',r2,shift,responses,flush=True)
    return out

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--L',type=int,default=40)
    parser.add_argument('--R',type=float,default=1.)
    parser.add_argument('--r2',type=float,nargs='+',default=[9,16])
    parser.add_argument('--output',default=str(HERE))
    parser.add_argument('--axis-only',action='store_true')
    parser.add_argument('--order',type=int,choices=[2,4],default=2)
    args=parser.parse_args()
    output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    prior.L=args.L;prior.R=args.R;prior.WIDTH=6*args.R
    if args.axis_only:prior.SHIFT=[(0,0,0),(1,0,0),(2,0,0)]
    if args.order==4:
        from improved_source import ImprovedSource
        source=ImprovedSource(prior.L,prior.R)
    else:
        source=prior.Source()
        # These large arrays do not enter the geometric kernel calculation.
        del source.P,source.H,source.p,source.h
    records=[]
    for r2 in args.r2:
        rec=compute(source,r2);records.append(rec)
        (output/f'BALL_r2_{r2}.json').write_text(json.dumps(rec,indent=2)+'\n')
    result=dict(status='CONNECTION_COMPLETION_AND_DETERMINISTIC_CAR_PULLBACK_COMPONENT_TEST',
                source_L=prior.L,source_R=prior.R,lapse=1,flat_source_unchanged=True,
                source_difference_order=args.order,
                proper_volume_constraint_used_only_geometrically=True,
                moving_CAR_algebra_is_proved_literal_geometric_ball=False,
                no_entropy_fit=True,balls=records)
    if args.order==4:
        result['flat_source_unchanged']=False
        result['source_change']='Standard fourth-order differences; new eta and modular checks required'
        result['positive_W_identity_error']=source.polynomial_identity_error
        result['minimum_W']=source.minimum_W
    (output/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
