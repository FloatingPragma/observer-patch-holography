"""One global coframe from actual entropy response; exact neutral mixtures."""
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
L=3;NC=L**3;NM=4*NC
coords=np.array(list(itertools.product(range(L),repeat=3)))
k=2*np.pi*(coords+.5)/L
F=np.exp(1j*coords@k.T)/math.sqrt(NC)
lap=2*np.sum(1-np.cos(k),axis=1)
Lap=(F*lap)@F.conj().T
Lap2=Lap@Lap
Kmom=sum(np.sin(k[:,j])[:,None,None]*ALPHA[j] for j in range(3))
KIN=np.einsum('xk,kab,yk->xayb',F,Kmom,F.conj()).reshape(NM,NM)


def source(sigma):
    v=np.exp(-sigma);vv=np.repeat(v,4)
    b=v[:,None]*Lap*v[None,:]
    return .5*(vv[:,None]*KIN+KIN*vv[None,:])+np.kron(b@b,BETA)


def vacuum_and_pair(H):
    e,u=np.linalg.eigh(H);P=(u*(e<0))@u.conj().T
    basis_p=np.zeros(NM);basis_p[0]=1
    basis_h=np.zeros(NM);basis_h[1]=1
    p=(np.eye(NM)-P)@basis_p;p/=np.linalg.norm(p)
    h=P@basis_h;h/=np.linalg.norm(h)
    Q=P+np.outer(p,p.conj())-np.outer(h,h.conj())
    return P,Q,float(np.vdot(p,H@p).real-np.vdot(h,H@h).real),e,u


def region(cells):return (4*np.array(cells)[:,None]+np.arange(4)).ravel()
def neighbor(c,a,sgn=1):
    x=coords[c].copy();x[a]=(x[a]+sgn)%L
    return int(x@np.array([L*L,L,1]))


def gaussian_entropy_modular(C):
    ev,u=np.linalg.eigh(C);ev=np.clip(ev,1e-14,1-1e-14)
    entropy=float(-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev)))
    mod=(u*np.log((1-ev)/ev))@u.conj().T
    return entropy,mod


comb_cache={}
def gaussian_density_blocks(C):
    """Exact number-sector density blocks via exterior powers, <=8 modes."""
    ev,u=np.linalg.eigh(C);ev=np.clip(ev,0,1);n=len(ev)
    blocks=[]
    for count in range(n+1):
        key=(n,count)
        if key not in comb_cache:
            comb_cache[key]=np.array(list(itertools.combinations(range(n),count)),dtype=int).reshape(-1,count) if count else np.empty((1,0),int)
        inds=comb_cache[key]
        if count==0:
            ext=np.ones((1,1),complex)
        else:
            ext=np.linalg.det(u[inds[:,None,:,None],inds[None,:,None,:]])
        occ=np.zeros((len(inds),n),bool)
        if count:occ[np.arange(len(inds))[:,None],inds]=True
        prob=np.prod(np.where(occ,ev[None,:],1-ev[None,:]),axis=1)
        blocks.append((ext*prob)@ext.conj().T)
    return blocks


def mixed_entropy(C0,C1,eps):
    b0=gaussian_density_blocks(C0);b1=gaussian_density_blocks(C1)
    values=np.concatenate([np.linalg.eigvalsh((1-eps)*a+eps*b) for a,b in zip(b0,b1)])
    values=np.maximum(values,0)
    nz=values>0
    return float(-np.sum(values[nz]*np.log(values[nz]))),float(values.sum())


def bernoulli_entropy(p):
    if p<0 or p>1:raise ValueError('selector probability outside physical domain')
    return -sum(x*math.log(x) for x in (p,1-p) if x>0)


def main():
    H=source(np.zeros(NC));P,Q,pair_energy,e,u=vacuum_and_pair(H)
    dC=Q-P
    singles={};pairs={}
    for c in range(NC):
        ids=region([c]);S,K=gaussian_entropy_modular(P[np.ix_(ids,ids)])
        singles[c]=(ids,S,K,float(np.trace(K@dC[np.ix_(ids,ids)]).real))
        for a in range(3):
            d=neighbor(c,a);ids2=region([c,d]);S2,K2=gaussian_entropy_modular(P[np.ix_(ids2,ids2)])
            pairs[c,a]=(ids2,S2,K2,float(np.trace(K2@dC[np.ix_(ids2,ids2)]).real))
    # All actual source-state derivatives, with no geometry equation supplied.
    de=e[:,None]-e[None,:];df=(e<0).astype(float)[:,None]-(e<0).astype(float)[None,:]
    divided=np.zeros_like(de);nz=df!=0;divided[nz]=df[nz]/de[nz]
    ds_single=np.zeros((NC,NC));ds_pair=np.zeros((3,NC,NC))
    for b in range(NC):
        mask=np.zeros(NC);mask[b]=1;M=np.diag(mask);mm=np.repeat(mask,4)
        dH=-.5*(mm[:,None]*KIN+KIN*mm[None,:])-np.kron(Lap2@M+2*Lap@M@Lap+M@Lap2,BETA)
        dP=u@(divided*(u.conj().T@dH@u))@u.conj().T
        for c,(ids,S,K,J) in singles.items():ds_single[c,b]=np.trace(K@dP[np.ix_(ids,ids)]).real
        for (c,a),(ids,S,K,J) in pairs.items():ds_pair[a,c,b]=np.trace(K@dP[np.ix_(ids,ids)]).real
    rows=[];rhs=[];metadata=[]
    for a in range(3):
        for c in range(NC):
            d=neighbor(c,a)
            ids,S0,K0,J0=singles[c];ids2,S1,K1,J1=pairs[c,a]
            row=.5*(ds_single[c]+ds_pair[a,c])
            row=row.copy();row[c]-=3*(S1-S0);row[d]-=1.5*(S1-S0)
            rows.append(row);rhs.append(.5*(J0+J1));metadata.append(dict(axis=a,centre_cell=c,neighbor=d))
    matrix=np.array(rows);J=np.array(rhs)
    train=matrix[:NC];jt=J[:NC]
    sigma,_,rank,sv=np.linalg.lstsq(train,-jt,rcond=1e-11)
    residual=matrix@sigma+J
    joint,_,jrank,jsv=np.linalg.lstsq(matrix,-J,rcond=1e-11)
    # Expected-volume selector uniquely fixed by branch volumes, no entropy fit.
    eps0=.01/max(1,float(np.max(abs(sigma))))
    validations=[]
    for eps in (eps0,eps0/2,eps0/4):
        sg=eps*sigma;HP=source(sg);PP,QQ,ep,_,_=vacuum_and_pair(HP)
        volumes=np.exp(3*sg)
        s0vals={};s1vals={};trace_error=0
        for c,(ids,S,K,J0) in singles.items():
            actual,tr=mixed_entropy(PP[np.ix_(ids,ids)],QQ[np.ix_(ids,ids)],eps)
            s0vals[c]=actual;trace_error=max(trace_error,abs(tr-1))
        for (c,a),(ids,S,K,J1) in pairs.items():
            actual,tr=mixed_entropy(PP[np.ix_(ids,ids)],QQ[np.ix_(ids,ids)],eps)
            s1vals[c,a]=actual;trace_error=max(trace_error,abs(tr-1))
        changes=[];selectors=[];reverse_changes=[]
        for meta in metadata:
            c,a,d=meta['centre_cell'],meta['axis'],meta['neighbor']
            w=(1.5-volumes[c])/volumes[d]
            actual=bernoulli_entropy(w)+(1-w)*s0vals[c]+w*s1vals[c,a]
            reference=math.log(2)+.5*(singles[c][1]+pairs[c,a][1])
            changes.append(actual-reference);selectors.append(w)
            if a==0:
                wr=(1.5-volumes[d])/volumes[c]
                actualr=bernoulli_entropy(wr)+(1-wr)*s0vals[d]+wr*s1vals[c,a]
                refr=math.log(2)+.5*(singles[d][1]+pairs[c,a][1])
                reverse_changes.append(actualr-refr)
        validations.append(dict(epsilon=eps,source_state_number_difference=float(np.trace(QQ-PP).real),
                 mixed_density_trace_error=trace_error,selector_min=min(selectors),selector_max=max(selectors),
                 all_region_entropy_changes=changes,trained_max_abs=max(abs(x) for x in changes[:NC]),
                 heldout_y_max_abs=max(abs(x) for x in changes[NC:2*NC]),
                 heldout_z_max_abs=max(abs(x) for x in changes[2*NC:]),
                 reversed_x_max_abs=max(abs(x) for x in reverse_changes),
                 external_pair_preparation_flag_entropy=bernoulli_entropy(eps),
                 geometry_vacuum_energy_shift=float(np.trace(HP@PP-H@P).real),
                 actual_pair_energy=ep))
        print('validated epsilon',eps,flush=True)
    out=dict(architecture='Declared entropy-based preparation with expected-volume mixed-region readout',
             literal_branchwise_fixed_volume=False,autonomous_native_law=False,Newton_G_used=False,
             L=L,one_particle_modes=NM,neutral_pair_energy=pair_energy,
             initial_pair_number_difference=float(np.trace(Q-P).real),
             global_sigma=sigma.tolist(),sigma_mean=float(np.mean(sigma)),
             source_defined_training_matrix=train.tolist(),matter_response=jt.tolist(),
             train_rank=int(rank),train_singular_values=sv.tolist(),
             trained_linear_residual_max=float(np.max(abs(residual[:NC]))),
             heldout_linear_residual_max=float(np.max(abs(residual[NC:]))),
             all_orientation_least_squares_residual_max=float(np.max(abs(matrix@joint+J))),
             all_orientation_rank=int(jrank),zero_mode_response_norm=float(np.linalg.norm(train@np.ones(NC))),
             validations=validations,region_metadata=metadata,
             matched_Einstein_stationarity_proved=False)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    np.savez_compressed(HERE/'SOURCE_ARRAYS.npz',coordinates=coords,H_reference=H,P_reference=P,
                        P_pair=Q,K_all=matrix,J_all=J,sigma=sigma,sigma_all_orientations=joint)
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_defined_training_matrix','validations','global_sigma','region_metadata')},indent=2))


if __name__=='__main__':main()
