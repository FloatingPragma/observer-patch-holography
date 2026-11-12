"""Own-ensemble loop ratios; no imported perturbative running curve."""
from pathlib import Path
from functools import lru_cache
import time
import numpy as np
from model.colour_deconfinement import su3_numba as K

RECTANGLES=((1,1),(2,1),(2,2),(3,2),(3,3))


@lru_cache(maxsize=1)
def fundamental_casimir():
    import sympy as s
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    raw=representation()[1]
    t=[x[:3,:3] for x in raw[:7]]+[(raw[6][:3,:3]+2*raw[7][:3,:3])/s.sqrt(3)]
    assert s.Matrix(8,8,lambda i,j:s.simplify(s.trace(t[i]*t[j])))==s.eye(8)/2
    casimir=sum((x*x for x in t),s.zeros(3)).applyfunc(s.simplify)
    assert casimir==s.Rational(4,3)*s.eye(3)
    return casimir[0,0]


def shifted(index,up,axis,count):
    out=index
    for _ in range(count):out=up[out,axis]
    return out


def loops(U,up):
    """Average both orientations of every plane; all configurations retained."""
    n=len(U);index=np.arange(n);lines={}
    for mu in range(4):
        p=np.broadcast_to(np.eye(3,dtype=complex),(n,3,3)).copy()
        for r in range(1,4):
            p=p@U[shifted(index,up,mu,r-1),mu]
            lines[mu,r]=p
    values=[]
    for r,t in RECTANGLES:
        total=0.
        for mu in range(4):
            for nu in range(4):
                if mu==nu:continue
                a=lines[mu,r]@lines[nu,t][shifted(index,up,mu,r)]
                b=lines[nu,t]@lines[mu,r][shifted(index,up,nu,t)]
                total+=np.einsum('nij,nij->',a,b.conj()).real/(3*n)
        values.append(total/12)
    return np.array(values)


def tree_loop(side,r,t):
    """Exact finite Fourier sum of the Gaussian loop, times C_F/2."""
    if not r or not t:return 0.
    k=2*np.pi*np.indices((side,)*4).reshape(4,-1).T/side
    hat=4*np.sin(k/2)**2;den=hat.sum(axis=1);den[0]=1
    ex=np.exp(1j*k[:,0]);ey=np.exp(1j*k[:,1])
    sx=sum(ex**i for i in range(r));sy=sum(ey**i for i in range(t))
    jx=sx*(1-ey**t);jy=sy*(ex**r-1)
    value=(abs(jx)**2+abs(jy)**2)/den;value[0]=0
    return float(fundamental_casimir())/2*float(value.mean())


def tree_coefficients(side):
    return np.array([tree_loop(side,r,r)+tree_loop(side,r-1,r-1)-
                     2*tree_loop(side,r,r-1) for r in (1,2,3)])


def creutz(w):
    w=np.asarray(w,dtype=float)
    if w.shape!=(5,) or not np.all(np.isfinite(w)) or np.any(w<=0):
        raise ValueError('five finite positive ensemble loop means are required')
    ratios=np.array([w[0],w[2]*w[0]/w[1]**2,w[4]*w[2]/w[3]**2])
    if np.any(ratios<=0):raise ValueError('nonpositive loop ratio; no real logarithmic coupling')
    return -np.log(ratios)


def jackknife(samples,side,block):
    """Nonlinear estimator with common block deletions across all scales."""
    samples=np.asarray(samples);nb=len(samples)//block
    if nb<4 or nb*block!=len(samples):raise ValueError('at least four complete equal blocks required')
    blocks=samples.reshape(nb,block,5).mean(axis=1)
    k=tree_coefficients(side);mean=blocks.mean(axis=0)
    estimates=np.array([creutz((nb*mean-row)/(nb-1))/k for row in blocks])
    center=estimates.mean(axis=0);delta=estimates-center
    covariance=(nb-1)/nb*delta.T@delta
    full=creutz(mean)/k
    differences=np.diff(full)
    transform=np.array([[-1,1,0],[0,-1,1]])
    dcov=transform@covariance@transform.T
    return {'block_samples':block,'blocks':nb,'g_squared':full.tolist(),
            'standard_errors':np.sqrt(np.diag(covariance)).tolist(),
            'covariance':covariance.tolist(),'successive_scale_differences':differences.tolist(),
            'difference_standard_errors':np.sqrt(np.diag(dcov)).tolist(),
            'loop_means':mean.tolist(),'creutz_ratios':creutz(mean).tolist(),
            'tree_geometry_coefficients':k.tolist()}


def generate(path,side=6,beta=6.,seed=101,start='cold',warmup=200,samples=400,stride=2):
    if not K.AVAILABLE:raise RuntimeError('Production requires the existing numba source kernel')
    up,dn=K.neighbours(*(side,)*4);K.seed_rng(seed)
    U=K.cold_start(side**4) if start=='cold' else K.hot_start(side**4)
    begin=time.monotonic();K.compound_sweeps(U,up,dn,beta,warmup,2)
    records=np.empty((samples,5));plaquette=np.empty(samples)
    for i in range(samples):
        K.compound_sweeps(U,up,dn,beta,stride,2)
        if i%25==0:K.reunitarize(U)
        records[i]=loops(U,up);plaquette[i]=sum(K.plaquette(U,up))/2
    if not np.allclose(records[:,0],plaquette,atol=2e-13):raise AssertionError('loop/independent plaquette disagreement')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,loops=records,plaquette=plaquette,side=side,beta=beta,
        seed=seed,start=start,warmup=warmup,stride=stride,samples=samples,
        overrelaxations=2,elapsed_seconds=time.monotonic()-begin,
        final_links=U)
    return {'path':str(path),'seconds':time.monotonic()-begin,'mean':records.mean(axis=0).tolist()}


def gauge_control(U,side,seed=812):
    from model.colour_deconfinement.su3 import random_su3
    up,_=K.neighbours(*(side,)*4);g=random_su3(side**4,np.random.default_rng(seed))
    transformed=np.empty_like(U)
    for mu in range(4):transformed[:,mu]=g@U[:,mu]@np.swapaxes(g[up[:,mu]].conj(),-1,-2)
    return float(np.max(abs(loops(U,up)-loops(transformed,up))))
