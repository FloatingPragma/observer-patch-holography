import numpy as np
from kernels import sibling
bound=sibling('actual_neutral_mixture_bound','neutral_mixture_entropy_bound/bound.py')

def response_with_bound(region,source,matter,centre):
    coords=region['cells']+np.asarray(centre);ix=coords%source['L']
    phase=np.exp(2j*np.pi*source['twist']*coords.sum(-1)/source['L'])
    local=[(matter[name][ix[:,0],ix[:,1],ix[:,2]]*phase[:,None,None]).reshape(4*len(coords),4) for name in ['particle','hole']]
    X=np.hstack(local);M=np.zeros((8,8),complex);G0=np.zeros_like(M);G1=np.zeros_like(M);G2=np.zeros_like(M);J=0.
    signs=np.r_[np.ones(4),-np.ones(4)]
    for name,d in [('G1',2),('G2',2),('H',4)]:
        e,U=region['eigensystems'][name]
        if np.min(e)<=0 or np.max(e)>=1:raise ValueError('Actual reference not faithful; no floor inserted in inverse.')
        for component in range(d):
            reduced=np.zeros((region['sizes'][name],8),complex)
            for o in region['orbits']:
                reduced[o['columns'][name]]=o['sectors'][name]['bases'][component].conj().T@X[o['rows']]
            Y=U.conj().T@reduced
            M+=Y.conj().T@Y
            G0+=Y.conj().T@(Y/(e*(1-e))[:,None])
            G1+=Y.conj().T@(Y/(1-e)[:,None])
            G2+=Y.conj().T@(Y*(e/(1-e))[:,None])
            J+=float(np.sum(np.log((1-e)/e)[:,None]*abs(Y)**2*signs[None,:]))
    z=bound.log_Z_from_grams(M,G0,G1,G2,signs)
    if not z['numerical_gate_pass']:
        return dict(centre=list(centre),actual_mixture_entropy_tangent=J,low_rank=z,finite_epsilon_intervals=[],reason='Low-rank determinant numerical gate failed; no entropy interval claimed.')
    intervals=[bound.finite_entropy_interval(J,eps,z['log_Z']) for eps in [1e-2,1e-3,1e-4]]
    if z['log_Z']<700 and J>0:
        chi=np.expm1(max(0.,z['log_Z']))
        recommended=min(.01,.005*J/max(chi,1e-300))
        controlled=bound.finite_entropy_interval(J,recommended,z['log_Z'])
    else:controlled=None
    return dict(centre=list(centre),actual_mixture_entropy_tangent=J,low_rank=z,finite_epsilon_intervals=intervals,half_percent_controlled_interval=controlled,external_preparation_flag='Retained in complete source record, excluded from this regional algebra; both branches exactly neutral.')

def matter_fourier_kernel(matter,source,shift):
    L=source['L'];grid=np.stack(np.meshgrid(*([np.arange(L)]*3),indexing='ij'),-1)
    shifted=(grid+np.asarray(shift))%L
    result=np.zeros((L,L,L,4,4),complex)
    for name,sign in [('particle',1),('hole',-1)]:
        A=np.fft.fftn(matter[name],axes=(0,1,2))/np.sqrt(L**3)
        Aq=A[shifted[...,0],shifted[...,1],shifted[...,2]]
        result+=sign*(Aq@A.conj().swapaxes(-1,-2))
    return np.fft.ifftn(result,axes=(0,1,2))
