"""Fixed-topology coframe diagnostic with a proper Wilson regulator.

Declared candidate H = v alpha.sin(k) + beta r v^2 sum(1-cos(k)),
v=exp(-sigma).  In physical momenta p=v*k/h, its low-momentum symbol
is alpha.p + (r*h/2) beta |p|^2, so r*h is an unchanged proper
regulator length.  No empirical constant is used and no native selection
or Einstein stationarity is asserted.  Transverse integration is quadrature.
"""
import json
import math
import os
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
import numpy as np

HERE=Path(__file__).resolve().parent
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)


def area_entropy(r,sigma,Nx,nquad,order=2):
    # Four-component Clifford source reduces at each transverse momentum
    # to two conjugate Pauli sources. Each has equal regional entropy.
    # A half-period slab has two boundaries: the factors cancel exactly.
    v=math.exp(-sigma)
    k=2*np.pi*(np.arange(Nx)+.5)/Nx
    h=Nx//2
    dif=np.arange(h)[:,None]-np.arange(h)[None,:]
    idx=dif%Nx
    phase=np.exp(1j*np.pi*dif/Nx)
    ab,wt=np.polynomial.legendre.leggauss(nquad)
    transverse=(ab+1)*np.pi/2
    wt=wt/2 # integral over each positive quadrant, normalized by pi
    total=0.
    for j,ky in enumerate(transverse):
        for l,kz in enumerate(transverse):
            if l<j:
                continue
            perp=math.sqrt(math.sin(ky)**2+math.sin(kz)**2)
            lap=2*(3-np.cos(k)-math.cos(ky)-math.cos(kz))
            if order==2:
                regulator=r*v*v*lap/2
            elif order==4:
                regulator=r**3*v**4*lap**2
            else:
                raise ValueError('order must be 2 or 4')
            vector=np.stack([v*np.sin(k),np.full(Nx,v*perp),regulator],axis=1)
            energy=np.linalg.norm(vector,axis=1)
            h2=vector[:,0,None,None]*X+vector[:,1,None,None]*Y+vector[:,2,None,None]*Z
            p=(np.eye(2)-h2/energy[:,None,None])/2
            kernel=np.fft.ifft(p,axis=0)
            corr=(kernel[idx]*phase[:,:,None,None]).transpose(0,2,1,3).reshape(2*h,2*h)
            ev=np.clip(np.linalg.eigvalsh(corr),1e-15,1-1e-15)
            entropy=-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev))
            total+=wt[j]*wt[l]*entropy*(1 if j==l else 2)
    return float(total)


def main():
    import sys
    order=int(sys.argv[1]) if len(sys.argv)>1 else 2
    records=[]
    for r,Nx,nquad in [(1,32,10),(2,64,12),(4,128,16)]:
        vals=[area_entropy(r,s,Nx,nquad,order) for s in (-.01,0,.01)]
        derivative=(vals[2]-vals[0])/.02
        rec=dict(regulator_order=order,r=r,Nx=Nx,transverse_quadrature=nquad,entropy_per_coordinate_face=vals[1],
                 entropy_minus=vals[0],entropy_plus=vals[2],dS_dsigma=derivative,
                 dlogS_dsigma=derivative/vals[1],geometric_area_exponent=2,
                 relative_area_scaling_error=abs(derivative/(2*vals[1])-1))
        records.append(rec)
        print(json.dumps(rec),flush=True)
    out=dict(scope='Declared fixed-topology Wilson coframe candidate; no Einstein realization',records=records,
             all_area_responses_nonzero=all(abs(x['dS_dsigma'])>1e-4 for x in records),
             same_sign_as_proper_area=all(x['dS_dsigma']>0 for x in records),
             selected_from_observed_gravity=False,
             still_missing=['fixed-proper-volume local-curvature entropy contrast',
                            'quantitative modular/stress error on the same refined source',
                            'native law selecting regulator and geometry state',
                            'generalized-entropy stationarity with error o(ell^4)',
                            'same-source physical clock, positive universal eta and G'])
    (HERE/f'REGULATOR_ORDER{order}_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
