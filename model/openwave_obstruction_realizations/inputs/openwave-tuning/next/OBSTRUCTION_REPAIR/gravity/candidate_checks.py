"""Same-source checks for the declared quartic proper-regulator candidate."""
import json
import math
import os
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
import numpy as np
from regulator_candidate import X,Y,Z,area_entropy

HERE=Path(__file__).resolve().parent


def modular_low_mode(R,L,sigma=0):
    v=math.exp(-sigma)
    k=2*np.pi*(np.arange(L)+.5)/L
    h=L//2
    H=v*np.sin(k)[:,None,None]*X+R**3*v**4*(2*(1-np.cos(k)))[:,None,None]**2*Z
    E,U=np.linalg.eigh(H)
    P=(U*(E<0)[:,None,:])@np.swapaxes(U.conj(),-1,-2)
    dif=np.arange(h)[:,None]-np.arange(h)[None,:]
    kernel=np.fft.ifft(P,axis=0)
    C=(kernel[dif%L]*np.exp(1j*np.pi*dif/L)[:,:,None,None]).transpose(0,2,1,3).reshape(2*h,2*h)
    c,u=np.linalg.eigh(C)
    floor_records=[]
    for floor in (1e-12,1e-13,1e-14):
        cf=np.clip(c,floor,1-floor)
        mod=(u*np.log((1-cf)/cf))@u.conj().T
        phase=np.exp(1j*k[0]*np.arange(h))/math.sqrt(L)
        waves=np.einsum('x,ab->xab',phase,U[0]).reshape(2*h,2)
        slopes=np.real(np.diag(waves.conj().T@mod@waves))
        ratio=(slopes[1]-slopes[0])/(E[0,1]-E[0,0])
        expected=L/(np.pi*v)
        floor_records.append(dict(eigenvalue_floor=floor,bulk_coefficient=ratio,
                                  expected_cft_coefficient=expected,ratio=ratio/expected))
    return dict(R=R,L=L,sigma=sigma,transverse_momentum=[0,0],floor_checks=floor_records)


def reduction_check():
    L=16;R=2;sigma=.1;v=math.exp(-sigma);ky=.43;kz=.71
    k=2*np.pi*(np.arange(L)+.5)/L
    alpha=[np.kron(p,X) for p in (X,Y,Z)]
    beta=np.kron(np.eye(2),Z)
    lap=2*(3-np.cos(k)-math.cos(ky)-math.cos(kz))
    vectors=[v*np.sin(k),np.full(L,v*math.sin(ky)),np.full(L,v*math.sin(kz))]
    h4=sum(p[:,None,None]*a for p,a in zip(vectors,alpha))+R**3*v**4*lap[:,None,None]**2*beta
    h2=v*np.sin(k)[:,None,None]*X+v*math.hypot(math.sin(ky),math.sin(kz))*Y+R**3*v**4*lap[:,None,None]**2*Z
    entropies=[]
    for source in (h2,h4):
        E,U=np.linalg.eigh(source); nb=source.shape[-1]
        P=(U*(E<0)[:,None,:])@np.swapaxes(U.conj(),-1,-2)
        dif=np.arange(L//2)[:,None]-np.arange(L//2)[None,:]
        ker=np.fft.ifft(P,axis=0)
        C=(ker[dif%L]*np.exp(1j*np.pi*dif/L)[:,:,None,None]).transpose(0,2,1,3).reshape(nb*L//2,nb*L//2)
        ev=np.clip(np.linalg.eigvalsh(C),1e-15,1-1e-15)
        entropies.append(float(-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev))))
    return dict(S2=entropies[0],S4=entropies[1],residual=abs(entropies[1]-2*entropies[0]))


def main():
    out=dict(scope='Diagnostic candidate checks, not Einstein realization')
    out['four_component_reduction']=reduction_check()
    out['modular_low_modes']=[modular_low_mode(R,L,s) for R,L,s in [(1,64,0),(2,128,0),(4,128,0),(4,256,0),(4,256,.1)]]
    # This is a quadrature refinement of the R=4 entry, chosen after its
    # initial 1.03-percent discrepancy. It is not a held-out physics test.
    vals=[area_entropy(4,s,128,24,4) for s in (-.01,0,.01)]
    out['area_quadrature_refinement']=dict(R=4,Nx=128,nquad=24,entropy_per_face=vals[1],
                dlogS_dsigma=(vals[2]-vals[0])/(.02*vals[1]),after_initial_diagnostic=True)
    out['checks']={
        'four_component_reduction':out['four_component_reduction']['residual']<1e-10,
        'low_mode_error_decreases':abs(out['modular_low_modes'][3]['floor_checks'][1]['ratio']-1)<abs(out['modular_low_modes'][2]['floor_checks'][1]['ratio']-1),
        'low_mode_subpercent_on_large_interval':abs(out['modular_low_modes'][3]['floor_checks'][1]['ratio']-1)<.01,
        'coframe_clock_consistent_subpercent':abs(out['modular_low_modes'][4]['floor_checks'][1]['ratio']-1)<.01,
        'modular_floor_stable':max(abs(r['floor_checks'][0]['ratio']-r['floor_checks'][-1]['ratio']) for r in out['modular_low_modes'])<.001,
        'area_exponent_within_three_percent':abs(out['area_quadrature_refinement']['dlogS_dsigma']/2-1)<.03,
    }
    out['checks']={key:bool(value) for key,value in out['checks'].items()}
    out['all_checks_pass']=all(out['checks'].values())
    (HERE/'CANDIDATE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    if not out['all_checks_pass']:
        raise SystemExit(1)


if __name__=='__main__':main()
