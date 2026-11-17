"""Independent planar normalization and low-mode modular checks for order4."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from regulator_candidate import X,Y,Z
R=1.5

def symbol(k):
    t=np.sin(k/2)**2
    return np.sin(k)*(1+2*t/3),4*t*(1+t/3)

def regional(H,L):
    e,u=np.linalg.eigh(H);p=(u*(e<0)[:,None,:])@np.swapaxes(u.conj(),-1,-2)
    delta=np.arange(L//2)[:,None]-np.arange(L//2)[None,:]
    kernel=np.fft.ifft(p,axis=0)
    C=(kernel[delta%L]*np.exp(1j*np.pi*delta/L)[:,:,None,None]).transpose(0,2,1,3).reshape(L,L)
    return C,e,u

def area(Nx,nquad,sigma=0):
    k=2*np.pi*(np.arange(Nx)+.5)/Nx;dx,lx=symbol(k);v=math.exp(-sigma)
    ab,wt=np.polynomial.legendre.leggauss(nquad);kt=(ab+1)*np.pi/2;wt=wt/2
    total=0.
    for j,ky in enumerate(kt):
        dy,ly=symbol(ky)
        for h,kz in enumerate(kt[j:],j):
            dz,lz=symbol(kz)
            H=v*dx[:,None,None]*X+v*math.hypot(dy,dz)*Y+R**3*v**4*(lx+ly+lz)[:,None,None]**2*Z
            C,_,_=regional(H,Nx);ev=np.clip(np.linalg.eigvalsh(C),1e-15,1-1e-15)
            S=-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev))
            total+=wt[j]*wt[h]*S*(1 if j==h else 2)
    return float(total)

def modular(L,sigma=0):
    v=math.exp(-sigma);k=2*np.pi*(np.arange(L)+.5)/L;d,l=symbol(k)
    H=v*d[:,None,None]*X+R**3*v**4*l[:,None,None]**2*Z
    C,E,U=regional(H,L);ev,u=np.linalg.eigh(C)
    phase=np.exp(1j*k[0]*np.arange(L//2))/math.sqrt(L)
    waves=np.einsum('x,ab->xab',phase,U[0]).reshape(L,2)
    weights=abs(u.conj().T@waves)**2
    checks=[]
    for floor in (1e-11,1e-13,1e-15):
        e=np.clip(ev,floor,1-floor)
        vals=np.log((1-e)/e)@weights
        coefficient=(vals[1]-vals[0])/(E[0,1]-E[0,0])
        checks.append(dict(floor=floor,coefficient=float(coefficient),ratio=float(coefficient/(L/(np.pi*v)))))
    return dict(L=L,sigma=sigma,expected=L/(np.pi*v),floor_checks=checks)

def main():
    records=[]
    for Nx,nq in [(64,16),(96,24)]:
        vals=[area(Nx,nq,sigma) for sigma in (-.01,0,.01)]
        rec=dict(Nx=Nx,nquad=nq,eta_per_coordinate_face=vals[1],uniform_area_exponent=(vals[2]-vals[0])/(.02*vals[1]))
        records.append(rec);print('area',rec,flush=True)
    modes=[modular(L,sigma) for L,sigma in [(128,0),(256,0),(256,.1)]]
    out=dict(status='IMPROVED_SOURCE_PLANAR_AND_LOW_MODE_COMPONENT_TEST',source_difference_order=4,R=R,
             area=records,eta_relative_resolution_difference=abs(records[1]['eta_per_coordinate_face']/records[0]['eta_per_coordinate_face']-1),
             modular=modes,full_modular_operator_identity_proved=False,stationarity=False,physical_G=False)
    (HERE/'IMPROVED_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
