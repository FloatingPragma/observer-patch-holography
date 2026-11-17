"""Propagation of the ORIGINAL finite-entropy regional CAR algebra."""
import json,time
import numpy as np
from scipy.linalg import svdvals
from source import ALPHA,I4,symbol

start=time.process_time()
rows=[]
R=1.; box=64.; time_value=4.
for count in [256,512]:
    h=box/count
    x=(np.arange(count)-count//2)*h
    p=2*np.pi*np.fft.fftfreq(count,d=h)
    momenta=np.zeros((count,3));momenta[:,0]=p
    H,_,E=symbol(momenta,R)
    sinc=np.divide(np.sin(E*time_value),E,out=np.full_like(E,time_value),where=E>0)
    U=np.cos(E*time_value)[:,None,None]*I4-1j*sinc[:,None,None]*H
    Udirac=np.cos(p*time_value)[:,None,None]*I4-1j*np.sin(p*time_value)[:,None,None]*ALPHA[0]
    flat_kernel=np.fft.ifft(U,axis=0)
    dirac_kernel=np.fft.ifft(Udirac,axis=0)
    for ell in [2.,4.,8.]:
        inside=np.flatnonzero(abs(x)<=ell)
        for margin in sorted(set([1.,2.,4.,ell])):
            outside=np.flatnonzero(abs(x)>ell+time_value+margin)
            delta=(outside[:,None]-inside[None,:])%count
            M=flat_kernel[delta].transpose(0,2,1,3).reshape(4*len(outside),4*len(inside))
            D=dirac_kernel[delta].transpose(0,2,1,3).reshape(M.shape)
            rows.append(dict(count=count,h=h,box=box,R=R,time=time_value,ell=ell,margin=margin,number_input_modes=4*len(inside),number_outside_modes=4*len(outside),flat_operator_norm=float(svdvals(M,check_finite=False)[0]),flat_mean_squared_leakage=float(np.linalg.norm(M)**2/(4*len(inside))),dirac_frobenius=float(np.linalg.norm(D))))
    print(json.dumps(dict(count=count,cpu_seconds=time.process_time()-start)),flush=True)

out=dict(status='FINITE_FOURIER_PROPAGATION_LEAKAGE',convention='Original chi_B regions, periodic box; all four spinor components, no band projection.',rows=rows,cpu_seconds=time.process_time()-start,strict_microcausality=False,static_isotony='Exact chi_A chi_B=chi_A for A subset B; also exact under any one common unitary.',cross_observer_geometric_isotony='Not established by static nesting or the Poincare Lie algebra.')
with open('LEAKAGE.json','w') as f:json.dump(out,f,indent=2)
