"""Independent continuum kernel outside the light cone in the ky=kz=0 sector.

Subtract the known Abel-regularized high-momentum sine transform before
quadrature. This removes the finite Fourier Brillouin edge from the test.
"""
import json
import numpy as np
from scipy.integrate import quad

def c(a):
    value,error=quad(lambda p:1/np.sqrt(1+p**6),0,np.inf,weight='cos',wvar=abs(a),epsabs=2e-11,limit=250,limlst=250)
    return value,error

def sb(a):
    def remainder(p):
        # Stable evaluation of p^3/sqrt(1+p^6)-1.
        z=p**3
        return -1/(np.sqrt(1+z*z)*(np.sqrt(1+z*z)+z))
    value,error=quad(remainder,0,np.inf,weight='sin',wvar=abs(a),epsabs=2e-11,limit=250,limlst=250)
    return np.sign(a)*value,error

rows=[]
t=4.
for gap in [1.,2.,4.,8.,12.,16.]:
    x=t+gap
    cm,em=c(x-t);cp,ep=c(x+t)
    sm,esm=sb(t-x);sp,esp=sb(t+x)
    a=(cm-cp)/(2*np.pi)
    b=(1/(t-x)+1/(t+x)+sm+sp)/(2*np.pi)
    # U(x,t)=a alpha_x - i b beta. Its operator norm is |a|+|b|.
    rows.append(dict(time=t,spacelike_gap=gap,kernel_alpha=a,kernel_beta_imag=-b,kernel_operator_norm=abs(a)+abs(b),quadrature_error_bound=(em+ep+esm+esp)/(2*np.pi)))
out=dict(status='CONTINUUM_NONZERO_SPACELIKE_PROPAGATOR',source='R=1 positive flattening, ky=kz=0 exact continuum sector',rows=rows,strict_microcausality=False,scope='Pointwise kernel, not a full 3D regional norm. Finite Fourier edge and worst-case UV operator norms require separate treatment.')
with open('CONTINUUM_LEAKAGE.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
