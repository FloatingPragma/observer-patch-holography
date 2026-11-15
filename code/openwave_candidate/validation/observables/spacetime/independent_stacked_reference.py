"""Independent cell-gauge Bloch reference, without importing the demo source."""
from pathlib import Path
import sys
import math
import numpy as np
BUNDLE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BUNDLE/'candidate'))
from oph_openwave.geometry import CARRIER_SEAMS,RING_PORTS

A=np.array([[math.sqrt(3),math.sqrt(3)/2,0.],[0.,1.5,0.],[0.,0.,1.]])
D=np.array([math.sqrt(3)/2,.5,0.])


def incidence(q,vertical=1.):
    """ifft convention: cell wave exp(+i q.n); horizontal A(R)->B(R-t)."""
    q=np.asarray(q,float)
    b=np.zeros((24,87),complex)
    for sub in (0,1):
        for e,(a,c) in enumerate(CARRIER_SEAMS):
            b[12*sub+a,30*sub+e]=1;b[12*sub+c,30*sub+e]=-1
    shifts=((0,0),(1,0),(0,1))
    for e,(port,t) in enumerate(zip(RING_PORTS,shifts)):
        b[port,60+e]=1;b[12+port,60+e]=-np.exp(1j*np.dot(q[:2],t))
    b[:,63:]=vertical*(1-np.exp(-1j*q[2]))*np.eye(24)
    return b


def positive_acoustic(q,mass=.23,vertical=1.):
    b=incidence(q,vertical);values,vectors=np.linalg.eigh(b@b.conj().T)
    lam=float(max(0,values[0]));u=vectors[:,0]
    phase=np.vdot(np.ones(24),u)
    if abs(phase)>1e-12:u*=np.conj(phase)/abs(phase)
    energy=math.sqrt(mass**2+lam)
    spinor=np.r_[math.sqrt((energy+mass)/(2*energy))*u,
                 b.conj().T@u/math.sqrt(2*energy*(energy+mass))]
    return energy,spinor,b


def hessian_squared_energy(vertical=1.,step=.002):
    f=lambda q:positive_acoustic(q,vertical=vertical)[0]**2
    zero=np.zeros(3);hessian=np.zeros((3,3));f0=f(zero)
    for i in range(3):
        a=np.eye(3)[i]*step
        hessian[i,i]=(f(a)+f(-a)-2*f0)/step**2
        for j in range(i):
            b=np.eye(3)[j]*step
            hessian[i,j]=hessian[j,i]=(f(a+b)-f(a-b)-f(-a+b)+f(-a-b))/(4*step**2)
    return hessian
