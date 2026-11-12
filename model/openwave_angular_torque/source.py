"""Projected oscillator differential operators and actual current derivatives."""
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.special import eval_hermite, factorial, roots_hermite
from scipy.linalg import eigh
from model.openwave_intrinsic_magnetic.solver import BETA,GAMMA,SIGMA,native_reduction


def cross(i,j,k):
    return ((i,j,k) in ((0,1,2),(1,2,0),(2,0,1)))-((i,j,k) in ((0,2,1),(2,1,0),(1,0,2)))


@lru_cache(maxsize=12)
def oscillator(shell=2,length=1.,radius=1.2):
    if shell<2 or length<=0 or radius<=0:raise ValueError('at least two complete shells and positive lengths required')
    basis=tuple(p for p in product(range(shell+1),repeat=3) if sum(p)<=shell)
    where={p:i for i,p in enumerate(basis)};n=len(basis);ann=[]
    for axis in range(3):
        a=np.zeros((n,n),complex)
        for j,p in enumerate(basis):
            if p[axis]:
                q=list(p);q[axis]-=1;a[where[tuple(q)],j]=np.sqrt(p[axis])
        ann.append(a)
    x=tuple(length*(a+a.conj().T)/np.sqrt(2) for a in ann)
    p=tuple(-1j*(a-a.conj().T)/(length*np.sqrt(2)) for a in ann)
    # These are P_N (-i r cross grad) P_N. Products P_N r P_N p P_N
    # have a different boundary term and cannot define the angular generator.
    orbital=tuple(sum((-1j*cross(i,j,k)*ann[j].conj().T@ann[k] for j,k in product(range(3),repeat=2)),np.zeros((n,n),complex))for i in range(3))
    number=np.diag([sum(q) for q in basis]).astype(complex)
    # Exact Gaussian-polynomial integrals through Gauss-Hermite quadrature.
    # With the changed weight, degree <= 2*shell+1 needs shell+1 nodes.
    nodes,weights=roots_hermite(shell+3);c=1+length**2/(2*radius**2)
    axis_moments=[]
    for power in (0,1):
        out=np.zeros((shell+1,shell+1))
        for i,j in product(range(shell+1),repeat=2):
            norm=np.sqrt(2.**(i+j)*factorial(i)*factorial(j)*np.pi)
            out[i,j]=sum(weights*eval_hermite(i,nodes/np.sqrt(c))*eval_hermite(j,nodes/np.sqrt(c))*(length*nodes/np.sqrt(c))**power)/(np.sqrt(c)*norm)
        axis_moments.append(out)
    radial_x=[]
    for axis in range(3):
        radial_x.append(np.array([[np.prod([axis_moments[int(d==axis)][u[d],v[d]] for d in range(3)]) for v in basis]for u in basis]))
    # A_b(r)=exp(-r^2/(2R^2)) (b cross r)/2, projected before coupling.
    potential=tuple(tuple(sum((cross(i,j,k)*radial_x[k]/2 for k in range(3)),np.zeros((n,n))) for i in range(3))for j in range(3))
    uniform=tuple(tuple(sum((cross(i,j,k)*x[k]/2 for k in range(3)),np.zeros((n,n),complex))for i in range(3))for j in range(3))
    full_J=tuple(np.kron(np.eye(4),orbital[i])+np.kron(np.kron(np.eye(2),SIGMA[i]/2),np.eye(n))for i in range(3))
    return {'basis':basis,'x':x,'p':p,'L':orbital,'N':number,'J':full_J,'potential':potential,'uniform':uniform,
        'electric_metric':np.pi**1.5*radius**5/4,'magnetic_metric':5*np.pi**1.5*radius**3/8,'spatial_dimension':n}


def operators(data,mass=8.,omega=.2,charge=-1.):
    n=data['spatial_dimension'];h=np.kron(BETA,mass*np.eye(n))
    for i in range(3):h+=np.kron(GAMMA[i],data['p'][i])
    derivatives=[];external=[]
    for j in range(3):
        derivatives.append(-charge*sum((np.kron(GAMMA[i],data['potential'][j][i])for i in range(3)),np.zeros_like(h)))
        external.append(-charge*sum((np.kron(GAMMA[i],data['uniform'][j][i])for i in range(3)),np.zeros_like(h)))
    return h,tuple(derivatives),tuple(external),omega*np.kron(np.eye(4),data['N'])


def kernel(h):
    values,vectors=eigh(h);absolute=np.abs(values)
    if min(absolute)<1e-12:raise ValueError('a nonzero signed gap is required')
    return (vectors*absolute)@vectors.conj().T,values,vectors


def evaluate(psi,b,ops,external=(0.,0.,0.)):
    h0,derivatives,probe,trap=ops
    h=h0+sum((v*x for v,x in zip(b,derivatives)),np.zeros_like(h0))+sum((v*x for v,x in zip(external,probe)),np.zeros_like(h0))
    k,w,u=kernel(h);v=u.conj().T@psi
    divided=(w[:,None]+w[None,:])/(np.abs(w[:,None])+np.abs(w[None,:]))
    def gradient(d):return float(np.vdot(v,(divided*(u.conj().T@d@u))@v).real)
    return k+trap,np.array([gradient(d)for d in derivatives]),np.array([gradient(d)for d in probe])


def angular_expectation(psi,data):return np.array([np.vdot(psi,j@psi).real for j in data['J']])


def field_energy(b,pi,data,coupling):
    return coupling*np.dot(pi,pi)/(2*data['electric_metric'])+data['magnetic_metric']*np.dot(b,b)/(2*coupling)


def rhs(t,state,data,ops,coupling):
    n=ops[0].shape[0];psi=state[:n];b=state[n:n+3].real;pi=state[n+3:n+6].real
    h,gradient,_=evaluate(psi,b,ops)
    return np.r_[-1j*h@psi,coupling*pi/data['electric_metric'],-data['magnetic_metric']*b/coupling-gradient]
