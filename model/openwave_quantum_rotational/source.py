"""Complete-shell matter and total-occupation Maxwell operators."""
from functools import lru_cache
from itertools import product
from fractions import Fraction
import numpy as np
import sympy as sp
from scipy.linalg import eigh
from scipy.special import gamma,gammainc
from model.openwave_intrinsic_magnetic.solver import SIGMA,native_reduction
from model.openwave_angular_torque.source import cross

MASS=16
ETA=Fraction(1,8)
RADIUS_SQUARED=Fraction(5,2)


def occupations(cap):
    return tuple(sorted((p for p in product(range(cap+1),repeat=3)if sum(p)<=cap),key=lambda p:(sum(p),p)))


def ladders(cap,exact=False):
    basis=occupations(cap);where={p:i for i,p in enumerate(basis)};n=len(basis)
    matrices=[]
    for axis in range(3):
        a=sp.zeros(n)if exact else np.zeros((n,n),complex)
        for j,p in enumerate(basis):
            if p[axis]:
                q=list(p);q[axis]-=1
                a[where[tuple(q)],j]=sp.sqrt(p[axis])if exact else np.sqrt(p[axis])
        matrices.append(a)
    return basis,tuple(matrices)


def exact_axis_moment(i,j,power):
    x=sp.Symbol('x');c=sp.Rational(6,5)
    poly=sp.Poly(sp.hermite(i,x)*sp.hermite(j,x)*x**power,x)
    integral=0
    for (degree,),coefficient in poly.terms():
        if degree%2==0:
            k=degree//2
            integral+=coefficient*sp.factorial2(2*k-1)/(2**k*c**(sp.Rational(1,2)+k))
    return sp.simplify(integral/sp.sqrt(2**(i+j)*sp.factorial(i)*sp.factorial(j)))


@lru_cache(maxsize=1)
def matter_exact():
    basis,a=ladders(2,True);n=len(basis);zero=sp.zeros(n)
    orbital=tuple(sum((-sp.I*cross(i,j,k)*a[j].H*a[k]for j,k in product(range(3),repeat=2)),zero)for i in range(3))
    p=tuple(-sp.I*(v-v.H)/sp.sqrt(2)for v in a)
    x=tuple((v+v.H)/sp.sqrt(2)for v in a)
    radial=[]
    for axis in range(3):
        radial.append(sp.Matrix([[sp.prod(exact_axis_moment(u[d],v[d],int(d==axis))for d in range(3))for v in basis]for u in basis]))
    sigma=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))
    j=tuple(sp.kronecker_product(sp.eye(2),orbital[i])+sp.kronecker_product(sigma[i]/2,sp.eye(n))for i in range(3))
    b0=sum((sp.kronecker_product(sigma[i],p[i])for i in range(3)),sp.zeros(2*n))
    # q=-1 and A_b=exp(-r²/5)(b cross r)/2.
    derivative=tuple(sum((sp.kronecker_product(sigma[i],sum((cross(i,k,d)*radial[d]/2 for d in range(3)),zero))for i in range(3)),sp.zeros(2*n))for k in range(3))
    external=tuple(sum((sp.kronecker_product(sigma[i],sum((cross(i,k,d)*x[d]/2 for d in range(3)),zero))for i in range(3)),sp.zeros(2*n))for k in range(3))
    number=sp.diag(*[sum(p)for p in basis])
    r2=number+sp.Rational(3,2)*sp.eye(n)+sum(((v*v+v.H*v.H)/2 for v in a),zero)
    return {'basis':basis,'B':b0,'derivative':derivative,'external':external,'J':j,'L':orbital,'N':number,'r2':r2,'sigma':sigma}


def array(m):return np.array(m.evalf(18),complex)


def photon(cap):
    if cap not in(1,2,3,4):raise ValueError('declared complete field caps are1 through4')
    basis,a=ladders(cap);n=len(basis);zero=np.zeros((n,n),complex)
    j=tuple(sum((-1j*cross(i,k,l)*a[k].conj().T@a[l]for k,l in product(range(3),repeat=2)),zero)for i in range(3))
    number=np.diag([sum(p)for p in basis]).astype(complex)
    q=tuple(float(ETA)*(v+v.conj().T)for v in a)
    pi=tuple((v-v.conj().T)/(2j*float(ETA))for v in a)
    top=np.diag([int(sum(p)==cap)for p in basis])
    return {'basis':basis,'a':a,'J':j,'N':number,'q':q,'pi':pi,'top':top}


@lru_cache(maxsize=4)
def source(cap):
    exact=matter_exact();f=photon(cap);d=len(f['basis']);n=20
    b=np.kron(array(exact['B']),np.eye(d))
    for derivative,q in zip(exact['derivative'],f['q']):b+=np.kron(array(derivative),q)
    eigen,u=eigh(b);energies=np.sqrt(MASS**2+eigen**2)
    # Stable positive square-root difference, retaining both signed B lines.
    kminus=(u*(eigen**2/(energies+MASS)))@u.conj().T
    nm=np.kron(np.kron(np.eye(2),array(exact['N'])),np.eye(d))
    nf=np.kron(np.eye(n),f['N'])
    hminus=kminus+nm+nf;h=hminus+MASS*np.eye(n*d)
    jm=tuple(np.kron(array(j),np.eye(d))for j in exact['J'])
    jf=tuple(np.kron(np.eye(n),j)for j in f['J']);j=tuple(x+y for x,y in zip(jm,jf))
    j2=sum((v@v for v in j),np.zeros_like(h))
    w,v=eigh(hminus);ground=v[:,:2];z,rotation=eigh(ground.conj().T@j[2]@ground)
    psi=ground@rotation[:,np.argmax(z)];seed=np.zeros(len(psi),complex);seed[0]=1
    overlap=np.vdot(seed,psi);psi*=np.exp(-1j*np.angle(overlap))
    assert w[2]>.99 and abs(w[1]-w[0])<1e-12
    ext=np.kron(array(exact['external'][2]),np.eye(d));rhs=b@ext+ext@b
    derivative=u@((u.conj().T@rhs@u)/(energies[:,None]+energies[None,:]))@u.conj().T
    return {'cap':cap,'B':b,'K_minus_mass':kminus,'H':h,'H_minus_mass':hminus,'Nm':nm,'Nf':nf,'Jm':jm,'Jf':jf,'J':j,'J2':j2,
            'psi':psi,'seed':seed,'E':float(MASS+w[0]),'all_energies':MASS+w,'current_derivative':derivative,'photon':f}


def radial_effect(radius):
    data=matter_exact();basis=data['basis'];n=len(basis);out=np.zeros((n,n));x=sp.Symbol('x')
    for i,u in enumerate(basis):
        for j,v in enumerate(basis):
            polynomials=[sp.Poly(sp.hermite(u[d],x)*sp.hermite(v[d],x),x).terms()for d in range(3)]
            value=0.
            for terms in product(*polynomials):
                powers=[t[0][0]for t in terms]
                if any(p%2 for p in powers):continue
                half=[p//2 for p in powers];coefficient=float(sp.prod(t[1]for t in terms))
                value+=coefficient*np.prod([gamma(p+.5)for p in half])*gammainc(sum(half)+1.5,radius**2)
            norm=np.pi**1.5*np.sqrt(np.prod([2**(u[d]+v[d])*float(sp.factorial(u[d])*sp.factorial(v[d]))for d in range(3)]))
            out[i,j]=value/norm
    return out
