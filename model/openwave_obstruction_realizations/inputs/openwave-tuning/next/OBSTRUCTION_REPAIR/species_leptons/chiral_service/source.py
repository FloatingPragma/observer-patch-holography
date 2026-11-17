"""Complete domain-wall Hamiltonian and separate Euclidean GW control."""
from itertools import product
from pathlib import Path
import importlib.util
import numpy as np

spec=importlib.util.spec_from_file_location('fractional_input',Path(__file__).parent/'inputs/fractional_source.py')
fractional=importlib.util.module_from_spec(spec);spec.loader.exec_module(fractional)
BETA=fractional.BETA;ALPHA=fractional.ALPHA;G5=fractional.GAMMA5
A4=np.kron(fractional.SIGMA[1],np.eye(2))
PP=(np.eye(4)+G5)/2;PM=(np.eye(4)-G5)/2
T=-(BETA+1j*A4)/2
H_CLOCK=1/8

def boundary_map(S,internal=1):
    E=np.zeros((4*S,4),complex);E[:4]=PP;E[-4:]=PM
    return np.kron(E,np.eye(internal))

def boundary_mass(S,coefficient,internal=1):
    M=np.zeros((4*S,4*S),complex)
    M[:4,-4:]=coefficient*PP@BETA@PM
    M[-4:,:4]=M[:4,-4:].conj().T
    return np.kron(M,np.eye(internal))

def fifth(S,internal=1):
    H=np.zeros((4*S,4*S),complex)
    for s in range(S-1):
        H[4*s:4*s+4,4*s+4:4*s+8]=T
        H[4*s+4:4*s+8,4*s:4*s+4]=T.conj().T
    return np.kron(H,np.eye(internal))

def axial_wall(S,block=4):
    assert S%2==0
    return np.kron(np.diag([1]*(S//2)+[-1]*(S//2)),np.eye(block))

def bloch(S,k,mass=0,onsite_shift=0):
    spatial=sum((a*np.sin(x) for a,x in zip(ALPHA,k)),np.zeros((4,4),complex))
    spatial+=BETA*(sum(1-np.cos(k))+onsite_shift)
    return np.kron(np.eye(S),spatial)+fifth(S)+boundary_mass(S,H_CLOCK*mass)

def one_particle(S,links,mass=0,side=2):
    spatial=fractional.one_particle(0,H_CLOCK,links,side=side)
    n=len(spatial);cells=n//24
    H=np.kron(np.eye(S),spatial)
    for s in range(S-1):
        H[s*n:(s+1)*n,(s+1)*n:(s+2)*n]=np.kron(np.eye(cells),np.kron(T,np.eye(6)))
        H[(s+1)*n:(s+2)*n,s*n:(s+1)*n]=H[s*n:(s+1)*n,(s+1)*n:(s+2)*n].conj().T
    M=H_CLOCK*mass*np.kron(np.eye(cells),np.kron(PP@BETA@PM,np.eye(6)))
    H[:n,-n:]+=M;H[-n:,:n]+=M.conj().T
    return H

def spatial_overlap(theta):
    """One spatial U1 holonomy, both spin copies and both energy signs."""
    gamma=-1j*BETA@ALPHA[0]
    V=-np.cos(theta)*np.eye(4)+1j*gamma*np.sin(theta)
    h=BETA@(np.eye(4)+V)
    q=G5@(np.eye(4)-V)/2
    return h,q

def u3_matter(W):
    """u,d,e representations in the physical electron unit."""
    R=np.zeros((7,7),complex);R[:3,:3]=W/np.linalg.det(W)
    R[3:6,3:6]=W;R[6,6]=np.linalg.det(W)
    return R

def euclidean_wilson(links,side=2):
    sites=list(product(range(side),repeat=4));pos={x:i for i,x in enumerate(sites)}
    d=links.shape[-1];block=4*d
    D=4*np.eye(block*len(sites),dtype=complex)
    gamma=[-1j*BETA@a for a in ALPHA]+[BETA]
    for x,p in enumerate(sites):
        xx=slice(x*block,(x+1)*block)
        for mu in range(4):
            q=list(p);q[mu]=(q[mu]+1)%side;y=pos[tuple(q)]
            yy=slice(y*block,(y+1)*block);U=links[x,mu]
            D[xx,yy]+=-np.kron(np.eye(4)-gamma[mu],U)/2
            D[yy,xx]+=-np.kron(np.eye(4)+gamma[mu],U.conj().T)/2
    return D,np.kron(np.eye(len(sites)),np.kron(G5,np.eye(d)))

def overlap(DW,chirality):
    H=chirality@(DW-np.eye(len(DW)))
    ev,U=np.linalg.eigh(H)
    sign=(U*np.sign(ev))@U.conj().T
    D=np.eye(len(DW))+chirality@sign
    return D,-sign,float(np.min(abs(ev)))

def overlap_derivative(DW,chirality,dDW):
    """Exact Frechet formula for the sign on a fixed nonzero-gap sector."""
    ev,V=np.linalg.eigh(chirality@(DW-np.eye(len(DW))))
    numerator=np.sign(ev)[:,None]-np.sign(ev)[None,:]
    denominator=ev[:,None]-ev[None,:]
    divided=np.zeros_like(denominator)
    mask=numerator!=0
    divided[mask]=numerator[mask]/denominator[mask]
    return chirality@V@(divided*(V.conj().T@chirality@dDW@V))@V.conj().T
