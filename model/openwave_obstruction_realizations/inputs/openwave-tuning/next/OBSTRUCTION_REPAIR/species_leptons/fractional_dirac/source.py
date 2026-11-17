"""Declared24-component fractional-charge Wilson/Dirac quark kinematics.

Index=(Dirac upper/lower,physical spin,flavour u/d,colour). Link W is U(3):
R_u(W)=det(W)^-1 W, R_d(W)=W. Both bands and all CAR modes are source data.
"""
from itertools import product
import numpy as np
from scipy import sparse

I2=np.eye(2,dtype=complex)
SIGMA=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
BETA=np.kron(SIGMA[2],I2)
ALPHA=np.array([np.kron(SIGMA[0],a)for a in SIGMA])
GAMMA5=np.kron(SIGMA[0],I2)
SPIN=np.array([np.kron(I2,a)/2 for a in SIGMA])
PL=(np.eye(4)-GAMMA5)/2
TAUP=np.array([[0,1],[0,0]],complex)
Q=np.kron(np.eye(4),np.kron(np.diag([2/3,-1/3]),np.eye(3)))

def reps(W):
    out=np.zeros((6,6),complex)
    out[:3,:3]=W/np.linalg.det(W);out[3:,3:]=W
    return out

def local_gauge(W):return np.kron(np.eye(4),reps(W))

def current_matrices():
    internal=np.kron(TAUP,np.eye(3))
    vector=[np.kron(a,internal)for a in [np.eye(4),*ALPHA]]
    axial=[np.kron(a,internal)for a in [GAMMA5,*[2*x for x in SPIN]]]
    left=[(v-a)/2 for v,a in zip(vector,axial)]
    return vector,axial,left

def edges(side=2):
    cells=list(product(range(side),repeat=3));index={x:i for i,x in enumerate(cells)}
    out=[]
    for i,x in enumerate(cells):
        for axis in range(3):
            y=list(x);y[axis]=(y[axis]+1)%side
            out.append((i,index[tuple(y)],axis))
    return cells,out

def one_particle(mass,h,links,side=2,r=1):
    """Cell-clock matrix. Links are all positive oriented periodic edges.

    Each oriented edge contributes r beta/2 at both endpoints and
    -(r beta+i alpha_axis) R(W)/2 from source row x to target column y.
    Thus the free Bloch symbol is beta[h m+r sum(1-cos k)]+sum alpha sin k.
    """
    cells,es=edges(side);assert len(links)==len(es)
    dim=24*len(cells);H=np.zeros((dim,dim),complex)
    beta=np.kron(BETA,np.eye(6))
    for x in range(len(cells)):H[24*x:24*x+24,24*x:24*x+24]+=h*mass*beta
    for (x,y,a),W in zip(es,links):
        xx=slice(24*x,24*x+24);yy=slice(24*y,24*y+24)
        H[xx,xx]+=r*beta/2;H[yy,yy]+=r*beta/2
        K=-np.kron(r*BETA+1j*ALPHA[a],reps(W))/2
        H[xx,yy]+=K;H[yy,xx]+=K.conj().T
    return H

def open_small(mass,h,r=1):
    """Two cells, physicalspin up, onecolour, bothDirac coordinates/flavours.
    Identity link along z; complete8-mode invariant one-particle control.
    """
    be=np.kron(SIGMA[2],I2);al=np.kron(SIGMA[0],I2)
    H=np.zeros((8,8),complex)
    H[:4,:4]=(h*mass+r/2)*be;H[4:,4:]=(h*mass+r/2)*be
    H[:4,4:]=-(r*be+1j*al)/2;H[4:,:4]=H[:4,4:].conj().T
    return H

def second_quantize(h):
    """Literal number-conserving CAR operator on every occupation."""
    n=len(h);dim=1<<n;rows=[];cols=[];values=[]
    nz=list(zip(*np.nonzero(h)))
    for w in range(dim):
        for a,b in nz:
            a,b=int(a),int(b)
            if not(w>>b)&1:continue
            z=w^(1<<b);sgn=(-1)**((w&((1<<b)-1)).bit_count())
            if(z>>a)&1:continue
            sgn*=(-1)**((z&((1<<a)-1)).bit_count());v=z|(1<<a)
            rows.append(v);cols.append(w);values.append(sgn*h[a,b])
    return sparse.csr_matrix((values,(rows,cols)),shape=(dim,dim))
