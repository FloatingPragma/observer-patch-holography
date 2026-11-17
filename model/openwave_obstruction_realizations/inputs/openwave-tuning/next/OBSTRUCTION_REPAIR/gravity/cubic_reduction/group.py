"""Exact projective cubic G1 matrix units, retaining every multiplicity."""
import itertools
import numpy as np
from scipy.spatial.transform import Rotation

X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.])

def rotations():
    out=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1,1),repeat=3):
            R=np.diag(signs)@np.eye(3,dtype=int)[list(perm)]
            if round(np.linalg.det(R))!=1:continue
            q=Rotation.from_matrix(R).as_quat()
            spin=q[3]*np.eye(2)-1j*(q[0]*X+q[1]*Y+q[2]*Z)
            out.append((R,spin,np.kron(spin,np.eye(2))))
    return out

def orbit_bases(cells):
    group=rotations();index={tuple(x):i for i,x in enumerate(cells)}
    remaining=set(index);result=[]
    while remaining:
        seed=np.array(min(remaining));orbit=sorted({tuple(R@seed) for R,_,_ in group})
        remaining.difference_update(orbit);local_index={x:i for i,x in enumerate(orbit)}
        dim=4*len(orbit);P11=np.zeros((dim,dim),complex);P21=np.zeros_like(P11)
        for R,spin,S in group:
            for x,j in local_index.items():
                i=local_index[tuple(R@np.array(x))]
                P11[4*i:4*i+4,4*j:4*j+4]+=(2/24)*spin[0,0].conj()*S
                P21[4*i:4*i+4,4*j:4*j+4]+=(2/24)*spin[1,0].conj()*S
        ev,u=np.linalg.eigh((P11+P11.conj().T)/2)
        B1=u[:,ev>.5];B2=P21@B1
        global_cells=np.array([index[x] for x in orbit])
        rows=(4*global_cells[:,None]+np.arange(4)).ravel()
        result.append(dict(cells=np.array(orbit),rows=rows,B1=B1,B2=B2,
                     P11_idempotence=float(np.max(abs(P11@P11-P11))),
                     component2_gram_error=float(np.max(abs(B2.conj().T@B2-np.eye(B2.shape[1])))) if B2.shape[1] else 0.))
    return result

def full_small_bases(cells):
    orbits=orbit_bases(cells);rank=sum(o['B1'].shape[1] for o in orbits)
    B1=np.zeros((4*len(cells),rank),complex);B2=np.zeros_like(B1);at=0
    for o in orbits:
        n=o['B1'].shape[1]
        B1[np.ix_(o['rows'],np.arange(at,at+n))]=o['B1']
        B2[np.ix_(o['rows'],np.arange(at,at+n))]=o['B2'];at+=n
    return B1,B2,orbits
