import importlib.util,itertools
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('old_G1_group',HERE.parent/'cubic_reduction/group.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def permutation_parity(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def sign_character(R):
    diagonals=[(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1)]
    index={x:i for i,x in enumerate(diagonals)}
    perm=[]
    for x in diagonals:
        y=R@np.array(x)
        if y[0]<0:y=-y
        perm.append(index[tuple(y)])
    return permutation_parity(perm)

def symmetric_cube(S):
    B=np.zeros((8,4),complex)
    for col in range(4):
        rows=[i for i in range(8) if i.bit_count()==col]
        B[rows,col]=1/np.sqrt(len(rows))
    return B.conj().T@np.kron(np.kron(S,S),S)@B

def rotations():
    return [(R,S,T,{'G1':S,'G2':sign_character(R)*S,'H':symmetric_cube(S)}) for R,S,T in old.rotations()]

def orbit_bases(cells):
    group=rotations();index={tuple(x):i for i,x in enumerate(cells)}
    remaining=set(index);result=[]
    while remaining:
        seed=np.array(min(remaining));orbit=sorted({tuple(R@seed) for R,_,_,_ in group})
        remaining.difference_update(orbit);local_index={x:i for i,x in enumerate(orbit)}
        dim=4*len(orbit);sectors={}
        for name,d in [('G1',2),('G2',2),('H',4)]:
            units=[np.zeros((dim,dim),complex) for _ in range(d)]
            for R,_,T,reps in group:
                rep=reps[name]
                for x,j in local_index.items():
                    i=local_index[tuple(R@np.array(x))]
                    for component in range(d):
                        units[component][4*i:4*i+4,4*j:4*j+4]+=(d/24)*rep[component,0].conj()*T
            P11=units[0]
            ev,u=np.linalg.eigh((P11+P11.conj().T)/2)
            B1=u[:,ev>.5]
            bases=[B1]+[units[j]@B1 for j in range(1,d)]
            sectors[name]=dict(d=d,bases=bases,multiplicity=B1.shape[1],idempotence=float(np.linalg.norm(P11@P11-P11)))
        assembled=np.hstack([b for sector in sectors.values() for b in sector['bases']])
        complete_error=float(np.linalg.norm(assembled@assembled.conj().T-np.eye(dim)))
        global_cells=np.array([index[x] for x in orbit]);rows=(4*global_cells[:,None]+np.arange(4)).ravel()
        result.append(dict(cells=np.array(orbit),rows=rows,sectors=sectors,complete_error=complete_error))
    return result

def full_small_bases(cells):
    orbits=orbit_bases(cells);out={}
    for name,d in [('G1',2),('G2',2),('H',4)]:
        rank=sum(o['sectors'][name]['multiplicity'] for o in orbits)
        bases=[np.zeros((4*len(cells),rank),complex) for _ in range(d)]
        at=0
        for o in orbits:
            n=o['sectors'][name]['multiplicity']
            for j in range(d):bases[j][np.ix_(o['rows'],np.arange(at,at+n))]=o['sectors'][name]['bases'][j]
            at+=n
        out[name]=bases
    return out,orbits
