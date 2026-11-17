"""Unfitted finite local auxiliary source, explicitly distinct from two-cell vacuum.

Quark modes are3*alpha+colour with alpha=(up u,up d,down u,down d).
Full fixed-Nq=3 CAR basis220; auxiliary0 plus(a=0..7,i=0..2), dimension25.
A vector is a(25,220) array, auxiliary index1+3*a+i.
"""
from itertools import combinations
from pathlib import Path
from fractions import Fraction
import json
import numpy as np
from scipy import sparse
from scipy.linalg import cholesky,solve_triangular
from staggered_vacuum import creation,annihilation
from staggered_chain4 import sector
HERE=Path(__file__).resolve().parent

def one_body(words,h):
    ix={w:i for i,w in enumerate(words)};rows=[];cols=[];data=[]
    nz=list(zip(*np.nonzero(h)))
    for j,w in enumerate(words):
        for a,b in nz:
            z,s=annihilation(int(b),w)
            if z is None:continue
            y,t=creation(int(a),z)
            if y is not None:rows.append(ix[y]);cols.append(j);data.append(h[a,b]*s*t)
    return sparse.csr_matrix((data,(rows,cols)),shape=(len(words),)*2,dtype=complex)

def generators():
    T=[]
    for a,b in [(0,1),(0,2),(1,2)]:
        x=np.zeros((3,3),complex);x[a,b]=x[b,a]=.5;T.append(x)
        y=np.zeros((3,3),complex);y[a,b]=-.5j;y[b,a]=.5j;T.append(y)
    T.append(np.diag([.5,-.5,0]).astype(complex))
    T.append(np.diag([1,1,-2]).astype(complex)/(2*np.sqrt(3)))
    S=[np.array([[0,1],[1,0]],complex)/2,np.array([[0,-1j],[1j,0]],complex)/2,np.diag([.5,-.5]).astype(complex)]
    return T,S

def source():
    words=[sum(1<<i for i in p) for p in combinations(range(12),3)];words.sort();ix={w:i for i,w in enumerate(words)}
    T,S=generators();G=[];colour=[];spin=[]
    for a in range(8):
        colour.append(one_body(words,np.kron(np.eye(4),T[a])))
        for i in range(3):G.append(one_body(words,np.kron(np.kron(S[i],np.eye(2)),T[a])))
    for i in range(3):spin.append(one_body(words,np.kron(np.kron(S[i],np.eye(2)),np.eye(3))))
    tau=np.array([[0,1],[0,0]],complex)
    V=one_body(words,np.kron(np.kron(np.eye(2),tau),np.eye(3)))
    A=[one_body(words,np.kron(np.kron(2*S[i],tau),np.eye(3))) for i in range(3)]
    f=np.array([[[-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c]) for c in range(8)] for b in range(8)] for a in range(8)]).real
    adj=[-1j*f[a] for a in range(8)]
    eps=np.zeros((3,3,3));eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1;eps[0,2,1]=eps[2,1,0]=eps[1,0,2]=-1
    vector=[-1j*eps[i] for i in range(3)]
    c=sector(4,1);integer=np.zeros((220,20),dtype=np.int64)
    for j,w in enumerate(c['words']):integer[ix[w]]=c['R'][j]
    gram=integer.T@integer;L=cholesky(gram.astype(float),lower=True)
    R=solve_triangular(L,integer.T,lower=True).T
    S2=sum((J@J for J in spin),sparse.csr_matrix((220,220),dtype=complex))
    ss=R.T@S2@R;P_half=(3.75*np.eye(20)-ss)/3;P_threehalf=(ss-.75*np.eye(20))/3
    assert np.linalg.norm(R.T@R-np.eye(20))<1e-12
    assert np.linalg.norm(P_half@P_half-P_half)<1e-12
    assert abs(np.trace(P_half)-4)<1e-12
    return dict(words=words,colour=colour,spin=spin,G=G,V=V,A=A,adjoint_colour=adj,vector_spin=vector,
       integer_singlet_basis=integer,singlet_basis=R,half_projector=P_half,threehalf_projector=P_threehalf)

def hamiltonian_action(src,kappa,psi):
    v=np.asarray(psi).reshape(25,220);out=np.zeros_like(v,dtype=complex);Delta=3*float(kappa);g=1/8
    for j,G in enumerate(src['G']):out[0]+=g*(G@v[j+1]);out[j+1]=Delta*v[j+1]+g*(G@v[0])
    return out.reshape(np.asarray(psi).shape)

def lower_isometry(src,kappa):
    Delta=3*float(kappa);g=1/8;R=src['singlet_basis'];P=src['half_projector'];Q=src['threehalf_projector']
    result=np.zeros((25,220,20),dtype=complex);E=np.zeros((20,20),dtype=complex)
    for lam,projector in [(4,P),(2,Q)]:
        disc=np.sqrt(Delta*Delta+4*g*g*lam);c=np.sqrt((1+Delta/disc)/2);s=np.sqrt((1-Delta/disc)/2)
        result[0]+=c*R@projector;E+=(Delta-disc)/2*projector
        for j,G in enumerate(src['G']):result[j+1]-=s/np.sqrt(lam)*(G@R@projector)
    return result,E

def total_generator(src,quark,auxiliary,psi):
    v=np.asarray(psi).reshape(25,220);out=np.stack([quark@r for r in v]);out[1:]+=auxiliary@v[1:]
    return out.reshape(np.asarray(psi).shape)

def checks():
    src=source();R=src['singlet_basis'];results=[]
    for kappa in [.5,1,2]:
        D,E=lower_isometry(src,kappa);flat=D.reshape(5500,20)
        iso=float(np.linalg.norm(flat.conj().T@flat-np.eye(20)))
        residual=max(float(np.linalg.norm(hamiltonian_action(src,kappa,D[:,:,j])-(D@E)[:,:,j])) for j in range(20))
        gauss=0.;rotation=0.
        for a in range(8):
            aux=np.kron(src['adjoint_colour'][a],np.eye(3))
            gauss=max(gauss,max(float(np.linalg.norm(total_generator(src,src['colour'][a],aux,D[:,:,j]))) for j in range(20)))
        for i in range(3):
            aux=np.kron(np.eye(8),src['vector_spin'][i]);reduced=R.T@src['spin'][i]@R
            rotation=max(rotation,max(float(np.linalg.norm(total_generator(src,src['spin'][i],aux,D[:,:,j])-(D@reduced)[:,:,j])) for j in range(20)))
        assert max(iso,residual,gauss,rotation)<1e-11
        results.append(dict(kappa=kappa,isometry_residual=iso,H_state_residual=residual,total_Gauss_residual=gauss,total_spin_intertwiner_residual=rotation))
    out=dict(schema='oph.local_auxiliary.interface.v1',full_fixed_N3_quark_dimension=220,auxiliary_dimension=25,
        full_fixed_N3_carrier_dimension=5500,physical_Gauss_sector_dimension=80,lower_dressed_dimension=20,
        H='Delta Pi1 + g sum_ai (|ai><0|+|0><ai|) G_ai',Delta='3*kappa',g='1/8',
        source_is_distinct_from_two_cell_vacuum=True,common_Dirac_transfer=False,physical_axial_identification=False,checks=results)
    (HERE/'LOCAL_AUXILIARY_INTERFACE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':checks()
