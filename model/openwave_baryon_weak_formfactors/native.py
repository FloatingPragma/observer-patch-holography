"""Literal weak and spin currents on the full native three-quark sector."""
from functools import lru_cache
from itertools import product
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_baryon_ordering.source import balanced_baryon
from model.openwave_native_composites.spectrum import one_body
from model.common_source import calibrated_matter as internal
from model.coherent_defect_wave.physical_attachment.chiral import representation
PAULI=(np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))

@lru_cache(None)
def currents():
 a=balanced_baryon();old=a['old'];B=a['B'];modes=a['modes'];index={v:i for i,v in enumerate(modes)}
 subsets,_,_,_=representation();raising=np.array(internal.weak_raising_current(families=False)).astype(complex)
 t3=(raising@raising.conj().T-raising.conj().T@raising)/2
 qi=[[subsets.index((c,w))for c in range(3)]for w in(3,4)];local=[];full=[]
 for x in range(2):
  row=[];frow=[]
  for sigma in PAULI:
   operators=[];fullops=[]
   for weak in (raising,t3):
    A=np.zeros((24,24),complex)
    for s,z,w,v,c in product(range(2),range(2),range(2),range(2),range(3)):
     A[index[(0,x,s,w,c)],index[(0,x,z,v,c)]]=sigma[s,z]*weak[qi[w][c],qi[v][c]]
    J=one_body(A,a['basis']);j=(B.conj().T@J@B).toarray()
    assert np.max(abs(J@B-B@j))<3e-13
    operators.append(j);fullops.append(J)
   row.append(operators);frow.append(fullops)
  local.append(row);full.append(frow)
 spin_lower=np.zeros((24,24))
 for x,w,c in product(range(2),range(2),range(3)):spin_lower[index[(0,x,1,w,c)],index[(0,x,0,w,c)]]=1
 Sminus=(B.conj().T@one_body(spin_lower,a['basis'])@B).toarray()
 labels=old['labels'];where={v:i for i,v in enumerate(labels)};parity=np.zeros((120,120))
 for i,label in enumerate(labels):parity[where[tuple(sorted((1-x,s,w)for x,s,w in label))],i]=1
 # Second-quantized site exchange is checked before using its reduced form.
 lookup={b:i for i,b in enumerate(a['basis'])};rr=[];cc=[];vv=[]
 for col,mask in enumerate(a['basis']):
  occupied=[index[(species,1-x,s,w,c)]for i,(species,x,s,w,c)in enumerate(modes)if mask>>i&1]
  sign=(-1)**sum(occupied[i]>occupied[j]for i in range(3)for j in range(i+1,3))
  rr.append(lookup[sum(1<<i for i in occupied)]);cc.append(col);vv.append(sign)
 Pfull=sp.csr_matrix((vv,(rr,cc)),shape=(2024,2024));assert np.max(abs(Pfull@B-B@parity))<3e-13
 total=[sum(local[x][mu][0]for x in range(2))for mu in range(4)]
 return {'native':a,'local':local,'full':full,'total_raising':total,'spin_lower':Sminus,'parity':parity,'full_parity':Pfull}

def quartet(h,parity_sign=1):
 d=currents();old=d['native']['old'];q=np.diag(old['q']).real;z=np.diag(old['z']).real
 ids=np.flatnonzero((abs(q)<1e-10)&(abs(z-.5)<1e-10));P=d['parity'][np.ix_(ids,ids)]
 e,U=eigh(P);basis=U[:,abs(e-parity_sign)<1e-10];energies,v=eigh(basis.conj().T@h[np.ix_(ids,ids)]@basis)
 nup=np.zeros(120,complex);nup[ids]=basis@v[:,0]
 # Fix a deterministic phase, then generate all partners by literal currents.
 k=int(np.argmax(abs(nup)));nup*=np.conj(nup[k])/abs(nup[k])
 ndown=d['spin_lower']@nup;pup=d['total_raising'][0]@nup;pdown=d['spin_lower']@pup
 W=np.column_stack((nup,ndown,pup,pdown));assert np.max(abs(W.conj().T@W-np.eye(4)))<2e-11
 assert np.linalg.norm(old['spin']@W-.75*W)<2e-9
 for i,target in enumerate((0,0,1,1)):assert np.linalg.norm(old['q']@W[:,i]-target*W[:,i])<1e-11
 return W,energies

def formfactors(h):
 d=currents();W,energies=quartet(h);n=W[:,:2];p=W[:,2:];blocks=[p.conj().T@J@n for J in d['total_raising']]
 gv=np.trace(blocks[0]).real/2;ga=np.trace(PAULI[3]@blocks[3]).real/2
 assert abs(gv-1)<1e-11 and np.max(abs(blocks[0]-np.eye(2)))<1e-11
 for mu in range(1,4):assert np.max(abs(blocks[mu]-ga*PAULI[mu]))<1e-11
 leakage=[float(np.linalg.norm(J@n-p@b,2))for J,b in zip(d['total_raising'],blocks)]
 return {'gV':float(gv),'gA':float(ga),'quartet':W,'charge_spin_ground_energies':energies,'current_blocks':blocks,'current_leakage_operator_norm':leakage}
