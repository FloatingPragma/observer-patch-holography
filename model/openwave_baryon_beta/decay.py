"""No emitted-pair cap: all charge-admitted three-quark/lepton states."""
from functools import lru_cache
from itertools import combinations,product
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_calibrated_hadrons.source import sector,parameters,INPUTS
from model.openwave_baryon_weak_formfactors.native import currents,quartet,PAULI
from model.openwave_native_composites.spectrum import one_body
from model.common_source import calibrated_matter as internal

G=1/64

def index(w,x,s):return 4*w+2*x+s

def wedge(U,masks):
    occupied=[[i for i in range(8)if m>>i&1]for m in masks]
    return np.array([[np.linalg.det(U[np.ix_(a,b)])for b in occupied]for a in occupied])

@lru_cache(None)
def construction(coupling=G,kinetic=None):
    d=sector('baryon');c=currents();masks=[sum(1<<i for i in row)for row in combinations(range(8),4)]
    reference=masks.index(15);ne=np.array([(m>>4).bit_count()for m in masks])
    lz=np.array([sum((1-2*s)for w,x,s in product(range(2),repeat=3)if m>>index(w,x,s)&1)for m in masks])
    qh=np.rint(np.diag(d['q']).real).astype(int);zh=np.rint(2*np.diag(d['z']).real).astype(int)
    charge=(qh[:,None]-ne).ravel();spin=(zh[:,None]+lz).ravel()
    keep=np.flatnonzero((charge==0)&(spin==1));outside=np.flatnonzero((charge!=0)|(spin!=1))
    q=float((parameters()['delta']-1)/2)if kinetic is None else float(kinetic)
    assert q>0
    lep=np.zeros((8,8));K=q/2*np.array([[1.,-1.],[-1.,1.]])
    for w,s in product(range(2),repeat=2):
        ids=[index(w,x,s)for x in range(2)]
        lep[np.ix_(ids,ids)]=(np.eye(2)+K)if w else -K
    Hl=one_body(lep,masks)+2*q*sp.eye(70)
    raising=np.array(internal.weak_raising_current(families=False)).astype(complex)[np.ix_([10,11],[10,11])]
    third=(raising@raising.conj().T-raising.conj().T@raising)/2
    V=sp.csr_matrix((8400,8400),dtype=complex)
    for x in range(2):
        for mu,sigma in enumerate(PAULI):
            pair=[]
            for weak in (raising,third):
                A=np.zeros((8,8),complex)
                for w,v,s,z in product(range(2),repeat=4):A[index(w,x,s),index(v,x,z)]=sigma[s,z]*weak[w,v]
                pair.append(one_body(A,masks))
            lp,l3=pair;qp,q3=map(sp.csr_matrix,c['local'][x][mu])
            V+=(1 if mu==0 else -1)*(sp.kron(qp,lp.conj().T)+sp.kron(qp.conj().T,lp)+2*sp.kron(q3,l3))
    H0=(sp.kron(sp.csr_matrix(d['h']),sp.eye(70))+sp.kron(sp.eye(120),Hl)).tocsr();V=V.tocsr()
    assert np.max(abs(H0[keep][:,outside].data),initial=0)<1e-11
    assert np.max(abs(V[keep][:,outside].data),initial=0)<1e-12
    H0=H0[keep][:,keep];V=V[keep][:,keep];H=(H0+coupling*V).tocsr()
    W,_=quartet(d['h']);initial_full=np.kron(W[:,0],np.eye(70)[:,reference]);initial=initial_full[keep]
    assert abs(np.vdot(initial,initial)-1)<1e-12
    U=np.zeros((8,8))
    for w,s,p,x in product(range(2),repeat=4):U[index(w,p,s),index(w,x,s)]=(-1)**(p*x)/np.sqrt(2)
    F=wedge(U,masks);assert np.max(abs(F@F.T-np.eye(70)))<2e-14
    T=sp.kron(sp.eye(120),sp.csr_matrix(F)).tocsr()[keep][:,keep]
    assert np.max(abs((T@T.T-sp.eye(len(keep))).data),initial=0)<2e-14
    measured_h=np.asarray(abs(d['B']).power(2).T@d['H'].diagonal().real).ravel()
    measured=(np.repeat(measured_h,70)+np.tile(Hl.diagonal().real,120))[keep]+coupling*V.diagonal().real
    measured_momentum=(np.repeat(measured_h,70)+np.tile(np.diag(F@Hl@F.T).real,120))[keep]+coupling*(T@V@T.T).diagonal().real
    seed_index=int(np.argmax(abs(initial)**2));seed=np.eye(len(keep),dtype=complex)[:,seed_index]
    success=float(abs(initial[seed_index])**2);failure=seed-initial*np.vdot(initial,seed)
    assert abs(np.vdot(failure,failure).real+success-1)<1e-12
    return dict(H=H,H0=H0,V=V,initial=initial,keep=keep,masks=masks,ne=ne[keep%70],
                neutron=W[:,:2],proton=W[:,2:],native=d,T=T,kinetic=q,coupling=coupling,
                source_dimension=8400,charge_zero_dimension=int(np.sum(charge==0)),
                source_spin_dimension=len(keep),reference=reference,measured_energy=measured,
                measured_momentum_energy=measured_momentum,
                seed_index=seed_index,seed=seed,success=success,failure=failure)

def eigensystem(d):
    offset=float(np.vdot(d['initial'],d['H']@d['initial']).real)
    e,U=eigh(d['H'].toarray()-offset*np.eye(len(d['keep'])))
    return e,U,offset

def evolve(d,time,eigen=None,state=None):
    e,U,_=eigensystem(d)if eigen is None else eigen
    v=d['initial']if state is None else state
    return U@(np.exp(-1j*time*e)*(U.conj().T@v))

def records(d,time,eigen=None):
    psi=evolve(d,time,eigen);position=abs(psi)**2;momentum=abs(d['T']@psi)**2
    assert abs(position.sum()-1)<1e-11 and abs(momentum.sum()-1)<1e-11
    full=np.zeros(8400,complex);full[d['keep']]=psi;full=full.reshape(120,70)
    proton=float(np.linalg.norm(d['proton'].conj().T@full)**2)
    energy=float(np.vdot(psi,d['H']@psi).real);initial=float(np.vdot(d['initial'],d['H']@d['initial']).real)
    assert abs(energy-initial)<2e-8
    by_n=[float(position[d['ne']==n].sum())for n in range(3)]
    return {'source_parameter':float(time),'conditional_electron_count_probability':by_n,
            'conditional_proton_band_probability':proton,'preparation_success_probability':d['success'],
            'original_electron_count_probability':[d['success']*p for p in by_n],
            'original_failure_probability':1-d['success'],'source_energy':energy,
            'position_readout_work':float(position@d['measured_energy']-energy),
            'momentum_readout_work':float(momentum@d['measured_momentum_energy']-energy),
            'position_probability':position.tolist(),'momentum_probability':momentum.tolist()}
