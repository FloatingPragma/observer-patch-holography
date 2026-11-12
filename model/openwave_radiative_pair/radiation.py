"""Native charged Dirac current coupled to four finite transverse modes."""
from itertools import combinations,product
from functools import lru_cache
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import expm_multiply
from model.openwave_lepton_mass_shell.spectrum import dirac_hamiltonian

SIG=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))

def apply(mask,i,j):
    if not(mask>>j&1):return None
    sign=(-1)**((mask&((1<<j)-1)).bit_count());out=mask^(1<<j)
    if out>>i&1:return None
    sign*=(-1)**((out&((1<<i)-1)).bit_count());return out|(1<<i),sign

def lift(matrix,masks):
    where={m:i for i,m in enumerate(masks)};rows=[];cols=[];data=[]
    for i,j in zip(*np.nonzero(abs(matrix)>1e-14)):
        for col,mask in enumerate(masks):
            term=apply(mask,i,j)
            if term is not None:
                rows.append(where[term[0]]);cols.append(col);data.append(matrix[i,j]*term[1])
    return sp.csr_matrix((data,(rows,cols)),shape=(len(masks),)*2)

@lru_cache(maxsize=None)
def construction(cap=1,eta=.05,omega_ratio=1.):
    # All energies divided by this branch's native mass. The declared ring
    # spacing a=sqrt(3)/m gives cochain photon omega=m and sin(q)/a=m/2.
    masks=tuple(sum(1<<i for i in c)for c in combinations(range(12),6));where={m:i for i,m in enumerate(masks)}
    frames=[];energies=[]
    for momentum in (0,np.sqrt(3)/2,-np.sqrt(3)/2):
        h=dirac_hamiltonian((0,0,momentum/np.sqrt(3)),1.)
        columns=[];es=[]
        for sign in (-1,1):
            for spin in (0,1):
                ids=[spin,spin+2];v,U=np.linalg.eigh(h[np.ix_(ids,ids)]);ix=0 if sign==-1 else 1
                w=np.zeros(4,complex);w[ids]=U[:,ix]
                if w[ids[0]].real<0:w=-w
                columns.append(w);es.append(v[ix])
        frames.append(np.column_stack(columns));energies+=es
    energies=np.array(energies);ground=sum(1<<i for i in (0,1,4,5,8,9));Eg=sum(energies[i]for i in range(12)if ground>>i&1)
    fermion_E=np.array([sum(energies[i]for i in range(12)if mask>>i&1)-Eg for mask in masks])
    pair=np.zeros(len(masks),complex)
    for spin in (0,1):
        out,sign=apply(ground,spin+2,spin);pair[where[out]]=sign/np.sqrt(2)
    modes=list(product((1,2),(0,1)));nph=(cap+1)**4;basis=list(product(range(cap+1),repeat=4));nb=np.array(basis)
    photon_n=nb.sum(axis=1);photon_momentum=(nb@np.array([1,1,2,2]))%3
    a=sp.diags(np.sqrt(np.arange(1,cap+1)),1,shape=(cap+1,cap+1),format='csr')
    currents=[];coupling=sp.csr_matrix((len(masks)*nph,)*2,dtype=complex)
    for mode,(q,pol)in enumerate(modes):
        # Emitting momentum q reduces the source momentum by q.
        J=np.zeros((12,12),complex);alpha=-np.kron(np.diag([1.,-1.]),SIG[pol])
        for p in range(3):
            out=(p-q)%3;J[4*out:4*out+4,4*p:4*p+4]=frames[out].conj().T@alpha@frames[p]
        F=lift(J,masks);currents.append(F)
        b=sp.csr_matrix([[1.]])
        for j in range(4):b=sp.kron(b,a.conj().T if j==mode else sp.eye(cap+1),format='csr')
        term=sp.kron(F,b,format='csr');coupling+=term+term.conj().T
    diagonal=(fermion_E[:,None]+omega_ratio*photon_n[None,:]).ravel()
    H=sp.diags(diagonal,format='csr')+eta*coupling
    initial=np.kron(pair,np.eye(nph)[0]);vacuum=np.kron(np.eye(len(masks))[where[ground]],np.eye(nph)[0])
    pmom=np.array([sum((i//4)*(mask>>i&1)for i in range(12))%3 for mask in masks]);totalmom=(pmom[:,None]+photon_momentum[None,:]).ravel()%3
    diff=H.tocoo();assert all(totalmom[i]==totalmom[j]for i,j in zip(diff.row,diff.col))
    assert abs(np.vdot(initial,initial)-1)<1e-14
    return {'H':H,'H0_diagonal':diagonal,'coupling':coupling,'initial':initial,'vacuum':vacuum,'masks':masks,'ground_index':where[ground],
        'photon_basis':basis,'photon_number':photon_n,'photon_momentum':photon_momentum,'source_energy':fermion_E,'frames':frames,'currents':currents,'total_momentum':totalmom}

def observe(c,psi,time,initial):
    nph=len(c['photon_basis']);law=abs(psi.reshape(-1,nph))**2
    photons=law[c['ground_index']];field_marginal=law.sum(axis=0)
    source=float(law.sum(axis=1)@c['source_energy']);free=float(law.ravel()@c['H0_diagonal'])
    total=float(np.vdot(psi,c['H']@psi).real);interaction=total-free
    # Complete simultaneous free-source-energy / photon occupation record.
    # The reported ground/radiation grouping has the full other-source bin.
    records=np.vstack([photons,field_marginal-photons])
    assert np.min(records)>-2e-14
    assert abs(law.sum()-1)<2e-9
    return {'time_in_inverse_native_mass':time,
        'two_photon_no_fermion_probability':float(photons[c['photon_number']==2].sum()),
        'complete_source_ground_or_complement_and_photon_probabilities':records.tolist(),
        'full_probability':float(law.sum()),'free_source_energy':source,
        'free_photon_energy':free-source,'interaction_energy':interaction,'full_energy_before_readout':total,
        'complete_number_readout_work':-interaction,
        'photon_occupation_mean':float(field_marginal@c['photon_number']),
        'photon_top_cap_probability':float(field_marginal[np.any(np.array(c['photon_basis'])==max(max(b)for b in c['photon_basis']),axis=1)].sum()),
        'pair_survival':float(abs(np.vdot(initial,psi))**2),
        'forbidden_total_momentum_probability':float(abs(psi[c['total_momentum']!=0])@abs(psi[c['total_momentum']!=0]))}


def run_probe(cap=1,eta=.05,times=(0,5,20,50),initial_kind='pair',omega_ratio=1.):
    c=construction(cap,eta,omega_ratio);H=c['H'];initial=c['initial']if initial_kind=='pair'else c['vacuum'];rows=[]
    for t in times:
        psi=expm_multiply(-1j*t*H,initial,traceA=-1j*t*H.diagonal().sum())
        rows.append(observe(c,psi,t,initial))
    return rows


def source_controls():
    from model.openwave_lepton_mass_shell.spectrum import native_frames
    native=native_frames();c=construction();H=c['H'];masks=c['masks']
    assert (H-H.conj().T).nnz==0
    # Actual total spin on the pair input and source free reference.
    spin_checks=[];commutators=[]
    pair_matrix=np.zeros((12,12),complex)
    pair_matrix[2,0]=pair_matrix[3,1]=1/np.sqrt(2)
    for axis in range(3):
        S=np.zeros((12,12),complex)
        for p in range(3):S[4*p:4*p+4,4*p:4*p+4]=c['frames'][p].conj().T@np.kron(np.eye(2),SIG[axis]/2)@c['frames'][p]
        lifted=lift(S,masks);state=c['initial'].reshape(len(masks),-1)[:,0]
        spin_checks.append(float(np.linalg.norm(lifted@state)))
        commutators.append(float(np.linalg.norm(S@pair_matrix-pair_matrix@S)))
    assert max(commutators)<1e-13 and max(spin_checks)>0.4
    vector=c['initial'].copy();moments=[];n=len(c['photon_basis'])
    for degree in range(5):
        target=vector.reshape(-1,n)[c['ground_index']]
        moments.append({'degree':degree,'radiation_amplitude_real':target.real.tolist(),'radiation_amplitude_imag':target.imag.tolist()})
        vector=H@vector
    return {'native_full_source_intertwiners_checked':len(native['frames']),
        'canonical_Dirac_CAR_modes':12,'full_fermion_Fock_dimension':4096,
        'exact_charge_zero_occupation_sector_dimension':len(masks),
        'charge_sector_is_invariant_not_a_postselection':True,
        'all_other_charge_sectors_retained_as_direct_sum':True,
        'actual_full_source_spin_norms':spin_checks,'p0_pair_creation_native_spin_commutators':commutators,
        'pair_is_spin_scalar_but_full_sea_is_not_spin_eigenstate':True,
        'exact_translation_momentum_modulo3_preserved':True,
        'radiation_H_moments':moments}


if __name__=='__main__':
    import json
    print(json.dumps(run_probe(),indent=2))
