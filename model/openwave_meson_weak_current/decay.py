"""Full native semileptonic current on exact colour and charge sectors."""
from functools import lru_cache
from itertools import product
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
import sympy as sy
from model.openwave_native_composites.spectrum import sector,one_body,norm
from model.native_flavour_transport.transport import source,word_unitary,numeric
from model.common_source import calibrated_matter as internal
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_weak_decay.decay import mass_data

G=1/64
PAULI=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
BAR=(np.eye(2),*(-p for p in PAULI))
ETA=(1,-1,-1,-1)

def _lepton_index(weak,x,spin,f):return weak*12+x*6+spin*3+f

def _restricted_current(a,masks):
    lookup={b:i for i,b in enumerate(masks)};out=np.zeros((len(masks),len(masks)),complex)
    for j,mask in enumerate(masks):
        for p,q in zip(*np.nonzero(abs(a)>1e-14)):
            p,q=int(p),int(q)
            if not (mask>>q)&1:continue
            before=mask^(1<<q)
            if(before>>p)&1:continue
            after=before|(1<<p)
            if after not in lookup:continue
            sign=(-1)**((mask&((1<<q)-1)).bit_count()+(before&((1<<p)-1)).bit_count())
            out[lookup[after],j]+=sign*a[p,q]
    return out

@lru_cache(maxsize=6)
def construction(zero_word=False,contact_resonance=True):
    native=sector('meson');subsets,_,_,_=representation();raising=numeric(internal.weak_raising_current(families=False))
    qi=[[subsets.index((c,w))for c in range(3)]for w in(3,4)]
    active=[i for i,(q,a)in enumerate(native['labels'])if a[2]==0]
    oldlabels=[native['labels'][i]for i in active];h=native['h'][np.ix_(active,active)]
    assert len(active)==32
    mlight,mheavy,_=mass_data();mass=np.array([1.,1/mlight,mheavy/mlight]);q=(mass[1]-2)/2-(11*G/4 if contact_resonance else 0)
    ev,F=eigh(numeric(source()['y']));C=np.eye(3,dtype=complex)if zero_word else F.conj().T@numeric(word_unitary((2,3,10)))@F
    assert np.max(abs(C.conj().T@C-np.eye(3)))<2e-13
    hadron_labels=[(f,*qq[:2],qq[2],*aa[:2])for f in range(3)for qq,aa in oldlabels]
    Hh=sp.block_diag([sp.csr_matrix(h+(mf-1)*np.eye(32))for mf in mass],format='csr')
    full_diagonal=np.asarray(native['H'].diagonal()).real
    measured_hadron=np.asarray(abs(native['B'][:,active]).power(2).T@full_diagonal).ravel()
    measured_hadron=np.concatenate([measured_hadron+(mf-1)for mf in mass])
    qops=[]
    for x in range(2):
        row=[]
        for bar in BAR:
            aq=np.zeros((48,48),complex)
            for i,(species,cell,spin,w,c)in enumerate(native['modes']):
                if species or cell!=x:continue
                for j,(other,site,ss,ww,cc)in enumerate(native['modes']):
                    if other or site!=x:continue
                    aq[i,j]=bar[spin,ss]*raising[qi[w][c],qi[ww][cc]]
            full=one_body(aq,native['basis']);reduced=(native['B'].conj().T@full@native['B']).toarray()
            assert norm(full@native['B']-native['B']@sp.csr_matrix(reduced))<2e-13
            plus=reduced[np.ix_(active,active)]
            # The diagonal native current follows from the same full generator.
            a3=np.zeros((48,48),complex)
            T3=(raising@raising.conj().T-raising.conj().T@raising)/2
            for i,(species,cell,spin,w,c)in enumerate(native['modes']):
                if species or cell!=x:continue
                for j,(other,site,ss,ww,cc)in enumerate(native['modes']):
                    if other or site!=x:continue
                    a3[i,j]=bar[spin,ss]*T3[qi[w][c],qi[ww][cc]]
            zz=(native['B'].conj().T@one_body(a3,native['basis'])@native['B']).toarray()[np.ix_(active,active)]
            row.append((sp.csr_matrix(np.kron(C,plus)),sp.csr_matrix(np.kron(np.eye(3),zz))))
        qops.append(row)
    # Filled native neutrals and empty charged lepton modes; all three families.
    reference=(1<<12)-1;masks=[reference];lepton_labels=[None]
    for f in range(3):
        for xe,se,xn,sn in product(range(2),repeat=4):
            e=_lepton_index(1,xe,se,f);nu=_lepton_index(0,xn,sn,f)
            masks.append(reference^(1<<nu)|(1<<e));lepton_labels.append((f,xe,se,xn,sn))
    lep=np.zeros((24,24),complex);constant=0.
    for weak,spin,f in product(range(2),range(2),range(3)):
        ids=[_lepton_index(weak,x,spin,f)for x in range(2)]
        kinetic=q/2*np.array([[1.,-1.],[-1.,1.]])
        block=mass[f]*np.eye(2)+kinetic if weak else -kinetic
        lep[np.ix_(ids,ids)]=block
        if weak==0:constant+=q
    Hl=sp.csr_matrix(_restricted_current(lep,masks)+constant*np.eye(49))
    assert abs(Hl[0,0])<1e-12
    jl=raising[np.ix_([10,11],[10,11])];t3=(jl@jl.conj().T-jl.conj().T@jl)/2
    lops=[]
    for x in range(2):
        row=[]
        for bar in BAR:
            pair=[]
            for weakmatrix in(jl,t3):
                A=np.zeros((24,24),complex)
                for w,ww,s,ss,f in product(range(2),range(2),range(2),range(2),range(3)):
                    A[_lepton_index(w,x,s,f),_lepton_index(ww,x,ss,f)]=bar[s,ss]*weakmatrix[w,ww]
                pair.append(sp.csr_matrix(_restricted_current(A,masks)))
            row.append(pair)
        lops.append(row)
    fullsize=96*49;V=sp.csr_matrix((fullsize,fullsize),dtype=complex)
    for x in range(2):
        for mu in range(4):
            qp,qz=qops[x][mu];lp,lz=lops[x][mu]
            V+=ETA[mu]*(sp.kron(qp,lp.conj().T)+sp.kron(qp.conj().T,lp)+2*sp.kron(qz,lz))
    H0=(sp.kron(Hh,sp.eye(49))+sp.kron(sp.eye(96),Hl)).tocsr()
    keep=[];total_charge=[];spin2=[]
    for hi,(f,xq,sq,w,xa,sa)in enumerate(hadron_labels):
        for li,ell in enumerate(lepton_labels):
            Ne=int(ell is not None);charge=-w-Ne
            sz=(1-2*sq)+(1-2*sa)
            if ell is not None:sz+=(1-2*ell[2])-(1-2*ell[4])
            total_charge.append(charge);spin2.append(sz)
            if charge==-1 and sz==0:keep.append(hi*49+li)
    keep=np.array(keep);outside=np.setdiff1d(np.arange(fullsize),keep)
    assert len(keep)==888
    assert norm(H0[keep][:,outside])<2e-13 and norm(V[keep][:,outside])<2e-13
    H0=H0.tocsr()[keep][:,keep];V=V.tocsr()[keep][:,keep]
    assert norm(V-V.conj().T)<2e-13
    H=H0+G*V
    lookup={v:i for i,v in enumerate(keep)}
    hadron_initial=[i for i,(f,xq,sq,w,xa,sa)in enumerate(hadron_labels)if f==1 and w==1]
    # Preserve the orbital seed rule independently of projector roundoff.
    _,seed_frame=eigh(Hh.toarray()[np.ix_(hadron_initial,hadron_initial)])
    seed_hadron=hadron_initial[int(np.argmax(abs(seed_frame[:,0])**2))]
    # The scalar family energy changes eigenvalues, not the native eigenspace.
    native_initial=np.array(hadron_initial)%32
    eh,Uh=eigh(h[np.ix_(native_initial,native_initial)])
    eh=eh+(mass[1]-1)
    assert eh[1]-eh[0]>1e-5
    initial=np.zeros(len(keep),complex)
    for i,z in zip(hadron_initial,Uh[:,0]):
        if i*49 in lookup:initial[lookup[i*49]]=z
        else:assert abs(z)<1e-12
    assert abs(np.vdot(initial,initial)-1)<1e-12
    hadron_final=[i for i,(f,xq,sq,w,xa,sa)in enumerate(hadron_labels)if f==0 and w==0]
    ef,Uf=eigh(Hh.toarray()[np.ix_(hadron_final,hadron_final)])
    final_h=np.zeros(96,complex);final_h[hadron_final]=Uf[:,0]
    final_l=np.zeros(49,complex)
    for k,label in enumerate(lepton_labels):
        if label is None:continue
        f,xe,se,xn,sn=label
        if f==0 and se==sn:
            hole=_lepton_index(0,xn,sn,f)
            final_l[k]=(-1)**(xe+xn+hole+11)/(2*np.sqrt(2))
    target=np.kron(final_h,final_l)[keep]
    assert abs(np.vdot(target,target)-1)<1e-12
    vi=float(np.vdot(initial,V@initial).real);vf=float(np.vdot(target,V@target).real)
    assert abs(vi+3)<1e-12 and abs(vf-2.5)<1e-12
    assert abs(np.vdot(initial,H@initial).real-np.vdot(target,H@target).real)<1e-10 if contact_resonance else True
    seed=lookup[seed_hadron*49];success=float(abs(initial[seed])**2)
    rawseed=np.eye(len(keep),dtype=complex)[:,seed]
    failure=rawseed-initial*np.vdot(initial,rawseed)
    measurement_diagonal=np.repeat(measured_hadron,49)+np.tile(Hl.diagonal().real,96)
    measurement_diagonal=measurement_diagonal[keep]+G*V.diagonal().real
    return dict(H=H,H0=H0,V=V,initial=initial,target=target,parent_energy=float(eh[0]),mass=mass,q=q,C=C,
                keep=keep,hadron_labels=hadron_labels,lepton_labels=lepton_labels,masks=masks,
                spectral_seed=seed,success=success,failure=failure,seed=rawseed,
                measured_energy=measurement_diagonal,zero_word=zero_word,contact_resonance=contact_resonance,
                full_colour_invariance_error=norm(native['H']@native['B']-native['B']@sp.csr_matrix(native['h'])))

def spin_parity_reduction(d):
    spin=[sy.Matrix([[0,1],[1,0]])/2,sy.Matrix([[0,-sy.I],[sy.I,0]])/2,sy.diag(1,-1)/2]
    total=[]
    for axis in range(3):
        generators=[]
        for factor in range(4):
            factors=[sy.eye(2)for _ in range(4)];factors[factor]=-spin[axis].T if factor==3 else spin[axis]
            generators.append(sy.kronecker_product(*factors))
        total.append(sum(generators,sy.zeros(16)))
    S2=sum((a*a for a in total),sy.zeros(16));null=sy.GramSchmidt(S2.nullspace(),orthonormal=True)
    assert len(null)==2
    spin0=np.array([numeric(x).ravel()for x in null]);lookup={int(v):k for k,v in enumerate(d['keep'])}
    hadron={v:i for i,v in enumerate(d['hadron_labels'])};lepton={v:i for i,v in enumerate(d['lepton_labels'])}
    columns=[];keys=[]
    for f,xq,xa in product(range(3),range(2),range(2)):
        v=np.zeros(888,complex)
        for sq,sa in((0,1),(1,0)):
            hi=hadron[(f,xq,sq,1,xa,sa)];v[lookup[hi*49]]=(1 if sq==0 else -1)/np.sqrt(2)
        columns.append(v);keys.append(('parent',f,xq,xa))
    for f,ell,xq,xa,xe,xn,k in product(range(3),range(3),range(2),range(2),range(2),range(2),range(2)):
        v=np.zeros(888,complex)
        for sq,sa,se,sn in product(range(2),repeat=4):
            value=spin0[k,((2*sq+sa)*2+se)*2+sn]
            if abs(value)==0:continue
            hi=hadron[(f,xq,sq,0,xa,sa)];li=lepton[(ell,xe,se,xn,sn)]
            phase=(-1)**(_lepton_index(0,xn,sn,ell)+11)
            v[lookup[hi*49+li]]=phase*value
        columns.append(v);keys.append(('emission',f,ell,xq,xa,xe,xn,k))
    B0=np.stack(columns,axis=1);where={k:i for i,k in enumerate(keys)};even=[]
    for i,key in enumerate(keys):
        if key[0]=='parent':other=('parent',key[1],1-key[2],1-key[3])
        else:other=('emission',key[1],key[2],*(1-x for x in key[3:7]),key[7])
        j=where[other]
        if i<j:even.append((B0[:,i]+B0[:,j])/np.sqrt(2))
    B=np.stack(even,axis=1);assert B.shape==(888,150)
    h=B.conj().T@(d['H']@B)
    assert np.max(abs(B.conj().T@B-np.eye(150)))<2e-14
    residual=float(np.linalg.norm(d['H']@B-B@h));assert residual<5e-11
    assert np.linalg.norm(d['initial']-B@(B.conj().T@d['initial']))<2e-12
    assert np.linalg.norm(d['target']-B@(B.conj().T@d['target']))<2e-12
    return B,h,residual


def explore(zero_word=False,contact_resonance=True):
    d=construction(zero_word,contact_resonance);B,h,_=spin_parity_reduction(d)
    w,U=eigh(h);frame=U.conj().T@B.conj().T@d['initial'];rows=[]
    for t in(0.,32.,128.,256.,512.,720.,1024.):
        psi=B@U@(np.exp(-1j*t*w)*frame);p=abs(psi)**2
        emission=sum(p[k]for k,i in enumerate(d['keep'])if i%49!=0)
        light=sum(p[k]for k,i in enumerate(d['keep'])if i%49!=0 and d['hadron_labels'][i//49][0]==0)
        rows.append((t,float(emission),float(light),float(abs(np.vdot(d['target'],psi))**2),float(np.vdot(psi,d['H']@psi).real)))
    return d,rows


def momentum_readout(d):
    lookup={int(v):k for k,v in enumerate(d['keep'])};lepton={v:i for i,v in enumerate(d['lepton_labels'])}
    T=sp.lil_matrix((888,888),dtype=complex)
    for row,value in enumerate(d['keep']):
        hi,li=divmod(int(value),49);label=d['lepton_labels'][li]
        if label is None:T[row,row]=1;continue
        family,pe,se,pn,sn=label
        for xe,xn in product(range(2),repeat=2):
            old=lepton[(family,xe,se,xn,sn)];column=lookup[hi*49+old]
            phase=(-1)**(_lepton_index(0,xn,sn,family)+11)
            T[row,column]=phase*(-1)**(pe*xe+pn*xn)/2
    T=T.tocsr();assert norm(T@T.conj().T-sp.eye(888))<1e-14
    return T


def complete_records(d,time,eigen=None):
    B,h,residual=spin_parity_reduction(d)
    if d['zero_word']:
        # Exact invariant evolution, embedded in every original outcome below.
        family=np.array([d['hadron_labels'][int(v)//49][0]==1 for v in d['keep']])
        assert np.count_nonzero(d['H'][family][:,~family].data)==0
        assert np.count_nonzero(d['H'][~family][:,family].data)==0
        assert np.count_nonzero(d['initial'][~family])==0
        columns=np.flatnonzero(np.any(B[family]!=0,axis=0))
        assert np.count_nonzero(B[~family][:,columns])==0
        B=B[:,columns];h=h[np.ix_(columns,columns)]
        w,U=eigh(h)
    else:
        w,U=eigh(h)if eigen is None else eigen
    initial=B.conj().T@d['initial'];state=B@U@(np.exp(-1j*time*w)*(U.conj().T@initial))
    p=abs(state)**2;T=momentum_readout(d);out=T@state;pm=abs(out)**2;weight=d['success']
    assert abs(p.sum()-1)<3e-12 and abs(pm.sum()-1)<3e-12
    emitted=np.array([i%49!=0 for i in d['keep']]);by_family=np.zeros((3,2,2));by_hadron=np.zeros(3)
    for i,value in enumerate(d['keep']):
        hi,li=divmod(int(value),49);label=d['lepton_labels'][li]
        if label is not None:
            f,pe,se,pn,sn=label;by_family[f,pe,pn]+=float(pm[i]);by_hadron[d['hadron_labels'][hi][0]]+=float(p[i])
    energy=float(np.vdot(state,d['H']@state).real);initial_energy=float(np.vdot(d['initial'],d['H']@d['initial']).real)
    assert abs(energy-initial_energy)<2e-9
    # A fine native colour record breaks the retained coherent colour singlet.
    # The exact full occupation diagonal adds its current-square energy.
    colour_cost=d['measured_energy']-d['H'].diagonal().real
    momentum_energy=(T@d['H']@T.conj().T).diagonal().real+colour_cost
    work_position=float(np.dot(p,d['measured_energy'])-energy)
    work_momentum=float(np.dot(pm,momentum_energy)-energy)
    failure_energy=float(np.vdot(d['failure'],d['H']@d['failure']).real)
    preparation_work=weight*initial_energy+failure_energy-float(np.vdot(d['seed'],d['H']@d['seed']).real)
    assert abs(preparation_work)<1e-10
    return {'source_time':time,'preparation_success_probability':weight,'preparation_failure_probability':1-weight,
            'conditional_total_lepton_pair_probability':float(p[emitted].sum()),
            'unconditional_total_lepton_pair_probability':float(weight*p[emitted].sum()),
            'conditional_light_bound_spin_zero_receiver_probability':float(abs(np.vdot(d['target'],state))**2),
            'unconditional_light_bound_spin_zero_receiver_probability':float(weight*abs(np.vdot(d['target'],state))**2),
            'conditional_emission_by_quark_family':by_hadron.tolist(),
            'conditional_lepton_family_and_two_momenta':by_family.tolist(),
            'complete_position_colour_triplet_probability_each':(weight*p/3).tolist(),
            'complete_momentum_colour_triplet_probability_each':(weight*pm/3).tolist(),
            'colour_expansion_rule':'Each listed probability occurs for each of three matched native colour occupations. Position and momentum lists are separate complete instruments.',
            'complete_probability_sum':float(3*sum(weight*p/3)+1-weight),
            'orthogonal_full_Fock_complement_probability':0.,
            'complete_coarse_record_count_each_instrument':2666,
            'source_energy_before_readout':energy,'spectral_preparation_mean_work':preparation_work,
            'successful_position_readout_mean_work':work_position,'successful_momentum_readout_mean_work':work_momentum,
            'total_position_readout_mean_work':weight*work_position,'total_momentum_readout_mean_work':weight*work_momentum,
            'exact_spin_parity_reduction_numeric_residual':residual,
            'arithmetic_scope':'floating evaluation of the exact finite source and complete instruments; no outward probability-roundoff enclosure'}


def certificate():
    d=construction();B,h,residual=spin_parity_reduction(d);eigen=eigh(h)
    times=(0.,32.,128.,256.,512.,720.,1024.)
    rows=[complete_records(d,t,eigen)for t in times]
    assert rows[5]['conditional_total_lepton_pair_probability']>.99
    assert rows[5]['unconditional_total_lepton_pair_probability']>.24
    full_values,full_vectors=eigh(d['H'].toarray())
    full=full_vectors@(np.exp(-720j*full_values)*(full_vectors.conj().T@d['initial']))
    reduced=B@eigen[1]@(np.exp(-720j*eigen[0])*(eigen[1].conj().T@B.conj().T@d['initial']))
    disagreement=float(np.linalg.norm(full-reduced));assert disagreement<1e-8
    zero=construction(True);zero_row=complete_records(zero,720.)
    family=np.array([zero['hadron_labels'][int(v)//49][0] for v in zero['keep']])
    zero_off_family=zero['H'].toarray()*(family[:,None]!=family[None,:])
    assert np.count_nonzero(zero_off_family)==0
    assert zero_row['conditional_emission_by_quark_family'][0]<1e-16
    untuned=construction(False,False);untuned_row=complete_records(untuned,720.)
    form=np.vdot(d['target'],d['V']@d['initial']);expected=d['C'][0,1]/np.sqrt(2)
    assert abs(form-expected)<2e-12
    return {'schema':'oph.native_bound_meson_semileptonic_current.v1',
            'source':'unchanged positive colour/spin two-cell hadron action in each family, plus the complete local native quark-lepton charged and neutral weak-current contraction',
            'full_CAR_modes':180,'actual_colour_source_module':15,'families':3,'sites':2,'spins':2,
            'filled_reference':'all12 neutral lepton coordinates filled, charged leptons and other spectators vacant; zero neutral kinetic modes filled by declaration',
            'same_hadronic_source_couplings':{'U':1,'J':.25,'t':.125},
            'native_family_number_energies':d['mass'].tolist(),'weak_contact_G':G,
            'declared_family_frame_port_word':[2,3,10],'declared_port_phase':'-i',
            'actual_current_frame_real':d['C'].real.tolist(),'actual_current_frame_imag':d['C'].imag.tolist(),
            'current_frame_definition':'F^dagger V F for ordered real eigenframe of native Y and native port word V; assigned to the conserved analyzer family-energy labels',
            'free_gap_outgoing_frequency':float((d['mass'][1]-2)/2),'contact_resonant_outgoing_frequency':d['q'],
            'resonance_formula':'q=(m_second-2)/2-11G/4; source-derived equal expectation energies, not an empirical weak/kaon fit',
            'native_contact_expectations_parent_receiver':['-3','5/2'],
            'supplied_seed_energy_above_prepared_reference':float(np.vdot(d['seed'],d['H']@d['seed']).real),
            'exact_bound_to_bound_current_form_factor':'<receiver|V|parent>=C_01/sqrt(2); vector hadron form factor1, spatial spin-vector contribution0 between spin-zero mesons',
            'form_factor_real_imag':[float(form.real),float(form.imag)],
            'source_reductions':{'all_colour_singlet_hadron_dimension':96,'lepton_vacuum_and_pair_dimension':49,
              'charge_and_spin_z_dimension':888,'total_spin_zero_even_parity_dimension':150,
              'full_colour_reduction_residual':d['full_colour_invariance_error'],'spin_parity_reduction_residual':residual,
              'full888_vs_reduced150_evolution_difference':disagreement},
            'record_labels':[{'hadron_family_cell_spin_weak_anticell_antispin':list(d['hadron_labels'][int(v)//49]),
                              'lepton_family_cell_or_momentum_spin_neutral_cell_or_momentum_spin':None if int(v)%49==0 else list(d['lepton_labels'][int(v)%49])}
                              for v in d['keep']],
            'complete_records':rows,'zero_port_word_control':zero_row,'untuned_free_gap_box_control':untuned_row,
            'zero_port_exact_family_offblock_nonzero_count':0,
            'zero_port_light_emission_exact':'0 by exact family-block conservation; numerical eigensolver leakage is retained in the floating record',
            'zero_contact_control':'all lepton-pair probabilities exactly0; the bound parent evolves by phase only',
            'architecture_inputs':['positive hadron action','native family-label mapping and port word','weak contactG','filled neutral reference','two-site outgoing frequency chosen for source resonance','spectral preparation and complete position/momentum apparatus'],
            'empirical_hadron_mass_lifetime_or_mixing_inputs':[],
            'physical_kaon_decay_or_lifetime_identified':False,
            'scope':'A full finite semileptonic bound-meson conversion with complete failures, colour, lepton momentum and work. The source-resonant box and native pulse are adjustable architecture, not forced physical dynamics. No continuum decay rate, observed strange identity, local electroweak gauge completion or transfer to the separately deformed nuclear source.',
            'verdict':'PASS_NATIVE_BOUND_MESON_CURRENT_CONVERSION_WITH_COMPLETE_RECORDS'}
