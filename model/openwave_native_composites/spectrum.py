"""Exact reducing colour sectors of a positive native-current CAR source."""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, sqrt
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
import sympy as sy
from model.coherent_defect_wave.physical_attachment.chiral import representation

U=1.; J=.25; T=.125; MASS=1.

@lru_cache(None)
def native():
    subsets,_,g,_=representation()
    qi=[[subsets.index((c,w)) for c in range(3)] for w in (3,4)]
    ai=[[subsets.index(tuple(c for c in range(3) if c!=missing)) for missing in range(3)],
        [subsets.index(tuple(c for c in range(5) if c!=missing)) for missing in range(3)]]
    # Orthogonalize the two Cartan currents before forming the Casimir.
    colour=list(g[:7])+[(g[6]+2*g[7])/sy.sqrt(3)]
    qg=[[np.array(h.extract(ids,ids)).astype(complex) for h in colour] for ids in qi]
    ag=[[np.array(h.extract(ids,ids)).astype(complex) for h in colour] for ids in ai]
    joins=[]
    for f in range(2):
        # q^T M anti is invariant; solve the native-current equations exactly.
        cols=[]
        for i,j in product(range(3),repeat=2):
            e=sy.zeros(3); e[i,j]=1
            cols.append(sy.Matrix.vstack(*[(g[k].extract(qi[0],qi[0]).T*e+e*g[k].extract(ai[f],ai[f])).reshape(9,1) for k in range(8)]))
        ns=sy.Matrix.hstack(*cols).nullspace(); assert len(ns)==1
        m=ns[0].reshape(3,3); m=m/sy.sqrt((m.H*m).trace()/3)
        assert m.H*m==sy.eye(3)
        joins.append(np.array(m).astype(complex))
    qcharge=[[float((g[10]+g[11])[i,i]) for i in ids] for ids in qi]
    acharge=[[float((g[10]+g[11])[i,i]) for i in ids] for ids in ai]
    return qg,ag,joins,qcharge,acharge


def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def create_word(ids):
    if len(set(ids))!=len(ids):return None,0
    return sum(1<<i for i in ids),parity(ids)

def one_body(a,basis):
    lookup={b:i for i,b in enumerate(basis)}; rows=[];cols=[];values=[]
    terms=[(i,j,a[i,j]) for i,j in zip(*np.nonzero(abs(a)>1e-14))]
    for col,b in enumerate(basis):
        for i,j,v in terms:
            if not (b>>j)&1:continue
            b1=b^(1<<j)
            if (b1>>i)&1:continue
            target=b1|(1<<i)
            if target not in lookup:raise AssertionError('claimed reducing sector is not invariant')
            sign=(-1)**((b&((1<<j)-1)).bit_count()+(b1&((1<<i)-1)).bit_count())
            rows.append(lookup[target]);cols.append(col);values.append(sign*v)
    return sp.csr_matrix((values,(rows,cols)),shape=(len(basis),len(basis)),dtype=complex)

def norm(a):return float(np.max(abs(a.data))) if sp.issparse(a) and a.nnz else (0. if sp.issparse(a) else float(np.max(abs(a))))

@lru_cache(None)
def sector(kind):
    qg,ag,joins,qcharge,acharge=native()
    # Modes (species, cell, spin, weak/conjugate type, native colour).
    modes=list(product(range(1 if kind=='baryon' else 2),range(2),range(2),range(2),range(3)))
    index={v:i for i,v in enumerate(modes)}; d=len(modes)
    if kind=='baryon':basis=[sum(1<<i for i in ids) for ids in combinations(range(d),3)]; number=3
    else:basis=[(1<<i)|(1<<j) for i in range(24) for j in range(24,48)]; number=2
    lookup={b:i for i,b in enumerate(basis)}
    pauli=[np.array([[0,1],[1,0]],complex)/2,np.array([[0,-1j],[1j,0]])/2,np.diag([1.,-1.])/2]
    hop=np.zeros((d,d),complex); charge=np.zeros(d); sz=np.zeros(d)
    for i,(species,x,s,f,c) in enumerate(modes):
        hop[i,i]+=T;hop[i,index[(species,1-x,s,f,c)]]-=T
        charge[i]=(qcharge if species==0 else acharge)[f][c];sz[i]=(.5 if s==0 else -.5)
    kinetic=one_body(hop,basis); C=sp.csr_matrix((len(basis),len(basis)),dtype=complex); S=C.copy(); global_c=[]; global_s=[]
    for axis in range(8):
        gs=sp.csr_matrix(C.shape,dtype=complex)
        for x in range(2):
            a=np.zeros((d,d),complex)
            for species,ss,f,c,cc in product(range(1 if kind=='baryon' else 2),range(2),range(2),range(3),range(3)):
                a[index[(species,x,ss,f,c)],index[(species,x,ss,f,cc)]]=(qg if species==0 else ag)[f][axis][c,cc]
            dx=one_body(a,basis);C=C+dx@dx;gs=gs+dx
        global_c.append(gs)
    for axis in range(3):
        gs=sp.csr_matrix(C.shape,dtype=complex)
        for x in range(2):
            a=np.zeros((d,d),complex)
            for species,f,c,s,ss in product(range(1 if kind=='baryon' else 2),range(2),range(3),range(2),range(2)):
                a[index[(species,x,s,f,c)],index[(species,x,ss,f,c)]]=pauli[axis][s,ss]
            dx=one_body(a,basis);S=S+dx@dx;gs=gs+dx
        global_s.append(gs)
    H=MASS*number*sp.eye(len(basis))+kinetic+U*C+J*S
    orbitals=list(product(range(2),range(2),range(2)))
    br=[];bc=[];bv=[]
    labels=list(combinations_with_replacement(orbitals,3)) if kind=='baryon' else list(product(orbitals,repeat=2))
    for col,label in enumerate(labels):
        terms={}
        if kind=='baryon':
            for cp in permutations(range(3)):
                mask,sign=create_word([index[(0,*orb,c)] for orb,c in zip(label,cp)])
                if sign:terms[mask]=terms.get(mask,0)+parity(cp)*sign
        else:
            a,b=label
            for c,cc in product(range(3),repeat=2):
                mask,sign=create_word([index[(0,*a,c)],index[(1,*b,cc)]])
                terms[mask]=joins[b[2]][c,cc]*sign
        scale=sqrt(sum(abs(v)**2 for v in terms.values()))
        for mask,v in terms.items():
            if v:br.append(lookup[mask]);bc.append(col);bv.append(v/scale)
    B=sp.csr_matrix((bv,(br,bc)),shape=(len(basis),len(labels)),dtype=complex)
    assert norm(B.conj().T@B-sp.eye(len(labels)))<1e-12
    assert max(norm(g@B) for g in global_c)<1e-12
    h=(B.conj().T@H@B).toarray()
    assert norm(H@B-B@sp.csr_matrix(h))<1e-12
    spin=(B.conj().T@sum((g@g for g in global_s),sp.csr_matrix(H.shape))@B).toarray()
    charge_diag=np.array([sum(charge[i] for i in range(d) if (b>>i)&1) for b in basis])
    sz_diag=np.array([sum(sz[i] for i in range(d) if (b>>i)&1) for b in basis])
    q=(B.conj().T@sp.diags(charge_diag)@B).toarray(); z=(B.conj().T@sp.diags(sz_diag)@B).toarray()
    assert norm(h@spin-spin@h)<1e-12 and norm(h@q-q@h)<1e-12
    return dict(H=H,B=B,h=h,spin=spin,q=q,z=z,basis=basis,modes=modes,labels=labels,C=C,S=S,kinetic=kinetic)

def selected(kind,charge,spin):
    a=sector(kind); h=a['h'];s=a['spin'];q=a['q'];z=a['z']
    # Orthogonal simultaneous charge/spin/z projection, including every other outcome.
    filter_matrix=(q-charge*np.eye(len(q)))@(q-charge*np.eye(len(q)))+(s-spin*(spin+1)*np.eye(len(s)))@(s-spin*(spin+1)*np.eye(len(s)))+(z-spin*np.eye(len(z)))@(z-spin*np.eye(len(z)))
    w,v=eigh(filter_matrix);P=v[:,w<1e-8];assert P.shape[1]
    ev,u=eigh(P.conj().T@h@P)
    ground=u[:,abs(ev-ev[0])<1e-8]; projector=P@ground@ground.conj().T@P.conj().T
    diagonal=np.diag(projector).real
    seed=int(np.flatnonzero(diagonal>max(diagonal)-1e-9)[0]);success=float(diagonal[seed])
    psi=projector[:,seed]/sqrt(success);time=.7
    evals,evec=eigh(h);evolved=evec@(np.exp(-1j*time*evals)*(evec.conj().T@psi))
    assert abs(np.vdot(evolved,evolved)-1)<1e-12
    full=a['B']@evolved; probs=abs(full)**2;records={};coloc=0.;charge_cells=[0.,0.]
    for b,p in zip(a['basis'],probs):
        occupied=[a['modes'][i] for i in range(len(a['modes'])) if (b>>i)&1]
        cells=[sum(v[1]==x for v in occupied) for x in range(2)]
        key=','.join(map(str,cells));records[key]=records.get(key,0.)+float(p)
        if 0 in cells:coloc+=float(p)
        for species,x,ss,f,c in occupied:
            ch=(native()[3] if species==0 else native()[4])[f][c]
            charge_cells[x]+=float(p)*ch
    return {'charge':charge,'internal_spin':spin,'ground_energy':round(float(ev[0]),12),
      'selected_sector_dimension':P.shape[1],'ground_multiplicity_at_fixed_charge_and_spin_projection':ground.shape[1],
      'eigen_residual':float(np.linalg.norm(h@psi-ev[0]*psi)),
      'co_location_probability':round(coloc,12),'cell_occupation_records':{k:round(v,12) for k,v in sorted(records.items())},
      'cell_charge_expectations':[round(v,12) for v in charge_cells],
      'successful_full_occupation_record_probabilities':[round(float(p),12) for p in probs],
      'complete_unconditional_records':{'preparation_failure':round(1-success,12),
       'preparation_success_and_every_full_sector_occupation':[round(float(success*p),12) for p in probs],
       'outside_fixed_native_occupation_sector':'identity failure channel; prepared probability0'},
      'spectral_preparation':{'basis_seed':seed,'success_probability':round(success,12),'failure_probability':round(1-success,12),'complete_record_count':len(a['basis'])+1,'time':time,'phase_only_evolution':True}}

def certificate():
    rows=[]
    for kind,charges,spin in [('baryon',[0,1],.5),('baryon',[-1,0,1,2],1.5),('meson',[-1,0,1],0),('meson',[-1,0,1],1)]:
        for charge in charges:rows.append({'kind':kind,**selected(kind,charge,spin)})
    sizes={}
    for kind in ('baryon','meson'):
        a=sector(kind);ev=np.linalg.eigvalsh(a['h']);sizes[kind]={'full_fixed_occupation_dimension':len(a['basis']),'exact_colour_singlet_dimension':len(a['labels']),'lowest_energy':round(float(ev[0]),12),'highest_energy':round(float(ev[-1]),12),'full_sector_invariance_residual':norm(a['H']@a['B']-a['B']@sp.csr_matrix(a['h']))}
    return {'verdict':'PASS: bound native-current colour singlets; physical hadron identification remains conditional',
      'architecture':{'cells':2,'native_internal_coordinates':15,'spin_coordinates':2,'full_matter_CAR_modes':60,'full_fock_dimension':'2^60','m':MASS,'U':U,'J':J,'t':T,'kinetic':'positive graph Laplacian','gauge_scope':'global native colour, weak and hypercharge; no local Gauss link field','reference':'empty matter vacuum; opposite-charge and conjugate-colour idle apparatus references; zero reference energy','clock':'autonomous H time with externally specified preparation/readout'},
      'all_fock_floor':'H >= m_min N; vacuum unique on the matter factor',
      'sector_sizes':sizes,'states':rows,
      'uniform_graph_bounds':{'maximum_degree':1,'baryon_trial_upper':3*MASS+3*J/4+3*T,'two_plus_one_threshold_lower':3*MASS+8*U/3+3*J/4,'baryon_binding_gap_lower':8*U/3-3*T,'meson_trial_upper':2*MASS+2*T,'one_plus_one_threshold':2*MASS+8*U/3+3*J/2,'meson_binding_gap_lower':8*U/3+3*J/2-2*T},
      'physical_limitations':['declared positive current interaction and number-energy action, distinct from the Weyl/Higgs action','internal SU(2) spin is not a proof of relativistic total angular momentum','finite bound eigenstates do not establish physical hadron masses, lifetimes or continuum confinement','neutral meson singlet, vector and spin-three-halves states are retained','no strange decay, neutron-proton splitting or local gauge-field dynamics'],
      'controls':{'turn_off_U_and_J':'free sectors have zero binding below fragment thresholds','equal_light_masses':'exact native weak flavour symmetry forces proton-like/neutron-like degeneracy'}}
