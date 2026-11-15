"""Uniform Gauss-fiber reduction of the actual CAR shared-record instrument."""
from pathlib import Path
from itertools import combinations
import json,hashlib,math
import numpy as np
from fermionic_shared_records import FermionicSharedRecords,Construction,transport,BASE,FOCK


def reduction(model):
    """Connected complete sector: equal fibers over all fixed-number occupations."""
    occupations=tuple(sum(1<<z for z in ids) for ids in combinations(range(model.positions),model.number))
    index={mask:i for i,mask in enumerate(occupations)}
    counts={mask:sum(s[0]==mask for s in model.states) for mask in occupations}
    if len(set(counts.values()))!=1 or min(counts.values())<1:raise ValueError('complete equal fibers required')
    fiber=next(iter(counts.values()));v=np.zeros((model.n,len(occupations)))
    for row,(mask,_) in enumerate(model.states):v[row,index[mask]]=1/math.sqrt(fiber)
    c=model.construction;h=np.zeros((len(occupations),)*2)
    for col,mask in enumerate(occupations):
        h[col,col]=c.mass*sum((1 if z<c.ports else -1)*((mask>>z)&1) for z in range(model.positions))
        for e,ends in enumerate(c.seams):
            for side,x in enumerate(ends):
                for src,dst in ((x,c.ports+e),(c.ports+e,x)):
                    moved=transport(mask,src,dst,model.positions)
                    if moved is not None:h[index[moved[0]],col]+=(1 if side==0 else -1)*moved[1]
    operators={}
    for law in model.operators:
        operators[law]=[]
        for e,full_ops in enumerate(model.operators[law]):
            reduced=[]
            for meta,_ in full_ops:
                k=np.zeros((len(occupations),)*2)
                for source,mask in enumerate(occupations):
                    # Any representative has identical occupation target and CARphase.
                    representative=next(s for s in model.states if s[0]==mask)
                    target,j,d,phase=model.target(representative,e,meta['coin'])
                    value=j if meta['law']=='C' else d
                    if value==meta['record_value']:
                        k[index[target[0]],source]=phase/(2 if law=='MIX' else math.sqrt(2))
                reduced.append((meta,k))
            operators[law].append(tuple(reduced))
    return occupations,v,h,operators,fiber


def run():
    model=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    masks,v,h,operators,fiber=reduction(model)
    # Independently construct the Slater reference from the unchanged one-particle H.
    one=model.construction.hamiltonian;energies,orbitals=np.linalg.eigh(one);negative=orbitals[:,energies<0]
    slater=np.array([np.linalg.det(negative[[z for z in range(6) if mask>>z&1],:]) for mask in masks])
    ground=float(energies[energies<0].sum());rho=np.outer(slater,slater.conj());lift=v@rho@v.T
    herror=float(np.max(abs(model.hamiltonian@v-v@h)))
    kerror=max(float(np.max(abs(full@v-v@small))) for law in operators for fullops,smallops in zip(model.operators[law],operators[law]) for (_,full),(_,small) in zip(fullops,smallops))
    grounderror=float(np.linalg.norm(h@slater-ground*slater));spectrum=np.linalg.eigvalsh(h)-ground
    output={};posterror=0.
    for law in operators:
        final=np.zeros_like(rho);records=[]
        for e,ops in enumerate(operators[law]):
            for meta,k in ops:
                branch=k@rho@k.T/3;final+=branch
                records.append({**meta,'probability':float(np.trace(branch).real),'unnormalized_excitation_energy':float(np.trace((h-ground*np.eye(20))@branch).real)})
        full,ledger=model.commit(lift,law)
        posterror=max(posterror,float(np.max(abs(full-v@final@v.T))))
        output[law]={'complete_records':records,'excitation_energy_after':float(np.trace((h-ground*np.eye(20))@final).real),'instrument_work':float(np.trace(h@(final-rho)).real),'trace':float(np.trace(final).real),'ground_survival':float(np.vdot(slater,final@slater).real),'poststate_real':final.real.tolist(),'poststate_imag':final.imag.tolist()}
    paths=[Path(__file__),Path(__file__).with_name('fermionic_shared_records.py'),FOCK,BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/fiber_preparation.py']
    return {'schema':'oph.c9.fermionic_uniform_reference.v1','checks':{'isometry':bool(np.max(abs(v.T@v-np.eye(20)))<1e-12),'H_intertwiner':bool(herror<1e-12),'every_C_E_MIX_K_intertwiner':bool(kerror<1e-12),'actual_Slater_ground':bool(grounderror<1e-12),'positive_excitation_spectrum':bool(min(spectrum)>-1e-12),'all_complete_poststates_intertwine':bool(posterror<1e-12),'monitoring_generates_positive_work':bool(all(row['instrument_work']>1e-4 for row in output.values())),'global_holonomy_ground_transfer_rejected':bool(np.linalg.eigvalsh(model.hamiltonian)[0]<ground-1e-4)},'parameters':{'positions':6,'fixed_fermion_number':3,'modulus':3,'mass':.23,'preparation':'uniform modular Gauss fiber on each occupation, negative-filled Slater occupation superposition','native_phase':'direct CAR transport; same phase as full model','full_dimension':model.n,'reduced_dimension':20,'fiber_cardinality':fiber},'errors':{'H_intertwiner':herror,'K_intertwiner':kerror,'ground_eigenvector':grounderror,'complete_poststate':posterror},'ground_energy':ground,'ground_domain':'20dim invariant uniform-fiber code only, not all60holonomy states','full_H_minimum':float(np.linalg.eigvalsh(model.hamiltonian)[0]),'excitation_spectrum':spectrum.tolist(),'slater_amplitudes':slater.tolist(),'occupation_masks':list(masks),'outputs':output,'scope':'Same actual shared-record CAR construction in an explicitly supplied uniform-fiber Slater preparation. Passive coherent code-ground is stable; native monitoring generally excites it. Other record-holonomy sectors have a lower global energy; no global-ground transfer. No physical species, vacuum preparation, spin or mass identification.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
 r=run();Path(__file__).with_name('FERMIONIC_UNIFORM_REFERENCE.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
