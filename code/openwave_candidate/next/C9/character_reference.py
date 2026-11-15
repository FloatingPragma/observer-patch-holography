"""All modular character codes and actual global fixed-number reference."""
from pathlib import Path
import importlib.util,json,hashlib
import numpy as np
from fermionic_uniform_reference import reduction,FermionicSharedRecords,Construction,transport,BASE
FABLE=BASE/'validation/bridges/providers/fable_candidate.py'
spec=importlib.util.spec_from_file_location('executed_fable_characters',FABLE)
fable=importlib.util.module_from_spec(spec);spec.loader.exec_module(fable)


def blocks(model):
    masks,v0,h0,k0,fiber=reduction(model);c=model.construction
    phases=fable.holonomy_sectors(c.ports,c.seams,model.modulus)
    index={mask:i for i,mask in enumerate(masks)};out=[]
    for theta in phases:
        angles=theta.ravel();v=v0.astype(complex).copy()
        for row,(_,f) in enumerate(model.states):v[row]*=np.exp(1j*angles@f)
        single=fable.Gauged(c.ports,c.seams,c.mass,phase=theta).hamiltonian()
        h=np.zeros(h0.shape,complex)
        for col,mask in enumerate(masks):
            h[col,col]=sum(single[z,z] for z in range(model.positions) if mask>>z&1)
            for src in range(model.positions):
                for dst in range(model.positions):
                    if src==dst or abs(single[dst,src])<1e-14:continue
                    move=transport(mask,src,dst,model.positions)
                    if move is not None:h[index[move[0]],col]+=single[dst,src]*move[1]
        operators={}
        for law in k0:
            operators[law]=[]
            for seam,ops in enumerate(k0[law]):
                adjusted=[]
                for meta,k in ops:
                    d=meta['record_value'];j=d if meta['law']=='C' else (d-(2*meta['coin']-1)*(abs(d)%2))//2
                    adjusted.append((meta,k*np.exp(-1j*j*(theta[seam,0]-theta[seam,1]))))
                operators[law].append(tuple(adjusted))
        out.append((theta,v,h,operators))
    return masks,out


def run():
    model=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    masks,sectors=blocks(model);w=np.concatenate([v for _,v,_,_ in sectors],axis=1)
    herr=kerr=commutator=0.;sector_reports=[]
    for char,(theta,v,h,ops) in enumerate(sectors):
        herr=max(herr,float(np.max(abs(model.hamiltonian@v-v@h))))
        p=v@v.conj().T;commutator=max(commutator,float(np.max(abs(model.hamiltonian@p-p@model.hamiltonian))))
        for law in ops:
            for fullops,smallops in zip(model.operators[law],ops[law]):
                for (_,k),(_,small) in zip(fullops,smallops):
                    kerr=max(kerr,float(np.max(abs(k@v-v@small))))
                    commutator=max(commutator,float(np.max(abs(k@p-p@k))))
        sector_reports.append({'character_index':char,'halfhop_angles':theta.tolist(),'spectrum':np.linalg.eigvalsh(h).tolist()})
    eigenvalues,eigenvectors=np.linalg.eigh(model.hamiltonian);ground=float(eigenvalues[0]);selected=abs(eigenvalues-ground)<1e-10
    multiplicity=int(sum(selected));pground=eigenvectors[:,selected]@eigenvectors[:,selected].T;rho=pground/multiplicity
    weights=[float(np.trace((v@v.conj().T)@rho).real) for _,v,_,_ in sectors]
    gap=float(min(eigenvalues[eigenvalues>ground+1e-10])-ground);outputs={}
    for law in model.operators:
        final,ledger=model.commit(rho,law)
        outputs[law]={**ledger,'ground_survival':float(np.trace(pground@final).real),'poststate_real':final.real.tolist(),'poststate_imag':final.imag.tolist(),'trace':float(np.trace(final).real),'excitation_energy':float(np.trace((model.hamiltonian-ground*np.eye(model.n))@final).real)}
    files=[Path(__file__),Path(__file__).with_name('fermionic_uniform_reference.py'),Path(__file__).with_name('fermionic_shared_records.py'),FABLE]
    return {'schema':'oph.c9.character_reference.v1','checks':{'all_character_direct_sum_unitary':bool(np.max(abs(w.conj().T@w-np.eye(model.n)))<1e-12),'all_H_character_intertwiners':bool(herr<1e-12),'all_native_character_intertwiners':bool(kerr<1e-12),'H_and_every_K_preserve_character_projectors':bool(commutator<1e-12),'full_global_reference_passive':bool(np.linalg.norm(model.hamiltonian@rho-ground*rho)<1e-12),'all_global_excitation_energies_nonnegative':bool(min(eigenvalues-ground)>-1e-12),'global_ground_degeneracy_retained':bool(multiplicity==2),'all_native_monitoring_work_positive':bool(all(o['mean_instrument_work']>1e-4 for o in outputs.values()))},'parameters':{'graph':'triangle','mass':.23,'fixed_CAR_number':3,'modulus':3,'full_dimension':60,'characters':3,'occupation_dimension_per_character':20,'character_isometry':'exp(+i sum theta_h F_h)/sqrt3 at fixedoccupation','character_coordinates':'existing Fable spanningforest/chord labels; sign convention matches Gauged','global_reference':'normalized projector onto complete lowest-eigenvalue multiplet; no arbitrary unique-vector selection'},'errors':{'H_intertwiner':herr,'native_K_intertwiner':kerr,'character_commutator':commutator},'sector_spectra':sector_reports,'global_ground_energy':ground,'global_ground_multiplicity':multiplicity,'global_ground_character_weights':weights,'fixed_number_gap_above_ground_multiplet':gap,'global_ground_density_real':rho.real.tolist(),'global_ground_density_imag':rho.imag.tolist(),'outputs':outputs,'scope':'Actual full60state fixednumber sharedrecord source, all modularcharacters retained. FiniteCAR statistics/phase/ideal preparation supplied. Global passive ground generally excited by native monitoring. No unique physical vacuum, species, spin, continuum or calibrated mass claim.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}

if __name__=='__main__':
 r=run();Path(__file__).with_name('CHARACTER_REFERENCE.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
