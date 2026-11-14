"""Independent direct Fourier decomposition of full shared CAR record sector."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'candidate'))
from fermionic_shared_records import FermionicSharedRecords
from oph_openwave.construction import Construction
from audit_car_native_algebra import phase


def run(modulus=3):
    c=Construction(3,((0,1),(1,2),(2,0)),mass=.23);m=FermionicSharedRecords(c,7,modulus)
    masks=sorted({mask for mask,_ in m.states});idx={mask:i for i,mask in enumerate(masks)};columns=[];blocks=[];kerror=0.;herror=0.
    for chi in range(modulus):
        v=np.zeros((m.n,20),complex)
        for row,(mask,f) in enumerate(m.states):v[row,idx[mask]]=np.exp(2j*np.pi*chi*f[5]/modulus)/np.sqrt(modulus)
        columns.append(v);h=np.zeros((20,20),complex);reduced_ops=[]
        for col,mask in enumerate(masks):
            u=[(mask>>z)&1 for z in range(6)];h[col,col]=.23*sum((1 if z<3 else -1)*u[z] for z in range(6))
            for e,(a,b) in enumerate(c.seams):
                for side,x in enumerate((a,b)):
                    for src,dst,delta in ((x,3+e,1),(3+e,x,-1)):
                        if u[src] and not u[dst]:
                            angle=-2j*np.pi*chi*delta/modulus if 2*e+side==5 else 0
                            h[idx[mask^(1<<src)^(1<<dst)],col]+=(1 if side==0 else -1)*phase(u,src,dst)*np.exp(angle)
        herror=max(herror,float(np.max(abs(m.hamiltonian@v-v@h))))
        for law in ('C','E','MIX'):
            for e,operations in enumerate(m.operators[law]):
                for meta,k in operations:
                    small=np.zeros((20,20),complex)
                    for col,mask in enumerate(masks):
                        u=[(mask>>z)&1 for z in range(6)];a,b=c.seams[e];d=u[a]-u[b];j=(d-(2*meta['coin']-1)*abs(d))//2
                        if (j if meta['law']=='C' else d)!=meta['record_value']:continue
                        target=mask if j==0 else mask^(1<<a)^(1<<b)
                        shift=-j if e==2 else 0
                        small[idx[target],col]=(phase(u,a,b) if j else 1)*np.exp(-2j*np.pi*chi*shift/modulus)/(2 if law=='MIX' else np.sqrt(2))
                    kerror=max(kerror,float(np.max(abs(k@v-v@small))))
                    if law=='MIX':reduced_ops.append((meta,small))
        eigen,vectors=np.linalg.eigh(h);psi=vectors[:,0];rho=np.outer(psi,psi.conj());out=np.zeros_like(rho);ledger=[]
        for meta,k in reduced_ops:
            branch=k@rho@k.conj().T/3;out+=branch
            ledger.append({**meta,'probability':float(np.trace(branch).real),'subnormalized_energy':float(np.trace(h@branch).real)})
        blocks.append({'character':chi,'spectrum':eigen.tolist(),'ground_energy':float(eigen[0]),'complete_native_ground_outcomes':ledger,'ground_instrument_energy_after':float(np.trace(h@out).real),'ground_instrument_work':float(np.trace(h@(out-rho)).real)})
    unitary=np.concatenate(columns,axis=1);isometry=float(np.max(abs(unitary.conj().T@unitary-np.eye(m.n))))
    actual=np.linalg.eigvalsh(m.hamiltonian);union=np.sort(np.concatenate([b['spectrum'] for b in blocks]));ground=float(actual[0]);multiplicity=int(sum(abs(actual-ground)<1e-10));gap=float(actual[multiplicity]-ground)
    producer=HERE/('CHARACTER_REFERENCE.json' if modulus==3 else 'CHARACTER_MINIMA.json');saved=json.loads(producer.read_text())
    fullvalues,fullvectors=np.linalg.eigh(m.hamiltonian)
    selected=fullvectors[:,abs(fullvalues-ground)<1e-10]
    globalrho=selected@selected.conj().T/multiplicity
    savedrho=(np.array(saved['global_ground_density_real'])+1j*np.array(saved['global_ground_density_imag'])) if modulus==3 else globalrho
    reference_error=float(np.max(abs(globalrho-savedrho)))
    checks={'producer_global_reference_or_energy_reproduced':(reference_error<1e-12 if modulus==3 else abs(ground-saved['global_ground'])<1e-12),'producer_gap_reproduced':(abs(gap-saved['fixed_number_gap_above_ground_multiplet'])<1e-12 if modulus==3 else abs(gap-saved['global_fixed_number_gap'])<1e-12),'complete_character_unitary':isometry<1e-12,'all_H_intertwiners':herror<1e-12,'every_C_E_MIX_K_intertwiner':kerror<1e-12,'entire_full_spectrum_reproduced':bool(np.allclose(actual,union)),
      'expected_global_ground_multiplicity':multiplicity==(2 if modulus==3 else 1),'positive_gap_above_full_ground_multiplet':gap>0,'nonuniform_characters_contain_global_ground':any(abs(b['ground_energy']-ground)<1e-12 for b in blocks[1:]),'uniform_ground_global_transfer_rejected':blocks[0]['ground_energy']-ground>.1}
    if not all(checks.values()):raise ValueError(checks)
    return {'schema':'oph.c9.car_character_independent_review.v1','checks':checks,'producer_global_reference_error':reference_error,'unitary_error':isometry,'H_error':herror,'all_K_error':kerror,'modulus':modulus,'minimizing_characters':[b['character'] for b in blocks if abs(b['ground_energy']-ground)<1e-10],'global_ground_energy':ground,'global_ground_multiplicity':multiplicity,'gap_above_ground_multiplet':gap,'character_blocks':blocks,
      'scope':f'All{modulus}holonomycharacters of actual{m.n}statefixedN3source, no characterselected asphysicalbackground. Nativeinstrument preservescharacters but depositswork. Positiveglobalreference from entiregroundmultiplet, no physicalspecies/spin/vacuumpreparation claim.',
      'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/('character_reference.py' if modulus==3 else 'character_minima.py'),producer,HERE/'fermionic_shared_records.py',BASE/'candidate/oph_openwave/construction.py')},'full_criterion_met':False}

if __name__=='__main__':
 r=run();HERE.joinpath('CAR_CHARACTERS_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'ground':r['global_ground_energy'],'multiplicity':r['global_ground_multiplicity'],'gap':r['gap_above_ground_multiplet']}))
