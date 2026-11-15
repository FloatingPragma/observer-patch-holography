"""Independent uniform-fiber CAR reconstruction, determinant reference and work."""
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1]
sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'candidate'))
from fermionic_shared_records import FermionicSharedRecords
from oph_openwave.construction import Construction
from audit_car_native_algebra import phase


def run(modulus=3):
    c=Construction(3,((0,1),(1,2),(2,0)),mass=.23);full=FermionicSharedRecords(c,7,modulus)
    masks=sorted({mask for mask,_ in full.states});index={mask:i for i,mask in enumerate(masks)}
    v=np.zeros((full.n,len(masks)))
    for mask in masks:
        ids=[i for i,(m,_) in enumerate(full.states) if m==mask];v[ids,index[mask]]=1/np.sqrt(len(ids))
    h=np.zeros((20,20));ops=[]
    for j,mask in enumerate(masks):
        u=[(mask>>z)&1 for z in range(6)]
        h[j,j]=.23*sum((1 if z<3 else -1)*u[z] for z in range(6))
        for e,(a,b) in enumerate(c.seams):
            for side,x in enumerate((a,b)):
                for src,dst in ((x,3+e),(3+e,x)):
                    if u[src] and not u[dst]:h[index[mask^(1<<src)^(1<<dst)],j]+=(1 if side==0 else -1)*phase(u,src,dst)
    residuals=[];work=[]
    for seam,sourceops in enumerate(full.operators['MIX']):
        a,b=c.seams[seam]
        for meta,k in sourceops:
            reduced=np.zeros((20,20))
            for j,mask in enumerate(masks):
                u=[(mask>>z)&1 for z in range(6)];d=u[a]-u[b];coin=meta['coin'];current=(d-(2*coin-1)*abs(d))//2
                if (current if meta['law']=='C' else d)!=meta['record_value']:continue
                dst=mask if not current else mask^(1<<a)^(1<<b)
                reduced[index[dst],j]=.5*(phase(u,a,b) if current else 1)
            residuals.append(float(np.max(abs(k@v-v@reduced))));ops.append((meta,reduced))
    for law in ('C','E'):
        for e,sourceops in enumerate(full.operators[law]):
            for meta,k in sourceops:
                small=next(reduced for saved,reduced in ops if saved==meta) * np.sqrt(2)
                residuals.append(float(np.max(abs(k@v-v@small))))
    values,columns=np.linalg.eigh(c.hamiltonian);negative=np.flatnonzero(values<0);reference_energy=float(values[negative].sum())
    slater=np.array([np.linalg.det(columns[[i for i in range(6) if mask>>i&1]][:,negative]) for mask in masks])
    rho=np.outer(slater,slater);after=np.zeros_like(rho)
    for meta,k in ops:
        branch=k@rho@k.T/3;after+=branch
        work.append({**meta,'probability':float(np.trace(branch)),
                     'subnormalized_energy':float(np.trace(h@branch)),
                     'subnormalized_excitation_energy':float(np.trace((h-reference_energy*np.eye(20))@branch))})
    eig=np.linalg.eigvalsh(h);expected=sorted(sum(values[list(ids)]) for ids in itertools.combinations(range(6),3))
    full_ground=float(np.linalg.eigvalsh(full.hamiltonian)[0])
    full_reference_minimum=full_ground-reference_energy
    herror=float(np.max(abs(full.hamiltonian@v-v@h)));kerror=max(residuals)
    reference_error=float(np.linalg.norm(h@slater-reference_energy*slater))
    actualwork=float(np.trace(h@after)-reference_energy)
    fullafter,_=full.commit(v@rho@v.T);outputerror=float(np.max(abs(fullafter-v@after@v.T)))
    producer=HERE/'FERMIONIC_UNIFORM_REFERENCE.json'
    saved=json.loads(producer.read_text())
    order=[masks.index(mask) for mask in saved['occupation_masks']]
    savedpost=np.array(saved['outputs']['MIX']['poststate_real'])+1j*np.array(saved['outputs']['MIX']['poststate_imag'])
    producer_post_error=float(np.max(abs(after[np.ix_(order,order)]-savedpost))) if modulus==3 else None
    checks={'producer_poststate_independently_reproduced':producer_post_error is None or producer_post_error<1e-12,'uniform_isometry':bool(np.allclose(v.T@v,np.eye(20))),
       'same_H_intertwines':herror<1e-12,'every_C_E_MIX_K_intertwines':kerror<1e-12,
       'determinant_reference_normalized':abs(np.linalg.norm(slater)-1)<1e-12,
       'filled_reference_is_uniform_code_ground':reference_error<1e-12 and abs(eig[0]-reference_energy)<1e-12,
       'every_fixed_number_energy_reproduced':bool(np.allclose(eig,expected)),
       'global_ground_transfer_rejected':full_reference_minimum < -.01,
       'positive_excitation_reference_within_uniform_code':eig[0]>=reference_energy-1e-12,
       'complete_native_work':abs(sum(x['subnormalized_excitation_energy'] for x in work)-actualwork)<1e-12,
       'native_monitoring_heats_reference':actualwork>1e-3,'full_and_reduced_outputs_agree':outputerror<1e-12}
    checks={k:bool(value) for k,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    paths=[Path(__file__),HERE/'fermionic_uniform_reference.py',producer,HERE/'fermionic_shared_records.py',BASE/'candidate/oph_openwave/construction.py']
    return {'schema':'oph.c9.car_uniform_fiber_independent_review.v1','checks':checks,'full_dimension':full.n,'reduced_dimension':20,
        'fiber_cardinality':modulus,'producer_poststate_error':producer_post_error,'H_error':herror,'every_MIX_K_error':kerror,'full_output_error':outputerror,
        'full_global_ground_energy':full_ground,'full_reference_shifted_minimum':full_reference_minimum,'reference_energy':reference_energy,'reference_residual':reference_error,'fixed_number_excitation_gap':float(eig[1]-eig[0]),
        'native_reference_work':actualwork,'complete_native_reference_outcomes':work,
        'scope':'Declared fixedN3 CAR/uniformGauss-fiber preparation; exact same nativeinstrument. Positive reference ONLYwithin invariantuniformcode: fullsource has lowerenergyoutsidecode. Monitoringwork retained, reference notstationary. No globalvacuum or physicalspecies/spin/sea preparation inference.',
        'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'full_criterion_met':False}

if __name__=='__main__':
    r=run();HERE.joinpath('CAR_UNIFORM_FIBER_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'work':r['native_reference_work'],'gap':r['fixed_number_excitation_gap']}))
