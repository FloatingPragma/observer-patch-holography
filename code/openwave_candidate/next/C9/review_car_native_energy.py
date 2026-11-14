"""Independent complete native-event then global-energy eigenspace readout."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'candidate'))
from fermionic_shared_records import FermionicSharedRecords
from oph_openwave.construction import Construction


def run():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),mass=.23),7,6)
    values,vectors=np.linalg.eigh(m.hamiltonian);psi=vectors[:,0];groups=[]
    for i,x in enumerate(values):
        if groups and abs(x-groups[-1][0])<1e-10:groups[-1][1].append(i)
        else:groups.append([float(x),[i]])
    source=HERE/'CHARACTER_MINIMA.json';saved=json.loads(source.read_text());outputs={};error=0.
    for law in ('C','E','MIX'):
        joint=[]
        for e,operations in enumerate(m.operators[law]):
            for meta,k in operations:
                amplitudes=vectors.conj().T@(k@psi/np.sqrt(3))
                joint.append({**meta,'energy_probabilities':[float(sum(abs(amplitudes[ids])**2)) for _,ids in groups]})
        probabilities=np.sum([r['energy_probabilities'] for r in joint],axis=0)
        first=min(x-values[0] for (x,_),p in zip(groups,probabilities) if x-values[0]>1e-10 and p>1e-12)
        expected=[g['probability'] for g in saved['native_ground_energy_readout'][law]['all_energy_groups']]
        error=max(error,float(np.max(abs(probabilities-expected))))
        outputs[law]={'complete_joint_outcomes':joint,'complete_energy_marginal':probabilities.tolist(),
           'first_numerically_supported_excitation':float(first),
           'mean_work':float(sum((x-values[0])*p for (x,_),p in zip(groups,probabilities)))}
    checks={'every_C_E_MIX_complete_joint_law':all(abs(sum(v['complete_energy_marginal'])-1)<1e-12 for v in outputs.values()),'producer_marginals_reproduced':error<1e-12,'native_accessible_gap_distinct_from_globalgap':all(abs(v['first_numerically_supported_excitation']-2.052218311973655)<1e-12 for v in outputs.values())}
    if not all(checks.values()):raise ValueError(checks)
    return {'schema':'oph.c9.car_native_energy_independent_review.v1','checks':checks,'producer_marginal_error':error,'outputs':outputs,
      'supplied':{'eigenvalue_grouping_tolerance':1e-10,'numerical_support_threshold':1e-12},
      'scope':'Full120state globalenergy eigenspace PVM after eachnativeevent, with allzero/degenerate groups. Actualcharacterpreservation explains inaccessiblelowergloballevels; no physicalmass identification.',
      'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'fermionic_shared_records.py',HERE/'character_minima.py',source)},'full_criterion_met':False}

if __name__=='__main__':
 r=run();HERE.joinpath('CAR_NATIVE_ENERGY_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
