"""Replay same-action full-Fock reference, source recovery and entropy."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F

import numpy as np

from . import modular as m, controls
from ..defect_field_clock import clock
from ..defect_field_bridge.dynamics_packet import compare_json,compare_arrays

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUTPUT=HERE/'results.json'
ARRAYS=HERE/'results.npz'


def plain(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [plain(v) for v in value]
    return value


def clock_records(histories,arrays):
    sequences=[(a,b) for a in (-1,0,1) for b in (-1,0,1)]
    runs={s:clock.history(s,F(1)) for s in sequences}
    arrays['clock_outputs']=np.array([[runs[s]['outputs'][x,y] for x in (0,1) for y in (0,1)] for s in sequences])
    result={}
    for tag,row in histories.items():
        entropy=0.;classical=0.;mass=0.
        for index,field in enumerate(row['outputs']):
            first,second=divmod(index,32)
            sequence=(m.joint.current(first%16),m.joint.current(second%16))
            p=row['probabilities'][index];sf=m.entropy_from_amplitude(field)
            for state in runs[sequence]['outputs'].values():
                q=float(np.vdot(state,state).real);prob=p*q
                entropy+=q*sf+p*m.entropy_from_amplitude(state)
                if prob>0:classical-=prob*np.log(prob)
                mass+=prob
        result[tag]={'joint_native_and_clock_records':4096,
                     'joint_probability_mass':float(mass),
                     'continuing_component_field_and_clock_work_CQ_entropy':float(entropy),
                     'continuing_component_including_purifiers_CQ_entropy':float(classical),
                     'first_export_abort_entropy_included':False,
                     'subsystems':'These entropies concern the unnormalized continuing component only. Work entropy traces declared purifiers for this observable; stored continuing quantum factors retain them and are conditionally pure. First-export abort weights/Grams and reconstructible prefix outputs remain separate.',
                     'clock_advanced_once_per_physical_interval':True}
    return result


def build():
    data=m.reference();regional=m.regional_certificate();prepared=m.prepare_reference();recovery=m.source_recovery()
    arrays={'H':data['H'],'reference':data['rho'],'modular_K':data['K'],
            'source_purification':prepared['columns'],
            'preparation_failure_weights':prepared['run'].failure_weights,
            'preparation_failure_gram':prepared['run'].failure_gram,
            'recovery_isometry':recovery['run'].success,'recovery_target':recovery['target'],
            'recovery_failure_weights':recovery['run'].failure_weights,
            'recovery_failure_gram':recovery['run'].failure_gram,
            'regional_reference':regional['outer'].target_reference,
            'regional_K':regional['regional_K'],'bare_regional_H':regional['bare_regional_H']}
    histories={};reports={}
    for law in ('C','E'):
        for policy in ('coherent','flagged'):
            tag=policy+'_'+law;row=m.history_entropy(law,policy);histories[tag]=row
            reports[tag]={k:v for k,v in row.items() if not isinstance(v,np.ndarray)}
            for k,v in row.items():
                if isinstance(v,np.ndarray):arrays[tag+'_'+k]=v
    clock_report=clock_records(histories,arrays)
    from ..defect_field_energy import energy
    files=list(HERE.glob('*.py'))+[HERE/'README.md',controls.PATH,Path(energy.__file__),
          ROOT/'model/defect_field_bridge/joint.py',ROOT/'model/defect_field_bridge/joint_results.json',
          ROOT/'model/defect_field_interactions/interaction.py',ROOT/'model/defect_field_clock/clock.py',
          ROOT/'model/field_thermal/thermal.py',ROOT/'model/reference_recovery/recovery.py',
          ROOT/'model/reference_observation/observation.py',ROOT/'model/modular_flow/flow.py',
          ROOT/'model/common_source/mass_response/history.py',ROOT/'model/observer/source_scalar/compiler.py',
          ROOT/'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
          ROOT/'Lean/Thermodynamics/EinsteinPremiseLink.lean']
    pins=m.joint.bindings()
    pins.update({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    return {'schema':'oph.same_action_full_fock_modular.v1',
            'source_sha256':dict(sorted(pins.items())),
            'domain':'Full16-dimensional four-mode Fock algebra of the actual two-site joint source factors; all regional number sectors remain.',
            'experiment_relation':'A new declared Gibbs preparation uses the same complete field instrument; the original localized-input joint witness is unchanged.',
            'full_clock_reference_stationarity_claimed':False,
            'beta_GeV_inverse':str(m.BETA),'inverse_temperature_is_declared_design_input':True,
            'reference_minimum_eigenvalue':data['minimum_eigenvalue'],
            'modular_identity':'-log(rho)=beta H+log(Z) I; sigma_s(A)=rho^(is) A rho^(-is) has Hamiltonian time t=-beta*s',
            'physical_clock':'The existing Omega_ref=1/4 GeV clock measures co-driven field action time. A field Gibbs inverse temperature is not a second clock or a species-dependent time conversion.',
            'reference_preparation':{'source_resources':prepared['resources'],
              'compiler_resource_certificate':prepared['preparation_certificate'],
              'complete_vector_error':prepared['complete_error'],
              'initialized_qubits':8,'preparation_calls':1,
              'retained_purifier_dimension':16,
              'purifier_charge_convention':'dual charge -Q^T; (Q tensor I + I tensor -Q^T)|sqrt(rho)> = 0',
              'proper_code_failure_family':'Every prefix and complement row remains defined by the frozen source word and execute_complete; all weights and Grams are retained.'},
            'native_Petz_recovery':{'source_resources':recovery['resources'],
              'complete_isometry_error':recovery['complete_error'],
              'full_input_dimension':4,'output_field_dimension':16,'retained_environment_dimension':4,
              'input_qubits':2,'fresh_initialized_qubits':4,'total_work_qubits':6,
              'generic_state_recovery_claimed':False,'regional_retraction_claimed':False},
            'regional_observations':{k:v for k,v in regional.items() if not isinstance(v,np.ndarray) and k not in ('outer','inner','direct')},
            'exact_boundary_obstruction':m.exact_boundary_obstruction(),
            'entropy_identity':'S(sigma)-S(rho)=beta Tr[H(sigma-rho)]-D(sigma||rho) for normalized states; numerical trace drift adds (Tr(sigma)-1)log(Z).',
            'entropy_component':'Subnormalized continuing source component, with trace corrections; no first-export abort entropy is added to these values.',
            'first_export_abort_entropy_included':False,
            'first_export_failures_zero_in_exact_proper_code':True,
            'remainder_identity':'The entropy Bregman remainder is D(rho+epsilon*delta||rho)-epsilon Tr(delta), including nonzero numerical trace drift.',
            'histories':reports,'clock_CQ_records':clock_report,
            'source_controls':'Frozen emitted native words are regenerated only by --freeze-controls; every replay independently checks their complete target and current source action.',
            'source_error_admission':m.ADMISSION,
            'no_output_normalization_or_eigenvalue_clipping':True,
            'remaining_geometric_obligations':['No geometric diamond kernel or bulk/edge/central split has been identified.','No Lorentz stress/null tomography, fixed-volume area response, shrinking-ball family or physical Newton normalization is supplied.','The selected two-site factors are not an invariant restriction of the full ring.','Gibbs reference stationarity under repaired driven histories is tested, not assumed; a global Gibbs entropy identity does not imply generalized-entropy stationarity.'],
            'EinsteinRegisterData_constructed':False,'canonical_physical_premises_discharged':0},arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--freeze-controls',action='store_true');group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.freeze_controls:
        controls.freeze();print('Frozen emitted Gibbs/Petz source controls.');return
    report,arrays=build();report=plain(report)
    if args.write:
        OUTPUT.write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ARRAYS,**arrays)
    elif args.check:
        compare_json(report,json.loads(OUTPUT.read_text()))
        with np.load(ARRAYS,allow_pickle=False) as saved:compare_arrays(arrays,{k:saved[k] for k in saved.files})
    print('Same-action faithful reference, native Petz service and quantum/CQ entropy verified.')


if __name__=='__main__':main()
