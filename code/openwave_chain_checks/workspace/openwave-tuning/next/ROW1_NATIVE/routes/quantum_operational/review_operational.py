"""Coarse quantum charge PVM: source records, executed processor, native pointer."""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
import operational_readout as op
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];C10=ROOT/'openwave-tuning/next/C10';SPATIAL=HERE.parent/'quantum_spatial'
sys.path.insert(0,str(C10))
from joint_source import JointNativeSource


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def continue_injections(model,rho,u,law):
    """Independent sparse partial-injection contraction of actual native Ks."""
    coherent=u@rho@u.conj().T;out=np.zeros_like(coherent);rows=[]
    for operators in model.field.operators_by_law[law]:
        for meta,k in operators:
            target,source=np.nonzero(k);values=k[target,source]
            if len(set(target))!=len(target) or len(set(source))!=len(source):raise ValueError('not a partial injection')
            branch=np.zeros_like(coherent)
            branch[np.ix_(target,target)]=coherent[np.ix_(source,source)]*values[:,None]*values.conj()[None,:]/3
            out+=branch
            rows.append({**meta,'joint_probability':.29*float(np.trace(branch).real),'unnormalized_source_energy':.29*float(np.einsum('ij,ji->',model.hamiltonian,branch).real)})
    result=.71*coherent+.29*out
    rows.append({'status':'no_commit','joint_probability':.71*float(np.trace(rho).real),'unnormalized_source_energy':.71*float(np.einsum('ij,ji->',model.hamiltonian,coherent).real)})
    if abs(sum(r['joint_probability'] for r in rows)-np.trace(rho).real)>1e-11:raise ValueError('continuation not complete')
    return result,rows


def run():
    files=[Path(__file__),HERE/'operational_readout.py',SPATIAL/'quantum_spatial.py',SPATIAL/'QUANTUM_SPATIAL.json',SPATIAL/'quantum_spatial_states.npz',C10/'joint_source.py',C10/'charged_reservoir.py',C10/'modular_field_adapter.py',C10/'joint_source_states.npz',ROOT/'reverse-engineering-reality/model/observer/source_scalar/compiler.py',ROOT/'reverse-engineering-reality/model/reference_recovery/recovery.py',ROOT/'reverse-engineering-reality/model/native_interface/pointer.py',ROOT/'reverse-engineering-reality/model/native_clock/configuration.py',HERE.parent/'auditor_bridge/native_pointer_bridge.py',HERE.parent/'auditor_bridge/NATIVE_POINTER_BRIDGE.json']
    files += [ROOT/'reverse-engineering-reality'/op.compiler.RESPONSE,ROOT/'reverse-engineering-reality'/op.compiler.CONTROL_RECEIPT,ROOT/'reverse-engineering-reality/model/observer/native.py',ROOT/'reverse-engineering-reality/model/observer/source_discrimination/circuit.py',ROOT/'openwave-tuning/next/C9/fermionic_shared_records.py',ROOT/'openwave-tuning/candidate/oph_openwave/fiber_preparation.py',ROOT/'openwave-tuning/validation/bridges/providers/fock.py']
    before={str(p.relative_to(ROOT)):digest(p) for p in files}
    model=JointNativeSource();archive=np.load(SPATIAL/'quantum_spatial_states.npz');report=json.loads((SPATIAL/'QUANTUM_SPATIAL.json').read_text())
    frozen=np.load(C10/'joint_source_states.npz');assert np.array_equal(model.hamiltonian,frozen['hamiltonian'])
    records=np.array([f for _,f in model.states]);boundaries=frozen['boundary_B']
    decoded=op.decoded_occupations(records,boundaries)
    actual=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in model.states])
    if not np.array_equal(decoded,actual):raise ValueError('record-only decoder not actual quantum occupation')
    region=tuple(report['parameters']['region']);q=decoded[:,region].sum(axis=1);rho=archive['input_state'];indices=op.encoded_indices(model.states,region)
    ideal=op.coarse_branches(rho,q)
    generator_error=max(float(np.max(abs(ideal[value]-archive[f'coarse_Q{value}']))) for value in range(4))
    oracle,expected,admitted,oracle_error=op.compiled_oracle(3)
    two,_,two_admitted,two_error=op.compiled_oracle(2)
    instrument,errors,pointer_run,pointer_receipt=op.native_pointer_instrument(oracle)
    arrays={'decision_success':oracle.success,'decision_failure_effect':oracle.failure_gram,'pointer_failure_effect':pointer_run.failure_gram,'coherent_code_indices':indices}
    compiled_errors={};fine_count={}
    for law,data in instrument.items():
        compiled_errors[law]=max(float(np.max(abs(op.apply_local_channel(rho,data['branches'][value],indices)-ideal[value]))) for value in range(4))
        fine_count[law]=len(data['fine_native_records'])
        for value,ks in data['branches'].items():arrays[f'{law}_Q{value}_local_kraus']=np.stack(ks)
    # Reset primitives are an existing complete native record policy. Separate
    # independent apparatus preparation is an explicit hypothesis.
    bridge=json.loads((HERE.parent/'auditor_bridge/NATIVE_POINTER_BRIDGE.json').read_text())
    reset_per_bit=Fraction(bridge['joint_reset_success_probability']);reset=reset_per_bit**2
    local_rho=rho[np.ix_(indices.ravel(),indices.ravel())].reshape(8,48,8,48)
    reduced=np.einsum('arcr->ac',local_rho)
    complete_laws={}
    for law,data in instrument.items():
        source_abort=float(np.trace(data['failure_effects']['oracle']@reduced).real)
        first_abort=float(np.trace(data['failure_effects']['pointer_first']@reduced).real)
        second_abort=float(np.trace(data['failure_effects']['pointer_second']@reduced).real)
        success=[float(np.trace(op.apply_local_channel(rho,data['branches'][value],indices)).real) for value in range(4)]
        total=float(1-reset)+float(reset)*(sum(success)+source_abort+first_abort+second_abort)
        complete_laws[law]={'reset_timeout_exact':str(1-reset),'reset_success_exact':str(reset),'unnormalized_charge_probabilities_after_reset':[float(reset)*p for p in success],'unconditional_source_word_abort':float(reset)*source_abort,'unconditional_first_pointer_export_abort':float(reset)*first_abort,'unconditional_second_pointer_export_abort':float(reset)*second_abort,'total_probability':total}
        if abs(total-1)>1e-11:raise ValueError('compiled complete instrument lost outcomes')
    # The measured function never exposes fine F or B labels. Explicit
    # controls demonstrate why measuring them first changes the instrument.
    control=archive['fine_forget_control_state'];field_equal=np.all(records[:,None,:]==records[None,:,:],axis=2)
    control_coarse=sum(op.coarse_branches(control,q).values());fine=control*field_equal
    fine_distance=float(np.linalg.norm(control_coarse-fine,ord='nuc')/2)
    crosscoarse=archive['crossB_coarse_control_poststate'];crossjoint=archive['crossB_joint_control_poststate']
    cross_distance=float(np.linalg.norm(crosscoarse-crossjoint,ord='nuc')/2)
    # Independent continued-state contraction on generator's selected JOINT
    # branch, then an additional genuinely COARSE charge branch on same H.
    u=expm(-.37j*model.hamiltonian)
    chosen=report['selected_joint_continuation']['branch'];joint=archive[chosen['poststate_key']]
    joint_after,_=continue_injections(model,joint,u,'MIX')
    continuation_error=float(np.max(abs(joint_after-archive['continued_joint_poststate'])))
    selected=max(ideal,key=lambda value:float(np.trace(ideal[value]).real));branch=ideal[selected]
    continuation={}
    for law in ('C','E','MIX'):
        after,ledger=continue_injections(model,branch,u,law)
        energy_rows=model.energy_number_readout(after)
        before_energy=float(np.einsum('ij,ji->',model.hamiltonian,branch).real);after_energy=float(np.einsum('ij,ji->',model.hamiltonian,after).real)
        if abs(sum(row['probability'] for row in energy_rows)-np.trace(branch).real)>1e-11:raise ValueError('energy branch probability lost')
        if abs(sum(row['probability']*row['energy'] for row in energy_rows)-after_energy)>1e-10:raise ValueError('energy branch ledger differs')
        continuation[law]={'joint_probability':float(np.trace(branch).real),'complete_native_outcomes':ledger,'complete_final_number_energy_records':energy_rows,'source_energy_before':before_energy,'source_energy_after':after_energy,'native_monitoring_work':after_energy-before_energy}
        if law=='MIX':arrays['continued_coarse_charge_state']=after
    arrays['selected_coarse_charge_state']=branch
    np.savez_compressed(HERE/'operational_states.npz',**arrays)
    afterpins={str(p.relative_to(ROOT)):digest(p) for p in files}
    checks={'same_frozen_C10_H':bool(np.array_equal(model.hamiltonian,frozen['hamiltonian'])),'all_record_B_Pauli_decoders_exact':bool(np.array_equal(decoded,actual)),'whole_coarse_PVM_matches_generator_archive':generator_error<1e-12,'actual_bounded_source_decision_oracles':oracle_error<1e-10 and two_error<1e-10,'actual_C_E_pointer_preserves_all_matrix_units':all(e['coarse_instrument_matrix_unit_error']<1e-10 for e in errors.values()),'complete_export_and_reset_outcomes':all(abs(x['total_probability']-1)<1e-10 for x in complete_laws.values()),'full384_coarse_states_preserved':max(compiled_errors.values())<1e-10,'fine_F_and_extra_B_measurements_rejected':abs(fine_distance-.5)<1e-12 and abs(cross_distance-.5)<1e-12,'independent_sameH_native_continuation':continuation_error<1e-11,'all_C_E_MIX_energy_branch_outcomes_retained':len(continuation)==3,'sources_unchanged':before==afterpins}
    return {'schema':'oph.C10.operational_quantum_charge.v1','checks':checks,'source_sha256':before,'decision_oracles':{'three_site':{'native_M6_pulses':len(oracle.word),'work_bits':5,'input_columns':8,'matrix_error':oracle_error,'admitted_word':admitted},'two_site':{'native_M6_pulses':len(two.word),'work_bits':4,'input_columns':4,'matrix_error':two_error,'admitted_word':two_admitted}},'region':list(region),'native_pointer_template_errors':errors,'native_pointer_fine_record_count':fine_count,'native_pointer_fine_records':{law:data['fine_native_records'] for law,data in instrument.items()},'native_pointer_source_word':[{'kind':kind,'angle':angle,'targets':list(targets)} for kind,angle,targets in pointer_run.word],'finite_native_reset_contract':bridge['reset'],'compiled_coarse_fullstate_error':compiled_errors,'unconditional_apparatus_laws':complete_laws,'mutations':{'fine_F_forget_trace_distance':fine_distance,'extra_B_forget_crossnumber_trace_distance':cross_distance},'same_source_continuation':{'independent_selected_joint_max_error':continuation_error,'coarse_charge_outcome':selected,'coherent_H_parameter':.37,'native_probability':.29,'laws':continuation},'clock_and_resources':{'three_site_decision_source_pulses':114,'two_site_decision_source_pulses':41,'two_native_pointer_coupling_source_pulses':12,'native_terminal_pointer_attempts':2,'native_auxiliary_resets':4,'native_reset_attempt_cap':2*sum(x['certificate']['native_attempts_max'] for x in bridge['reset']),'source_during_apparatus':'paused under supplied addressed schedule','physical_clock_conversion':'not derived; M6 integrated pulse angle, native attempt count, and H hopparameter are separately typed'},'apparatus_ancestry':{'coherent_source_transfer':'Explicit finite Gauss-code reindexing |mask,F> to |occupation bits> tensor |cycle>; declared reference-preserving proper-code load/export and routing service. No fine F or B measurement is performed.','implemented':'Bounded114pulse (or41pulse two-site) M6 population-count oracle followed by two actual native C/E pointer instruments, with all prefix-export and finite-native-reset failures. Whole-degenerate charge class source coherences preserved.','fine_work':'Occupation/cycle register is a coherent reindexing, not an exposed measurement record. Source bits remain intact; only sum pointer outputs are read.','failure_continuation':'Processor/export failures are retained flagged stop spaces with their exact lazy prefix Kraus families; no artificial C10 H-energy or successful native continuation is assigned to an abort space.','energy_scope':'Complete same-H source/number/energy records on conclusive coarse branches; no physical total apparatus-energy or native preparation-cost claim.'},'archive':'operational_states.npz','archive_sha256':digest(HERE/'operational_states.npz'),'official_rows_changed':[]}

if __name__=='__main__':
 r=run()
 if not all(r['checks'].values()):raise ValueError(r['checks'])
 (HERE/'OPERATIONAL_QUANTUM_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'compiled_coarse_fullstate_error':r['compiled_coarse_fullstate_error'],'native_pointer_template_errors':r['native_pointer_template_errors']},indent=2))
