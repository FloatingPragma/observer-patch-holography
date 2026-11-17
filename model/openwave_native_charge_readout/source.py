"""Complete coherent charge service compiled through actual source primitives.

Exact reversible circuits retain blank-scratch admission and its failure.
Primitive source exponentials are executed; numerical error bounds are
conditional on the checked routine bounds, not physical hardware guarantees.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
import math
import numpy as np

from model.openwave_active_native_charge import source as charge
from model.openwave_coherent_native_architecture import source as architecture
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler
from model.exchange_sector import statistics_service as source_service
from model.common_source import record_code as matching_code,parameter_selection as selection

SOURCE_DIM=512
MEMORY_DIM=8
SIZE=SOURCE_DIM*MEMORY_DIM


def source_address(pair):
    return circuit.encode(tuple(int(v in pair)for v in range(12)))


def shift_permutation(seam,current,extra=True):
    """All4096 source/memory addresses, with explicit identity padding."""
    if seam not in charge.SEAMS or current not in (-1,0,1):
        raise ValueError('one native seam and binary current required')
    result=list(range(SIZE))
    if current==0:
        return tuple(result)
    origin,destination=seam if current==1 else tuple(reversed(seam))
    for after in charge.PAIRS:
        if destination not in after or origin in after:
            continue
        before=tuple(sorted((set(after)-{destination})|{origin}))
        v=charge.voltage(before,after,seam,extra)
        address=source_address(after)
        for memory in range(6):
            result[address+SOURCE_DIM*memory]=address+SOURCE_DIM*((memory+v)%6)
    if sorted(result)!=list(range(SIZE)):
        raise ArithmeticError('complete memory service is not a permutation')
    return tuple(result)


def clean_circuit_columns(program):
    """Execute every admitted input through the emitted reversible program.

    The complete instrument first tests whether scratch is blank. Dirty
    scratch has a retained identity failure branch, so it is not discarded
    or asserted to have the clean-input logical action.
    """
    bits=program['address_bits']
    out=architecture.execute_basis(program,np.arange(1<<bits,dtype=np.int64))
    if any(out>>bits) or len(np.unique(out))!=len(out):
        raise ArithmeticError('scratch is not restored or clean action is nonunitary')
    return tuple(map(int,out))


def _hash(x):
    return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


@lru_cache(maxsize=1)
def compiled_transport():
    matrices,errors,costs=architecture.execute_primitives()
    reports=[]
    # Every seam/current has its own exact full-address map. Runtime uses
    # voltage-equivalent distinct permutations once, without persisted results.
    groups={}
    for seam in charge.SEAMS:
        for current in (-1,0,1):
            perm=shift_permutation(seam,current)
            groups.setdefault(perm,[]).append((seam,current))
    for perm,uses in groups.items():
        program=architecture.compile_permutation(perm)
        actual=clean_circuit_columns(program)
        if actual!=perm:
            raise ArithmeticError('emitted source circuit differs on a full input column')
        inverse={**program,'gates':list(reversed(program['gates']))}
        restored=architecture.execute_basis(inverse,np.array(actual,dtype=np.int64))
        assert tuple(map(int,restored))==tuple(range(SIZE))
        counts=Counter(k for k,_ in program['gates'])
        # Product expansion bounds arbitrary-reference operator error whenever
        # each executed routine has operator error at most the fixed ceiling.
        epsilon=F(sum(counts.values()),10**10)
        assert epsilon<1
        reports.append({'seam_current_uses':[[list(s),j]for s,j in uses],
            'source_memory_input_columns':SIZE,'matrix_unit_identities_implied':SIZE**2,
            'restored_scratch_bits':program['scratch_bits'],
            'dirty_scratch_identity_failure_projector_rank':(1<<(program['address_bits']+program['scratch_bits']))-SIZE,
            'complete_scratch_admission':True,'full_source_and_memory_padding_identity':True,
            'actual_emitted_forward_inverse_identity_on_all_admitted_columns':True,
            'inverse_source_primitives':'reverse order of the same self-inverse flip/CCX source routines; numerical inverse can instead use compiler.inverse on each underlying word',
            'gate_counts':dict(counts),'native_pulses':sum(counts[k]*costs[k]['native_pulses']for k in counts),
            'proper_code_service_events':sum(counts[k]*costs[k]['proper_code_service_events']for k in counts),
            'operator_error_upper_if_each_routine_error_at_most_1e_minus10':str(epsilon/(1-epsilon)),
            'schedule_sha256':_hash(program),'executed_clean_columns_sha256':_hash(actual)})
    return {'distinct_complete_permutations':len(groups),'seam_current_maps':90,
        'executed_full_input_columns':SIZE*len(groups),'primitive_source_errors':errors,
        'primitive_resources':costs,'compiled_services':reports,
        'execution_method':'Actual full M6 primitive words plus exact integer execution of every emitted reversible gate on all admitted source/memory columns. Dirty scratch is a separate complete failure outcome.',
        'full_dense_21_qubit_unitary_numerically_executed':False,
        'physical_hardware_error_certificate':False,'persisted_scientific_result_cache':False}


def _apply_small(state,matrix,targets,bits):
    return compiler._apply(state,matrix,list(reversed(targets)),bits)


def executed_witness_unitary(extra=True):
    """Actual source matrices execute a full16-address controlled charge shift."""
    perm=[]
    for address in range(16):
        source=address&1;memory=address>>1
        shift=(5 if extra else 3)if source==0 else 0
        perm.append(source+2*((memory+shift)%6 if memory<6 else memory))
    program=architecture.compile_permutation(tuple(perm))
    matrices,errors,costs=architecture.execute_primitives()
    bits=program['address_bits']+program['scratch_bits']
    out=np.eye(1<<bits,dtype=complex)
    for kind,targets in program['gates']:
        out=_apply_small(out,matrices[kind],targets,bits)
    expected=np.zeros_like(out)
    truth=architecture.execute_basis(program)
    expected[truth,np.arange(len(truth))]=1
    error=float(np.linalg.norm(out-expected,'fro'))
    if error>2e-8:
        raise ArithmeticError('actual complete coherent witness execution failed')
    return out,{'complete_clean_and_dirty_basis_columns':len(out),
        'full_operator_matrix_units':len(out)**2,'raw_source_execution_error_frobenius':error,
        'primitive_gate_count':len(program['gates']),
        'native_pulses':sum(costs[k]['native_pulses']for k,_ in program['gates']),
        'polar_projection_used':False,'extra_voltage':extra}


def charge_character(k):
    if type(k)is not int or k not in range(6):
        raise ValueError('six discrete charge characters required')
    # Diagnostic complex amplitudes; exact expressions use sixth roots.
    value=np.zeros(8,complex)
    value[:6]=np.exp(-2j*math.pi*k*np.arange(6)/6)/math.sqrt(6)
    return value


@lru_cache(maxsize=6)
def prepared_character(k):
    target=charge_character(k)
    amplitudes=np.zeros(8);amplitudes[:6]=1/math.sqrt(6)
    word=compiler.compile_real_state(amplitudes)
    phases=[-math.pi*k*m/3 if m<6 else 0. for m in range(8)]
    word+=compiler.compile_diagonal(phases)
    executed=compiler.execute(word,3)
    error=float(np.linalg.norm(executed[:,0]-target))
    assert error<2e-10
    return executed[:,0],{'character':k,'complete_unitary_columns':8,
        'prepared_character_error':error,'word_sha256':_hash(word),
        'resources':compiler.resources(word),
        'blank_admission_failure_rank':7,'all_nonblank_inputs_retained_as_failure':True}


@lru_cache(maxsize=1)
def prepared_source_and_export():
    # Source code is the literal inherited nine-bit native encoding.
    pairs=((0,1),(0,6));outpairs=((1,5),(5,6))
    real=np.zeros(512)
    for pair in pairs:real[source_address(pair)]=1/math.sqrt(2)
    prep=compiler.compile_real_state(real)
    blank=np.zeros(512);blank[0]=1
    actual=compiler.execute(prep,9,blank)
    error=float(np.linalg.norm(actual-real));assert error<2e-9
    base=np.zeros(512)
    for pair in outpairs:base[source_address(pair)]=1/math.sqrt(2)
    plusword=compiler.compile_real_state(base)
    quarters=[0]*512;quarters[source_address(outpairs[1])]=-1
    yword=plusword+source_service.diagonal_word(tuple(quarters))
    words={'X':plusword,'Y':yword}
    outputs={}
    for axis,word in words.items():
        columns=np.zeros((512,2),complex)
        for j,pair in enumerate(outpairs):columns[source_address(pair),j]=1
        # Actual source readout inverse acts on both source code columns.
        outputs[axis]=compiler.execute(compiler.inverse(word),9,columns)
    return {'prep_actual':actual,'readout_columns':outputs,
        'report':{'source_dimension':512,'prepared_pairs':[list(x)for x in pairs],
        'readout_pairs':[list(x)for x in outpairs],'preparation_error':error,
        'preparation_resources':compiler.resources(prep),
        'quadrature_readout_resources':{a:compiler.resources(w)for a,w in words.items()},
        'complete_source_outside_code_failure_rank':510,
        'readout_within_code':'inverse prepared eigenstate word; computational zero is plus, all other addresses are minus; initial code-complement projection supplies the separate failure',
        'word_sha256':{'prepare':_hash(prep),**{a:_hash(w)for a,w in words.items()}}}}


def quadrature_records(extra=True):
    U,execution=executed_witness_unitary(extra)
    # Complete source quadratures use actual compiled one-qubit routines.
    h=[compiler.pulse('I',-math.pi/2,[0]),compiler.pulse('H',math.pi/2,[0])]
    xread=compiler.execute(h,1)
    # Y=i(|0><1|-|1><0|); Z(-pi/4) before H yields this orientation.
    yword=[compiler.pulse('Z',-math.pi/4,[0])]+h
    yread=compiler.execute(yword,1)
    expected_x=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
    assert np.linalg.norm(xread-expected_x)<1e-10
    cos=[F(1),F(1,2),F(-1,2),F(-1),F(-1,2),F(1,2)]
    sin=['0','-sqrt(3)/2','-sqrt(3)/2','0','sqrt(3)/2','sqrt(3)/2']
    source_execution=prepared_source_and_export()
    source_initial=np.array([source_execution['prep_actual'][source_address(pair)]for pair in ((0,1),(0,6))])
    native_rows=charge.native.raw_local_branches('C')
    moving_slots=[i for i,r in enumerate(native_rows)if r.weight and charge.target((0,1),(0,5),r)==(1,5)]
    hold_slots=[i for i,r in enumerate(native_rows)if r.weight and charge.target((0,1),(0,5),r)==(0,1)]
    assert len(moving_slots)==len(hold_slots)==1
    rows=[]
    for k in range(6):
        input=np.zeros(U.shape[0],complex)
        memory,character_report=prepared_character(k)
        full_prepared=np.outer(memory,source_execution['prep_actual'])
        initial_projector=np.zeros((8,512),bool)
        for pair in ((0,1),(0,6)):initial_projector[:6,source_address(pair)]=True
        preparation_norm=float(np.sum(abs(full_prepared)**2))
        initial_failure=float(np.sum(abs(full_prepared[~initial_projector])**2))
        initial=np.kron(memory[:6],source_initial)
        input[:12]=initial
        admitted_norm=float(np.sum(abs(input)**2))
        assert abs(admitted_norm+initial_failure-preparation_norm)<2e-14
        output=U@input
        output_failure=float(np.sum(abs(output[12:])**2))
        admitted_output=np.zeros_like(output);admitted_output[:12]=output[:12]
        records=[]
        for axis,read in [('X',xread),('Y',yread)]:
            measured=_apply_small(admitted_output.reshape(-1,1),read,[0],int(math.log2(len(output)))).reshape(-1)
            for sign in range(2):
                probability=float(np.sum(abs(measured[np.arange(len(measured))%2==sign])**2))/4
                records.append({'axis':axis,'sign':1 if sign==0 else -1,'joint_probability':probability})
        # All512 actual readout outputs are summed, on the explicitly admitted
        # six-memory / blank-scratch branch. No amplitude is renormalized.
        source_coefficients=output[:12].reshape(6,2).T
        full_readout=[]
        for axis in ('X','Y'):
            amplitudes=source_execution['readout_columns'][axis]@source_coefficients
            plus=float(np.sum(abs(amplitudes[0])**2))/4
            minus=float(np.sum(abs(amplitudes[1:])**2))/4
            full_readout.extend([plus,minus])
        assert max(abs(records[n]['joint_probability']-full_readout[n])for n in range(4))<2e-9
        for n in range(4):records[n]['joint_probability']=full_readout[n]
        holding=admitted_norm/2
        moving_failure=output_failure/2
        records.extend([
            {'axis':'other_native_raw_outcomes','sign':0,'joint_probability':holding},
            {'axis':'moving_scratch_or_memory_failure','sign':0,'joint_probability':moving_failure},
            {'axis':'initial_source_or_memory_admission_failure','sign':0,'joint_probability':initial_failure}])
        raw_records=[]
        for slot,branch in enumerate(native_rows):
            for outcome in range(5):
                probability=0.
                if slot in moving_slots:
                    probability=records[outcome]['joint_probability']if outcome<4 else moving_failure
                if slot in hold_slots and outcome==4:probability=holding
                raw_records.append({'native_slot':slot,'native_record':list(branch.record),
                    'quadrature_or_failure':('X+','X-','Y+','Y-','source_memory_or_scratch_failure')[outcome],
                    'joint_probability':probability})
        raw_total=sum(x['joint_probability']for x in raw_records)+initial_failure
        total=sum(x['joint_probability']for x in records)
        assert len(raw_records)==125 and abs(raw_total-total)<2e-14
        assert abs(total-preparation_norm)<2e-9
        for n in range(4):
            angle=(5 if extra else 3)*math.pi*k/3
            wanted=(1+(1 if n%2==0 else -1)*(math.cos(angle)if n<2 else math.sin(angle)))/8
            if abs(records[n]['joint_probability']-wanted)>2e-9:
                raise ArithmeticError(('quadrature source sign or phase mismatch',k,n,records[n],wanted))
        rows.append({'charge_character':k,'character_preparation':character_report,
            'complete_unconditional_records':records,'all_125_raw_native_readout_slots':raw_records,
            'initial_admission_failure_probability':initial_failure,
            'actual_initial_admitted_norm':admitted_norm,
            'moving_output_scratch_or_memory_failure_probability':moving_failure,
            'preparation_norm_defect':preparation_norm-1,
            'complete_protocol_norm_defect_relative_to_preparation':total-preparation_norm,
            'total_record_weight':total,'numerical_state_renormalized':False,
            'exact_X_joint_probabilities':list(map(str,[(1+cos[k])/8,(1-cos[k])/8]))if extra else [str(F(1+(-1)**k,8)),str(F(1-(-1)**k,8))],
            'exact_Y_formula':'(1 +/- '+sin[k]+')/8'if extra else '1/8 each'})
    # Label C=+1 for k0, -1 for k2; score X+/- as +/-2, other outcomes0.
    score=lambda row:sum(2*x['sign']*x['joint_probability']for x in row['complete_unconditional_records']if x['axis']=='X')
    covariance=(score(rows[0])-score(rows[2]))/2
    expected=F(3,8)if extra else F(0)
    assert abs(covariance-float(expected))<2e-10
    return {'execution':execution,'complete_six_character_experiments':rows,
        'unconditional_label_record_covariance':covariance,
        'unconditional_label_record_covariance_exact':str(expected),
        'X_readout_native_pulses':len(h),'Y_readout_native_pulses':len(yword),
        'quadrature_choice':'fair classical apparatus coin; both choices and both signs retained',
        'complete_source_code_failure_effect':'I-P_{a,b}; outside the prepared two-state source code is a retained failure',
        'coherent_character_preparation':'Actual full8-address source word from admitted blank memory; every other memory input has a retained preparation-failure outcome.',
        'conditioned_branch_probabilities_relabelled_unconditional':False}


def uncompensated_clock_leakage():
    """Exact matching-code populations and cross-label coherence transfer."""
    weights,timeout=selection.completion_weights()
    cos=(F(1),F(1,2),F(-1,2),F(-1),F(-1,2),F(1,2))
    sin=(F(0),F(1,2),F(1,2),F(0),F(-1,2),F(-1,2))
    rows=[];witness=None
    supports=[tuple(tuple(i for i,n in enumerate(state)if n)for state in support)for support in matching_code.SUPPORTS]
    for law in ('C','E'):
        leakage=[[F(0)for _ in range(6)]for _ in range(6)]
        transfer=[[[[timeout,F(0)]for _ in range(6)]for _ in range(6)]for _ in range(6)]
        for seam,weight in weights.items():
            for branch in charge.native.raw_local_branches(law):
                if not branch.weight:continue
                bylogical=[]
                for logical,support in enumerate(supports):
                    phases=[]
                    for pair in support:
                        after=charge.target(pair,seam,branch)
                        if after is not None:phases.append(charge.voltage(pair,after,seam))
                    bylogical.append(phases)
                    n=len(phases)
                    if not n:continue
                    for k in range(6):
                        modulus=sum((cos[(k*(x-y))%6]for x in phases for y in phases),F(0))
                        loss=weight*branch.weight*(F(n,6)-modulus/F(6*n))
                        assert loss>=0
                        leakage[k][logical]+=loss
                assert len(set(map(len,bylogical)))==1
                n=len(bylogical[0])
                if not n:continue
                factor=weight*branch.weight/F(6*n)
                for k in range(6):
                    for l in range(6):
                        for m in range(6):
                            transfer[k][l][m][0]+=factor*sum((cos[k*(x-y)%6]for x in bylogical[l]for y in bylogical[m]),F(0))
                            transfer[k][l][m][1]+=factor*sum((sin[k*(x-y)%6]for x in bylogical[l]for y in bylogical[m]),F(0))
        assert all(x==0 for row in leakage for x in row)
        assert all(transfer[0][l][m]==[F(1),F(0)]for l in range(6)for m in range(6))
        for k in range(1,6):
            for l in range(6):
                assert transfer[k][l][l]==[F(1),F(0)]
                for m in range(l+1,6):
                    real,imag=transfer[k][l][m]
                    if witness is None and real<1:
                        witness={'law':law,'charge_character':k,'logical_labels':[l,m],
                            'coherence_multiplier_real_exact':str(real),
                            'coherence_multiplier_imag_over_sqrt3_exact':str(imag),
                            'real_coherence_loss_exact':str(1-real),
                            'matching_code_population_leakage_exact':'0'}
        rows.append({'law':law,'complete_selector_timeout_exact':str(timeout),
            'leakage_probability_by_charge_and_logical_label_exact':[[str(x)for x in row]for row in leakage],
            'logical_coherence_multipliers_real_plus_i_sqrt3_exact':[[[[str(t)for t in value]for value in row]for row in block]for block in transfer]})
    assert witness is not None
    return {'source_domain':'actual six matching supports of the512-address continuation, same complete selector and source recovery',
        'complete_uncompensated_one_event_transfer':rows,'nonzero_charge_witness':witness,
        'uncompensated_clock_theorem_transferred':False,
        'compensated_policy':'On each committed raw record, execute the actual compiled memory-feedback inverse before the original source recovery and matching mixer.',
        'compensated_identity':'(R_r tensor I6) V_r^dagger V_r (K_r tensor I6) = (R_r K_r) tensor I6 on every raw branch and arbitrary reference.',
        'padding_and_initialization':'Feedback fixes source and memory padding; inverse restores scratch on all4096 admitted columns. Dirty scratch remains an identity failure branch. Exact identity is on admitted clean scratch, with full failure complement retained.',
        'clock_time_recovery_cost':'Two feedback service invocations per active native event; no physical energy or second is inferred.',
        'active_charge_record_timing':'The source quadrature covariance is measured before inverse feedback; a recovery that undoes feedback also undoes this phase signature.',
        'simultaneous_charge_dependent_clock_frequency_derived':False}


def run():
    transport=compiled_transport()
    source_prep=prepared_source_and_export()
    coupled=quadrature_records(True)
    control=quadrature_records(False)
    assert coupled['unconditional_label_record_covariance_exact']=='3/8'
    assert control['unconditional_label_record_covariance_exact']=='0'
    return {'verdict':'PASS_COMPILED_ACTIVE_CHARGE_TRANSPORT_AND_COMPLETE_SOURCE_RECORDS',
        'shared_construction':architecture.construction_descriptor(),
        'record_policy':'Commit original native law/coin/raw outcomes, then apply the declared current-dependent coherent memory feedback. Source interference is within one raw branch; no erasure of distinct committed histories is used.',
        'complete_transport_compiler':transport,'native_source_preparation_and_readout':source_prep['report'],
        'actual_coupled_source_records':coupled,
        'altered_voltage_control':control,'matching_clock_recovery_join':uncompensated_clock_leakage(),
        'observation_domain':{'prepared_source_pairs':[[0,1],[0,6]],'moving_output_pairs':[[1,5],[5,6]],
            'native_law':'C','seam':[0,5],'moving_raw_record_probability_exact':'1/2',
            'source_initial_coherence_supplied':True,'memory_character_coherence_supplied':True,
            'character_initialization_conserves_charge_without_reservoir':False,
            'charge_preserving_domain':'Post-initialization transport, source readout and compensated recovery; character-preparation pulses are explicitly outside this conservation domain.',
            'all_other_native_raw_outcomes_retained':True,
            'complete_readout_effects':'(P+X)/4,(P-X)/4,(P+Y)/4,(P-Y)/4,I-P',
            'service_work':'Actual pulse and proper-code event counts; no physical energy normalization.'},
        'scope':{'memory_is_idle_uncoupled_label':False,'physical_charge_unit_selected':False,
            'localized_particle_core_derived':False,'physical_quark_string_derived':False,
            'rest_mass_GeV':None,'OpenWave_criterion_awards':0}}
