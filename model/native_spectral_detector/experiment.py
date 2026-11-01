"""Complete classical preparation ensemble and eighteen-bit clock execution."""
from fractions import Fraction as F
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import math
import time

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.ground_spectral_source import preparation as old
from model.reference_recovery.recovery import execute_complete
from model.cube_eta_protection import cube
from model.spin_resolved_spectrum import spin
from model.observer.source_scalar import compiler
from . import preparation as p, circuits as c, execution as e


@lru_cache(maxsize=1)
def ideal_ground_ensemble():
    states=[]
    for k,source in enumerate(spin.ground_family()):
        v=np.zeros(65536,complex)
        for mask,(a,b) in source.items():v[mask]=(a+b*math.sqrt(3))/(36*math.sqrt(spin.NORMS[k]))
        states.append(v)
    basis=np.column_stack(states);columns=[]
    for record in range(15):
        azimuth,polar=divmod(record,3);x=(-math.sqrt(3/5),0.,math.sqrt(3/5))[polar]
        co,si=math.sqrt((1+x)/2),math.sqrt((1-x)/2);phi=2*math.pi*azimuth/5
        coefficient=np.array([math.sqrt(math.comb(4,k))*co**(4-k)*si**k*np.exp(-1j*phi*(2-k))
                              for k in range(5)])
        columns.append(math.sqrt(p.SELECTOR_WEIGHTS[record]/90)*(basis@coefficient))
    return cube.operators().frame.matrix().conj().T@np.column_stack(columns)


@lru_cache(maxsize=1)
def prepared():
    blank=np.zeros(256);blank[0]=1
    tree=execute_complete(old.tree_word(),8,blank)
    selector=p.selector_run()
    probabilities=np.abs(selector.success[:,0])**2
    initial=old.embed_blank_modes(tree.success)
    runs=[]
    for record in range(15):
        state=math.sqrt(float(probabilities[record]))*initial
        runs.append(e.complete(p.rotation_blocks(record)+old.frame_blocks(),16,state))
    columns=np.column_stack([run.success[:,0] for run in runs])
    tree_norm=float(np.vdot(tree.success,tree.success).real)
    abort=tree.failure_probability+tree_norm*(selector.failure_probability+float(probabilities[15]))
    abort+=math.fsum(run.failure_probability for run in runs)
    return tree,selector,tuple(runs),columns,abort


def ideal_charge(ground):
    masks=np.arange(65536,dtype=np.int64);sign=(1-2*(masks&1))[:,None]
    occupied=(masks&2)!=0
    add=np.zeros_like(ground);remove=np.zeros_like(ground)
    add[occupied]=sign[occupied]*ground[masks[occupied]^2]
    remove[~occupied]=sign[~occupied]*ground[masks[~occupied]^2]
    return np.column_stack((add,remove))


@lru_cache(maxsize=1)
def charged():
    _,_,_,columns,abort=prepared()
    state=np.zeros((131072,15),complex);state[:65536]=columns
    run=e.complete(c.charge_blocks(),17,state)
    field=np.column_stack((run.success[:65536],run.success[65536:]))
    return run,field,abort+run.failure_probability


def _clock_chunk(task):
    indices,field=task
    initial=np.zeros((c.M*65536,len(indices)),complex);initial[:65536]=field
    run=e.complete(c.clock_blocks(),16+c.CLOCK_BITS,initial)
    return (indices,run.success,run.failure_weights,run.normalization_error,
            run.rounded_negative_weight_count,run.rounded_negative_weight_total)


@lru_cache(maxsize=1)
def clock_run():
    _,field,_=charged()
    initial=np.zeros((c.M*65536,30),complex);initial[:65536]=field
    groups=[list(range(6*j,6*j+6)) for j in range(5)]
    tasks=[(group,np.ascontiguousarray(field[:,group])) for group in groups]
    # These are orthogonal classical preparation/charge labels. Linear
    # propagation separates by columns; failure probabilities add. Every
    # worker still uses every field/clock row and every source instrument.
    success=np.zeros_like(initial);table=None;errors=[];rounded=0;negative=0.
    with ProcessPoolExecutor(max_workers=5,mp_context=multiprocessing.get_context('spawn')) as pool:
        for indices,columns,weights,error,count,total in pool.map(_clock_chunk,tasks):
            success[:,indices]=columns
            table=np.array(weights) if table is None else table+weights
            errors.append(error);rounded+=count;negative+=total
    failures=float(table.sum())
    error=abs(float(np.vdot(success,success).real)+failures-float(np.vdot(initial,initial).real))
    if error>e.ADMISSION or sum(errors)>e.ADMISSION:raise ArithmeticError('parallel complete ensemble accounting failed')
    from model.reference_observation.observation import frozen
    return e.Run(18,c.clock_blocks(),frozen(initial),frozen(success),frozen(table),failures,error,rounded,negative)


@lru_cache(maxsize=1)
def ideal_clock():
    field=ideal_charge(ideal_ground_ensemble())
    h=cube.operators().H.matrix(cube.operators().scale)
    evolved=np.zeros((c.M,65536,30),complex)
    ground_energy=-math.sqrt(3)/4
    for j in range(c.M):
        duration=j*float(c.TAU)
        for number,columns in ((5,slice(0,15)),(3,slice(15,30))):
            masks=cube.weights(number)
            evolved[j,masks,columns]=expm_multiply(-1j*duration*h[masks][:,masks],field[masks,columns])*np.exp(1j*duration*ground_energy)/math.sqrt(c.M)
    fourier=np.exp(2j*np.pi*np.arange(c.M)[:,None]*np.arange(c.M)[None,:]/c.M)/math.sqrt(c.M)
    return np.tensordot(fourier,evolved,axes=(1,0)).reshape(c.M*65536,30)


def pointer_probabilities(state):
    if state.shape!=(c.M*65536,30):raise ValueError('all clock, field and preparation rows required')
    data=state.reshape(c.M,65536,2,15)
    # Branch columns are ordered [15 addition,15 removal], independently
    # of the field and clock tensor factors.
    binary=np.sum(np.abs(data)**2,axis=(1,3))
    rows=[]
    for tag in range(6):
        for raw in np.ndindex(*(6,)*c.CLOCK_BITS):
            valid=tag<2 and all(x<2 for x in raw)
            index=sum(value<<j for j,value in enumerate(raw)) if valid else 0
            rows.append({'charge_record':tag,'clock_record':list(raw),
                         'probability':float(binary[index,tag]) if valid else 0.})
    return rows


def joint_pointer_probabilities(state):
    if state.shape!=(c.M*65536,30):raise ValueError('all clock, field and preparation rows required')
    binary=np.sum(np.abs(state.reshape(c.M,65536,2,15))**2,axis=1)
    return [{'selector_record':record,'charge_record':tag,'clock_record':list(raw),
             'probability':float(binary[sum(value<<j for j,value in enumerate(raw)),tag,record])
                if tag<2 and all(x<2 for x in raw) else 0.}
            for record in range(15) for tag in range(6)
            for raw in np.ndindex(*(6,)*c.CLOCK_BITS)]


def selector_pointer_probabilities():
    run=p.selector_run();probabilities=np.abs(run.success[:,0])**2
    return [{'record':list(raw),
             'probability':float(probabilities[sum(value<<j for j,value in enumerate(raw))])
                if all(x<2 for x in raw) else 0.,
             'continues':all(x<2 for x in raw) and sum(value<<j for j,value in enumerate(raw))<15}
            for raw in np.ndindex(*(6,)*p.SELECTOR_BITS)]


def preparation_diagnostics():
    tree,selector,runs,columns,abort=prepared()
    ideal=ideal_ground_ensemble()
    error=math.sqrt(float(np.linalg.norm(columns-ideal)**2)+abort)
    target=p.certificate()
    proved=F(target['complete_mixed_ground_preparation_error'])
    if error>float(proved)+e.ADMISSION:raise ArithmeticError('mixed-ground execution exceeds proof')
    return {'all15_conditional_ground_words_executed':True,
            'complete_preparation_error_diagnostic':error,
            'source_first_failure_and_invalid_selector_mass':abort,
            'ideal_ensemble_gram_rank':5,
            'native_prepared_state_columns':15,
            'complete_probability':float(np.vdot(columns,columns).real)+abort,
            'selector_valid_record_probabilities':[float(abs(x)**2) for x in selector.success[:,0]],
            'invalid_selector_record':15,
            'native_ground_rank_projection_used':False}


def arrays():
    tree,selector,runs,columns,_=prepared();charge,_,_=charged();clock=clock_run()
    tree_norm=float(np.vdot(tree.success,tree.success).real)
    result={'ground_tree_failure_weights':tree.failure_weights,
            'selector_failure_weights':tree_norm*selector.failure_weights,
            'prepared_ground_ensemble':columns,
            'charge_instrument_failure_weights':charge.failure_weights,
            'phase_clock_failure_weights':clock.failure_weights,
            'phase_clock_complete_success_state':clock.success}
    result.update({f'conditional_ground_{j}_failure_weights':run.failure_weights for j,run in enumerate(runs)})
    return result


def report():
    from model.protected_particle_readout.experiment import source_record
    from . import bounds
    tree,selector,runs,_,preparation_abort=prepared();charge,_,abort=charged();clock=clock_run()
    tree_norm=float(np.vdot(tree.success,tree.success).real)
    selector_record=source_record('selector',selector)
    selector_record['first_failure_probability']*=tree_norm
    selector_record['input_branch_weight']=tree_norm
    records=[source_record('ground_tree',tree),selector_record]
    records.extend(e.record(f'conditional_ground_{j}',run) for j,run in enumerate(runs))
    records.extend([e.record('charge_instrument',charge),e.record('phase_clock',clock)])
    probabilities=pointer_probabilities(clock.success)
    joint=joint_pointer_probabilities(clock.success)
    ideal=ideal_clock();ideal_probabilities=pointer_probabilities(ideal)
    total_abort=abort+clock.failure_probability
    total=math.fsum(row['probability'] for row in probabilities)+total_abort
    error=math.sqrt(float(np.linalg.norm(clock.success-ideal)**2)+total_abort)
    variation=(math.fsum(abs(a['probability']-b['probability']) for a,b in zip(probabilities,ideal_probabilities,strict=True))+total_abort)/2
    proof=bounds.certificate();bound=F(proof['complete_state_and_total_variation_error'])
    if abs(total-1)>e.ADMISSION or error>float(bound)+e.ADMISSION or variation>float(bound)+e.ADMISSION:
        raise ArithmeticError('complete native clock exceeds its source theorem')
    intervals=bounds.ideal_probabilities()
    for row in intervals:
        raw=[(row['index']>>j)&1 for j in range(c.CLOCK_BITS)]
        actual=next(x['probability'] for x in probabilities if x['charge_record']==row['charge_record'] and x['clock_record']==raw)
        lo,hi=map(F,row['native_unconditioned_probability_interval'])
        if not float(lo)-e.ADMISSION<=actual<=float(hi)+e.ADMISSION:
            raise ArithmeticError('native probability outside exact spectral interval')
    return {'invariant_ground_preparation_executed':True,
            'complete_binary_phase_clock_source_executed':True,
            'field_bits':16,'clock_bits':c.CLOCK_BITS,'full_quantum_rows':c.M*65536,
            'classical_selector_and_charge_columns':30,
            'column_partition_workers':5,
            'column_partition_drops_no_coherences':'Columns are orthogonal retained classical selector/charge records; each complete quantum column is propagated.',
            'preparation':preparation_diagnostics(),
            'all_selector_pointer_outcomes':selector_pointer_probabilities(),
            'all_joint_selector_charge_clock_outcomes':joint,
            'all_charge_clock_marginal_outcomes':probabilities,
            'independent_original_H_marginal_outcomes':ideal_probabilities,
            'source_records':records,
            'invalid_selector_stop_probability':tree_norm*float(abs(selector.success[15,0])**2),
            'first_failure_and_invalid_selector_stop_mass':total_abort,
            'complete_probability':total,
            'complete_output_error_diagnostic':error,
            'total_variation_error_diagnostic':variation,
            'host_numerical_admission':e.ADMISSION,
            'host_diagnostics_prove_literal_control_accuracy':False,
            'particle_number_or_spectral_projection_used':False,
            'success_charge_or_bin_postselection_used':False,
            'all78_lines_separately_resolved_by_native_bins':False,
            'finite_shot_measurement_records_generated':False,
            'passive_repair_applied_during_active_clock':False,
            'reference_action':'Original zero-chemical-potential H; not the separately centered H-N/2 reference.',
            'physical_species_or_relativistic_mass_identified':False}
