"""Destructive spectral filters using existing coherent source services."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from model.native_coexistent_effects import menus
from model.native_effect_synthesis import synthesis as syn, intervals as iv
from model.native_transfer_control import control as bounds
from . import spectral

STAGES=66
WORD_PULSES=menus.WORD_PULSES


def vector(index):
    if type(index) is not int or not 0<=index<STAGES:
        raise ValueError('an original energy label from zero through65 required')
    columns,_=spectral.load()
    return tuple((x,F()) for x in columns[index])


def word(index):
    """At survival zero the first filter already has exact residual I-P."""
    return menus.phase_word(vector(index),F(1,2),iv.pi_grid())


def budget():
    item=bounds.word_bound(WORD_PULSES,23)
    rounded,gain=bounds.compose([item]*STAGES)
    ideal=6*STAGES*spectral.GRAM_GATE
    complete=rounded+ideal
    if complete>=F(1,12000) or gain>=F(1001,1000):
        raise ArithmeticError('the complete source energy instrument exceeds its envelope')
    return {'stages':STAGES,'native_words':STAGES,'native_pulses':STAGES*WORD_PULSES,
            'maximum_literal_nodes_per_word':WORD_PULSES,
            'total_streamed_word_nodes':STAGES*WORD_PULSES,
            'word_admission_ceiling':syn.interpreter.MAX_WORD,
            'stage_counter_bits':7,'pulse_counter_bits':16,
            'terminal_stage_flags':STAGES,'final_quantum_residual_labels':1,
            'complete_halfdiamond_strict_upper':'1/12000',
            'complete_gain_strict_upper':'1001/1000',
            'ideal_ray_comparison_halfdiamond_upper':str(ideal),
            'single_monolithic_repair_AST_admitted':False,
            'single_monolithic_coherent_word_admitted':False,
            'sequential_declared_source_services':True,
            'all_first_processor_failures_and_padding_flags_retained':True,
            'nondestructive_Luders_postmeasurement_state_claimed':False}


def resolution(margin):
    if type(margin) is not F or margin<=0:
        raise ValueError('a positive exact counted-energy margin required')
    epsilon=60*(spectral.RESIDUAL_GATE+2*spectral.GRAM_GATE)/spectral.FLOOR
    return {'counted_energy_margin':str(margin),
            'nearby_spectral_bad_label_probability_upper':str(min(F(1),(epsilon/margin)**2)),
            'ideal_compiled_bad_label_probability_upper':str(min(F(1),(epsilon/margin)**2+6*STAGES*spectral.GRAM_GATE)),
            'stored_unnormalized_bad_label_score_upper':str((epsilon/margin)**2+F(2,12000)),
            'input_domain':'any normalized eigenstate of the actual fair transfer generator',
            'failure_and_final_residual_flags_count_as_bad_labels':True}


def ideal_probabilities(columns=None):
    """Diagnostic sequential instrument, with all terminal reference blocks."""
    columns=np.eye(STAGES,dtype=complex) if columns is None else menus.source_columns(columns)
    current=columns.copy();outputs=[]
    for index in range(STAGES):
        ray=np.array([float(x) for x,_ in vector(index)])
        ray/=np.linalg.norm(ray)
        row=ray@current
        outputs.append(np.outer(row,row.conj()))
        current-=np.outer(ray,row)
    return tuple(outputs),current


def sampling_budget(error=F(1,100),failure=F(1,100),*,preparation_rejection=False):
    """Independent normalized ideal trials of the declared coarse histogram."""
    if (type(error) is not F or error<=0 or type(failure) is not F
            or not 0<failure<1 or type(preparation_rejection) is not bool):
        raise ValueError('positive exact error, failure probability and Boolean preparation flag required')
    bins=68+int(preparation_rejection)
    needed=F(bins,4)/(failure*error*error)
    trials=(needed.numerator+needed.denominator-1)//needed.denominator
    return {'coarse_histogram_labels':bins,'independent_ideal_repetitions':trials,
            'total_variation_error':str(error),'failure_probability_upper':str(failure),
            'all_processor_failure_labels_coarsened_only_for_this_histogram':True,
            'full_archival_word_distribution_sampling_bound':False,
            'normalized_hardware_trial_law_or_independent_fresh_repetitions_derived':False}


def prepared_reference_histogram():
    """Exact complete energy effects on the native accepted q=2 reference.

    Integer matrix prefixes keep one common denominator. No floating
    eigensolver, polar basis, clipping or outcome renormalization is used.
    """
    from model.native_noncentral_reference import reference
    columns,_=spectral.load()
    integer_columns=tuple(tuple(int(x*spectral.GRID) for x in column) for column in columns)
    masks=reference.p.masks(2)
    weights=dict(zip(masks,reference.occupation_weights(2),strict=True))
    states=syn.interpreter.local.basis()
    ordered_weights=tuple(weights[sum(x<<i for i,x in enumerate(state))] for state in states)
    partition=sum(ordered_weights)
    prefix=[[int(i==j) for j in range(STAGES)] for i in range(STAGES)]
    denominator=1;probabilities=[]
    for column in integer_columns:
        norm=sum(x*x for x in column)
        row=[sum(column[i]*prefix[i][j] for i in range(STAGES)) for j in range(STAGES)]
        numerator=sum(weight*x*x for weight,x in zip(ordered_weights,row,strict=True))
        probabilities.append(F(numerator,partition*norm*denominator*denominator))
        prefix=[[norm*prefix[i][j]-column[i]*row[j] for j in range(STAGES)] for i in range(STAGES)]
        denominator*=norm
    residual=F(sum(ordered_weights[j]*sum(prefix[i][j]**2 for i in range(STAGES))
                   for j in range(STAGES)),partition*denominator*denominator)
    if min(probabilities)<0 or residual<0 or sum(probabilities)+residual!=1:
        raise ArithmeticError('the exact sequential energy instrument lost its residual or normalization')
    acceptance=reference.filter_budget(2)['acceptance']
    complete=(1-acceptance,)+tuple(acceptance*x for x in probabilities)+(acceptance*residual,)
    if sum(complete)!=1:raise ArithmeticError('the full preparation/energy law lost rejection')
    program=syn.interpreter.program(biased=False,commands=reference.p.commands(1),
                                    capacity=1,overflow='export_oldest')
    return {'native_repair_preparation_attempts':71,'repair_program_digest':syn.interpreter.digest(program),
            'repair_resources':program['declared_resources'],'reference_weights_in_source_order':ordered_weights,
            'reference_partition':partition,'preparation_acceptance':str(acceptance),
            'complete_outcome_order':'preparation rejection;66 energy flags;final quantum residual',
            'complete_exact_ideal_probabilities':[str(x) for x in complete],
            'complete_ideal_law_sums_to_one':True,
            'maximum_coherent_pulses_on_accepted_branch':STAGES*WORD_PULSES,
            'unconditional_stored_halfdiamond_upper':str(acceptance*F(1,12000)),
            'full_preparation_records_retained_by_source':True,
            'mixing_record_values_marginalized_only_in_reported_law':True,
            'preparation_rejection_removed_from_reported_law':False,
            'conditional_source_factorization_from_original_reference':True,
            'physical_pulse_to_attempt_duration_map':False}


@lru_cache(maxsize=11)
def fair_kernel(edge):
    """All25 actual fair E slots, independent of weighted detector weights."""
    source=syn.interpreter
    command={'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'E','repeat':1}
    program=source.program(biased=False,commands=[command],capacity=1,overflow='export_oldest')
    marker={(i,i,i,i):F(1) for i in range(STAGES)}
    rows=[]
    for out in source._outcomes(program,command,marker):
        pairs=[]
        for (y,i,z,j),value in out.density.items():
            if y!=z or i!=j or value!=F(1,2):
                raise ArithmeticError('the fair source E chart lost its partial-injection weight')
            pairs.append((i,y))
        rows.append((out.record,tuple(sorted(pairs))))
    if len(rows)!=25:
        raise ArithmeticError('the complete fair E alphabet changed')
    return tuple(rows)


def preparation_branch(slots):
    """Decode the fair eleven-E preparation; every zero raw word survives."""
    detector=bounds.detector
    if type(slots) is not tuple or len(slots)!=11 or any(type(s) is not int or not 0<=s<25 for s in slots):
        raise ValueError('all eleven original raw addresses required')
    mapping={i:i for i in range(STAGES)};records=[]
    for edge,slot in zip(detector.tree(),slots,strict=True):
        record,pairs=fair_kernel(edge)[slot]; records.append(record)
        step=dict(pairs);mapping={i:step[j] for i,j in mapping.items() if j in step}
    if len(mapping)>1:raise ArithmeticError('the fair E tree did not isolate one initial configuration')
    if mapping:
        initial,terminal=detector.decode(tuple(records))
        states=syn.interpreter.local.basis()
        if mapping!={states.index(initial):states.index(terminal)}:
            raise ArithmeticError('the source record decoder and actual fair branch disagree')
    return F(1,2**11),mapping,tuple(records)


def reset_word(terminal,target=0):
    if any(type(x) is not int or not 0<=x<STAGES for x in (terminal,target)):
        raise ValueError('original source configuration indices required')
    difference=syn.tomo.codes()[terminal]^syn.tomo.codes()[target]
    return syn.tomo.route_word(tuple(('X',bit) for bit in range(9) if difference>>bit&1))


def preparation_budget():
    return {'native_fixed_E_attempts':11,'repair_AST_nodes':11,
            'original_raw_slots_per_attempt':25,'formal_raw_words':25**11,
            'conditional_reset_native_pulses_upper':18,
            'retained_repair_FIFO_slots':1,'retained_repair_archive_records':10,
            'exact_prepared_target_source_index':0,
            'maximum_full_protocol_native_pulses':STAGES*WORD_PULSES+18,
            'native_attempts_and_service_pulses_are_distinct':True,
            'physical_pulse_to_attempt_conversion_derived':False}


def failure_branch(stage,pulse,complement,columns):
    """One complete first-processor-failure query with every success prefix."""
    if type(stage) is not int or not 0<=stage<STAGES:raise ValueError('an original stage required')
    current=menus.source_columns(columns)
    for index in range(stage):current,_=menus.successful_word(word(index),current)
    return syn.processor_failure(word(stage),pulse,complement,current)


def padding_branch(stage,address,columns):
    """One internal row behind the stage's observed destructive energy flag."""
    if type(stage) is not int or not 0<=stage<STAGES:raise ValueError('an original stage required')
    if type(address) is not int or not 0<=address<512 or address in syn.tomo.codes():
        raise ValueError('an original source-export padding address required')
    current=menus.source_columns(columns)
    for index in range(stage):current,_=menus.successful_word(word(index),current)
    return syn.work_row(word(stage),address)[list(syn.tomo.codes())]@current
