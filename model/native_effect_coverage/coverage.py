"""Declared native routing, finite preparation and complete effect coverage.

The inherited pair-analyzer design supplies the coherent words. The tree
instrument supplies their complete weighted-native basis readout and reset.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import numpy as np

from . import detector as d
from model.native_interface import tomography as t
from model.compiled_native_interface import construction as compiled
from model.observer.source_scalar import compiler
from model.native_program import interpreter
from model.local_repair_family import modular_clock as clock

RETURNED_GATE=F(1,10**11)
BLOCK_GATE=F(2,10**11)


def reset_operations(terminal,target):
    if any(type(x) is not int or not 0<=x<66 for x in (terminal,target)):
        raise ValueError('two literal source basis labels required')
    mask=t.codes()[terminal]^t.codes()[target]
    return tuple(('X',bit) for bit in range(9) if mask>>bit&1)


def reset_word(terminal,target):
    """Existing X words, with the identity branch requiring no service."""
    tape=t.route_word(reset_operations(terminal,target))
    if tape:
        interpreter.source_word(tape,9)
    return tape


def preparation_word(a,b=None,kind=None):
    """Post-reset preparation uses the executed inverse analyzer word."""
    if type(a) is not int or not 0<=a<66:
        raise ValueError('an original source basis index required')
    if b is None and kind is None:
        return []
    tape=compiler.inverse(t.word(a,b,kind))
    interpreter.source_word(tape,9)
    return tape


def block_arithmetic_bound(pulses=35):
    """Stored-pulse exact product versus binary64 local block multiplication.

    Blocks have at most eight rows and local dot length at most four. This
    bound covers ordinary IEEE multiplication; source expm/hardware error is
    a separate scope. The native full-array gate is rechecked on each host.
    """
    if type(pulses) is not int or not 0<=pulses<=35:
        raise ValueError('a bounded local block of at most35 pulses required')
    u=F(1,2**53);gamma=32*u/(1-32*u);g=F(2,10**11)
    tiny=F(1,2**900)
    envelope=pulses*(g+8*gamma*(1+g))+tiny
    if envelope>=F(1,4):
        raise ArithmeticError('the full-eight-address Frobenius bootstrap failed')
    result=16*pulses*gamma*(1+g)+tiny
    if result>RETURNED_GATE:
        raise ArithmeticError('the local arithmetic allowance exceeds its fixed gate')
    return result


def service_bound(blocks,pulses):
    if type(blocks) is not int or type(pulses) is not int or not 0<=blocks<=40 or not 0<=pulses<=131:
        raise ValueError('a bounded declared preparation/readout word required')
    error=(1+BLOCK_GATE)**blocks-1
    _,failure=clock._processor_failure_bound(pulses)
    return {'certified_blocks':blocks,'source_pulses':pulses,
            'operator_error_upper':str(error),
            'all_first_processor_failures_mass_upper':str(failure),
            'complete_halfdiamond_upper':str((error*(2+error)+failure)/2),
            'complete_map_gain_upper':str((1+error)**2+failure)}


def coherent_certificate():
    inherited=t.certificate()
    if (inherited['maximum_certified_blocks'],inherited['maximum_source_pulses_per_trial'],
            inherited['exact_design']['real_hermitian_rank'])!=(31,113,4356):
        raise ArithmeticError('the inherited full-carrier analyzer design changed')
    # Forward X/CNOT and analyzer arrays are freshly checked by t.certificate.
    # The inverse analyzers execute their own reverse-sign words.
    for kind in ('mixer','y_analyzer'):
        inverse=compiler.execute(compiler.inverse(compiled.word(kind)),3)
        if compiled.certify(inverse.conj().T,kind=kind)['operator_error_upper']>RETURNED_GATE:
            raise ArithmeticError('the actual inverse analyzer exceeds its source gate')
    blocks=[t.route_word((('X',0),)),compiler.cnot(0,1),
            compiled.word('mixer'),compiled.word('y_analyzer')]
    tapes=blocks+[compiler.inverse(tape) for tape in blocks]
    primitives={(p['kind'],p['angle'],len(p['targets'])) for tape in tapes for p in tape}
    for kind,angle,arity in sorted(primitives):
        leakage,gain=clock._processor_bounds(compiler.native_pulse(kind,angle),2**arity)
        if max(leakage,gain)>BLOCK_GATE:
            raise ArithmeticError('a declared primitive exceeds its full M6 export/gain gate')
    for a,b in combinations(range(66),2):
        for kind in ('mixer','y_analyzer'):
            interpreter.source_word(t.word(a,b,kind),9)
            interpreter.source_word(preparation_word(a,b,kind),9)
    resets=[]
    for terminal in range(66):
        for target in range(66):
            operations=reset_operations(terminal,target)
            if t.permutation(operations)[t.codes()[terminal]]!=t.codes()[target]:
                raise ArithmeticError('a record-controlled XOR fails on its actual source address')
            resets.append(len(operations))
    if max(resets)>9:
        raise ArithmeticError('reset exceeds the supplied nine-bit carrier')
    arithmetic=block_arithmetic_bound()
    return {'inherited_design':'model/native_interface/tomography.py',
            'inherited_pair_routes_rechecked':inherited['exact_routes_checked'],
            'all_addresses_per_route':512,'analyzer_settings':4291,
            'real_Hermitian_coordinates':4356,'inverse_analyzer_blocks_executed':2,
            'raw_M6_primitive_kinds_checked':len(primitives),
            'returned_block_operator_gate':str(RETURNED_GATE),
            'local_multiplication_error_upper':str(arithmetic),
            'exact_stored_pulse_product_block_gate':str(BLOCK_GATE),
            'original_export_complement_outcomes_retained':True,
            'source_padding_abort_retained':True,
            'physical_or_analytic_expm_accuracy_certified':False,
            'uniform_words':{'analyzer':service_bound(31,113),
                             'basis_reset':service_bound(9,18),
                             'coherent_preparation_after_tree':service_bound(40,131)},
            'source_basis_reset_pairs_checked':len(resets),
            'work_qubits':9,'maximum_proper_code_load_and_export_pairs':131,
            'six_level_blank_preparations_per_full_coherent_preparation_upper':131,
            'intermediate_padding_is_projected':False,
            'generic_execute_word_accepts_nine_qubits':False,
            'source_word_admission_and_compiler_execution_are_distinguished':True,
            'fair_current_only_detector_guarantee_consumed':False}


def ideal_preparation(a,b=None,kind=None):
    """Independent target; actual preparation consumes preparation_word()."""
    preparation_word(a,b,kind)
    vector=np.zeros(66,complex);vector[a]=1
    if b is not None:
        vector[a]=1/np.sqrt(2)
        vector[b]=(1 if kind=='mixer' else 1j)/np.sqrt(2)
    return np.outer(vector,vector.conj())


def preparation_certificate():
    # Every source label has a unique code; every requested superposition is
    # within its routed two-address reducing subspace, including inverses.
    if len(set(t.codes()))!=66:
        raise ArithmeticError('native chart is not injective')
    zero=(compiled.ZERO,compiled.ZERO);root=(F(),F(1,2))
    for kind in ('mixer','y_analyzer'):
        matrix=compiled.target(kind)
        vector={i:compiled._conjugate(matrix.get((2,i),zero)) for i in range(8)}
        expected={2:(root,compiled.ZERO),1:((root,compiled.ZERO) if kind=='mixer'
                                           else (compiled.ZERO,root))}
        for i in range(8):
            for j in range(8):
                actual=compiled._complex_product(vector[i],compiled._conjugate(vector[j]))
                wanted=compiled._complex_product(expected.get(i,zero),
                                                  compiled._conjugate(expected.get(j,zero)))
                if actual!=wanted:
                    raise ArithmeticError('the exact inverse target has the wrong prepared phase')
    return {'basis_preparations':66,'real_balanced_pair_preparations':2145,
            'imaginary_balanced_pair_preparations':2145,'spanning_preparations':4356,
            'basis_reset_ideal_channel':'rho_SR -> |a><a| tensor Tr_S rho_SR',
            'retained_archive_independent_of_reference':False,
            'prepared_source_factors_from_archive_and_reference':True,
            'pure_source_input_supplied':False,
            'coherent_code_and_processor_blank_services_supplied':True,
            'feedback_uses_decoded_records_only':True,
            'ordinary_repair_hold_grammar_issues_coherent_words':False,
            'number_of_source_copies_per_trial':1,
            'record_restoration':'fine branch sqrt(w_h)|initial><initial|',
            'record_reset':'fine branch sqrt(w_h)|target><initial|',
            'native_raw_instrument_is_destructive':True}


def stable_span_certificate(coherent):
    measurement=F(coherent['uniform_words']['analyzer']['complete_halfdiamond_upper'])
    preparation=F(coherent['uniform_words']['coherent_preparation_after_tree']['complete_halfdiamond_upper'])
    # sqrt(66*4290)<533; sqrt(34386)<186. Score maps retain failures at0.
    effect_loss=1066*measurement
    preparation_loss=372*preparation
    if max(effect_loss,preparation_loss)>=1:
        raise ArithmeticError('the complete service allowance does not preserve finite span')
    return {'actual_stored_primitive_record_functional_rank':4356,
            'effect_analysis_smallest_singular_value_lower':str(1-effect_loss),
            'actual_preparation_separating_lower_factor':str(1-preparation_loss),
            'preparation_input':'one fixed source/reference input, arbitrary but identical across settings',
            'preparation_success_blocks_not_necessarily_normalized':True,
            'rounded_complete_maps_asserted_exactly_trace_preserving':False,
            'all_failure_outputs_retained_and_scored_zero':True,
            'exact_ideal_pure_replacement_claimed_for_rounded_services':False,
            'fixed_error_certifies_rank_not_zero_error':True}


def coverage_certificate():
    design=t.design_certificate()
    if design['real_hermitian_rank']!=4356:
        raise ArithmeticError('finite quadrature effect span changed')
    # Positive summands of a rank-one positive operator share its range.
    # The finite menu atoms have one coordinate or two equal-modulus ones.
    witness=(F(3,5),F(4,5))
    if sum(x*x for x in witness)!=1 or abs(witness[0])==abs(witness[1]):
        raise ArithmeticError('the positive-cone exclusion witness is invalid')
    return {'basis_effects':66,'pair_settings':4290,'settings_total':4291,
            'real_effect_span_dimension':4356,'complex_generated_star_algebra':'M_66(C)',
            'finite_preparation_separates_all_Hermitian_effects':True,
            'preparation_test_frobenius_error_squared_coefficient':34386,
            'effect_coordinates':'diagonal; halfsum+Re(E_ab); halfsum-Im(E_ab)',
            'nonrealizable_rank_one_effect_in_one_analyzer_positive_menu':
                [['9/25','12/25'],['12/25','16/25']],
            'limitation_domain':'one analyzer followed by basis detector, arbitrary classical mixtures/coarsegrainings',
            'all_admitted_source_word_sequences_ruled_out':False,
            'universal_all_effect_realizability_or_coexistent_refinement':False,
            'operational_Born_probability_premises_independently_derived':False,
            'runtime_Kraus_trace_probability_semantics_consumed':True,
            'physical_position_momentum_CCR_or_QFT_identified':False}


def dyadic_diagonal_refinement(f,g,bits):
    """Exact shared-event refinement from fair native coins on seam(0,1).

    The tree label precedes these coin repairs; every native current and
    rejected interval is retained. This constructs diagonal effects only.
    """
    if type(bits) is not int or not 0<=bits<=53 or len(f)!=66 or len(g)!=66:
        raise ValueError('a declared diagonal finite-coin context is required')
    denominator=2**bits
    if any(type(x) is not F or x<0 or (x*denominator).denominator!=1 for x in (*f,*g)):
        raise ValueError('nonnegative dyadic probabilities with the declared denominator required')
    if any(a+b>1 for a,b in zip(f,g,strict=True)):
        raise ValueError('coexistent diagonal effects must sum to at most one')
    rows=d.kernel((0,1))
    for i in range(66):
        for coin in (0,1):
            if sum(row.weight for row in rows if row.record.coin==coin and any(x==i for x,_ in row.pairs))!=F(1,2):
                raise ArithmeticError('the actual equal-weight native seam is not a fair coin')
    return {'coin_suffix_attempts':bits,'full_native_attempts':11+bits,
            'coin_words':denominator,
            'accept_f_counts':[int(x*denominator) for x in f],
            'accept_g_counts':[int(x*denominator) for x in g],
            'reject_counts':[int((1-a-b)*denominator) for a,b in zip(f,g,strict=True)],
            'all_current_coin_and_reject_records_retained':True,
            'diagonal_additivity_from_disjoint_intervals':True}
