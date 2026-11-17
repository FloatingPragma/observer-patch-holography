"""Replay quantitative preparation, occupation and source-history bounds."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..observer.source_scalar import compiler
from ..occupation_control import intervals as iv
from ..scalar_observations.bounds import sqrt_upper
from .construction import product_incident_clock_bound, SOURCE_PATHS as CONSTRUCTION_SOURCE_PATHS
from .thermal_history import finite_gibbs_seed, history, preparation_bracket
from .thermal_spectrum import certify_spectrum
from .thermal_forms import fixed_code_form
from .thermal_source import run as source_run

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys(CONSTRUCTION_SOURCE_PATHS+tuple(
    'model/interacting_band_extension/'+name for name in (
    'thermal_forms.py', 'thermal_spectrum.py', 'thermal_history.py',
    'thermal_source.py', 'thermal_controller.py', 'thermal_certificate.py'))+(
    'model/occupation_control/intervals.py',)))

# Replays certify these declared budgets; eigensolver-dependent residuals are
# not serialized as if their rounded last digits were mathematical invariants.
# Every ceiling/floor is checked against a fresh outward enclosure, and the
# declared upper budgets (not the smaller measured values) compose below.
CANONICAL_SEED_ERROR_CEILING = F(1, 10**9)
TARGET_GIBBS_ERROR_CEILING = F(1, 10**9)
PATH_ARITHMETIC_CEILING = F(1, 10**10)
BRACKET_ENVELOPES = {
    'hot': (F('5.365084368'), F('5.36524283'), F('0.0054348')),
    'cold': (F('0.001840676370'), F('0.001840680344'), F('0.0014685')),
}
OCCUPATION_CEILINGS = {-1: F('0.0008034'), 1: F('0.0009881')}
# The twelve-bit tree has J/H/X/Z counts 40962/81924/8190/4095.
# Their angle bounds pi, pi/2, pi/4, pi give 88066.5*pi < 1937463/7.
# This template bound covers backend-dependent signed atan2/Walsh choices;
# the actual emitted binary angles are independently summed exactly below.
PREPARATION_PROCESSOR_NORM_CEILING = F(277000)


def _exact_budget(value, bound):
    if (type(value) not in (int, F) or type(bound) not in (int, F)
            or value < 0 or bound <= 0):
        raise ValueError('nonnegative exact value and positive exact budget required')
    return F(value), F(bound)


def certified_upper(value, ceiling, label):
    """Publish a fixed ceiling only after the fresh enclosure proves it."""
    value, ceiling = _exact_budget(value, ceiling)
    if value > ceiling:
        raise ArithmeticError(f'{label} exceeds its declared certificate ceiling')
    return str(ceiling)


def certified_lower(value, floor, label):
    """Publish a fixed lower endpoint only after the fresh enclosure proves it."""
    value, floor = _exact_budget(value, floor)
    if value < floor:
        raise ArithmeticError(f'{label} falls below its declared certificate floor')
    return str(floor)


def _certify_preparation_resources(executed, target):
    """Certify the actual amplitude-tree word's cost with exact angle sums.

    Recompilation uses the very target returned by the executed preparation.
    Check the whole resource record before replacing its platform-dependent
    floating sum with an exactly verified, declared upper budget.
    """
    word = compiler.compile_real_state(target.T.reshape(-1))
    resources = executed['preparation']['source_resources']
    compare_receipt(resources, compiler.resources(word))
    exact_cost = sum((abs(F(event['angle'])) for event in word), F())
    resources['integrated_processor_norm_upper'] = certified_upper(
        exact_cost, PREPARATION_PROCESSOR_NORM_CEILING, 'preparation processor norm')


def upper(value, digits=9):
    """Portable outward display; all acceptance comparisons precede rounding."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('a finite real bound is required')
    value = F(value)
    return str(F(math.ceil(value*10**digits),10**digits))


def compact(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: compact(item) for key,item in value.items()}
    if type(value) in (tuple,list):
        return [compact(item) for item in value]
    return value


def bracket_summary(row, protocol):
    lower, ceiling, preparation = BRACKET_ENVELOPES[protocol]
    return {key: row[key] for key in ('occupation','lower_form_shell','inverse_temperature','coupling')} | {
        'infinite_partition_lower': certified_lower(row['infinite_partition_lower'],lower,protocol+' partition'),
        'infinite_partition_upper': certified_upper(row['infinite_partition_upper'],ceiling,protocol+' partition'),
        'canonical_preparation_trace_distance_upper': certified_upper(
            row['canonical_preparation_trace_distance_upper'],preparation,protocol+' preparation'),
        'harmonic_complement_upper': upper(row['harmonic_complement_upper'],15)}


def run():
    # The hotter original benchmark is retained separately. Its preparation
    # certificate does not claim a complete source execution at d=24.
    hot = preparation_bracket(24,48)
    old_lower = certify_spectrum(fixed_code_form(4),bits=48)['partition_function_interval'][0]
    old_allowance = sqrt_upper(1-old_lower/hot['infinite_partition_upper'])
    cold = preparation_bracket(8,16,F(4))
    hot_summary, cold_summary = bracket_summary(hot,'hot'), bracket_summary(cold,'cold')
    seed, seed_budget = finite_gibbs_seed()
    seed_error = F(certified_upper(seed_budget['canonical_seed_vector_error_upper'],
                                  CANONICAL_SEED_ERROR_CEILING,'canonical seed'))
    rows=[]
    ideal_max=F()
    for sign in (-1,1):
        outputs, budget = history(seed,sign)
        # Both exact complete-prefix isometries act on the normalized finite
        # canonical seed. Their vector difference is at most two seed errors,
        # two finite path defects and the escaped-row defect. Every stopped
        # complete prefix uses a subset of these nonnegative contributions.
        finite_residual = F(certified_upper(F(budget['finite_residual']),PATH_ARITHMETIC_CEILING,'finite path'))
        endpoint_jump = F(certified_upper(F(budget['endpoint_jump']),PATH_ARITHMETIC_CEILING,'path rounding'))
        occupation = F(certified_upper(F(budget['occupation_residual']),OCCUPATION_CEILINGS[sign],'escaped rows'))
        vector = 2*seed_error+2*(finite_residual+endpoint_jump)+occupation
        bound = min(F(1), F(cold_summary['canonical_preparation_trace_distance_upper'])+vector)
        if not 0 < bound < F(1,400):
            raise ArithmeticError('the complete interacting occupation bound is not below 1/400')
        ideal_max=max(ideal_max,bound)
        rows.append({'native_sign':sign,
            'record_alphabet':[[0,0],[0,1],[1,0],[1,1]],
            'probabilities_numerical':[float(np.linalg.norm(value)**2) for value in outputs],
            'complete_infinite_output_trace_distance_upper':upper(bound),
            'every_completed_zero_one_or_two_read_prefix_trace_distance_upper':upper(bound),
            'candidate_vector_error_upper':certified_upper(F(budget['vector_error_upper']),occupation,'candidate vector'),
            'finite_path_residual_upper':str(finite_residual),
            'occupation_residual_upper':str(occupation),
            'path_rounding_upper':str(endpoint_jump),
            'retained_work_dimension':64,'retained_reference_dimension':64})
    executed, source_target=source_run(include_target=True)
    _certify_preparation_resources(executed, source_target)
    split=executed['target_scaled_splitting_error_upper']
    target_seed=F(iv.frobenius_upper(*iv.complex_sub(iv.complex_exact(source_target),iv.complex_exact(seed))))
    target_vector=target_seed+seed_budget['canonical_seed_vector_error_upper']
    target_error=F(certified_upper(
        target_vector*(1+F(iv.frobenius_upper(*iv.complex_exact(source_target))))/2,
        TARGET_GIBBS_ERROR_CEILING,'source target to finite Gibbs'))
    source_finite=max([split+executed['source_numerical_guard']]+[
        row['finite_source_budget_under_numerical_guards']['target_scaled_splitting_error_upper']
        +row['finite_source_budget_under_numerical_guards']['source_numerical_guard']
        for row in executed['fixed_controller_histories']])
    guarded=ideal_max+source_finite+target_error
    if guarded >= F(1,200):
        raise ArithmeticError('the source/occupation bound under numerical guards exceeds 1/200')
    # This is a resource/accuracy schedule, not execution of half a million
    # services or a measurement of physical service fidelity.
    clock=product_incident_clock_bound(8,occupation=8)
    clock_error=clock['nominal_work_law_with_retained_clock_tape_error_upper']
    if ideal_max+clock_error >= F(3,1000):
        raise ArithmeticError('the stated clock/occupation schedule exceeds 3/1000')
    clock_keys=('level','incident_rate_per_normalized_attempt_unit',
        'incident_ticks_per_read_interval','maximum_interacting_ticks',
        'maximum_guarded_local_word_calls','complete_preparation_half_trace_error_required',
        'per_word_complete_instrument_half_diamond_error_required',
        'same_H_product_formula_error_upper','timeout_probability_upper',
        'model_time_and_native_arrival_time_discrepancy_upper')
    return compact({'schema':'oph.interacting_band_extension.thermal_certificate.v1',
        'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
        'original_hot_protocol':hot_summary|{
            'original_d4_preparation_allowance_upper':certified_upper(old_allowance,F('0.566265'),'old preparation'),
            'd24_source_history_executed':False},
        'cold_protocol':cold_summary|{
            'temperature_selected_from_empirical_fit':False,
            'canonical_seed_vector_error_upper':str(seed_error),
            'histories':rows},
        'source_execution':executed,
        'source_numerical_target_to_exact_finite_Gibbs_error_upper':upper(target_error,12),
        'complete_source_to_infinite_output_error_upper_under_numerical_guards':upper(guarded),
        'complete_stopped_controller_to_infinite_prefix_output_error_upper_under_numerical_guards':upper(guarded),
        'incident_clock_schedule':{key:clock[key] for key in clock_keys}|{
            'modular_increment_per_incident_at_beta_four':-F(1,4*clock['incident_rate_per_normalized_attempt_unit']),
            'complete_nominal_work_output_with_same_clock_tape_error_upper_under_required_services':upper(ideal_max+clock_error),
            'full_refined_clock_schedule_executed':False,
            'physical_service_accuracy_certified':False},
        'composition':{
            'same_frequency_two_code_and_interacting_H':True,
            'reference_retained_in_common_Hermite_basis':True,
            'first_read_disturbance_and_feedback_retained':True,
            'all_outcomes_retained_without_renormalization':True,
            'common_native_continuation_CPTP_contraction':True,
            'independent_classical_native_tapes_average_complete_prefix_bounds':True,
            'unfinished_mid_service_states_covered_by_prefix_bound':False,
            'source_clock_word_and_occupation_limits_interchangeable_without_budget':False},
        'boundary':{
            'finite_occupation_preparation_and_history_accuracy_certified':True,
            'spatial_two_mode_projection_removed':False,
            'coherent_services_and_Gibbs_parameters_derived_from_committed_repair':False,
            'all_alternative_finite_theorem_premises_conjoined':False,
            'physical_seconds_or_universal_empirical_agreement_established':False}})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path)
    mode.add_argument('--check',type=Path)
    args=parser.parse_args()
    result=run()
    if args.check:
        compare_receipt(strict_json.load(args.check),result)
        print('Interacting thermal occupation and source-history certificate reproduced')
    elif args.write:
        with args.write.open('xb') as stream:
            stream.write(canonical(result).encode())
        print(f'Wrote {args.write}')
    else:
        print(canonical(result),end='')


if __name__=='__main__':
    main()
