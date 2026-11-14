"""Join the native event clock and one-history circulation without a new law.

The parameter is recorded nonzero native movements, not physical seconds.
Exact predictions concern diagonal native preparations and the declared uniform
seam/fair-coin scheduler. Finite stopped checks retain their no-event tail.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp

from ..native_circulation import construction as circulation
from ..native_clock import clock
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/event_circulation/__init__.py', 'model/event_circulation/__main__.py',
    'model/event_circulation/construction.py',
    *circulation.SOURCE_PATHS, *clock.SOURCE_PATHS,
)))


@lru_cache(maxsize=1)
def embedded():
    """Actual marked jump chain, using all native branch addresses."""
    kernel, triangle, weights, branches = circulation.source()
    degree = tuple(sum(i == x and current != 0 for i, _, _, _, current, _ in branches)
                   for x in range(66))
    if set(degree) != {8, 10} or sum(degree) != 600:
        raise ArithmeticError('native moving degrees changed')
    palm = tuple(F(d, 600) for d in degree)
    marked = [[[0]*66 for _ in range(66)] for _ in range(3)]
    for i, j, _, _, current, score in branches:
        if current:
            marked[score+1][i][j] += 40//degree[i]
        elif i != j:
            raise ArithmeticError('a no-current operation moved the diagonal state')
    marked = tuple(tuple(tuple(row) for row in matrix) for matrix in marked)
    transition = tuple(tuple(F(sum(m[i][j] for m in marked), 40)
                             for j in range(66)) for i in range(66))
    if any(sum(row) != 1 for row in transition):
        raise ArithmeticError('jump chain lost probability')
    if clock._row_times(palm, transition) != palm:
        raise ArithmeticError('Palm law is not event stationary')
    if any(palm[i]*transition[i][j] != palm[j]*transition[j][i]
           for i in range(66) for j in range(66)):
        raise ArithmeticError('event chain lost detailed balance')
    # Check irreducibility, needed for the stated long-event limit.
    reached = {0}
    while True:
        future = reached | {j for i in reached for j, p in enumerate(transition[i]) if p}
        if future == reached:
            break
        reached = future
    if len(reached) != 66:
        raise ArithmeticError('jump chain is not irreducible')
    return degree, palm, marked, transition


def geometry():
    kernel, triangle, weights, branches = circulation.source()
    degree, palm, _, transition = embedded()
    previous = circulation.exact_geometry()
    p = sp.Matrix([[sp.Rational(v) for v in row] for row in previous['cycle_projection']])
    b = sp.zeros(12, 30)
    for e, (a, z) in enumerate(kernel.seams):
        b[a, e] = -1; b[z, e] = 1
    means = [[F(0)]*30 for _ in range(66)]
    seconds = [F(0)]*30
    for i, _, e, _, current, _ in branches:
        if current:
            means[i][e] += F(current, degree[i])
            seconds[e] += palm[i]*F(current*current, degree[i])
    for i, state in enumerate(kernel.states):
        expected = -(b.T*sp.Matrix(state))/degree[i]
        if sp.Matrix(means[i]) != expected or p*expected != sp.zeros(30, 1):
            raise ArithmeticError('event-conditioned cycle drift is nonzero')
    if set(seconds) != {F(1, 30)}:
        raise ArithmeticError('event covariance is not source-derived I30/30')
    # A single actual movement at a specified attempt prepares Palm, with
    # failed preparation retained. Waiting for the first movement is different.
    prepared = [F(0)]*66
    for i, j, _, _, current, _ in branches:
        if current:
            prepared[j] += F(1, 66*60)
    if sum(prepared) != F(5, 33) or tuple(v/F(5, 33) for v in prepared) != palm:
        raise ArithmeticError('one-attempt herald does not prepare the claimed law')
    uniform = (F(1, 66),)*66
    first_event = clock._row_times(uniform, transition)
    if first_event == palm:
        raise ArithmeticError('first-event negative control stopped distinguishing Palm')
    return {
        'source_triangle': list(triangle), 'moving_branch_counts': list(degree),
        'event_stationary_configuration_law': list(map(str, palm)),
        'adjacent_configuration_event_mass': '2/5',
        'nonadjacent_configuration_event_mass': '3/5',
        'all_66_event_rows_normalized_reversible_and_irreducible': True,
        'projected_conditional_mean_zero_for_every_configuration': True,
        'stationary_edge_second_moment_per_event': '1/30',
        'cycle_covariance_per_event': 'P/30', 'cycle_rank': 19,
        'single_attempt_Palm_preparation_success': '5/33',
        'single_attempt_preparation_failure': '28/33',
        'first_event_from_uniform_is_Palm': False,
        'TV_first_event_from_uniform_to_Palm': str(sum(abs(a-z) for a,z in zip(first_event,palm))/2),
    }


def count_laws(horizons=(1, 4, 16, 64, 256)):
    """Exact joint counter/terminal-state recursion in movement-event count."""
    horizons = tuple(circulation.horizon(n) for n in horizons)
    if not horizons or len(set(horizons)) != len(horizons):
        raise ValueError('distinct nonempty event horizons required')
    degree, _, marked, _ = embedded()
    sparse = [[(i,j,value) for i,row in enumerate(matrix) for j,value in enumerate(row) if value]
              for matrix in marked]
    rows = [list(degree)]
    denominator = 600
    saved = {}
    for n in range(max(horizons)+1):
        if n in horizons:
            saved[n] = (tuple(tuple(row) for row in rows), denominator)
        if n == max(horizons):
            break
        future = [[0]*66 for _ in range(len(rows)+2)]
        for g, transitions in enumerate(sparse):
            for s, row in enumerate(rows):
                target = future[s+g]
                for i, j, value in transitions:
                    target[j] += row[i]*value
        rows = future
        denominator *= 40
        if sum(map(sum, rows)) != denominator:
            raise ArithmeticError('event counter lost complete history mass')
    return saved


def holding_control():
    """State-dependent identity holding cancels exactly from the marked jumps."""
    degree, palm, marked, transition = embedded()
    # Nonuniform holding is a stronger control than a global time rescaling.
    activity = tuple(F(1,2) if i % 2 else F(3,4) for i in range(66))
    _, _, _, branches = circulation.source()
    counts = [[[0]*66 for _ in range(66)] for _ in range(3)]
    for i,j,_,_,current,score in branches:
        if current:
            counts[score+1][i][j] += 1
    for i in range(66):
        escape = activity[i]*F(degree[i],60)
        for score in (-1,0,1):
            for j in range(66):
                count = counts[score+1][i][j]
                if activity[i]*F(count,60)/escape != F(marked[score+1][i][j],40):
                    raise ArithmeticError('identity holding changed an embedded marked row')
    return {
        'positive_state_dependent_activity': list(map(str,activity)),
        'all_marked_event_rows_identical': True,
        'general_identity': '(a(x) Q1_g(x,y))/(a(x) h(x)) = Q1_g(x,y)/h(x)',
        'stationary_law_in_attempt_parameter_changes': True,
        'all_event_counter_laws_for_same_initial_state_preserved': True,
        'waiting_costs_and_attempt_records_preserved': False,
        'work_evolving_during_holds_covered': False,
        'changed_relative_seam_weights_covered': False,
    }


def attempt_instrument():
    """Keep no-move separate from a moving event with triangle score zero."""
    _, _, _, branches = circulation.source()
    groups = {}
    for i,j,e,coin,current,g in branches:
        groups.setdefault((e,coin,current), np.zeros((66,66),complex))[j,i] = 1/math.sqrt(60)
    def outcome(key):
        e, _, current = key
        return 0 if current == 0 else 2+circulation.source()[2][e]*current
    return runtime.Instrument(
        kraus=tuple(tuple(groups[key] for key in sorted(groups) if outcome(key)==k)
                    for k in range(4)), visible=True,
        resources=(('native_attempts',1),('movement_and_circulation_reads',1)),
        binding={'outcomes':['no movement','movement score -1','movement score 0','movement score +1'],
                 'native_candidate':'C','no_movement_and_zero_circulation_are_distinct':True})


def runtime_check():
    """Actual herald plus three attempt prefixes; preparation failures retained."""
    degree, palm, _, _ = embedded()
    operation = attempt_instrument()
    herald, _ = clock.transfer.instruments(circulation.source()[0])
    zero = np.zeros((66,66))
    inactive = runtime.Instrument(kraus=((np.eye(66),),(zero,),(zero,),(zero,)),visible=True,
        binding={'outcome_zero':'inactive slot after failed preparation',
                 'other_outcomes':'retained zero-probability padding'})
    program = runtime.Program('event_circulation_three_attempt_prefixes',66,
        np.eye(66)/66,
        (runtime.Event('Palm_movement_herald',instruments=(herald,)),)+
        tuple(runtime.Event(f'attempt_{i+1}',instruments=(inactive,operation),control=0)
              for i in range(3)),
        binding={'initial_law':'uniform stationary diagonal native state',
                 'Palm_preparation':'actual first native attempt conditioned on movement',
                 'parameter':'native attempts until first nonzero movement'})
    histories = [runtime.execute(program,stop_after=1)]
    for n in range(1,4):
        histories.append(runtime.execute(program,stop_after=n+1,checkpoint=histories[-1]))
    _,_,_,branches = circulation.source()
    matrices = [[[0]*66 for _ in range(66)] for _ in range(4)]
    for i,j,_,_,current,g in branches:
        matrices[0 if current==0 else g+2][i][j] += 1
    residual = 0.
    for branch in histories[-1].branches:
        flag,*word = branch.record
        if flag:
            counts = list(degree)
            for k in word:
                counts = [sum(counts[i]*matrices[k][i][j] for i in range(66)) for j in range(66)]
            expected = np.diag([float(F(value,3960*60**3)) for value in counts])
            if dict(branch.resources).get('native_attempts')!=4:
                raise ArithmeticError('successful preparation attempt cost changed')
        else:
            expected = (np.diag([float(F(60-d,3960)) for d in degree])
                        if word==[0,0,0] else np.zeros((66,66)))
            if dict(branch.resources).get('native_attempts')!=1:
                raise ArithmeticError('failed preparation consumed a native continuation')
        residual = max(residual,float(np.linalg.norm(branch.density-expected,'fro')))
    retained = []
    stopped_error = 0.
    for n in range(1,4):
        for k in range(1,4):
            word = (0,)*(n-1)+(k,)
            counts = list(degree)
            for outcome in word:
                counts = [sum(counts[i]*matrices[outcome][i][j] for i in range(66)) for j in range(66)]
            population = tuple(F(v,600*60**n) for v in counts)
            branch = next(b for b in histories[n].branches if b.record==(1,*word))
            stopped_error = max(stopped_error,float(np.linalg.norm(
                branch.density-np.diag([float(F(5,33)*p) for p in population]),'fro')))
            retained.append({'attempt_of_first_event':n,'triangle_score':k-2,
                             'conditional_probability':str(sum(population)),
                             'unconditional_probability':str(F(5,33)*sum(population))})
    tail = sum(p*(1-F(d,60))**3 for p,d in zip(palm,degree))
    no_event = next(b for b in histories[3].branches if b.record==(1,0,0,0))
    tail_state = np.diag([float(F(5,33)*p*(1-F(d,60))**3) for p,d in zip(palm,degree)])
    stopped_error = max(stopped_error,float(np.linalg.norm(no_event.density-tail_state,'fro')))
    if sum(F(row['unconditional_probability']) for row in retained)+F(5,33)*tail+F(28,33) != 1:
        raise ArithmeticError('stopped experiment lost timeout mass')
    prefix = histories[1]
    restored = runtime.checkpoint_from_json(program,runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program,checkpoint=restored)
    return {'executed_complete_herald_and_three_attempt_records':len(histories[3].branches),
            'preparation_failure_probability':'28/33','preparation_success_probability':'5/33',
            'inactive_slots_after_failed_preparation_have_explicit_zero_padding_records':True,
            'full_terminal_density_dimension':66,'finite_cap_attempts':3,
            'first_event_records':retained,'conditional_no_event_timeout_probability':str(tail),
            'unconditional_no_event_timeout_probability':str(F(5,33)*tail),
            'complete_three_attempt_output_residual_diagnostic':circulation.bounded(residual,'event runtime'),
            'complete_stopped_state_residual_diagnostic':circulation.bounded(stopped_error,'event prefixes'),
            'checkpoint_complete_distance_diagnostic':circulation.bounded(
                runtime.complete_distance(continued,histories[-1]),'checkpoint'),
            'physical_halt_device_implemented':False,
            'stopped_law_from_executed_prefixes':True,
            'roundoff_interval_certificate':False}


def run():
    facts = geometry()
    summaries = []
    for n,(rows,den) in sorted(count_laws().items()):
        probabilities = [F(sum(row),den) for row in rows]
        moments = [sum(F(s-n)**k*p for s,p in enumerate(probabilities)) for k in range(1,5)]
        if moments[0]!=0 or moments[1]!=F(n,10) or moments[2]!=0:
            raise ArithmeticError('event martingale moments changed')
        summaries.append({'movement_events':n,'counter_min':-n,'counter_max':n,
            'complete_counter_probabilities':list(map(str,probabilities)),
            'mean':str(moments[0]),'variance':str(moments[1]),'fourth_moment':str(moments[3]),
            'excess_kurtosis':str(moments[3]/moments[1]**2-3),
            'terminal_counter_joint_counts_sha256':hashlib.sha256(canonical(rows).encode()).hexdigest(),
            'all_counter_and_terminal_configuration_mass_retained':True})
    _,palm,marked,_ = embedded()
    second_product = sum(palm[i]*F(marked[0][i][j]+marked[2][i][j],40)
                         *F(sum(marked[0][j])+sum(marked[2][j]),40)
                         for i in range(66) for j in range(66))
    return {'schema':'oph.event_circulation.v1','all_pass':True,
        'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'exact_event_geometry':facts,'event_counter_predictions':summaries,
        'holding_invariance':holding_control(),'complete_execution':runtime_check(),
        'exact_first_event_timeout_laws':[
            {'attempt_cap':n,'Palm_start_no_event_probability':str(
                F(2,5)*F(13,15)**n+F(3,5)*F(5,6)**n),
             'all_initial_diagonal_states_upper':str(F(13,15)**n)}
            for n in (3,16,64,128)],
        'dependence_control':{'successive_signed_covariance':'0',
            'successive_squared_product':str(second_product),'independent_squared_product':'1/100',
            'increments_independent':False},
        'limit':{'parameter':'rescaled number of recorded nonzero native movements',
            'stationary_triangle_variance_after_N_events':'N/10',
            'cycle_covariance':'P/30','finite_partition_Brownian_increment_action_coefficient':'15',
            'attempt_parameter_action_coefficient':'99',
            'mean_stationary_events_per_attempt':'5/33',
            'finite_chain_ergodic_bounded_martingale_FCLT_imported':True,
            'new_Lean_proof_or_finite_N_Gaussian_error_bound':False},
        'scope':{'source_carrier_native_rule_uniform_seams_fair_coins_and_record_access_supplied':True,
            'Palm_preparation_required_for_exact_stationary_variance':True,
            'identity_holding_invariant_event_observables':True,
            'source_derived_physical_seconds_or_relativistic_field_action':False,
            'cycle_rank_is_spatial_dimension':False,'new_empirical_observation':False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path); mode.add_argument('--check',type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check),result)
        print('Event circulation receipt reproduced')
    elif args.write:
        with args.write.open('x',encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result),end='')
