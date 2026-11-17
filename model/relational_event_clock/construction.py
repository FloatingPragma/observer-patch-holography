"""Complete joint records and interventions indexed by local current events."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np
import sympy as sp

from . import kernel as k
from ..clocked_quantum_readout import construction as prepared
from ..relative_modular_dynamics import construction as modular
from ..native_interface import construction as interface
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*k.SOURCE_PATHS, *modular.SOURCE_PATHS,
    'model/relational_event_clock/__init__.py',
    'model/relational_event_clock/__main__.py',
    'model/relational_event_clock/construction.py',
    'model/local_tick/experiment.py',
    'model/local_tick/SPEC.md',
    'tools/strict_json.py',
))))


def _bound(value, name):
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 2e-11:
        raise ArithmeticError(f'{name} lost its numerical comparison')
    return value


def _sum(states):
    result = defaultdict(F)
    for state in states:
        for key, value in k.sparse(state).items():
            result[key] += value
    return k.sparse(dict(result))


def _probabilities(outputs):
    return {record: k.trace(state) for record, state in outputs.items()}


def joint_states(*, intervention=False, bias=False):
    """Complete two-event record/terminal-state law from the same Palm input."""
    if type(intervention) is not bool or type(bias) is not bool:
        raise ValueError('Boolean protocol and scheduler controls required')
    _, _, palm = k.source_palm()
    outputs = {}
    for first in range(1, len(k.marks())):
        state = k.event_apply('C', palm, mark=first, bias=bias)
        stages = tuple(enumerate(k.exact_intervention('C', state))) if intervention else ((-1, state),)
        for outcome, continued in stages:
            for second in range(1, len(k.marks())):
                outputs[first, outcome, second] = k.event_apply('C', continued, mark=second, bias=bias)
    if sum(map(k.trace, outputs.values()), F()) != 1:
        raise ArithmeticError('complete continuing local-event law lost probability')
    return outputs


def joint_report():
    reports, laws = [], []
    for intervention in (False, True):
        outputs = joint_states(intervention=intervention)
        law = _probabilities(outputs)
        laws.append(law)
        same_neighbour = sum((value for (a, _, b), value in law.items()
                              if k.marks()[a][0] == k.marks()[b][0]), F())
        opposite_sign = sum((value for (a, _, b), value in law.items()
                             if k.marks()[a][1] == -k.marks()[b][1]), F())
        digest = hashlib.sha256()
        rows = []
        for record, state in outputs.items():
            density = [str(state.get((i, i), F())) for i in range(66)]
            if any(i != j for i, j in state):
                raise ArithmeticError('a diagonal native intervention produced coherence')
            digest.update(canonical([list(record), density]).encode())
            rows.append({'record': list(record), 'conditional_probability': str(law[record]),
                         'unconditional_prepared_probability': str(F(5, 198)*law[record])})
        reports.append({'event_triggered_native_intervention': intervention,
            'intervention_seam': list(k.intervention_seam()) if intervention else None,
            'complete_joint_records': rows,
            'complete_native_terminal_diagonal_sha256': digest.hexdigest(),
            'same_neighbour_probability': str(same_neighbour),
            'opposite_clock_current_probability': str(opposite_sign),
            'preparation_success': '5/198', 'preparation_failure': '193/198',
            'conditioned_reset_between_events': False})
    pooled = {(a, b): sum((p for (x, _, y), p in laws[1].items() if (x, y) == (a, b)), F())
              for a in range(1, 11) for b in range(1, 11)}
    baseline = {(a, b): p for (a, _, b), p in laws[0].items()}
    control = {(a, b): F(1, 50) if k.marks()[a][1] == -k.marks()[b][1] else F()
               for a, b in baseline}
    intervention_tv = sum((abs(p-pooled[key]) for key, p in baseline.items()), F())/2
    if intervention_tv <= 0:
        raise ArithmeticError('the designated boundary intervention no longer changes the event law')
    return {'protocols': reports,
            'event_triggered_intervention_changes_two_event_TV':
                str(intervention_tv),
            'independent_neighbours_with_required_sign_alternation_TV':
                str(sum((abs(p-control[key]) for key, p in baseline.items()), F())/2),
            'independent_neighbour_same_edge_probability': '1/5',
            'records': 'first/second marks index signed neighbouring seam; intervention outcome is -1=absent or0/1/2=current-1/0/+1'}


def controls():
    if k.quantum.mixing.seams()[0] != (0, 1):
        raise ArithmeticError('the declared remote seam-weight control changed')
    _, _, palm = k.source_palm()
    ordinary = _probabilities(joint_states())
    biased = _probabilities(joint_states(bias=True))
    tv = sum((abs(ordinary[key]-biased[key]) for key in ordinary), F())/2
    if tv <= 0:
        raise ArithmeticError('relative seam scheduler control lost its event-law difference')
    diagonal_counts = k.block_counts('C')
    q = sp.Matrix(diagonal_counts[0])/60
    r = sum((sp.Matrix(m)/60 for m in diagonal_counts[1:]), sp.zeros(66))
    jump = sp.Matrix(k.event_block('C'))
    if jump != r*(sp.eye(66)-q).inv(method='DM'):
        raise ArithmeticError('source event resolvent order changed')
    phase = tuple(1 if v == 0 else -1 if v == 1 else 0 for v in range(12))
    # Algebraic identity on the exact diagonal and both phase sectors. The
    # CP-series proof covers all sectors; no state-dependent hidden measure
    # is used to create the inserted identity hold.
    sectors_checked = 0
    for law in ('C', 'E'):
        for delta in (k.ZERO, phase, tuple(-x for x in phase)):
            raw = k.block_counts(law, delta)
            size = len(raw[0]); original_q = sp.Matrix(raw[0])/60
            original_r = sum((sp.Matrix(a)/60 for a in raw[1:]), sp.zeros(size))
            for activity in (F(1, 2), F(2, 3)):
                lazy_q = (1-sp.Rational(activity))*sp.eye(size)+sp.Rational(activity)*original_q
                lazy_r = sp.Rational(activity)*original_r
                if sp.eye(size)-lazy_q != sp.Rational(activity)*(sp.eye(size)-original_q):
                    raise ArithmeticError('identity holding failed its resolvent factorization')
                if lazy_r/sp.Rational(activity) != original_r:
                    raise ArithmeticError('identity holding changed a relative source operation')
            sectors_checked += size
    hazard = [sum(matrix[j][i] for matrix in diagonal_counts[1:] for j in range(66))
              for i in range(66)]
    zero = [i for i, value in enumerate(hazard) if value == 0]
    if len(zero) != 15:
        raise ArithmeticError('the local zero-hazard domain changed')
    stopped_columns = [sum(jump[j, i] for j in range(66)) for i in zero]
    if set(stopped_columns) != {sp.Integer(1)}:
        raise ArithmeticError('remote motion no longer permits eventual local ticks')
    # A diagonal reciprocal hazard is undefined on 15 valid initial states.
    # Even on the Palm support it drops remote state changes before the tick.
    deleted = sp.zeros(66)
    for i, rate in enumerate(hazard):
        if rate:
            for j in range(66):
                deleted[j, i] = sum(matrix[j][i] for matrix in diagonal_counts[1:])/sp.Rational(rate)
    original_first = k.event_apply('C', k.diagonal((F(1,66),)*66))
    palm_tv = sum((abs(original_first.get((i,i), F())-palm.get((i,i), F()))
                   for i in range(66)), F())/2
    if palm_tv <= 0:
        raise ArithmeticError('first-hit preparation no longer distinguishes the local Palm law')
    return {'identity_holding': {
        'checked_activity': ['1/2', '2/3'],
        'checked_diagonal_and_phase_matrix_units': sectors_checked,
        'exact_map_identity': 'R_a(I-Q_a)^(-1)=R(I-Q)^(-1), Q_a=(1-a)Id+aQ, R_a=aR',
        'complete_event_instrument_and_event_triggered_protocol_invariant': True,
        'attempt_timestamps_costs_and_finite_cap_success_invariant': False,
        'independent_work_evolving_during_holds_covered': False},
        'changed_remote_seam_weight': {'biased_seam': [0, 1],
            'relative_weight': 2, 'other_seam_relative_weights': 1,
            'same_initial_local_Palm_state': True,
            'complete_two_event_record_TV': str(tv)},
        'instantaneous_hazard_normalization_control': {
            'zero_instantaneous_hazard_states': len(zero),
            'all_these_states_tick_eventually': True,
            'one_state_index': zero[0],
            'one_state_occupied_ports': [i for i,x in enumerate(k.quantum.STATES[zero[0]]) if x],
            'remote_no_tick_maps_cannot_be_erased': True,
            'full_resolvent_differs_from_diagonal_hazard_normalization': jump != deleted},
        'first_tick_from_uniform_is_local_Palm': False,
        'first_tick_from_uniform_to_local_Palm_TV': str(palm_tv)}


def palm_report():
    failed, success, palm = k.source_palm()
    no = sp.Matrix(k.block_counts('C')[0])/60
    waits = (sp.eye(66)-no).inv(method='DM')
    expected = sum(F(waits[j, i])*palm.get((i, i), F()) for i in range(66) for j in range(66))
    if expected != F(198, 5):
        raise ArithmeticError('local Palm intertick attempt cost changed')
    tails = k.tail_certificate()
    if next(F(row['maximum_timeout_probability']) for row in tails if row['attempt_cap'] == 2) != F(1799, 1800):
        raise ArithmeticError('the stated two-attempt geometric clock certificate changed')
    return {'clock_port': k.CLOCK_PORT, 'clock_marks': [list(mark) for mark in k.marks()],
        'one_attempt_success': str(k.trace(success)), 'one_attempt_failure': str(k.trace(failed)),
        'Palm_diagonal': [str(palm.get((i, i), F())) for i in range(66)],
        'Palm_support_dimension': len(palm),
        'Palm_is_faithful_on_full_M66': False,
        'Palm_mean_attempts_per_local_tick': str(expected),
        'selected_event_chain_stationary': True,
        'local_Palm_faithful_reference_on_support': True,
        'event_clock_identified_with_reference_modular_flow': False,
        'all_input_reference_no_tick_tail_bounds': tails,
        'geometric_tail': 'after2m attempts: probability(no local tick)<=(1799/1800)^m',
        'finite_cap_instrument_comparison': 'discarding the unresolved state is forbidden; its trace bounds the missing eventual branch mass',
        'quantum_event_map': 'sum_n>=0 R Q^n, where Q includes real remote native repairs'}


def _preparation():
    source = k.source('C').instrument
    herald = runtime.Instrument(kraus=(source.kraus[0], tuple(m for group in source.kraus[1:] for m in group)),
        visible=True, resources=(('native_attempts', 1), ('local_clock_reads', 1)),
        binding={'success': 'actual local movement at specified herald attempt',
                 'clock_port': k.CLOCK_PORT, 'failed_preparation_retained': True})
    program = runtime.Program('local_clock_Palm_herald', 66, np.eye(66)/66,
        (runtime.Event('native_local_movement_herald', instruments=(herald,)),))
    history = runtime.execute(program)
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(history))
    if runtime.complete_distance(history, restored) != 0:
        raise ArithmeticError('the actual local Palm preparation checkpoint changed')
    return history


def finite_joint(intervention=False, cap=2):
    """Execute two stopped blocks on retained full66 matrices, with all tails."""
    if type(intervention) is not bool or type(cap) is not int or not 1 <= cap <= 4:
        raise ValueError('a Boolean intervention and cap from1 through4 required')
    history = _preparation()
    initial = {b.record[0]: b.density for b in history.branches}
    failed, success, _ = k.source_palm()
    operation, extra = k.source('C').instrument, k.fixed_intervention('C')
    rows, residual, mass = [], 0., F()
    digest = hashlib.sha256()
    def retain(record, actual, target, attempts):
        nonlocal residual, mass
        residual = max(residual, float(np.linalg.norm(actual-k.dense(target), 'fro')))
        value = k.trace(target); mass += value
        rows.append({'record': list(record), 'probability': str(value), 'native_attempts': attempts})
        digest.update(canonical([list(record), [[i,j,str(v)] for (i,j),v in sorted(target.items())]]).encode())
    retain((0,0,0,-1,0,0), initial[0], failed, 1)
    waiting, exact_waiting = initial[1], success
    for n in range(1, cap+1):
        actual_parts = runtime._apply(operation, waiting, 1)
        exact_parts = k.exact_step('C', exact_waiting)
        waiting, exact_waiting = actual_parts[0], exact_parts[0]
        for first in range(1, 11):
            if intervention:
                actual_middle = runtime._apply(extra, actual_parts[first], 1)
                exact_middle = k.exact_intervention('C', exact_parts[first])
                middles = tuple((j, actual_middle[j], exact_middle[j]) for j in range(3))
            else:
                middles = ((-1, actual_parts[first], exact_parts[first]),)
            for mark, current, exact_current in middles:
                for m in range(1, cap+1):
                    actual_again = runtime._apply(operation, current, 1)
                    exact_again = k.exact_step('C', exact_current)
                    current, exact_current = actual_again[0], exact_again[0]
                    for second in range(1, 11):
                        retain((1,n,first,mark,m,second), actual_again[second], exact_again[second],
                               1+n+int(intervention)+m)
                retain((1,n,first,mark,cap+1,0), current, exact_current,
                       1+n+int(intervention)+cap)
    retain((1,cap+1,0,-1,0,0), waiting, exact_waiting, 1+cap)
    if mass != 1 or len(rows) != 2+10*cap*(3 if intervention else 1)*(10*cap+1):
        raise ArithmeticError('complete local-clock experiment dropped a record or probability')
    return {'intervention': intervention, 'attempt_cap_per_event': cap,
        'complete_records': rows, 'full_native_terminal_state_dimension': 66,
        'exact_complete_terminal_state_sha256': digest.hexdigest(),
        'full_state_maximum_frobenius_residual_diagnostic': _bound(residual, 'local event execution'),
        'execution': 'original native Kraus instruments via primitive runtime on continuing unnormalized matrices',
        'preparation_checkpoint_serialized_and_replayed': True,
        'reset_between_local_clock_events': False,
        'reference_hardware_or_arithmetic_error_certified': False}


def coherent_initial(law):
    states = modular.prepared_states(law)
    if any(value[1] for state in states for value in state.values()):
        raise ArithmeticError('existing native preparation left the rational-real protocol')
    return tuple({key:value[0] for key,value in state.items()} for state in states)


def coherent_predictions():
    reports = []
    for law in ('C', 'E'):
        initial = coherent_initial(law)
        rows = []
        for flag in (1, 2):
            state = initial[flag]
            for tick in (1, 2):
                state = k.event_apply(law, state)
                outcomes = prepared._pointer_exact(prepared._balanced_exact(state), law)
                success = k.trace(state)
                match = k.trace(outcomes[0 if flag == 1 else 1])
                rows.append({'herald': flag, 'local_ticks': tick,
                             'herald_weight': str(success),
                             'unconditional_match_probability': str(match),
                             'match_given_herald': str(match/success)})
        reports.append({'native_candidate': law, 'predictions': rows,
                        'failed_preparation_probability': str(k.trace(initial[0])),
                        'preparation_and_analysis_mixer': 'existing ideal66 balanced service surrogate',
                        'event_parameter': 'number of actual port2 occupation-changing events'})
    if any(F(a['match_given_herald']) == F(b['match_given_herald'])
           for a,b in zip(reports[0]['predictions'],reports[1]['predictions'],strict=True)):
        raise ArithmeticError('the C/E local-event coherent discriminator disappeared')
    return {'phase_predictions': reports,
        'exact_source_sectors_used': '66 diagonal matrix units plus the two10-dimensional signed(0,1) coherence sectors',
        'native_C_and_E_share_all_diagonal_clock_records': True,
        'native_C_and_E_local_event_phase_probabilities_differ': True,
        'modular_flow_identification_claimed': False,
        'compiled_mixer_comparison': 'existing compiled_native_interface supplies the separate gate implementation bound; this packet executes ideal66 service targets'}


def local_modular_corner():
    """The faithful local-Palm corner and an exact noncovariance witness."""
    palm = k.source_palm()[2]
    support = {i for i,j in palm if i == j}
    census = Counter(palm.values())
    if census != {F(1,100):30, F(1,50):10, F(1,25):5, F(1,20):6}:
        raise ArithmeticError('the local-Palm modular weight census changed')
    # Completed event maps land inside this corner on every native input.
    # The no-tick map does not; its unresolved state is retained on full M66.
    no_tick_leaves_corner = False
    for law in ('C','E'):
        for mark,_,_,_,mapping in k.source(law).groups:
            for i,j in mapping:
                if mark and j not in support:
                    raise ArithmeticError('completed local tick left its Palm support')
                if not mark and i in support and j not in support:
                    no_tick_leaves_corner = True
    def address(ports):
        return k.quantum.INDEX[tuple(int(i in ports) for i in range(12))]
    initial = address((0,2)), address((2,3))
    final = address((0,4)), address((3,4))
    before = palm[initial[0],initial[0]]/palm[initial[1],initial[1]]
    after = palm[final[0],final[0]]/palm[final[1],final[1]]
    if (before,after) != (F(5,4),F(1,2)) or not no_tick_leaves_corner:
        raise ArithmeticError('the local corner modular domain/control changed')
    rows=[]
    for law in ('C','E'):
        input_unit={initial:F(1)}
        coefficient=k.event_apply(law,input_unit).get(final,F())
        one_step=k.exact_step(law,input_unit)[k.marks().index((4,-1))].get(final,F())
        if coefficient <= 0 or one_step != F(1,60):
            raise ArithmeticError('the local-event modular noncovariance witness disappeared')
        rows.append({'candidate':law,'input_matrix_unit':list(initial),
            'output_matrix_unit':list(final),'eventual_channel_coefficient':str(coefficient),
            'one_attempt_recorded_coefficient':str(one_step),
            'input_modular_weight_ratio':str(before),'output_modular_weight_ratio':str(after)})
    return {'corner_dimension':len(support),
        'faithful_corner_reference_weight_census':{str(w):count for w,count in sorted(census.items())},
        'full_native_M66_reference_faithful':False,
        'completed_event_channel_maps_entire_native_algebra_into_corner':True,
        'no_tick_map_can_leave_corner':True,
        'modular_observable_flow':'alpha_s(X)=sigma^(is) X sigma^(-is) on P M66 P',
        'parameter':'reference-state modular parameter s, distinct from local event count',
        'covariance_witnesses':rows,
        'stopped_channel_commutes_with_corner_modular_flow':False,
        'modular_flow_executed_or_physically_identified':False}


def coherent_execution(law, cap=3):
    """Complete source herald/local tick/native pointer execution, incl failures."""
    law = k.quantum.candidate(law)
    if type(cap) is not int or not 1 <= cap <= 16:
        raise ValueError('bounded coherent local-clock cap required')
    history = runtime.execute(prepared.preparation_program(law))
    initial = {b.record[0]: b.density for b in history.branches}
    exact_initial = coherent_initial(law)
    operation = k.source(law).instrument
    balanced = runtime.Instrument(kraus=((interface.mixer(),),))
    _, pointer = interface.herald_and_pointer(law)
    rows, residual, mass = [], 0., F()
    def retain(record, actual, target):
        nonlocal residual, mass
        residual = max(residual, float(np.linalg.norm(actual-k.dense(target), 'fro')))
        value = k.trace(target); mass += value
        flag,n,mark,bit = record
        rows.append({'record': list(record), 'probability': str(value),
            'native_attempts': 1+(min(n,cap) if flag else 0)+int(bit >= 0),
            'balanced_mixer_calls': int(flag != 0)+int(bit >= 0)})
    retain((0,0,0,-1), initial[0], exact_initial[0])
    for flag in (1, 2):
        waiting, exact_waiting = initial[flag], exact_initial[flag]
        for n in range(1, cap+1):
            actual_parts = runtime._apply(operation, waiting, 1)
            exact_parts = k.exact_step(law, exact_waiting)
            waiting, exact_waiting = actual_parts[0], exact_parts[0]
            for mark in range(1, 11):
                analyzed = runtime._apply(balanced, actual_parts[mark], 1)[0]
                actual = runtime._apply(pointer, analyzed, 1)
                target = prepared._pointer_exact(prepared._balanced_exact(exact_parts[mark]), law)
                for bit in (0, 1):
                    index = bit if flag == 1 else 1-bit
                    retain((flag,n,mark,bit), actual[index], target[index])
        retain((flag,cap+1,0,-1), waiting, exact_waiting)
    if mass != 1 or len(rows) != 40*cap+3:
        raise ArithmeticError('local-event phase instrument dropped a complete outcome')
    return {'candidate': law, 'attempt_cap': cap, 'complete_records': rows,
        'full_native_terminal_density_comparison': True,
        'maximum_frobenius_residual_diagnostic': _bound(residual, 'coherent local clock'),
        'direct_source_Kraus_execution_not_resolvent_output_substitution': True,
        'physical_service_or_composite_roundoff_bound_claimed': False}


def run():
    return {'schema': 'oph.relational-event-clock.v1',
        'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}),
        'local_clock': palm_report(), 'continuing_two_event_laws': joint_report(),
        'invariance_and_controls': controls(),
        'complete_finite_executions': [finite_joint(False), finite_joint(True)],
        'coherent_discriminator': coherent_predictions(),
        'local_Palm_modular_corner': local_modular_corner(),
        'coherent_finite_executions': [coherent_execution(law) for law in ('C', 'E')],
        'scope': {'new_primitive_or_field_action': False,
            'existing_native_source_and_phase_interface_reused': True,
            'clock_ticks_are_local_current_events_not_incident_attempts': True,
            'physical_seconds_or_preferred_global_clock_derived': False,
            'local_event_flow_is_modular_automorphism': False,
            'source_seam_weights_and_detector_design_declared': True,
            'new_empirical_measurement_comparison': False}}
