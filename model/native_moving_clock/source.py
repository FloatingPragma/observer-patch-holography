"""A reversible tracking apparatus for the actual charged two-load source.

The moving code is transported by the complete coherent event. Its clock
apparatus reverses that event coherently, executes the native corrected clock,
and reapplies the transport. This is an explicitly controlled apparatus;
neither a local autonomous clock nor a stationary particle pole is inferred.
"""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np

from model.common_source import parameter_selection as selection
from model.openwave_native_coherent_clock import source as clock
from model.native_charged_propagation import source as motion


@lru_cache(None)
def geometry():
    g = motion.carrier()
    geo = selection.code_geometry()
    support = tuple(i for row in geo['supports'] for i in row)
    positions = tuple(g.index[tuple(i for i, v in enumerate(geo['states'][j]) if v)]
                      for j in support)
    if len(set(positions)) != 36:
        raise ArithmeticError('complete matching-support embedding required')
    embedding = np.zeros((g.states, 36), complex)
    embedding[list(positions), np.arange(36)] = 1
    logical = np.zeros((36, 6))
    for j in range(6):
        logical[6*j:6*j+6, j] = 1/np.sqrt(6)
    cost = selection.diagonal_cost(2)['values']
    energy = np.zeros(g.states)
    for j, state in enumerate(geo['states']):
        energy[g.index[tuple(i for i, v in enumerate(state) if v)]] = float(cost[j])
    return g, embedding, logical, positions, energy


def reflect(arc):
    return 2*arc.mean(axis=1, keepdims=True)-arc


def step(arc, character, inverse=False):
    if type(character) is not int or character not in range(6) or type(inverse) is not bool:
        raise ValueError('one of six characters and an explicit direction required')
    g = geometry()[0]
    return reflect(motion.apply_S(g, arc, character)) if inverse else motion.apply_S(g, reflect(arc), character)


def transport(arc, character, steps):
    if type(steps) is not int:
        raise ValueError('integer event count required')
    out = np.asarray(arc, complex).copy()
    for _ in range(abs(steps)):
        out = step(out, character, inverse=steps < 0)
    return out


@lru_cache(maxsize=32)
def frame(character, steps):
    g, embedding, _, _, _ = geometry()
    initial = np.repeat(embedding[:, None, :], g.events, axis=1)/np.sqrt(g.events)
    return transport(initial, character, steps)


def source_probabilities(columns, rho):
    """Trace the entire event register; retain every two-load source outcome."""
    return np.einsum('xer,rs,xes->x', columns, rho, columns.conj(), optimize=True).real


@lru_cache(maxsize=72)
def physical_effects(character, steps, spectator, axis):
    """Execute the inverse tracking word, then the complete native quadrature.

    The native source instrument is tensored with event identity. It has two
    rank-one source outcomes and the complete source complement; no event
    state is projected or normalized during inverse transport.
    """
    _, embedding, _, _, _ = geometry()
    returned = transport(frame(character, steps), character, -steps)
    result = []
    for ray in clock.pair_vectors(spectator, axis):
        amplitudes = np.einsum('x,xer->er', (embedding@ray).conj(), returned)
        result.append(amplitudes.conj().T@amplitudes)
    flat = returned.reshape(-1, 36)
    result.append(flat.conj().T@flat-result[0]-result[1])
    return result


def frame_controls():
    g, embedding, logical, positions, _ = geometry()
    rows = []
    for character in range(6):
        for n in (0, 1, 4, 8):
            columns = frame(character, n)
            flattened = columns.reshape(-1, 36)
            inverse_error = float(np.linalg.norm(transport(columns, character, -n)-frame(character, 0)))
            gram_error = float(np.max(abs(flattened.conj().T@flattened-np.eye(36))))
            logical_arc = columns@logical
            admitted = motion.matching_code_projection(g, logical_arc)
            code_prob = np.sum(abs(admitted)**2, axis=(0, 1))
            support_prob = np.sum(abs(logical_arc[list(positions)])**2, axis=(0, 1))
            if inverse_error > 2e-11 or gram_error > 2e-12:
                raise ArithmeticError('transported code lost a reference column')
            rows.append({'character': character, 'coherent_events': n,
                'complete_source_event_dimension': g.states*g.events,
                'matching_support_columns': 36, 'matrix_unit_Gram_checks': 36**2,
                'inverse_transport_error_diagnostic': inverse_error,
                'Gram_error_diagnostic': gram_error,
                'fixed_six_code_population_diagnostic': code_prob.tolist(),
                'fixed_36_support_population_diagnostic': support_prob.tolist(),
                'all_other_source_event_outcomes_retained': True})
    if not any(min(r['fixed_36_support_population_diagnostic']) < .99 for r in rows if r['coherent_events'] == 1):
        raise ArithmeticError('expected fixed-support leakage control disappeared')
    return rows


def repeated_records():
    g, _, _, _, energy = geometry()
    records = []
    # These settings are selected from the existing finite clock protocol.
    # No measured frequency, particle mass, or fit enters the experiment.
    time, steps = F(64), 4096
    errors = {law: F(clock.rc.error_bound(law,time,steps)['trace_norm_upper_exact'])
              for law in ('C','E')}
    h = clock.rc.logical_hamiltonian().diagonal()
    for character in (0, 1, 3):
        for law in ('C', 'E'):
            for spectator in range(3):
                for axis in ('X', 'Y'):
                    n1, n2 = 1, 4
                    psi = clock.pair_vectors(spectator, 'X')[0]
                    initial = np.outer(psi, psi.conj())
                    first = clock.evolve(initial, law, time, steps)
                    source_effects = clock.effects(spectator, axis)
                    moving_effects = [physical_effects(character, n, spectator, axis) for n in (n1, n2)]
                    effect_error = max(float(np.max(abs(p-q))) for effects in moving_effects for p,q in zip(effects, source_effects, strict=True))
                    if effect_error > 2e-12:
                        raise ArithmeticError('actual tracked quadrature differs from the source instrument')
                    probabilities = np.zeros((3,3)); after_first = np.zeros((36,36), complex)
                    for i, p in enumerate(source_effects):
                        branch = p@first@p
                        after_first += branch
                        second = clock.evolve(branch, law, time, steps)
                        for j, q in enumerate(moving_effects[1]):
                            probabilities[i,j] = np.trace(q@second).real
                    if abs(probabilities.sum()-1)>2e-10 or probabilities.min() < -2e-12:
                        raise ArithmeticError('complete two-time record probabilities failed')
                    second_unread = clock.evolve(after_first, law, time, steps)
                    p1 = source_probabilities(frame(character, n1), first)
                    p2 = source_probabilities(frame(character, n2), second_unread)
                    if max(abs(p1.sum()-1), abs(p2.sum()-1))>2e-10:
                        raise ArithmeticError('complete moving position record lost probability')
                    rest1 = source_probabilities(frame(character,0), first)
                    rest2 = source_probabilities(frame(character,0), second_unread)
                    source_work = [float(energy@(p1-rest1)), float(energy@(p2-rest2))]
                    error = errors[law]
                    gap = F(h[2*spectator+1]-h[2*spectator])
                    intervals = clock.interval_joint(gap,time,axis,min(F(2),2*error))
                    if not all(float(intervals[i][j][0])-2e-11 <= probabilities[i,j] <= float(intervals[i][j][1])+2e-11 for i in range(3) for j in range(3)):
                        raise ArithmeticError('inherited clock bound does not cover tracked records')
                    records.append({'character':character,'law':law,'spectator':spectator,'axis':axis,
                        'primary_manifest_preparation_and_law':character==1 and law=='C' and spectator==0,
                        'frame_indices':[n1,n2], 'time_per_block_exact':str(time), 'native_steps_per_block':steps,
                        'complete_two_time_probabilities_diagnostic':probabilities.tolist(),
                        'complete_two_time_intervals_exact':[[list(map(str,z)) for z in row]for row in intervals],
                        'source_position_probabilities_diagnostic':[p1.tolist(),p2.tolist()],
                        'native_clock_signed_gap_exact':str(gap),
                        'two_block_trace_norm_bound_exact':str(min(F(2),2*error)),
                        'tracking_effect_error_diagnostic':effect_error,
                        'reencoding_work_in_declared_H2_observable_diagnostic':source_work,
                        'transport_event_calls':n2+2*n1+2*n2,
                        'additional_primitive_implementation_error_budget_required':True})
    return records


def full_program_probe():
    """Use the executable padded program and its stored native reflection."""
    from model.native_observer_program import source as program
    from model.native_observer_program.contract import manifest,descriptor_sha256
    p=manifest();initial=program.prepare(p)
    baseline=program.observe_fixed_clock(initial,'X',p)
    rows=[]
    for n in (1,4):
        moved=initial.copy()
        for _ in range(n):moved=program.spatial_step(moved,stored=True)
        restored=moved.copy()
        for _ in range(n):restored=program.spatial_step(restored,inverse=True,stored=True)
        difference=float(np.linalg.norm(restored-initial))
        observed=program.observe_fixed_clock(restored,'X',p)
        record_error=float(np.max(abs(np.asarray(observed['joint_probabilities'])-baseline['joint_probabilities'])))
        if difference>2e-10 or record_error>2e-10:
            raise ArithmeticError('complete stored tracking word failed')
        rows.append({'frame_index':n,'full_state_shape':list(moved.shape),
            'actual_stored_forward_inverse_vector_error':difference,
            'complete_native_clock_records':observed,
            'same_primary_record_error_diagnostic':record_error,
            'moving_complete_position_records':program.terminal_records(moved),
            'tracking_forward_work_H2':program.source_energy(moved)-program.source_energy(initial),
            'tracking_inverse_work_H2':program.source_energy(restored)-program.source_energy(moved),
            'state_reset_or_event_postselection_used':False})
    return {'program_sha256':descriptor_sha256(),'executed_carrier':[512,8,64],
        'primary_law':'C','primary_spectator':0,'primary_character':1,
        'baseline':baseline,'tracked_records':rows,
        'other_repeated_record_rows':'Controls over the existing C/E and spectator family; the primary manifest rows are marked explicitly.',
        'applied_stored_native_reflection':True,'retained_full_source_memory_event_padding':True}


def run():
    controls = frame_controls()
    records = repeated_records()
    return {'schema':'oph.native_transported_clock.v1',
        'verdict':'PASS_COMPLETE_NATIVE_TRACKED_CODE_AND_CLOCK_RECORDS',
        'construction':'actual two-load charged coherent event, complete native corrected clock, reversible tracking apparatus',
        'intertwiner':'C_n=W^n A V; W C_n=C_(n+1); clock_n=Ad(W^n) clock_0 Ad(W^-n)',
        'frames':controls,'repeated_records':records,'single_program_execution':full_program_probe(),
        'exact_reference_safe_error_transfer':'Unitary conjugation preserves the inherited complete source/reference trace-norm bound. All apparatus errors add by complete-channel contraction.',
        'apparatus_error_bound':'two source-clock block bounds +14*epsilon_W +readout/preparation implementation errors for frames1 and4',
        'scope':{'fixed_matching_code_invariant_under_motion':False,
                 'inverse_tracking_is_a_source_reset':False,
                 'record_postselection_or_renormalization_used':False,
                 'local_autonomous_clock_proved':False,
                 'stationary_joint_generator_or_rest_pole_proved':False,
                 'physical_work_energy_or_GeV_scale_derived':False,
                 'literal_OpenWave_row3_complete':False,
                 'new_Nature_validation':False}}
