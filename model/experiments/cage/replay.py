"""Recomputation of the cage against its receipts.

    python -m model.experiments.cage              recompute every statement, read the stored counts, compare with results.json
    python -m model.experiments.cage --record     draw the counts from the record laws, write counts.npz and results.json

`counts.npz` holds the counts of every experiment: what the observers of the construction hold. A replay
recomputes the record laws from the law and the prepared reference, checks that the stored counts are draws
from them, and recomputes every reading from the stored counts. `results.json` holds the statements and the
readings. Numbers under `residuals` are differences that are zero up to rounding; a replay holds each one under
its bound and does not compare it with the stored value.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.scalar_channel.carrier import require
from . import construction as cn
from . import instruments as ins
from . import registers

HERE = Path(__file__).resolve().parent
COUNTS = HERE / 'counts.npz'
RESULTS = HERE / 'results.json'
SCHEMA = 'oph.scalar_channel.cage.v1'
CANDIDATE = 'A-native-caged-exciton-v1'
SEED = 20261119
FAMILIES = ('clock_prepared', 'clock_no_preparation', 'clock_cage_stopped', 'clock_sharp', 'clock_changed_law',
            'reference', 'dial', 'persistence', 'persistence_changed_law', 'persistence_no_converter')
BOUNDS = dict(                                # every residual and the bound it is held under
    pair_on_one_cell=1e-12, pair_on_20_cells=1e-12, excitation_raises_by_the_level_difference=1e-12,
    excitation_is_kept=1e-12, coherent_pulse_outside_e=1e-12, sharp_pulse_outside_the_upper_mode=1e-12,
    swing_against_the_closed_form=1e-12, circulation_level=1e-12, circulation_and_pair_in_the_reference=1e-10,
    clock_law_on_151_modes_against_four_modes=1e-10, clock_law_on_sixteen_states_against_four_modes=1e-12,
    reference_law_on_151_modes_against_four_modes=1e-10, dial_law_on_151_modes_against_three_modes=1e-10,
    dial_difference_against_the_survival_of_one_object=1e-12, persistence_law_on_two_modes_against_150_modes=1e-10,
    excitation_spread_over_the_waits=1e-10, law_keeps_the_pair=1e-12)


def simulate():
    """Everything on the simulator's side: statements about the construction and the law of every record."""
    out = dict(construction={}, object={}, residuals={}, laws={})
    car = cn.cell()
    E, W = car['C'].shape
    D = cn.law(car)
    profile = cn.ring_profile(car, cn.RING)
    Phi, h, residual = cn.modes(car, profile)
    e, chi = Phi[:, 0], cn.circulation(car, 0)
    J = cn.changed_law_term(car)
    supports = cn.small_supports(car)
    ring_seams = sorted(tuple(sorted(tuple(sorted((r[k], r[(k + 1) % 4]))) for k in range(4))) for r in cn.rings(car))
    levels = np.round(np.linalg.eigvalsh(cn.HOP ** 2 * (car['C'] @ car['C'].T).astype(float)), 9)
    out['residuals']['pair_on_one_cell'] = residual
    out['construction']['one_cell'] = dict(
        largest_difference_in_the_incidence_identity=cn.line_graph_identity(car), **cn.flat_branch(car),
        profiles_on_one_two_three_four_seams=[len(supports[k]) for k in (1, 2, 3, 4)],
        every_profile_on_four_seams_is_a_ring=sorted(tuple(sorted(s)) for s in supports[4]) == ring_seams,
        ring=list(cn.RING), seams_of_the_pair=int(np.count_nonzero(Phi[:E, 0])), wedges_of_the_pair=int(np.count_nonzero(np.abs(Phi[E:, 1]) > 1e-12)),
        ports_of_the_pair=sorted({int(b) for w, (a, b, c) in enumerate(car['wedges']) if abs(Phi[E + w, 1]) > 1e-12}),
        levels_of_K=[[float(v), int(n)] for v, n in zip(*np.unique(levels, return_counts=True))],
        coupling=float(h[0, 1]), level_of_the_cage=float(math.hypot(h[0, 0], h[0, 1])))
    federations = []
    for level in (0, 1, 2):
        fed = cn.federation(level)
        exact = cn.exact_rings(fed)
        per_cell = np.bincount([r[0] // 12 for r in exact], minlength=fed['ports'] // 12)
        kept = cn.kept_profiles(fed)
        six_ring, off = cn.leaking_ring(fed)
        federations.append(dict(
            cells=fed['ports'] // 12, ports=fed['ports'], seams=len(fed['edges']), wedges=len(fed['wedges']),
            ports_with_six_seams=int((fed['degree'] == 6).sum()), largest_difference_in_the_incidence_identity=cn.line_graph_identity(fed),
            rings_on_ports_with_five_seams=len(exact), fewest_and_most_such_rings_in_a_cell=[int(per_cell.min()), int(per_cell.max())],
            profiles_that_keep_the_level=int(sum(kept)), cells_with_six_and_with_seven_such_profiles=[kept.count(6), kept.count(7)],
            ring_through_a_port_with_six_seams=list(six_ring), largest_entry_of_its_defect=off))
        if level == 0:
            twenty = fed
    out['construction']['federations'] = federations
    # 20 cells: the pair of a ring on ports with five seams, and what a ring through a port with six seams keeps
    inner = cn.exact_rings(twenty)[0]
    out['residuals']['pair_on_20_cells'] = cn.pair(twenty, cn.ring_profile(twenty, inner))[3]
    times = np.arange(256) * 0.4
    out['construction']['weight_left_in_the_pair_on_20_cells'] = dict(
        ring_on_ports_with_five_seams=_weight_left_on_a_federation(twenty, inner, times),
        ring_through_a_port_with_six_seams=_weight_left_on_a_federation(twenty, cn.leaking_ring(twenty)[0], times),
        parameter_range=[float(times[0]), float(times[-1])])

    # the object
    raises, kept_x = ins.excitation_identities()
    outside_e, outside_upper = ins.pulses()
    ts = np.linspace(0, 6, 601)
    computed, closed = ins.swing(ts)
    out['residuals'].update(excitation_raises_by_the_level_difference=raises, excitation_is_kept=kept_x,
                            coherent_pulse_outside_e=outside_e, sharp_pulse_outside_the_upper_mode=outside_upper,
                            swing_against_the_closed_form=float(np.abs(computed - closed).max()),
                            law_keeps_the_pair=float(np.linalg.norm(D @ Phi @ Phi.T - Phi @ Phi.T @ D, 2)))
    occupation, distance = ins.reference(256)
    G, from_ideal = ins.pair_in_the_reference(Phi, occupation)
    three = np.column_stack([Phi, chi])
    S3 = three.T @ occupation @ three
    out['residuals'].update(circulation_level=float(np.abs(D @ chi + chi).max()),
                            circulation_and_pair_in_the_reference=float(np.abs(S3 - block_diag(G, 1.0)).max()))
    lower = ins.rotation()[0]
    out['object'] = dict(
        level_difference=2 * cn.ENERGY, coherent_pulse=dict(potential=-2 * cn.ENERGY, time=cn.T_COHERENT),
        sharp_pulse=dict(potential=-2 * cn.ENERGY ** 2 / cn.REST, time=cn.T_SHARP),
        work_of_the_coherent_pulse_on_the_ideal_reference=cn.ENERGY + cn.REST,
        weight_of_e_on_the_excitation=(1 + cn.REST / cn.ENERGY) / 2,
        reference_from_the_filled_lower_branch=distance, pair_from_the_ideal_lower_mode=from_ideal,
        sectors_in_the_reference=ins.sectors(G), after_the_coherent_pulse=ins.after_the_coherent_pulse(G),
        changed_law=dict(norm_of_its_commutator_with_the_pair=float(np.linalg.norm(J @ Phi @ Phi.T - Phi @ Phi.T @ J, 2)),
                         weight_of_the_lower_mode_left_in_the_pair_at_1_4_16_64=[float(x) for x in ins.weight_left(D, Phi, Phi @ lower, [1, 4, 16, 64], term=J)],
                         reference_keeps_its_level=bool(np.abs((D + J) @ chi + chi).max() < 1e-12)))

    # the laws of the records
    n = E + W + 1
    H4, S4 = block_diag(h, -cn.REST, 0.0), block_diag(S3, 0.0)
    e4, r4 = np.eye(4)[0], np.eye(4)[2]
    H, S = block_diag(D, 0.0), block_diag(occupation, 0.0)
    e_all, r_all = np.r_[e, 0.0], np.r_[chi, 0.0]
    waits = ins.schedules()
    coherent, sharp, none = ins.PREPARATIONS
    stopped = H4.copy()
    stopped[:2, :2] = 0
    laws = out['laws']
    laws['clock_prepared'] = np.array([ins.joint_law(H4, e4, r4, 3, S4, w, coherent[1], coherent[2]) for w in waits])
    laws['clock_no_preparation'] = np.array([ins.joint_law(H4, e4, r4, 3, S4, w, none[1], none[2]) for w in waits])
    laws['clock_cage_stopped'] = np.array([ins.joint_law(H4, e4, r4, 3, S4, w, coherent[1], coherent[2], wait_law=stopped) for w in waits])
    laws['clock_sharp'] = ins.joint_law(H4, e4, r4, 3, S4, waits[0], sharp[1], sharp[2])
    laws['clock_changed_law'] = ins.joint_law(block_diag(D + J, 0.0), e_all, r_all, n - 1, S, waits[0], coherent[1], coherent[2])
    out['residuals']['clock_law_on_151_modes_against_four_modes'] = float(np.abs(
        ins.joint_law(H, e_all, r_all, n - 1, S, waits[0], coherent[1], coherent[2]) - laws['clock_prepared'][0]).max())
    out['residuals']['clock_law_on_sixteen_states_against_four_modes'] = float(np.abs(
        ins.fock_joint_law(S4, waits[1][:32], coherent[1], coherent[2]) - laws['clock_prepared'][1][:32]).max())
    own_waits = np.arange(ins.REFERENCE_QUERIES) * cn.STEP
    own = ins.reference_law(H4, r4, 3, S4, own_waits)
    out['residuals']['reference_law_on_151_modes_against_four_modes'] = float(np.abs(ins.reference_law(H, r_all, n - 1, S, own_waits[::16]) - own[:, :, ::16]).max())
    laws['reference'] = np.stack([1 - own, own], axis=-1)
    dial = ins.dial_law(G)
    laws['dial'] = ins.pointers(dial)
    worst_dial, worst_survival = 0.0, 0.0
    for i, q in enumerate(ins.DIAL_WAITS):
        for j in range(0, len(ins.DIALS), 20):
            Hd = block_diag(D, float(ins.DIALS[j]))
            Hd[:-1, -1] = Hd[-1, :-1] = ins.DIAL_STRENGTH * e
            values, vectors = np.linalg.eigh(Hd)
            row = (vectors[-1] * np.exp(-1j * values * q * cn.STEP)) @ vectors.T
            for k in (0, 1):
                worst_dial = max(worst_dial, abs(float((row[:-1] @ occupation @ row[:-1].conj()).real) + k * abs(row[-1]) ** 2 - dial[i, k, j]))
            worst_survival = max(worst_survival, abs(dial[i, 1, j] - dial[i, 0, j] - abs(row[-1]) ** 2))
    out['residuals'].update(dial_law_on_151_modes_against_three_modes=worst_dial, dial_difference_against_the_survival_of_one_object=worst_survival)
    kept_law, excitation = ins.persistence_law(D, e, Phi, occupation)
    changed_law, _ = ins.persistence_law(D, e, Phi, occupation, term=J)
    without, _ = ins.persistence_law(D, e, Phi, occupation, converter=False)
    laws['persistence'], laws['persistence_changed_law'], laws['persistence_no_converter'] = ins.pointers(kept_law), ins.pointers(changed_law), ins.pointers(without)
    two = np.einsum('ij,pqj->pqi', ins.effects(cn.SIGMA), np.array(
        [[_pair_before_the_analyzer(h, G, potential, time, q) for q in ins.WAITS] for _, potential, time in ins.PREPARATIONS]))
    out['residuals'].update(persistence_law_on_two_modes_against_150_modes=float(np.abs(two - kept_law).max()),
                            excitation_spread_over_the_waits=float(np.ptp(excitation, axis=1).max()))
    out['object']['excitation_in_the_state'] = {name: float(excitation[i, 0]) for i, (name, _, _) in enumerate(ins.PREPARATIONS)}
    out['object']['largest_difference_between_the_record_laws_of_the_two_laws'] = float(np.abs(kept_law - changed_law).max())
    out['waits'] = waits
    return out


def _weight_left_on_a_federation(fed, ring, times):
    """The smallest weight that an object placed on the profile of a ring keeps in its pair, on a federation."""
    E, W = fed['C'].shape
    e, ep, _, _ = cn.pair(fed, cn.ring_profile(fed, ring))
    C = fed['C'].astype(float)
    D = sp.block_array([[cn.REST * sp.eye_array(E), cn.HOP * C], [cn.HOP * C.T, -cn.REST * sp.eye_array(W)]], format='csr')
    start = np.r_[e, np.zeros(W)].astype(complex)
    step = float(times[1] - times[0])
    require(bool(np.allclose(np.diff(times), step)), 'the times of the federation check are evenly spaced')
    psi = expm_multiply(-1j * D, start, start=float(times[0]), stop=float(times[-1]), num=len(times), endpoint=True,
                        traceA=-1j * (E - W) * cn.REST)
    return float((np.abs(psi[:, :E] @ e) ** 2 + np.abs(psi[:, E:] @ ep) ** 2).min())


def _pair_before_the_analyzer(h, G, potential, time, steps):
    """(G_ee, G_e'e', Re G_ee', Im G_ee') of the pair after a preparation and a wait, on the two modes alone."""
    from scipy.linalg import expm
    U = expm(-1j * steps * cn.STEP * h) @ expm(-1j * cn.STEP * h) @ expm(-1j * time * (h + potential * np.diag([1.0, 0.0]))) @ expm(-1j * cn.STEP * h)
    S = U @ G @ U.conj().T
    return [S[0, 0].real, S[1, 1].real, S[0, 1].real, S[0, 1].imag]


def draw_counts(laws):
    """Counts for every family from one generator with a fixed seed, in the order of FAMILIES."""
    rng = np.random.Generator(np.random.PCG64(SEED))
    shots = dict(clock=ins.SHOTS['clock'], reference=ins.SHOTS['reference'], dial=ins.SHOTS['dial'], persistence=ins.SHOTS['persistence'])
    out = {}
    for name in FAMILIES:
        law = laws[name]
        law = law.reshape(law.shape[:-2] + (4,)) if name.startswith('clock') else law
        out[name] = ins.draw(rng, law, shots[name.split('_')[0]]).reshape(laws[name].shape)
    return out


def observe(counts, sim):
    """Every reading from counts and settings, and how the stored counts sit in the record laws."""
    out, gates = {}, {}
    laws, waits = sim['laws'], sim['waits']
    shots = ins.SHOTS
    flat = lambda a: np.asarray(a).reshape(np.shape(a)[:-2] + (4,))
    fits = {}
    for name in FAMILIES:
        if name == 'persistence_no_converter':
            continue
        kind = name.split('_')[0]
        a, law = (flat(counts[name]), flat(laws[name])) if kind == 'clock' else (counts[name], laws[name])
        fits[name] = ins.counts_against_a_law(a, law, shots[kind])
    out['counts_against_the_record_laws'] = fits
    gates['the_stored_counts_are_draws_from_the_record_laws'] = all(
        f['upper_tail'] > 1e-6 and f['lower_tail'] > 1e-6 and f['largest_cell_in_standard_errors'] < 6 and f['counts_on_outcomes_of_probability_zero'] == 0
        for f in fits.values())

    # the swing of a cage
    clock = {}
    for k, schedule in enumerate(('steps_of_equal_length', 'steps_of_unequal_length')):
        phase, record = ins.read_clock(counts['clock_prepared'][k], shots['clock'])
        line = ins.read_line(phase, record)
        _, quiet = ins.read_clock(counts['clock_no_preparation'][k], shots['clock'])
        _, stopped = ins.read_clock(counts['clock_cage_stopped'][k], shots['clock'])
        clock[schedule] = dict(line=line, smallest_and_largest_step_of_the_phase=[float(np.diff(phase).min()), float(np.diff(phase).max())],
                               amplitude_without_preparation=ins.read_line(phase, quiet)['amplitude'],
                               amplitude_with_the_cage_stopped=ins.read_line(phase, stopped)['amplitude'],
                               largest_difference_between_the_phase_and_the_parameter=float(np.abs(phase - waits[k]).max()))
    phase, record = ins.read_clock(counts['clock_sharp'], shots['clock'])
    clock['sharp_excitation'] = dict(range_of_the_seam_record=float(np.ptp(record)), amplitude_of_the_best_line=ins.read_line(phase, record)['amplitude'])
    phase, record = ins.read_clock(counts['clock_changed_law'], shots['clock'])
    clock['changed_law'] = dict(best_line=ins.read_line(phase, record), smallest_and_largest_step_of_the_phase=[float(np.diff(phase).min()), float(np.diff(phase).max())])
    out['clock'] = clock
    lines = [clock[s]['line'] for s in ('steps_of_equal_length', 'steps_of_unequal_length')]
    gates['one_line_carries_the_seam_record'] = all(l['share_of_the_variance'] > 0.9999 for l in lines)
    gates['the_line_sits_at_the_level_difference_of_the_law'] = all(abs(l['frequency'] - 2 * cn.ENERGY) < 5e-4 for l in lines)
    gates['no_line_without_preparation'] = all(clock[s]['amplitude_without_preparation'] < 0.01 for s in ('steps_of_equal_length', 'steps_of_unequal_length'))
    gates['no_line_with_the_cage_stopped'] = all(clock[s]['amplitude_with_the_cage_stopped'] < 0.001 for s in ('steps_of_equal_length', 'steps_of_unequal_length'))
    gates['no_swing_of_a_sharp_excitation'] = clock['sharp_excitation']['range_of_the_seam_record'] < 0.01
    gates['no_single_line_under_the_changed_law'] = clock['changed_law']['best_line']['share_of_the_variance'] < 0.9

    # the fuel dial
    own = ins.read_reference(counts['reference'], shots['reference'])
    phases = [float(own[q]) for q in ins.DIAL_WAITS]
    dial = ins.read_dial(counts['dial'], shots['dial'], phases)
    out['dial'] = dict(reference_phase_per_step=float(own[-1] / (len(own) - 1)), reference_phase_at_the_two_waits=phases, **dial,
                       records_with_disagreeing_pointers=int(counts['dial'][..., 1:3].sum()))
    gates['the_line_and_the_dial_agree'] = all(abs(l['frequency'] - dial['level_difference']) < 5e-4 for l in lines)

    # a cage kept
    read, condition = ins.read_sectors(counts['persistence'], shots['persistence'], dial['coupling'])
    excitation = {name: [s['excitation'] for s in read[i]] for i, (name, _, _) in enumerate(ins.PREPARATIONS)}
    frequencies = counts['persistence'][..., 3] / shots['persistence']
    changed = counts['persistence_changed_law'][..., 3] / shots['persistence']
    law, law_changed = laws['persistence'][..., 3], laws['persistence_changed_law'][..., 3]
    out['persistence'] = dict(
        waits_in_steps=list(ins.WAITS), condition_number_of_the_four_effects=condition, excitation_from_counts=excitation,
        empty_and_doubly_occupied_from_counts={name: [[s['empty'], s['two']] for s in read[i]] for i, (name, _, _) in enumerate(ins.PREPARATIONS)},
        counts_of_the_law_from_the_law=float(np.abs(frequencies - law).max()),
        counts_of_the_changed_law_from_the_law=float(np.abs(changed - law).max()),
        counts_of_the_changed_law_from_the_changed_law=float(np.abs(changed - law_changed).max()),
        records_of_a_full_fuel_with_the_converter_off=int(counts['persistence_no_converter'][..., 1:].sum()),
        records_with_disagreeing_pointers=int(counts['persistence'][..., 1:3].sum() + counts['persistence_changed_law'][..., 1:3].sum()))
    gates['the_sharp_excitation_is_read_at_every_wait'] = min(excitation['sharp']) >= 0.98
    gates['the_coherent_state_holds_its_share_at_every_wait'] = 0.65 <= min(excitation['coherent']) and max(excitation['coherent']) <= 0.76
    gates['no_excitation_is_read_without_preparation'] = max(excitation['none']) <= 0.02
    gates['the_changed_law_is_told_apart_by_its_counts'] = (out['persistence']['counts_of_the_law_from_the_law'] < 1e-3
                                                            and out['persistence']['counts_of_the_changed_law_from_the_law'] > 0.1
                                                            and out['persistence']['counts_of_the_changed_law_from_the_changed_law'] < 1e-3)
    gates['no_record_with_the_converter_off'] = out['persistence']['records_of_a_full_fuel_with_the_converter_off'] == 0
    return out, gates


def compute(counts=None):
    """All statements and readings. Without counts they are drawn from the record laws with the fixed seed."""
    sim = simulate()
    drawn = counts is None
    if drawn:
        counts = draw_counts(sim['laws'])
    read, gates = observe(counts, sim)
    for name, value in sim['residuals'].items():
        gates['residual_' + name] = bool(value <= BOUNDS[name])
    require(set(sim['residuals']) == set(BOUNDS), 'every residual has a bound')
    result = dict(schema=SCHEMA, candidate=CANDIDATE,
                  protocol=dict(seed=SEED, shots=dict(ins.SHOTS), step=cn.STEP, strength_of_the_converter=cn.STRENGTH,
                                strength_of_the_dial_converter=ins.DIAL_STRENGTH, dials=[float(ins.DIALS[0]), float(ins.DIALS[-1]), len(ins.DIALS)],
                                waits_of_the_dial=list(ins.DIAL_WAITS), waits_of_the_clock=ins.QUERIES, phase_settings=ins.SETTINGS,
                                bounds_on_a_step_of_the_reference=list(ins.PHASE_STEP_BOUNDS), families=list(FAMILIES)),
                  construction=sim['construction'], object=sim['object'], registers=registers.run(), **read,
                  residuals=sim['residuals'], gates=gates)
    return _plain(result), counts


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        require(math.isfinite(float(value)), 'a result is no finite number')
        return float(value)
    return value


def compare(saved, fresh, path='results'):
    """Statements are equal; numbers agree to rounding and to the reach of a fit. Returns the list of differences."""
    if isinstance(saved, dict) and isinstance(fresh, dict):
        if set(saved) != set(fresh):
            return ['%s: keys differ: %s' % (path, sorted(set(saved) ^ set(fresh)))]
        return [d for k in saved for d in compare(saved[k], fresh[k], path + '.' + k)]
    if isinstance(saved, list) and isinstance(fresh, list):
        if len(saved) != len(fresh):
            return ['%s: lengths differ' % path]
        return [d for i, (a, b) in enumerate(zip(saved, fresh)) for d in compare(a, b, '%s[%d]' % (path, i))]
    if isinstance(saved, float) and isinstance(fresh, float):
        return [] if abs(saved - fresh) <= 1e-6 * abs(saved) + 1e-8 else ['%s: %r against %r' % (path, saved, fresh)]
    return [] if type(saved) is type(fresh) and saved == fresh else ['%s: %r against %r' % (path, saved, fresh)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_counts(path=COUNTS):
    with np.load(path, allow_pickle=False) as z:
        require(sorted(z.files) == sorted(FAMILIES), 'the counts file holds other families')
        return {name: z[name].astype(np.int64) for name in FAMILIES}


def verify(results=RESULTS, counts=COUNTS):
    """Recompute from the stored counts and compare with the stored results. Returns the list of differences."""
    saved = json.loads(Path(results).read_text(encoding='utf-8'))
    if saved.get('counts_sha256') != digest(counts):
        return ['the counts file is not the one the results name']
    fresh, _ = compute(load_counts(counts))
    failed = ['gate reads false: ' + name for name, ok in fresh['gates'].items() if not ok]
    body = {k: v for k, v in saved.items() if k not in ('residuals', 'counts_sha256')}
    return failed + compare(body, {k: v for k, v in fresh.items() if k != 'residuals'})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--record', action='store_true', help='draw the counts and write counts.npz and results.json')
    args = parser.parse_args(argv)
    if args.record:
        result, counts = compute()
        require(all(result['gates'].values()), 'a gate reads false: ' + ', '.join(k for k, v in result['gates'].items() if not v))
        np.savez_compressed(COUNTS, **counts)
        result['counts_sha256'] = digest(COUNTS)
        RESULTS.write_text(json.dumps(result, indent=1, allow_nan=False) + '\n', encoding='utf-8')
        print(json.dumps(dict(recorded=True, gates=len(result['gates'])), indent=1))
        return 0
    differences = verify()
    for line in differences:
        print(line, file=sys.stderr)
    print(json.dumps(dict(verified=not differences, candidate=CANDIDATE)))
    return 1 if differences else 0
