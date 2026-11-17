"""Actual native clock maps, full source archives, and calibration freedom."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, permutations
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
from scipy import linalg as la

from model.clocked_common_observer import clock
from model.lepton_cycle import instrument as lepton
from model.lepton_cycle import composition
from model.lepton_mass_compiler import analyzer
from model.native_flavor_orbit import intervals as iv
from model.particle_common_contract.provenance import import_closure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
THETA, KICK = F(2, 9), F(1, 7)
CURRENT = np.array([-1, 0, 1])
MARKS = ('idle', 'zero', 'minus', 'plus')


def encoded(array):
    a = np.asarray(array)
    return {'shape': list(a.shape), 'real': a.real.ravel().tolist(),
            'imaginary': a.imag.ravel().tolist()}


@lru_cache(maxsize=1)
def source_blocks():
    records, _ = lepton.source_data()
    masks = (np.array([x not in (2, 3) for x in records]),
             np.array([x == 2 for x in records]), np.array([x == 3 for x in records]))
    counts = [int(x.sum()) for x in masks]
    if counts != [50, 5, 5]:
        raise ArithmeticError('actual native source partition differs')
    basis = np.column_stack([m/np.sqrt(n) for m, n in zip(masks, counts, strict=True)])
    lift = np.kron(basis, basis)
    uniform = basis.T @ np.full(60, 1/np.sqrt(60))
    standing = np.kron(uniform, uniform)
    p = np.outer(uniform, uniform)
    returns = (np.kron(p, p), np.kron(np.eye(3)-p, p),
               np.kron(p, np.eye(3)-p), np.kron(np.eye(3)-p, np.eye(3)-p))
    bank, errors = {}, {}
    for kind, parts in (('weak', lepton.weak_instrument(float(KICK))),
                        ('mass', lepton.mass_instrument(float(THETA), float(KICK)))):
        bank[kind] = {}
        for name in lepton.LABELS:
            original = parts[name].reshape(3600, 3)
            reduced = lift.T @ original
            errors[kind+'/'+name] = float(la.norm(lift @ reduced-original))
            bank[kind][name] = reduced
        diagonal = sum(np.sum(abs(v)**2, axis=0) for v in bank[kind].values())
        if max(abs(diagonal-1)) > 3e-12:
            raise ArithmeticError('complete native return instrument differs')
    if max(errors.values()) > 3e-12 or la.norm(lift.T@lift-np.eye(9)) > 3e-12:
        raise ArithmeticError('native category isometry failed')
    # Pull the actual full analyzer projector through the same isometry.
    x, effects, amplitudes = lepton.analyzer_data(float(THETA))
    pulled = x.reshape(3600, 3)/amplitudes[None, :]
    pulled = lift.T @ pulled
    return {'lift': lift, 'standing': standing, 'returns': returns, 'bank': bank,
            'effects': effects, 'mass_projector_vectors': pulled, 'errors': errors}


def archive_columns(outcomes):
    values = np.ones((1, 3), complex)
    for kind, label in zip(('weak', 'mass'), outcomes):
        factor = source_blocks()['bank'][kind][label]
        values = np.einsum('ak,bk->abk', values, factor).reshape(-1, 3)
    return values


def clock_words(cap):
    if type(cap) is not int or not 2 <= cap <= 4:
        raise ValueError('a bounded cap between two and four is required')
    result = []
    def visit(word, events):
        if events == 2 or len(word) == cap:
            result.append(word)
            return
        for mark in MARKS:
            visit(word+(mark,), events+(mark in ('minus', 'plus')))
    visit((), 0)
    return result


@lru_cache(maxsize=256)
def base_clock_state(law, active_word):
    if not active_word:
        return clock.relational.exact_preparation(1)
    before = base_clock_state(law, active_word[:-1])
    source = clock.source(law)
    mark = active_word[-1]
    if mark == 'zero':
        return {key: v*F(source.no_current_counts[key[0]][key[1]], 60)
                for key, v in before.items() if source.no_current_counts[key[0]][key[1]]}
    if mark not in ('minus', 'plus'):
        raise ValueError('invalid native mark')
    return clock._moving(before, source.moving_maps, F(1))[-1 if mark == 'minus' else 1]


@lru_cache(maxsize=256)
def numerical_clock_state(law, active_word):
    if not active_word:
        return clock.quantum.phase_apparatus(1)[0]
    before = numerical_clock_state(law, active_word[:-1])
    index = {'minus': 0, 'zero': 1, 'plus': 2}[active_word[-1]]
    return clock.runtime._apply(clock.source(law).instrument, before, 1)[index]


def stopped_clock(law='C', cap=3, activity=F(1)):
    if type(activity) is not F or not 0 < activity <= 1:
        raise ValueError('exact activity in (0,1] required')
    rows = []
    error_sum = 0.
    for word in clock_words(cap):
        active = tuple(x for x in word if x != 'idle')
        weight = activity**len(active)*(1-activity)**word.count('idle')
        state = {key: weight*v for key, v in base_clock_state(law, active).items()}
        matrix = np.zeros((66, 66), complex)
        for key, value in state.items():
            matrix[key] = float(value)
        executed = float(weight)*numerical_clock_state(law, active)
        error = float(la.norm(matrix-executed, 'fro'))
        if error > 3e-14:
            raise ArithmeticError('actual native execution disagrees with rational history')
        error_sum += error
        mass = sum((v for (i, j), v in state.items() if i == j), F())
        rows.append({'word': list(word), 'active_attempts': len(active),
                     'event_count': sum(x in ('minus', 'plus') for x in word),
                     'probability_exact': str(mass),
                     'source_state': [[i, j, str(v)] for (i, j), v in sorted(state.items())],
                     'numeric_frobenius_error': error})
    if sum(F(r['probability_exact']) for r in rows) != 1:
        raise ArithmeticError('a stopped native history was lost')
    return rows, error_sum


def joint_run(law='C', activity=F(1), initial_idles=0, work_policy='event', cap=3):
    if work_policy not in ('event', 'slot') or initial_idles not in (0, 1):
        raise ValueError('undeclared work clock or idle prefix')
    native, error = stopped_clock(law, cap, activity)
    fourier = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    rows, detector = [], np.zeros(3)
    timeout = F()
    mean_slots = mean_attempts = mean_events = F()
    for native_index, entry in enumerate(native):
        cycles = entry['event_count'];probability = F(entry['probability_exact'])
        slots = len(entry['word'])+initial_idles
        mean_slots += probability*slots
        mean_attempts += probability*entry['active_attempts']
        mean_events += probability*cycles
        if cycles < 2:
            timeout += probability
        phase_time = float(KICK)*slots if work_policy == 'slot' else 0.
        phase = np.exp(-1j*phase_time*CURRENT)
        for outcomes in product(lepton.LABELS, repeat=cycles):
            archive = archive_columns(outcomes)*phase[None, :]
            output = archive @ fourier.conj()/np.sqrt(3)
            values = np.sum(abs(output)**2, axis=0)
            for d, value in enumerate(values):
                original = float(probability)*float(value)
                detector[d] += original
                rows.append({'clock_history_index': native_index,
                             'lepton_returns': list(outcomes), 'detector': d,
                             'probability': original, 'controller_slots': slots,
                             'native_attempts': entry['active_attempts'],
                             'completed_cycles': cycles, 'phase_time': phase_time,
                             'timeout': cycles < 2})
    if abs(sum(r['probability'] for r in rows)-1) > 2e-12:
        raise ArithmeticError('joint timeline lost original weight')
    delta = F(2, 15)*activity
    tail_bound = (1-delta)**cap + cap*delta*(1-delta)**(cap-1)
    if timeout > tail_bound:
        raise ArithmeticError('two-event all-input tail bound failed')
    return {'native_law': law, 'activity_exact': str(activity), 'cap': cap,
            'initial_idle_slots': initial_idles, 'work_policy': work_policy,
            'clock_histories': native, 'complete_joint_records': rows,
            'complete_record_count': len(rows), 'detector_probabilities': detector.tolist(),
            'timeout_probability_exact': str(timeout),
            'all_input_two_event_timeout_bound_exact': str(tail_bound),
            'mean_controller_slots_exact': str(mean_slots),
            'mean_native_attempts_exact': str(mean_attempts),
            'mean_completed_cycles_exact': str(mean_events),
            'native_numeric_frobenius_error_sum': error}


def no_reset_control():
    blocks = source_blocks();v = blocks['mass_projector_vectors']
    differences = [];rows = []
    for first in lepton.LABELS:
        old = blocks['bank']['weak'][first]
        evolved = old+np.expm1(-1j*float(KICK))*v*np.sum(v.conj()*old, axis=0)[None, :]
        for label, projector in zip(lepton.LABELS, blocks['returns'], strict=True):
            actual = projector @ evolved
            fresh = archive_columns((first, label))
            old_probability = float(np.sum(abs(actual)**2)/3)
            fresh_probability = float(np.sum(abs(fresh)**2)/3)
            differences.append(abs(old_probability-fresh_probability))
            rows.append({'returns': [first, label], 'without_fresh_pair': old_probability,
                         'with_fresh_pair': fresh_probability,
                         'no_reset_full_source_columns': encoded(actual)})
    tv = sum(differences)/2
    if tv < 1e-7 or abs(sum(r['without_fresh_pair'] for r in rows)-1) > 2e-12:
        raise ArithmeticError('correlated-source reset control disappeared')
    return {'complete_sixteen_return_laws': rows, 'return_total_variation': tv,
            'same_initial_declared_source_and_weak_cycle': True,
            'unrefreshed_source_correlations_retained': True,
            'fresh_reset_claimed_autonomous': False}


def generator_and_scale():
    effect_rows = analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']
    effects = [iv.I(F(lo), F(hi)) for lo, hi in effect_rows]
    curvatures = []
    for order in permutations(range(3)):
        value = effects[order[0]]-2*effects[order[1]]+effects[order[2]]
        if value.lo <= 0 <= value.hi:
            raise ArithmeticError('affine current/energy obstruction not separated')
        curvatures.append({'ordering': list(order), 'second_difference_interval': value.strings()})
    weak = lepton.unconditional_weak_multiplier(float(KICK))
    mass = lepton.unconditional_mass_multiplier(float(THETA), float(KICK))
    if np.max(abs(weak-np.ones((3, 3)))) < 1e-6 or abs(weak[0, 1]) >= 1:
        raise ArithmeticError('finite-round decoherence control failed')
    scaling = []
    e = source_blocks()['effects']
    for factor in (F(1), F(2), F(7, 3)):
        energy = float(factor)*e;time = float(KICK/factor)
        error = float(max(abs(np.exp(-1j*energy*time)-np.exp(-1j*e*float(KICK)))))
        if error > 5e-16:
            raise ArithmeticError('dimensionful rescaling changed an admitted word')
        scaling.append({'unit_factor': str(factor), 'time_factor': str(1/factor),
                        'native_counts_changed': False, 'phase_identity_error': error})
    return {'limiting_weak_generator': 'diag(-2/9,0,2/9)',
            'limiting_mass_generator': 'diag(e_0,e_1,e_2)',
            'source_effect_intervals': effect_rows,
            'all_six_affine_dictionary_obstructions': curvatures,
            'distinct_generators_do_not_prevent_a_shared_timeline': True,
            'finite_weak_multiplier': encoded(weak), 'finite_mass_multiplier': encoded(mass),
            'finite_round_channels_are_unitary': False,
            'weak_half_diamond_bound': composition.error_budget(1, KICK, cycle='weak'),
            'mass_half_diamond_bound': composition.error_budget(1, KICK, cycle='mass'),
            'dimensionful_rescaling_controls': scaling,
            'exact_rescaling_identity': 'exp(-i (c H)(t/c))=exp(-i H t); all dimensionless service parameters and original records are fixed',
            'absolute_GeV_or_SI_scale_derived': False,
            'electron_anchor_role': 'one empirical energy input removes this units freedom; it does not identify the interacting particle or native attempt duration'}


def bindings():
    paths = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    paths.update(('model/lepton_cycle/PROOF.md', 'model/lepton_service_compiler/PROOF.md',
                  'model/lepton_service_compiler/results.json',
                  'model/clocked_common_observer/README.md', 'model/lepton_mass_compiler/PROOF.md',
                  'Lean/Dynamics/ChoiCPTP.lean'))
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def certificate():
    pins = bindings()
    settings = [('C', F(1), 0, 'event'), ('C', F(1), 1, 'event'),
                ('C', F(1), 0, 'slot'), ('C', F(1), 1, 'slot'),
                ('C', F(1, 2), 0, 'event'), ('E', F(1), 0, 'event')]
    runs = [joint_run(*setting) for setting in settings]
    baseline = np.array([x['probability'] for x in runs[0]['complete_joint_records']])
    delayed = np.array([x['probability'] for x in runs[1]['complete_joint_records']])
    if not np.array_equal(baseline, delayed):
        raise ArithmeticError('event-only coupling changed under deterministic idle prefix')
    drift = sum(abs(np.array(runs[2]['detector_probabilities'])-runs[3]['detector_probabilities']))/2
    if drift < 1e-5:
        raise ArithmeticError('slot-driven work failed to observe inserted idle')
    blocks = source_blocks()
    result = {'schema': 'oph.lepton_clock_energy.v1', 'source_bindings': pins,
              'all_checks_pass': True,
              'source_contract': {'native_clock_dimension': 66, 'original_address_pair_dimension': 3600,
                  'family_dimension': 3, 'source': 'declared product of the actual native C/E clock carrier and fresh coherent lepton cells, with event-triggered services',
                  'family_wire_map': 'weak current j=k-1 and mass character k share one declared three-coordinate register',
                  'source_record_counts': [50, 5, 5], 'theta_exact': str(THETA), 'kick_exact': str(KICK),
                  'source_primitive_is_not_a_single_passive_repair': True,
                  'quantum_reference_extension': 'every branch tensors with identity on arbitrary untouched references'},
              'exact_category_transfer': {'source_dimension': 9, 'lift_errors': blocks['errors'],
                  'standing_source': encoded(blocks['standing']),
                  'full_archive_factors': {kind: {label: encoded(v) for label, v in row.items()}
                                           for kind, row in blocks['bank'].items()},
                  'full_original_source_reconstruction': 'B tensor B on each retained archive; B has the actual native record-type supports and entries 1/sqrt(50),1/sqrt(5),1/sqrt(5)'},
              'timeline_runs': runs, 'idle_prefix_event_law_identical': True,
              'idle_prefix_slot_work_detector_TV': float(drift),
              'correlated_source_control': no_reset_control(),
              'clock_energy_result': generator_and_scale(),
              'resources': {'maximum_native_attempts': 3, 'maximum_lepton_cycles': 2,
                  'maximum_address_pairs_including_final_fresh_pair': 3,
                  'per_cycle_coherent_forward_ticks': 2, 'per_cycle_inverse_ticks': 2,
                  'per_cycle_return_tests': 2, 'per_cycle_archive_swaps': 2,
                  'retained_clock_record_alphabet': list(MARKS),
                  'archive_storage': 'two complete address pairs and two two-bit return records; compression is a separate service',
                  'supplied_services': ['fresh blanks', 'coherent preparation and kick', 'routing', 'return tests', 'SWAP/archive', 'final three-outcome Fourier detector'],
                  'service_overhead_identified_with_native_attempts': False},
              'scope': {'complete_same_declared_protocol_timeline': True,
                  'waiting_idles_timeouts_and_all_returns_retained': True,
                  'full_quantum_archives_and_native_state_retained': True,
                  'numerical_execution_is_interval_certified': False,
                  'physical_energy_clock_conversion_selected': False,
                  'physical_particle_masses_or_couplings_derived': False}}
    if bindings() != pins:
        raise ArithmeticError('source changed during timeline execution')
    return result


def equal(a, b):
    if isinstance(b, dict):
        return isinstance(a, dict) and set(a) == set(b) and all(equal(a[k], v) for k, v in b.items())
    if isinstance(b, list):
        return isinstance(a, list) and len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    if isinstance(b, float):
        return type(a) in (int, float) and np.isfinite(a) and abs(a-b) <= 2e-13+2e-9*abs(b)
    return type(a) is type(b) and a == b


def main():
    parser = argparse.ArgumentParser(__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true');group.add_argument('--check', action='store_true')
    args = parser.parse_args();result = certificate();path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()), result):
        raise ArithmeticError('complete clock and lepton receipt differs')
    print('Verified complete native-event lepton timeline and explicit energy-scale freedom.')
