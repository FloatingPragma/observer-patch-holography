"""Read calibrated running action gaps from native quantum output records.

This is a source-execution sanity check of a supplied calibrated theory.
The packet supplies the action; it does not supply values to the inverse
decoder. Complete factorized field/pointer/clock outputs supply its scores.
No new laboratory counts, pole-mass prediction or parameter derivation is
claimed. Physical service errors remain separate calibration obligations.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from . import calibrated_matter as matter
from . import calibrated_source_frame as frame_source
from .chiral_candidate import gaussian_source as g, history
from .chiral_candidate import dispersion_identifiability as interval_math
from .mass_response import history as instruments, model
from ..calibrated_measurement.bounds import root_interval
from .halfline import pi_interval
from ..scalar_observations.experimental_readout.stats import concentration_certificate


ROOT = Path(__file__).resolve().parents[2]
BASE_DENOMINATOR = 512
SELECTOR_BITS = 5
IDEAL_ACCEPTED = F(15, 16)
NUMERIC_SCORE_ADMISSION = F(1, 10**8)
RECORD_RADIUS = F(1, 200)
PHYSICAL_SCORE_ALLOWANCE = F(1, 5000)
SAMPLER_SCORE_ALLOWANCE = F(1, 10**7)
WEAK_FRAME_SCORE_ALLOWANCE = F(1, 10**7)
TRIALS_PER_SETTING = 2_000_000
SAMPLER_SEED = 20261026
RECORD_LABELS = tuple(f'{channel}_{family}' for channel in ('up', 'down', 'charged_lepton') for family in range(3))+('clock',)


def _outward(box, bits=40):
    denominator = 2**bits
    lo, hi = box
    return (F((lo*denominator).__floor__(), denominator),
            F((hi*denominator).__ceil__(), denominator))


def catalogue():
    return tuple({'label': f'{channel}_{family}', 'channel': channel, 'family': family,
                  'mass_GeV': mass, 'multiplicity': 2 if channel == 'charged_lepton' else 6}
                 for channel in matter.CHANNELS
                 for family, mass in enumerate(matter.mass_coordinates(channel)))


def _norm(value):
    """Exact squared norm of the original returned binary components."""
    return sum((F(float(z.real))**2+F(float(z.imag))**2 for z in np.asarray(value).ravel()), F())


def _reduced(box):
    lo, hi = box
    pl, ph = pi_interval()
    middle = (lo+hi)/2
    winding = round(float(middle/(pl+ph)))
    candidates = [x-2*winding*p for x in box for p in (pl, ph)]
    reduced = min(candidates), max(candidates)
    if reduced[0] < -4 or reduced[1] > 4:
        raise ValueError('phase reduction is unresolved')
    return float(sum(reduced)/2), reduced, winding


@lru_cache(maxsize=4096)
def _trig(box, cosine=False):
    return interval_math._trig(_outward(box, 64), cosine=cosine)


def _pair_boxes(box):
    cosine, sine = _trig(box, True), _trig(box)
    zero, one = (F(), F()), (F(1), F(1))
    result = [[(zero, zero) for _ in range(4)] for _ in range(4)]
    result[0][0] = result[3][3] = cosine, zero
    result[1][1] = result[2][2] = one, zero
    result[0][3] = result[3][0] = zero, (-sine[1], -sine[0])
    return result


def _operator_error(actual, boxes):
    squared = F()
    for i, row in enumerate(boxes):
        for j, (real, imaginary) in enumerate(row):
            z = actual[i, j]
            squared += max(abs(F(float(z.real))-x) for x in real)**2
            squared += max(abs(F(float(z.imag))-x) for x in imaginary)**2
    return root_interval(squared)[1]


def _projection_error(outputs, box):
    target = _pair_boxes(box)
    squared = F()
    for bit, matrix in enumerate(outputs):
        rows = [[target[i][j] if (i & 1) == bit else ((F(), F()), (F(), F()))
                 for j in range(4)] for i in range(4)]
        squared += _operator_error(matrix, rows)**2
    return root_interval(squared)[1]


@lru_cache(maxsize=4096)
def _source_pair(angle):
    reduced, box, winding = _reduced(angle)
    gate = g.local_gate('pairing', 1., reduced)
    error = _operator_error(gate.fock_unitary, _pair_boxes(box))
    if error > NUMERIC_SCORE_ADMISSION:
        raise ArithmeticError('source pair matrix exceeds its exact binary-array error admission')
    return gate, box, winding, error


def _field_read(field, offset, quadrature):
    control, _, _, _ = _source_pair((-offset/2, -offset/2))
    state = control.fock_unitary@field
    extra = F()
    costs = dict(control.resources)
    if quadrature == 'sine':
        pl, ph = pi_interval()
        detector, _, _, _ = _source_pair((pl/4, ph/4))
        state = detector.fock_unitary@state
        extra = pl/4, ph/4
        for key in g._ADDITIVE_RESOURCES:
            costs[key] += detector.resources[key]
    pointer_word, kraus, _ = g._occupation_pointer()
    pointer_cost = history.compiler.resources(pointer_word)
    outputs = tuple(model.apply_subsystem(state, k, [0], 2) for k in kraus)
    return outputs, costs, pointer_cost, extra


def _clock(clock_initial, phase, offset):
    base, _, _ = _reduced((phase, phase))
    control, _, _ = _reduced((-offset, -offset))
    source = history.compiler.execute([history.compiler.pulse('Z', -base/2, [i]) for i in (0, 1)], 2, clock_initial)
    source = history.compiler.execute([history.compiler.pulse('Z', -control/2, [i]) for i in (0, 1)], 2, source)
    branches = {(): source}
    costs = []
    for site, axis in ((0, 'X'), (1, 'Y')):
        out = {}
        for record, state in branches.items():
            read = instruments.pauli_read(state, site, axis=axis)
            costs.append(read.resources)
            out.update({(*record, bit): value for bit, value in enumerate(read.outputs)})
        branches = out
    probabilities = {key: _norm(value) for key, value in branches.items()}
    means = tuple(-sum(((1 if key[i] == 0 else -1)*value for key, value in probabilities.items()), F()) for i in (0, 1))
    return tuple(branches.values()), probabilities, means, costs


def _distance_to_interval(value, box):
    return max(abs(value-box[0]), abs(value-box[1]))


def _counts(value):
    if (type(value) is not np.ndarray or value.dtype != np.dtype('int64')
            or value.ndim != 5 or value.shape[1:] != (2, 32, 512, 4)
            or not 1 <= value.shape[0] <= 64 or np.any(value < 0)
            or np.any(value > TRIALS_PER_SETTING)
            or np.any(value.sum(axis=(2, 3, 4)) != TRIALS_PER_SETTING)):
        raise ValueError('complete original int64 level/setting/address/field/clock counts required')
    return value


def _count_means(counts):
    patterns = np.arange(512)
    result = {}
    for q, quadrature in enumerate(('cosine', 'sine')):
        accepted = counts[q, :30].sum(axis=0)
        for index, label in enumerate(RECORD_LABELS[:-1]):
            signs = 1-2*((patterns >> index) & 1)
            value = int(np.sum(accepted*signs[:, None]))
            result[label, quadrature] = F(value if q == 0 else -value, TRIALS_PER_SETTING)
        signs = np.asarray([-1 if ((record >> (1-q)) & 1) == 0 else 1 for record in range(4)])
        result['clock', quadrature] = F(int(np.sum(accepted*signs)), TRIALS_PER_SETTING)
    return result


def _score_radius():
    return RECORD_RADIUS+PHYSICAL_SCORE_ALLOWANCE+NUMERIC_SCORE_ADMISSION+SAMPLER_SCORE_ALLOWANCE+WEAK_FRAME_SCORE_ALLOWANCE


def decode_records(counts):
    """Blind inverse: only complete simulated counts and fixed protocol inputs.

    No action, calibration packet, mass catalogue, numerical state or target
    enters this replay. Adaptive offsets are reconstructed from earlier count
    intervals. Caller separately establishes source and sampling assumptions.
    """
    from .calibrated_phase_inference import phase_update
    counts = _counts(counts)
    priors = {label: (F(), F(1)) for label in RECORD_LABELS}
    offsets = []
    radius = _score_radius()
    for level, setting in enumerate(counts):
        offsets.append([str(sum(priors[label])*2**level/2) for label in RECORD_LABELS])
        means = _count_means(setting)
        for label in RECORD_LABELS:
            boxes = [_outward((means[label, q]-radius, means[label, q]+radius)) for q in ('cosine', 'sine')]
            update = phase_update(priors[label], level, *boxes)
            if not update['certified']:
                raise ValueError('unresolved observer records: '+update['reason'])
            priors[label] = update['phase_base_interval']
    clock = priors['clock']
    if clock[0] <= 0:
        raise ValueError('record clock does not resolve a positive common reference')
    intervals = {label: (priors[label][0]/(2*clock[1]), priors[label][1]/(2*clock[0]))
                 for label in RECORD_LABELS[:-1]}
    return {'phase_intervals': priors, 'mass_over_reference_intervals': intervals,
            'analyzer_offsets': np.asarray(offsets)}


def execute(*, levels=21, seed=SAMPLER_SEED, retained_counts=None):
    """Actual source arrays -> complete-record means -> blind interval decoder."""
    if type(levels) is not int or not 1 <= levels <= 64 or type(seed) is not int or seed < 0:
        raise ValueError('one through64 fixed dyadic levels and a nonnegative integer seed required')
    if retained_counts is not None:
        retained_counts = _counts(retained_counts)
        if len(retained_counts) != levels:
            raise ValueError('retained record levels differ from the protocol')
    from .calibrated_phase_inference import phase_update
    frame = frame_source.source_frame()
    rows = catalogue()
    if sum(row['multiplicity'] for row in rows) != 42:
        raise ArithmeticError('calibrated90-mode charged field lost a pair factor')
    masses = [interval_math.expression_interval(row['mass_GeV']) for row in rows]
    if any(lo <= 0 or hi >= 256 for lo, hi in masses):
        raise ValueError('the common fixed256GeV action window was exceeded')
    labels = [row['label'] for row in rows]+['clock']
    priors = {label: (F(), F(1)) for label in labels}
    clock_initial, calibration, clock_cost = history._clock_preparation()
    contrast = -(F(calibration['populations'][0])-F(calibration['populations'][1]))
    coin = history.compiler.execute([history.compiler.pulse('H', math.pi/2, [0])], 1)[:, 0]
    selector_read = instruments.pauli_read(coin, 0, axis='Z')
    coin_weights = tuple(_norm(output) for output in selector_read.outputs)
    selector = tuple(math.prod(coin_weights[(address >> bit) & 1] for bit in range(SELECTOR_BITS)) for address in range(32))
    accepted, failed = sum(selector[:30], F()), sum(selector[30:], F())
    arrays, field_arrays, success_arrays, failure_arrays = {}, [], [], []
    clock_arrays, clock_failure_arrays, offsets_records, count_records = [], [], [], []
    rng = np.random.Generator(np.random.PCG64(seed))
    max_bias = max_operator_error = F()
    resources = {'field_pair_native_pulses': 0, 'analyzer_native_pulses': 0, 'number_pointer_native_pulses': 0,
                 'full90_mode_weak_mass_frame_native_pulses': 6*levels*frame['resources']['native_pulses'],
                 'reference_clock_preparation_native_pulses': 2*levels*clock_cost['native_pulses'],
                 'clock_advance_and_compensation_native_pulses': 8*levels,
                 'clock_pointer_native_pulses': 0,
                 'fair_selector_native_pulses': 10*levels*(1+selector_read.resources['native_pulses'])}
    for level in range(levels):
        scale = F(2**level, BASE_DENOMINATOR)
        offsets = {label: sum(priors[label])*2**level/2 for label in labels}
        offsets_records.append([str(offsets[label]) for label in labels])
        field, field_boxes, unobserved = [], [], F(1)
        for row, mass in zip(rows, masses, strict=True):
            angle = mass[0]*scale, mass[1]*scale
            gate, box, _, error = _source_pair(angle)
            field.append(gate.fock_unitary); field_boxes.append(angle)
            max_operator_error = max(max_operator_error, error)
            unobserved *= _norm(gate.fock_unitary[:, 0])**(row['multiplicity']-1)
            resources['field_pair_native_pulses'] += 2*row['multiplicity']*gate.resources['native_pulses']
        field_arrays.append(np.stack(field))
        clocks, clock_prob, clock_means, clock_pointer_costs = _clock(clock_initial, scale, offsets['clock'])
        clocks_failed, clock_prob_failed, _, _ = _clock(clock_initial, F(), offsets['clock'])
        clock_arrays.append(np.stack(clocks)); clock_failure_arrays.append(np.stack(clocks_failed))
        resources['clock_pointer_native_pulses'] += 2*sum(cost['native_pulses'] for cost in clock_pointer_costs[:2])
        clock_mass, clock_mass_failed = sum(clock_prob.values(), F()), sum(clock_prob_failed.values(), F())
        means, sampled_means, all_settings, failed_settings, setting_counts = {}, {}, [], [], []
        for quadrature in ('cosine', 'sine'):
            matrices, failed_matrices, probabilities, failed_probabilities = [], [], [], []
            for index, row in enumerate(rows):
                outputs, cost, pointer_cost, extra = _field_read(field[index], offsets[row['label']], quadrature)
                failed_outputs, _, _, _ = _field_read(np.eye(4, dtype=complex), offsets[row['label']], quadrature)
                matrices.append(np.stack(outputs)); failed_matrices.append(np.stack(failed_outputs))
                probabilities.append(tuple(_norm(value[:, 0]) for value in outputs))
                failed_probabilities.append(tuple(_norm(value[:, 0]) for value in failed_outputs))
                target = (field_boxes[index][0]-offsets[row['label']]/2,
                          field_boxes[index][1]-offsets[row['label']]/2)
                if extra:
                    target = target[0]+extra[0], target[1]+extra[1]
                _, target_box, _ = _reduced(target)
                max_operator_error = max(max_operator_error, _projection_error(outputs, target_box))
                resources['analyzer_native_pulses'] += cost['native_pulses']
                resources['number_pointer_native_pulses'] += pointer_cost['native_pulses']
            totals = [sum(p, F()) for p in probabilities]
            field_mass = unobserved*math.prod(totals)
            failed_mass = math.prod(sum(p, F()) for p in failed_probabilities)
            if abs(accepted*field_mass*clock_mass+failed*failed_mass*clock_mass_failed-1) > F(1, 10**7):
                raise ArithmeticError('complete source selector/field/clock instrument lost probability')
            for index, row in enumerate(rows):
                sign = 1 if quadrature == 'cosine' else -1
                mean = sign*accepted*(probabilities[index][0]-probabilities[index][1])*unobserved*clock_mass
                mean *= math.prod(totals[k] for k in range(9) if k != index)
                means[row['label'], quadrature] = mean
            clock_index = int(quadrature == 'sine')
            means['clock', quadrature] = accepted*field_mass*clock_means[clock_index]
            # A declared simulator samples a normalized probability law. The
            # source matrices above remain unchanged; their raw-mass error is
            # bounded separately. No experimental counts are asserted.
            field_law = np.ones(512)
            failed_field_law = np.ones(512)
            patterns = np.arange(512)
            for index in range(9):
                bit = (patterns >> index) & 1
                p = probabilities[index]
                pf = failed_probabilities[index]
                field_law *= np.asarray([float(p[0]/sum(p)), float(p[1]/sum(p))])[bit]
                failed_field_law *= np.asarray([float(pf[0]/sum(pf)), float(pf[1]/sum(pf))])[bit]
            cq = np.asarray([float(value/clock_mass) for value in clock_prob.values()])
            cf = np.asarray([float(value/clock_mass_failed) for value in clock_prob_failed.values()])
            selector_mass = accepted+failed
            law = np.stack([float(weight/selector_mass)*(field_law[:, None]*cq if address < 30
                       else failed_field_law[:, None]*cf)
                       for address, weight in enumerate(selector)])
            law /= law.sum()
            counts = (rng.multinomial(TRIALS_PER_SETTING, law.ravel()).reshape(32, 512, 4)
                      if retained_counts is None else retained_counts[level, int(quadrature == 'sine')])
            setting_counts.append(counts)
            accepted_counts = counts[:30].sum(axis=0)
            for index, row in enumerate(rows):
                signs = 1-2*((patterns >> index) & 1)
                total = int(np.sum(accepted_counts*signs[:, None]))
                if quadrature == 'sine':
                    total = -total
                sampled_means[row['label'], quadrature] = F(total, TRIALS_PER_SETTING)
                normalized = (accepted/selector_mass)*(probabilities[index][0]-probabilities[index][1])/totals[index]
                if quadrature == 'sine':
                    normalized = -normalized
                if abs(normalized-means[row['label'], quadrature]) > F(1, 10**9):
                    raise ArithmeticError('simulator normalization exceeds numerical score budget')
            # branches insertion order is (Xbit,Ybit): index=2*Xbit+Ybit.
            clock_signs = np.asarray([-1 if ((record >> (1-clock_index)) & 1) == 0 else 1 for record in range(4)])
            sampled_means['clock', quadrature] = F(int(np.sum(accepted_counts*clock_signs)), TRIALS_PER_SETTING)
            normalized_clock = (accepted/selector_mass)*clock_means[clock_index]/clock_mass
            if abs(normalized_clock-means['clock', quadrature]) > F(1, 10**9):
                raise ArithmeticError('clock sampler normalization exceeds numerical score budget')
            all_settings.append(np.stack(matrices)); failed_settings.append(np.stack(failed_matrices))
        count_records.append(np.stack(setting_counts))
        success_arrays.append(np.stack(all_settings)); failure_arrays.append(np.stack(failed_settings))
        for index, label in enumerate(labels):
            if label == 'clock':
                phase = scale-offsets[label], scale-offsets[label]
                amplitude = IDEAL_ACCEPTED*contrast
            else:
                phase = (2*field_boxes[index][0]-offsets[label], 2*field_boxes[index][1]-offsets[label])
                amplitude = IDEAL_ACCEPTED
            _, phase_box, _ = _reduced(phase)
            boxes = []
            for quadrature in ('cosine', 'sine'):
                ideal = _trig(phase_box, cosine=quadrature == 'cosine')
                ideal = ideal[0]*amplitude, ideal[1]*amplitude
                mean = means[label, quadrature]
                bias = _distance_to_interval(mean, ideal)
                max_bias = max(max_bias, bias)
                if bias > NUMERIC_SCORE_ADMISSION:
                    raise ArithmeticError('exact returned-array record bias exceeds fixed numerical allowance')
                radius = _score_radius()
                sampled = sampled_means[label, quadrature]
                boxes.append(_outward((sampled-radius, sampled+radius)))
            update = phase_update(priors[label], level, boxes[0], boxes[1])
            if not update['certified']:
                raise ArithmeticError('observer phase is unresolved: '+str(update.get('reason')))
            priors[label] = update['phase_base_interval']
    clock_box = priors['clock']
    if clock_box[0] <= 0:
        raise ArithmeticError('common reference clock remained unresolved')
    recovered = []
    for row, mass in zip(rows, masses, strict=True):
        lower, upper = priors[row['label']]
        interval = lower/(2*clock_box[1]), upper/(2*clock_box[0])
        if not interval[0] <= mass[0] <= mass[1] <= interval[1]:
            raise ArithmeticError('record-derived running-action interval missed its supplied action')
        recovered.append({'channel': row['channel'], 'family': row['family'],
            'record_gap_over_common_reference_interval': list(map(str, interval)),
            'record_interval_numeric_GeV': list(map(float, interval)),
            'record_interval_width_GeV': float(interval[1]-interval[0]),
            'supplied_running_action_mass_GeV': str(row['mass_GeV']),
            'calibrated_action_inside_record_interval': True})
    arrays.update(field_unitaries=np.stack(field_arrays), success_read_matrices=np.stack(success_arrays),
        failed_read_matrices=np.stack(failure_arrays), clock_reference_outputs=np.stack(clock_arrays),
        failed_clock_reference_outputs=np.stack(clock_failure_arrays), selector_coin=coin,
        analyzer_offsets=np.asarray(offsets_records), CKM_fock_unitary=frame['CKM_fock_unitary'],
        singlet_phase_unitary=frame['singlet_phase_unitary'],
        one_body_basis_change90=frame['one_body_basis_change90'],
        selector_read_outputs=np.stack(selector_read.outputs), complete_simulated_record_counts=np.stack(count_records))
    down_angles = tuple(_reduced(field_boxes[index])[0] for index in (3, 4, 5))
    weak = frame_source.source_defined_down_sector(down_angles)
    arrays.update(weak_down_source_unitary=weak['weak_unitary'], weak_down_mass_source_unitary=weak['mass_pair_product'])
    frame_error = (F(1)+F(1, 5_000_000_000))**24-1
    three_frames = (1+frame_error)**3-1
    if 2*three_frames+three_frames**2 > WEAK_FRAME_SCORE_ALLOWANCE:
        raise ArithmeticError('weak-frame bounded-score propagation exceeds its separate allowance')
    blind = decode_records(arrays['complete_simulated_record_counts'])
    if blind['phase_intervals'] != priors or not np.array_equal(blind['analyzer_offsets'], arrays['analyzer_offsets']):
        raise ArithmeticError('independent packet-blind complete-record inverse disagrees')
    sampling = concentration_certificate(20*levels, TRIALS_PER_SETTING, math.nextafter(float(RECORD_RADIUS), 0.))
    if F(sampling['width_exact']) > RECORD_RADIUS or F(sampling['familywise_error_upper']) >= F(1, 100):
        raise ArithmeticError('conditional simultaneous sampling design is insufficient')
    resources['total_native_pulses_per_complete_setting_schedule'] = sum(resources.values())
    resources['proper_code_load_export_service_events'] = 2*resources['total_native_pulses_per_complete_setting_schedule']
    resources['hypothetical_repeated_native_campaign_pulses'] = TRIALS_PER_SETTING*resources['total_native_pulses_per_complete_setting_schedule']
    resources['hypothetical_repeated_native_campaign_service_events'] = TRIALS_PER_SETTING*resources['proper_code_load_export_service_events']
    resources['accounting_basis'] = 'sum of one-trial pulse costs over all settings; clock pointers count one retained branch path, not repeated host queries'
    resources['service_event_counts_are_charged_interface_requirements_not_executed_protection_histories'] = True
    resources['classical_committed_record_bits_per_attempt'] = 16
    return {'schema': 'oph.calibrated_spectrum_observer.v1',
        'classification': 'native implementation sanity check of a supplied empirically calibrated running action',
        'records': {'input_kind': 'fixed-seed simulated complete instrument records; no acquired laboratory counts',
            'sampler': 'numpy.Generator(PCG64).multinomial; normalized source Born law', 'sampler_seed': seed,
            'trials_per_setting': TRIALS_PER_SETTING,
            'all_attempt_simulated_trials': 2*levels*TRIALS_PER_SETTING,
            'simulator_normalization_is_separate_from_unchanged_source_arrays': True,
            'dyadic_levels': levels, 'joint_selector_field_clock_outcomes_per_setting': 65536,
            'complete_settings': 2*levels, 'explicit_selector_addresses': 32,
            'accepted_addresses': list(range(30)), 'timeout_addresses': [30, 31],
            'factorization': 'address weight times nine complete binary field factors times four clock/reference factors times unobserved field norms',
            'all_outcomes_count_in_means': True, 'success_postselection_or_probability_renormalization': False},
        'field': {'CAR_modes': 90, 'charged_rest_pair_factors': 42,
            'distinct_charged_running_gaps': 9, 'neutral_zero_modes_retained': 6,
            'renormalization_convention': matter.parameters()['common_scale'],
            'mass_basis_is_same_weak_action': 'W45.T K45 W45; currents transform by W45.H T W45',
            'source_frame': {'mapping': frame['mapping'], 'diagnostics': frame['diagnostics'],
                'complete_frame_native_pulses': frame['resources']['native_pulses'],
                'application': 'source-defined conjugation and mass-coordinate readout; conservative three complete frames charged per setting',
                'mass_coordinate_outputs_are_literal_composed90_mode_binary_arrays': False,
                'independent_source_weak_down_block_columns': weak['source_all_input_columns'],
                'weak_down_block_level': levels-1, 'weak_down_block_errors': weak['errors'],
                'weak_frame_operator_admission_per_CKM_or_sign_factor_assumed': '1/5000000000',
                'three_frame_bounded_score_propagation_upper': str(_outward((F(), 2*three_frames+three_frames**2))[1]),
                'weak_frame_admission_is_an_outward_analytic_or_physical_certificate': False},
            'spectral_compilation': 'exact commuting rest-pair rotations; integer2pi reduction retains the whole four-column Fock action',
            'per_particle_or_per_dataset_clock_fit': False},
        'clock': {'native_generator_gap': 'log(5)', 'one_global_energy_dictionary': 'E0=1GeV/log(5)',
            'source_modular_parameter_s': '2**level/(512*log(5))',
            'field_generator_gap': 'm_GeV*log(5)', 'clock_phase': 'log(5)*s=2**level/512',
            'large_modular_parameter_implementation': 'programmed exact commuting spectral phase reduction; not autonomous hardware waiting or a coherence-time demonstration',
            'common_reference_energy_GeV': '1',
            'actual_physical_clock_hardware_calibrated': False, 'base_phase_interval': list(map(str, clock_box))},
        'inversion': {'decoder_uses_calibration_masses': False, 'same_visibility_cancels_in_sine_cosine_ratio': True,
            'aliases_resolved_by_record_selected_analyzer_offsets': True,
            'input_score_radius_assumed': str(RECORD_RADIUS), 'physical_score_bias_allowance_assumed': str(PHYSICAL_SCORE_ALLOWANCE),
            'separate_sampler_normalization_and_float_score_allowance': str(SAMPLER_SCORE_ALLOWANCE),
            'sampler_float_budget_contract': 'binary64 bounded probability operations;65536 cells times at most64 absolute rounding allowances2**-52, plus normalization; total below1e-7',
            'separate_source_defined_weak_frame_score_allowance': str(WEAK_FRAME_SCORE_ALLOWANCE),
            'exact_binary_source_score_bias_upper': str(_outward((F(), max_bias))[1]),
            'exact_binary_pair_analyzer_pointer_all_column_operator_error_upper': str(_outward((F(), max_operator_error))[1]),
            'mean_rectangle_outward_rounding_bits': 40,
            'simultaneous_sampling_certificate': sampling,
            'sampling_premise': 'independent simulated trials conditional on past records and chosen offsets; ten scores within a trial may be correlated',
            'adaptive_offsets_and_inverse_use_simulated_counts_only': True,
            'hardware_or_finite_shot_coverage_authenticated_here': False},
        'observer_read_running_spectrum': recovered, 'source_resources': resources,
        'scope': {'masses_are_empirical_inputs_not_derived_predictions': True,
            'calibration_recovery_is_independent_empirical_confirmation': False,
            'outputs_are_running_action_parameters_not_physical_poles': True,
            'interacting_self_energies_confinement_or_scattering_executed': False,
            'native_CE_protection_and_physical_export_errors_certified_here': False,
            'numerical_record_certificate_applies_to_returned_binary_arrays': True}}, arrays


def source_bindings():
    names = set(matter.source_bindings())|set(frame_source.source_bindings())|{
        'model/common_source/calibrated_spectrum_observer.py', 'model/common_source/test_calibrated_spectrum_observer.py',
        'model/common_source/calibrated_phase_inference.py', 'model/common_source/test_calibrated_phase_inference.py',
        'model/common_source/chiral_candidate/gaussian_source.py', 'model/common_source/chiral_candidate/history.py',
        'model/common_source/chiral_candidate/dispersion_identifiability.py',
        'model/common_source/mass_response/history.py', 'model/common_source/mass_response/model.py',
        'model/native_modular_reference/phase_inference.py', 'model/calibrated_measurement/bounds.py',
        'model/native_modular_reference/reference_tomography.py',
        'model/scalar_observations/experimental_readout/stats.py',
        'model/common_source/halfline.py'}
    names.update(history.compiler.source_summary()['source_sha256'])
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _same_report(expected, actual, path='report'):
    if type(expected) is not type(actual):
        raise ValueError(path+' changed its original JSON type')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' keys differ')
        for key in expected:
            _same_report(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _same_report(left, right, path+f'[{i}]')
    elif type(expected) is float:
        if path.endswith('.numerical_operator_admission') and expected != actual:
            raise ValueError(path+' fixed admission differs')
        if not math.isfinite(actual) or not math.isclose(expected, actual, rel_tol=2e-11, abs_tol=2e-12):
            raise ValueError(path+' numerical diagnostic differs')
    elif expected != actual:
        raise ValueError(path+' differs')


def compare_execution(report, arrays, fresh_report, fresh_arrays):
    """Strict schemas and exact records; platform-tolerant fresh source arrays."""
    _same_report(fresh_report, report)
    if set(arrays) != set(fresh_arrays):
        raise ValueError('complete artifact array keys differ')
    for name, expected in fresh_arrays.items():
        actual = arrays[name]
        if (type(actual) is not np.ndarray or actual.dtype != expected.dtype
                or actual.shape != expected.shape):
            raise ValueError('artifact schema differs: '+name)
        if actual.dtype.kind in 'Ui':
            if not np.array_equal(actual, expected):
                raise ValueError('retained records/offsets differ: '+name)
        elif (not np.isfinite(actual).all()
              or not np.allclose(actual, expected, rtol=2e-11, atol=2e-12)):
            raise ValueError('fresh source output differs: '+name)


def check(path):
    report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != source_bindings():
        raise ValueError('retained source bindings differ from the current closed source set')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('retained complete quantum/count artifact custody differs')
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError('duplicate complete artifact array names are inadmissible')
        arrays = {key: archive[key] for key in archive.files}
    counts = _counts(arrays.get('complete_simulated_record_counts'))
    seed = report['records']['sampler_seed']
    fresh, fresh_arrays = execute(levels=len(counts), seed=seed, retained_counts=counts)
    scientific = {key: value for key, value in report.items() if key not in ('source_bindings', 'artifact')}
    compare_execution(scientific, arrays, fresh, fresh_arrays)
    return scientific


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--output', type=Path)
    choice.add_argument('--check', type=Path)
    parser.add_argument('--levels', type=int, default=21)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('calibrated observer complete-record/source replay: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('choose fresh JSON and quantum artifact paths')
    report, arrays = execute(levels=args.levels)
    report['source_bindings'] = source_bindings()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps([row['record_interval_numeric_GeV'] for row in report['observer_read_running_spectrum']]))


if __name__ == '__main__':
    main()
