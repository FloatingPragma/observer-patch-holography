"""Complete17-bit native experiment and independent original-H comparator."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..cube_eta_protection import cube, repair, source as preparation
from ..observer.source_scalar import compiler
from ..reference_observation.observation import frozen
from . import circuits as c, execute, bounds


PHASES = (F(), F(1, 2))


def repair_layer(state, law):
    records = []
    for seam in repair.SEAMS:
        fields = state.reshape(2, 65536, -1).transpose(1, 0, 2).reshape(65536, -1)
        for label, out in repair.apply(fields, law, seam):
            value = out.reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
            records.append({'law': law, 'native_seam': [repair.PORTS[j] for j in seam],
                            'record': list(label), 'probability': float(np.vdot(value, value).real)})
            if label == ('outside',): continuation = value
        state = continuation
    return state, records


@dataclass(frozen=True)
class Prefix:
    state: np.ndarray
    runs: tuple
    records: tuple


@lru_cache(maxsize=1)
def prefix():
    tree, copied, _ = preparation.runs()
    initial = np.zeros((1 << 17, 1), complex); initial[:65536] = copied.success
    pointer = execute.complete([c.block(c.pointer_prepare(), 3, (0, 1, 16))], 17, initial)
    blocks = c.field_binding(bounds.TIME/2, bounds.STEPS)
    first = execute.complete(blocks, 17, pointer.success)
    state, r1 = repair_layer(first.success, 'C')
    second = execute.complete(blocks, 17, state)
    state, r2 = repair_layer(second.success, 'E')
    return Prefix(frozen(state), (tree, copied, pointer, first, second), tuple(r1+r2))


def repair_output(law, seam_index, record):
    """Reconstruct any native or outside vector from its actual source prefix."""
    if law not in ('C', 'E') or type(seam_index) is not int or not 0 <= seam_index < len(repair.SEAMS):
        raise ValueError('existing repair layer and seam index required')
    stored = prefix(); state = stored.runs[3 if law == 'C' else 4].success
    for index, seam in enumerate(repair.SEAMS):
        fields = state.reshape(2, 65536, -1).transpose(1, 0, 2).reshape(65536, -1)
        outcomes = repair.apply(fields, law, seam)
        if index == seam_index:
            for label, value in outcomes:
                if label == record:
                    return value.reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
            raise ValueError('no such complete repair outcome')
        state = outcomes[-1][1].reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
    raise ArithmeticError('unreachable repair outcome')


def ideal_analyzed(duration, phase, *, scalar_shift=0.):
    """Direct original H, using its exact invariant N6 domain only here."""
    reference = preparation.ideal_reference()[cube.weights(6)]
    occupation = np.array([(int(mask) & 3).bit_count() for mask in cube.weights(6)])
    local_phase = np.exp(1j*np.pi*occupation/2)
    arms = np.column_stack((reference, local_phase*reference))/math.sqrt(2)
    h = cube.restrict(cube.operators().H.matrix(cube.operators().scale), 6)
    evolved = expm_multiply(-1j*float(duration)*h, arms)*np.exp(-1j*float(duration)*scalar_shift)
    evolved[:, 1] *= local_phase.conj()*np.exp(1j*np.pi*float(phase))
    output = np.zeros((1 << 17, 1), complex)
    output[cube.weights(6), 0] = (evolved[:, 0]+evolved[:, 1])/math.sqrt(2)
    output[65536+cube.weights(6), 0] = (-evolved[:, 0]+evolved[:, 1])/math.sqrt(2)
    return output


def pointer_outputs(state):
    if state.shape != (131072, 1): raise ValueError('full17-bit field/pointer vector required')
    kraus = compiler.transfer_kraus(2, coherent=False)
    if len(kraus) != 6 or not np.array_equal(sum(k.conj().T@k for k in kraus), np.eye(6)):
        raise ArithmeticError('complete native pointer export lost an outcome')
    return tuple(k[:, :2]@state.reshape(2, 65536) for k in kraus)


def source_record(stage, run):
    if isinstance(run, execute.Run):
        tape = tuple(event for instrument, _ in run.blocks for event in instrument.tape)
        rounded = run.rounded_negative_weight_count, run.rounded_negative_weight_total
        summation = run.failure_sum_error, run.failure_sum_admission
    else:
        tape = run.word; rounded = 0, 0.; summation = 0., 0.
    return {'stage': stage, 'pulses': len(tape),
            'native_abort_slots': sum(6-(1 << len(targets)) for _, _, targets in tape),
            'padding_slots': sum(2 for _, _, targets in tape if len(targets) == 2),
            'failure_weight_array': stage+'_failure_weights',
            'first_failure_probability': run.failure_probability,
            'rounded_negative_diagnostic_count': rounded[0],
            'rounded_negative_diagnostic_total': rounded[1],
            'failure_sum_error_diagnostic': summation[0],
            'failure_sum_admission': summation[1],
            'per_rounded_entry_absolute_bound': 1e-40}


@lru_cache(maxsize=1)
def outputs():
    stored = prefix()
    analyzers = tuple(execute.complete([c.block(c.pointer_analyzer(phase), 3, (0, 1, 16))], 17, stored.state)
                      for phase in PHASES)
    return tuple((analysis, pointer_outputs(analysis.success)) for analysis in analyzers)


def arrays():
    stored = prefix()
    result = {name+'_failure_weights': run.failure_weights for name, run in
              zip(('tree', 'copy', 'pointer', 'field_first', 'field_second'), stored.runs, strict=True)}
    for index, (run, _) in enumerate(outputs()):
        result[f'analyzer_{index}_failure_weights'] = run.failure_weights
        result[f'analyzer_{index}_full_success_state'] = run.success
    return result


def zero_time_control():
    copied = preparation.runs()[1]
    state = np.zeros((131072, 1), complex); state[:65536] = copied.success
    prepared = execute.complete([c.block(c.pointer_prepare(), 3, (0, 1, 16))], 17, state)
    result = []
    for phase in PHASES:
        analysis = execute.complete([c.block(c.pointer_analyzer(phase), 3, (0, 1, 16))], 17, prepared.success)
        probabilities = [float(np.vdot(x, x).real) for x in pointer_outputs(analysis.success)]
        result.append({'phase_over_pi': str(phase), 'all_six_pointer_probabilities': probabilities,
                       'pointer_source_first_failures': prepared.failure_probability+analysis.failure_probability})
    return result


def free_control():
    """Exact integer bipartite witness, with an independent numerical check."""
    masks = cube.weights(6)
    odd = sum(3 << (2*j) for j, coordinate in enumerate(cube.GEOMETRY.sites)
              if sum(coordinate) % 2)
    signs = np.array([(-1)**(int(mask) & odd).bit_count() for mask in masks], dtype=np.int64)
    gamma = sparse.diags(signs, format='csr', dtype=np.int64)
    doublons = np.array([sum((int(mask) >> (2*j)) & 3 == 3 for j in range(8))
                         for mask in masks], dtype=np.int64)
    ops = cube.operators()
    real = cube.restrict(ops.H.real, 6)-ops.scale*sparse.diags(doublons, format='csr', dtype=np.int64)
    imag = cube.restrict(ops.H.imag, 6)
    lookup = {int(mask): j for j, mask in enumerate(masks)}
    reference = np.zeros(len(masks), dtype=np.int64)
    for mask, value in enumerate(preparation.reference_values()):
        if value:
            reference[lookup[sum(3 << (2*j) for j in range(8) if mask >> j & 1)]] = value
    # On a pair-only input exp(i*pi*n_origin/2) is exactly +/-1.
    kicked = reference*np.array([(-1)**((int(mask) & 3) == 3) for mask in masks], dtype=np.int64)
    for part in (real, imag):
        if (gamma@part@gamma+part).nnz or np.any(part@reference):
            raise ArithmeticError('free bipartite reversal or stationary reference failed')
    if np.any(gamma@reference-reference) or np.any(gamma@kicked-kicked):
        raise ArithmeticError('pair-only probe is not bipartite even')
    h0 = (real+1j*imag)/ops.scale
    value = np.vdot(kicked, expm_multiply(-1j*float(bounds.TIME)*h0, kicked))/56
    if abs(value.imag) > 1e-13:
        raise ArithmeticError('free imaginary response contradicts the exact parity witness')
    return {'domain': 'original cube N6; same R, local F and complete protected repairs',
            'Gamma': '(-1)^(sum of occupation on odd bipartite sites)',
            'exact_Gamma_H0_Gamma_equals_minus_H0': True,
            'exact_H0_R_equals_zero': True, 'exact_Gamma_R_and_FR_equal_R_and_FR': True,
            'free_characteristic_real_for_every_time': True,
            'free_pi_over_two_pointer_zero_probability': '1/2',
            'free_characteristic_at_time_one_diagnostic': [float(value.real), float(value.imag)],
            'additional_native_U0_word_executed': False}


def report():
    stored = prefix(); certificate = bounds.source_certificate()
    prefix_mass = math.fsum(run.failure_probability for run in stored.runs)
    native_mass = math.fsum(r['probability'] for r in stored.records if r['record'][0] == 'native')
    interval = bounds.characteristic_interval(bounds.TIME)
    settings = []
    random = np.random.default_rng(20261101)
    for index, (phase, (analysis, outcome)) in enumerate(zip(PHASES, outputs(), strict=True)):
        target = pointer_outputs(ideal_analyzed(bounds.TIME, phase))
        probabilities = [float(np.vdot(value, value).real) for value in outcome]
        ideal_probabilities = [float(np.vdot(value, value).real) for value in target]
        abort = prefix_mass+native_mass+analysis.failure_probability
        error = math.sqrt(sum(float(np.linalg.norm(a-b)**2) for a, b in zip(outcome, target, strict=True))+abort)
        total = abort+math.fsum(probabilities)
        bound = F(certificate['complete_native_state_errors'][index])
        if error > float(bound)+execute.ADMISSION or abs(total-1) > execute.ADMISSION:
            raise ArithmeticError('complete native Ramsey run exceeds the proved source/field error')
        real_interval = tuple(map(F, interval['real' if index == 0 else 'imag']))
        ideal_interval = (tuple((1+x)/2 for x in real_interval) if index == 0
                          else ((1-real_interval[1])/2, (1-real_interval[0])/2))
        native_interval = max(F(), ideal_interval[0]-bound), min(F(1), ideal_interval[1]+bound)
        if not float(native_interval[0])-execute.ADMISSION <= probabilities[0] <= float(native_interval[1])+execute.ADMISSION:
            raise ArithmeticError('native probability fails independent exact spectral enclosure')
        count = int(random.binomial(1 << 24, probabilities[0]))
        sample_interval = (max(F(), native_interval[0]-F(1, 2048)-F(execute.ADMISSION)),
                           min(F(1), native_interval[1]+F(1, 2048)+F(execute.ADMISSION)))
        observed = F(count, 1 << 24)
        if not sample_interval[0] <= observed <= sample_interval[1]:
            raise ArithmeticError('seeded unconditional simulated frequency exceeds its source-plus-sampling enclosure')
        settings.append({'phase_over_pi': str(phase), 'pointer_records': ['zero', 'one', 'abort2', 'abort3', 'abort4', 'abort5'],
                         'pointer_probabilities': probabilities, 'independent_original_H_probabilities': ideal_probabilities,
                         'rigorous_actual_source_pointer_zero_interval': [str(x) for x in native_interval],
                         'complete_probability': total, 'complete_output_error_diagnostic': error,
                         'analyzer_source_record': source_record(f'analyzer_{index}', analysis),
                         'simulated_counts': {'trials': 1 << 24, 'pointer_zero': count,
                                              'all_other_complete_outcomes': (1 << 24)-count,
                                              'observed_zero_frequency': count/(1 << 24),
                                              'source_plus_sampling_interval': [str(x) for x in sample_interval]}})
    return {'field_duration': str(bounds.TIME), 'field_segments': 2, 'Strang_steps_per_segment': bounds.STEPS,
            'source_prefix_records': [source_record(name, run) for name, run in
                                      zip(('tree', 'copy', 'pointer', 'field_first', 'field_second'), stored.runs, strict=True)],
            'all130_repair_records': list(stored.records), 'terminal_native_repair_probability': native_mass,
            'pointer_settings': settings, 'exact_known_characteristic_interval': interval,
            'executed_zero_field_duration_control': zero_time_control(),
            'exact_contact_free_control': free_control(),
            'sample_convention': 'NumPy PCG64 seed20261101, two independent binomial draws; event is unconditional final pointer zero, complement includes every other complete output',
            'sample_error_per_setting': '1/2048', 'simultaneous_sampling_failure_probability_bound': str(4*F(3, 8)**8),
            'sample_trials_per_setting': 1 << 24, 'sample_count_is_simulated_not_laboratory_data': True,
            'independent_reset_repetitions_are_the_sampling_service': True,
            'resolved_all_eight_frequencies_from_these_two_settings': False,
            'observed_inelastic_cyclic_dimension': 8, 'elastic_reference_dimension': 1,
            'all1820_protected_states_identified_by_this_probe': False,
            'absolute_energy_zero_determined_by_Ramsey': False,
            'measured_frequency_convention': 'omega=E_H-3; H+cI and H-muN produce the same fixedN6 two-arm probabilities',
            'numerical_simulation_admission': execute.ADMISSION,
            'complete_exact_source_errors_are_separate_from_numerical_admission': True}
