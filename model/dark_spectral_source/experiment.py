"""Actual complete native preparation, field wait, repair and Ramsey export."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..common_spin_protection import spin
from ..certified_reference.source import complete_export_certificate
from ..full_fock_field.local_simulation import split_bound
from ..observer.source_scalar import compiler
from . import circuits as c, execute, exact


TIME = F(1, 2)
STEPS = 2


def pauli_action(state, label):
    addresses = np.arange(len(state), dtype=np.int64)
    flip = sum(1 << bit for bit, letter in enumerate(label) if letter in 'XY')
    sign = sum(1 << bit for bit, letter in enumerate(label) if letter in 'YZ')
    phase = (1j**label.count('Y'))*np.array([(-1)**(int(x) & sign).bit_count() for x in addresses])
    output = np.empty_like(state)
    output[addresses ^ flip] = phase[:, None]*state
    return output


def ideal_plan(state, plan):
    out = np.asarray(state, complex).copy()
    for label, fraction in plan:
        theta = math.pi*float(fraction)
        out = math.cos(theta)*out+1j*math.sin(theta)*pauli_action(out, label)
    return out


def pointer_rotation(state, kind, fraction):
    return ideal_plan(state, [(c.label((c.POINTER, kind), bits=c.BITS), fraction)])


def ideal_analyzer(state, phase):
    out = ideal_plan(state, c.inverse(c.controlled(c.preparation_plan())))
    out *= np.exp(1j*float(phase)*math.pi/2)
    out = pointer_rotation(out, 'Z', -phase/2)
    return pointer_rotation(out, 'Y', F(1, 4))


def ideal_wait(time):
    dark = expm_multiply(-1j*float(time)*spin.operators().H.matrix(spin.operators().scale), exact.target())
    out = np.zeros((8192, 1), complex)
    out[0, 0] = 1/math.sqrt(2)
    out[4096:, 0] = dark/math.sqrt(2)
    return out


def repair(state, law, seam):
    retained = []
    fields = state.reshape(2, 4096, -1).transpose(1, 0, 2).reshape(4096, -1)
    for label, operator in spin.repair_outcomes('hole', law, seam):
        output = (operator@fields).reshape(4096, 2, -1).transpose(1, 0, 2).reshape(state.shape)
        if label == ('outside',):
            continuation = output
        else:
            retained.append({'law': law, 'seam': list(seam), 'record': list(label),
                             'probability': float(np.vdot(output, output).real)})
    return continuation, retained


def source_record(stage, run):
    return {'stage': stage, 'pulses': len(run.word),
            'pair_code_event_indices': [j for j, (_, _, targets) in enumerate(run.word) if len(targets) == 2],
            'source_abort_slots': sum(6-(1 << len(targets)) for _, _, targets in run.word),
            'padding_slots_in_four_column_table': sum(2 for _, _, targets in run.word if len(targets) == 2),
            'all_first_export_failure_weights': run.failure_weights.tolist(),
            'first_export_failure_probability': run.failure_probability}


@dataclass(frozen=True)
class Prefix:
    state: np.ndarray
    runs: tuple
    repair_records: tuple


@lru_cache(maxsize=1)
def prefix():
    initial = np.zeros((8192, 1), complex); initial[0, 0] = 1
    prepared = execute.complete(c.preparation_word(), c.BITS, initial)
    word = c.field_word(TIME/2, STEPS)
    first = execute.complete(word, c.BITS, prepared.success)
    state, first_records = repair(first.success, 'E', (0, 1))
    second = execute.complete(word, c.BITS, state)
    state, second_records = repair(second.success, 'C', (1, 2))
    state, third_records = repair(state, 'E', (0, 2))
    return Prefix(state, (prepared, first, second), tuple(first_records+second_records+third_records))


def source_certificate():
    complete_export_certificate()
    preparation = c.certify_word(c.preparation_word(), c.preparation_template())
    segment = c.certify_word(c.field_word(TIME/2, STEPS), c.field_template(TIME/2, STEPS))
    analyzers = [c.certify_word(c.analyzer_word(phase), c.analyzer_template(phase)) for phase in (F(), F(1, 2))]
    split = 2*split_bound(spin.GEOMETRY.sides, TIME/2, STEPS, refined=True)
    source_totals = [preparation['complete_source_operator_error']+2*segment['complete_source_operator_error']+
                     analyzer['complete_source_operator_error'] for analyzer in analyzers]
    totals = [value+split for value in source_totals]
    if max(totals) >= F(1, 100):
        raise ArithmeticError('the complete source/field approximation does not resolve the declared fringe')
    return {'preparation': preparation, 'each_field_segment': segment, 'analyzers': analyzers,
            'Hamiltonian_duration': TIME, 'field_segments': 2, 'Strang_steps_each_segment': STEPS,
            'exact_two_segment_Strang_error_bound': split,
            'literal_source_only_error_bounds': source_totals,
            'complete_output_error_bounds': totals,
            'source_error_bounds_apply_to_all_inputs_and_untouched_references': True,
            'all_export_abort_complements_retained': True,
            'physical_device_precision_certified': False}


def report():
    stored = prefix()
    source = [source_record(name, run) for name, run in zip(('preparation', 'field_first', 'field_second'), stored.runs, strict=True)]
    prefix_fail = math.fsum(record['first_export_failure_probability'] for record in source)
    repair_fail = math.fsum(record['probability'] for record in stored.repair_records)
    ideal = ideal_wait(TIME)
    certificate = source_certificate()
    settings = []
    for index, phase in enumerate((F(), F(1, 2))):
        analysis = execute.complete(c.analyzer_word(phase), c.BITS, stored.state)
        outputs = execute.pointer_outputs(analysis.success)
        probabilities = [float(np.vdot(value, value).real) for value in outputs]
        target = ideal_analyzer(ideal, phase)
        target_outputs = execute.pointer_outputs(target)
        earlier = prefix_fail+repair_fail+analysis.failure_probability
        error = math.sqrt(sum(float(np.vdot(a-b, a-b).real)
                              for a, b in zip(outputs, target_outputs, strict=True))+earlier)
        total = earlier+math.fsum(probabilities)
        analytic = (1+(np.exp(1j*math.pi*float(phase))*exact.amplitude(float(TIME))).real)/2
        if (abs(total-1) > execute.ADMISSION or error > float(certificate['complete_output_error_bounds'][index])
                or abs(probabilities[0]-analytic) > error+execute.ADMISSION):
            raise ArithmeticError('complete native Ramsey history lost probability or its spectral target')
        settings.append({'phase_over_pi': str(phase), 'analyzer_source': source_record('analyzer', analysis),
                         'pointer_records': ['zero', 'one', 'abort2', 'abort3', 'abort4', 'abort5'],
                         'pointer_probabilities': probabilities, 'analytic_pointer_zero_probability': analytic,
                         'complete_probability': total, 'complete_output_reference_error_diagnostic': error,
                         'complete_field_outputs_retained_without_postselection': True})
    prepared = stored.runs[0]
    expected = ideal_wait(F())
    prep_error = math.sqrt(float(np.linalg.norm(prepared.success-expected)**2)+prepared.failure_probability)
    if prep_error > float(certificate['preparation']['complete_source_operator_error']):
        raise ArithmeticError('actual native preparation disagrees with the sparse proved target')
    return {'same_source_and_field': 'six-site c=1/8,U=1 full H; original C/E hole repair plus Q',
            'preparation_complete_output_error_diagnostic': prep_error,
            'common_prefix_source_records': source,
            'all40_native_repair_records': list(stored.repair_records),
            'repair_native_probability': repair_fail,
            'pointer_settings': settings,
            'nominal_ideal_characteristic_real': exact.amplitude(float(TIME)).real,
            'nominal_ideal_characteristic_imag': exact.amplitude(float(TIME)).imag,
            'success_state_renormalized': False,
            'spectral_fit_or_postselection_used': False,
            'coherent_number_phase_reference_supplied': True,
            'physical_charge_superselection_evaded_or_derived': False,
            'host_numerical_admission': execute.ADMISSION}


def json_rationals(value):
    if isinstance(value, dict):
        return {key: json_rationals(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_rationals(item) for item in value]
    return str(value) if type(value) is F else value
