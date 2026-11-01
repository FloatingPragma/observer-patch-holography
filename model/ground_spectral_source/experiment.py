"""Actual complete17-bit source execution; no number or spectrum projection."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.cube_eta_protection import cube, repair
from model.interacting_vacuum import charged
from model.observer.source_scalar import compiler
from model.reference_observation.observation import frozen
from model.protected_particle_readout import execute
from model.protected_particle_readout.experiment import source_record
from . import preparation as p, circuits as c, bounds


def repair_layer(state, law):
    records = []
    for seam in repair.SEAMS:
        fields = state.reshape(2, 65536, -1).transpose(1, 0, 2).reshape(65536, -1)
        for label, output in repair.apply(fields, law, seam):
            value = output.reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
            records.append({'law':law, 'seam':[repair.PORTS[j] for j in seam],
                            'record':list(label), 'probability':float(np.vdot(value, value).real)})
            if label == ('outside',): continuation = value
        state = continuation
    return state, records


@dataclass(frozen=True)
class Prefix:
    state: np.ndarray
    runs: tuple
    repair_records: tuple


@lru_cache(maxsize=1)
def prefix():
    tree, frame = p.runs()
    initial = np.zeros((131072, 1), complex); initial[:65536] = frame.success
    pointer = execute.complete(c.pointer_blocks(), 17, initial)
    blocks = c.field_blocks(bounds.TIME/2, bounds.STEPS)
    first = execute.complete(blocks, 17, pointer.success)
    state, record_c = repair_layer(first.success, 'C')
    second = execute.complete(blocks, 17, state)
    state, record_e = repair_layer(second.success, 'E')
    return Prefix(frozen(state), (tree, frame, pointer, first, second), tuple(record_c+record_e))


def repair_output(law, seam_index, record):
    """Every zero/native/outside branch is addressable from its actual prefix."""
    if law not in ('C', 'E') or type(seam_index) is not int or not 0 <= seam_index < len(repair.SEAMS):
        raise ValueError('declared C/E service and seam index required')
    stored = prefix(); state = stored.runs[3 if law == 'C' else 4].success
    for index, seam in enumerate(repair.SEAMS):
        fields = state.reshape(2, 65536, -1).transpose(1, 0, 2).reshape(65536, -1)
        outcomes = repair.apply(fields, law, seam)
        if index == seam_index:
            for label, value in outcomes:
                if label == record:
                    return value.reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
            raise ValueError('no such native record')
        state = outcomes[-1][1].reshape(65536, 2, -1).transpose(1, 0, 2).reshape(state.shape)
    raise ArithmeticError('unreachable repair branch')


def pointer_outputs(state):
    if state.shape != (131072, 1): raise ValueError('complete17-bit state required')
    maps = compiler.transfer_kraus(2, coherent=False)
    if len(maps) != 6 or not np.array_equal(sum(k.conj().T@k for k in maps), np.eye(6)):
        raise ArithmeticError('complete native pointer readout changed')
    return tuple(k[:, :2]@state.reshape(2, 65536) for k in maps)


@lru_cache(maxsize=1)
def outputs():
    stored = prefix()
    analyses = [execute.complete(c.pointer_blocks(phase), 17, stored.state) for phase in c.PHASES]
    return tuple((run, pointer_outputs(run.success)) for run in analyses)


def ideal_analyzed(duration, phase):
    """Independent original-H exponential; invariant sectors only in comparator."""
    if type(duration) is not F or type(phase) is not F or phase not in c.PHASES:
        raise ValueError('exact comparison time and phase required')
    ground = p.ideal_ground(); probe = c.majorana(ground)
    ops = cube.operators(); h = ops.H.matrix(ops.scale)
    first = np.zeros(65536, complex); second = np.zeros(65536, complex)
    for number, state, target in ((4, ground, first), (5, probe, second)):
        masks = cube.weights(number)
        target[masks] = expm_multiply(-1j*float(duration)*h[masks][:, masks], state[masks])
    second = c.majorana(second)*np.exp(1j*math.pi*float(phase))
    return np.concatenate(((first+second)/2, (-first+second)/2))[:, None]


def zero_time_control():
    _, frame = p.runs()
    state = np.zeros((131072, 1), complex); state[:65536] = frame.success
    prepared = execute.complete(c.pointer_blocks(), 17, state)
    rows = []
    for phase in c.PHASES:
        analyzed = execute.complete(c.pointer_blocks(phase), 17, prepared.success)
        probabilities = [float(np.vdot(x, x).real) for x in pointer_outputs(analyzed.success)]
        rows.append({'phase_over_pi':str(phase), 'all_six_pointer_probabilities':probabilities,
                     'pointer_first_failures':prepared.failure_probability+analyzed.failure_probability})
    return rows


def arrays():
    result = {name+'_failure_weights':run.failure_weights for name, run in
              zip(('ground_tree', 'ground_frame', 'charged_pointer', 'field_first', 'field_second'), prefix().runs, strict=True)}
    result['prepared_ground_success_state'] = prefix().runs[1].success
    for j, (run, _) in enumerate(outputs()):
        result[f'analyzer_{j}_failure_weights'] = run.failure_weights
        result[f'analyzer_{j}_full_success_state'] = run.success
    return result


def report():
    stored = prefix(); complete = bounds.certificate(); intervals = bounds.probability_intervals()
    source_abort = math.fsum(run.failure_probability for run in stored.runs)
    repair_abort = math.fsum(row['probability'] for row in stored.repair_records if row['record'][0] == 'native')
    settings = []
    for index, (phase, (run, branches)) in enumerate(zip(c.PHASES, outputs(), strict=True)):
        target = pointer_outputs(ideal_analyzed(bounds.TIME, phase))
        probabilities = [float(np.vdot(x, x).real) for x in branches]
        target_probabilities = [float(np.vdot(x, x).real) for x in target]
        abort = source_abort+repair_abort+run.failure_probability
        total = math.fsum(probabilities)+abort
        error = math.sqrt(math.fsum(float(np.linalg.norm(x-y)**2) for x, y in zip(branches, target, strict=True))+abort)
        proved = F(complete['complete_native_state_and_probability_errors'][index])
        low, high = map(F, intervals[index]['literal_complete_pointer_zero_interval'])
        if (error > float(proved)+execute.ADMISSION or abs(total-1) > execute.ADMISSION or
                not float(low)-execute.ADMISSION <= probabilities[0] <= float(high)+execute.ADMISSION):
            raise ArithmeticError('native ground Ramsey output exceeds exact source/spectrum bound')
        settings.append({'phase_over_pi':str(phase), 'all_six_pointer_probabilities':probabilities,
                         'independent_original_H_probabilities':target_probabilities,
                         'unconditioned_pointer_zero_interval':intervals[index]['literal_complete_pointer_zero_interval'],
                         'first_failure_and_repair_stop_mass':abort, 'complete_probability':total,
                         'complete_output_error_diagnostic':error,
                         'analyzer_source_records':source_record(f'analyzer_{index}', run)})
    tree, frame = stored.runs[:2]
    prep_error = math.sqrt(float(np.linalg.norm(frame.success[:, 0]-p.ideal_ground())**2)+tree.failure_probability+frame.failure_probability)
    protection = charged.repair_certificate()
    if len(stored.repair_records) != 130 or protection['complete_outcomes_including_zero_and_outside'] != 130:
        raise ArithmeticError('complete protected repair catalogue changed')
    return {'ground_preparation_executed':True, 'full17bit_charged_Ramsey_executed':True,
            'input_initialization':'16 field qubits and1 pointer blank; tree executes on8 modes with8 exact blank spectators',
            'native_qubit_parity_phase_reference_supplied':True,
            'source_prefix_records':[source_record(name, run) for name, run in
                                     zip(('ground_tree', 'ground_frame', 'charged_pointer', 'field_first', 'field_second'), stored.runs, strict=True)],
            'all130_repair_outcome_slots':list(stored.repair_records),
            'repair_native_first_stop_mass':repair_abort,
            'ideal_protected_domain':'N4 directsum N5; all native repair branches vanish, outside branch is identity',
            'complete_ground_preparation_error_diagnostic':prep_error,
            'zero_time_native_control':zero_time_control(), 'settings':settings,
            'host_numerical_admission':execute.ADMISSION,
            'host_diagnostics_admit_literal_control_theorem':False,
            'success_or_charge_postselection_used':False,
            'all78_lines_native_resolved':False,
            'ground_selected_by_passive_scheduler':False,
            'physical_species_or_relativistic_mass_identified':False}
