"""Retain the reduced full-field schedule and a complete reference control."""
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from . import budget
from ..band_clock import clock
from ..band_clock.experiment import retained_counts, source_bindings as parent_bindings

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def source_bindings():
    paths = set(parent_bindings())
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py'))
    paths.update(('model/band_clock_budget/README.md',
                  'model/pair_cyclic_control/spectral_bound.py',
                  'model/pair_cyclic_control/__init__.py',
                  'model/pair_cyclic_control/README.md',
                  'model/band_clock/results.json'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def run():
    bindings = source_bindings()
    space = clock.ClockField()
    state, _, _, _ = clock.packet_with_reference()
    counts = retained_counts()
    steps = budget.minimum_steps()
    decision, outputs = budget.execute(counts, state, space=space)
    if decision['status'] != 'accept':
        raise ArithmeticError('retained count decoder rejected the fixed rate window')
    reference = clock.pointers(space, clock.exact_state(space, state, np.log(25/8)))
    probabilities = [float(np.vdot(v, v).real) for v in outputs.values()]
    distance = clock.cq_distance(outputs.values(), reference)
    reference_error = float(np.linalg.norm(
        sum(v.conj().T@v for v in outputs.values())-state.conj().T@state))
    certificate = budget.error_budget(steps=steps)
    if distance > float(F(certificate['complete_average_halfdiamond_upper'])):
        raise ArithmeticError('complete source output violates the certified ceiling')
    if abs(sum(probabilities)-1) > 2e-9 or reference_error > 2e-9:
        raise ArithmeticError('source evolution loses complete mass or reference')
    duration = decision['Hamiltonian_model_duration']
    before = budget.resources(space, duration, 32)
    after = budget.resources(space, duration, steps)
    result = {
        'schema': 'oph.band_clock_budget.v1',
        'source_bindings': bindings,
        'full_N2_dimension': space.dimension,
        'retained_band_momenta': 64,
        'previous_steps': 32,
        'selected_steps': steps,
        'previous_uniform_error_budget': clock.error_budget(),
        'improved_uniform_error_budget': certificate,
        'one_fewer_step_bound': budget.error_budget(steps=steps-1),
        'minimal_within_this_bound_not_all_programs': True,
        'source_resources_before': before,
        'source_resources_after': after,
        'supplied_pulses_saved': before['native_pulses']-after['native_pulses'],
        'source_execution': {
            'complete_pointer_probabilities': probabilities,
            'complete_CQ_half_trace_vs_true_rate': distance,
            'untouched_reference_marginal_residual': reference_error,
            'all_four_pointer_outputs_retained': True,
            'successful_output_normalization': False,
        },
        'scope': {
            'source_law_action_reference_or_detector_changed': False,
            'full_spatial_field_replaced_by_cyclic_control': False,
            'record_acquisition_reexecuted': False,
            'packet_and_global_Gibbs_preparation_cost_included': False,
            'physical_clock_mass_or_continuum_claimed': False,
            'floating_or_hardware_error_included_in_ideal_bound': False,
        },
    }
    if bindings != source_bindings():
        raise ArithmeticError('source changed during complete schedule replay')
    return result
