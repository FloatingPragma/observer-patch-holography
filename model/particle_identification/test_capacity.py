"""Independent complete-carrier checks of the conditional local error bound."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from ..full_fock_field.field import FieldGeometry, FullFock
from ..full_fock_field.local_simulation import local_observer_schedule, split_bound
from .capacity import (certificate, coefficient_counterexample, duhamel_action_bound,
                       history_error_budget, local_history_budget)


def _action_bound(sides, time, speed, contact, other_time, other_speed, other_contact):
    return duhamel_action_bound(sides, time, speed, contact, simulated_duration=other_time,
                                simulated_speed=other_speed, simulated_interaction=other_contact)


@pytest.mark.parametrize('sides', [(1, 1, 1), (2, 1, 1), (1, 2, 1)])
def test_actual_full_car_generator_and_unitary_differences_obey_term_count_bound(sides):
    space = FullFock(2*np.prod(sides).item())
    first = FieldGeometry(sides, speed=.125, interaction=1.)
    second = FieldGeometry(sides, speed=.2, interaction=.75)
    a = space.hamiltonian(first).toarray()/3
    b = space.hamiltonian(second).toarray()*(-.4)
    bound = _action_bound(sides, F(1, 3), F(1, 8), 1, F(-2, 5), F(1, 5), F(3, 4))
    assert np.linalg.norm(a-b, 2) <= float(bound)+1e-14
    assert np.linalg.norm(expm(-1j*a)-expm(-1j*b), 2) <= float(bound)+1e-14


def test_clock_rescaling_cancels_exactly_and_zero_evolution_needs_no_action_error():
    assert _action_bound((10000, 5000, 2), F(1, 8), F(1, 8), 1,
                         F(1, 24), F(3, 8), 3) == 0
    assert _action_bound((2, 1, 1), 0, 1, 7, 0, 100, 3) == 0
    assert _action_bound((2, 1, 1), 1, 1, 2, 1, 3, 5) == 8


def _trace_norm(matrix):
    return np.abs(np.linalg.eigvalsh(matrix)).sum()


def _adaptive_history(initial, unitaries, flip):
    """Keep all two-bit branches and quantum field/reference output blocks."""
    branches = {(): initial}
    for step, unitary in enumerate(unitaries):
        output = {}
        for history, state in branches.items():
            value = unitary@state@unitary.conj().T
            # The second instrument depends on the first retained record.
            mode = history[-1] if step else 0
            strength = .7 if not history or history[-1] == 0 else .8
            occupancy = (np.arange(16) >> mode) & 1
            p0 = np.where(occupancy, strength, 1-strength)
            operators = [np.kron(np.diag(np.sqrt(p)), np.eye(2)) for p in (p0, 1-p0)]
            raw = [operator@value@operator for operator in operators]
            for record in (0, 1):
                output[history+(record,)] = (1-flip)*raw[record]+flip*raw[1-record]
        branches = output
    return branches


def test_complete_adaptive_history_includes_initial_correlation_and_instrument_errors():
    space = FullFock(4)
    h = space.hamiltonian(FieldGeometry((2, 1, 1))).toarray()
    other = space.hamiltonian(FieldGeometry((2, 1, 1), speed=.126, interaction=.999)).toarray()
    times, other_times = (F(1, 5), F(-1, 10)), (F(2001, 10000), F(-1001, 10000))
    target_units = [np.kron(expm(-1j*float(t)*h), np.eye(2)) for t in times]
    other_units = [np.kron(expm(-1j*float(t)*other), np.eye(2)) for t in other_times]
    vector = np.zeros(32, complex)
    vector[0], vector[3*2+1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    initial = np.outer(vector, vector.conj())
    initial_error, flip = F(1, 50), F(1, 100)
    perturbed = (1-float(initial_error))*initial+float(initial_error)*np.eye(32)/32
    target = _adaptive_history(initial, target_units, 0.)
    actual = _adaptive_history(perturbed, other_units, float(flip))
    assert set(actual) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    assert abs(sum(np.trace(block) for block in actual.values())-1) < 1e-13
    assert abs(sum(np.trace(block) for block in target.values())-1) < 1e-13
    quantum_distance = sum(_trace_norm(actual[key]-target[key]) for key in actual)/2
    record_distance = sum(abs(np.trace(actual[key]-target[key])) for key in actual)/2
    actions = [_action_bound((2, 1, 1), t, F(1, 8), 1, s, F(63, 500), F(999, 1000))
               for t, s in zip(times, other_times, strict=True)]
    budget = history_error_budget(initial_error, locality_errors=(0, 0), schedule_errors=(0, 0),
                                  action_errors=actions, instrument_errors=(flip, flip))
    assert 0 < record_distance <= quantum_distance <= float(budget['complete_record_error'])


def test_existing_nested_schedule_reused_and_actual_simulator_split_error_recomputed():
    times = (F(1, 8), F(-1, 16))
    plan = local_observer_schedule(times, F(1, 100), refined=True)
    budget = local_history_budget(times, F(1, 100), refined=True)
    assert budget['initial_restriction_radius'] == plan['initial_restriction_radius']
    assert budget['complete_record_error'] == plan['half_diamond_bound']
    first, last = budget['stages']
    assert first['input_radius'] > first['output_radius'] == last['input_radius']
    altered = local_history_budget((F(1, 8),), F(1, 100), simulated_durations=(F(1, 4),),
                                    simulated_speed=F(1, 4), simulated_interaction=2, refined=True)
    row, = altered['stages']
    assert row['schedule_operator_error'] == split_bound(
        row['sides'], F(1, 4), row['steps'], F(1, 4), 2, refined=True)
    assert row['duhamel_error'] > 0
    assert altered['complete_record_error'] == 1  # Planning epsilon is not an output guarantee.


def test_zero_hopping_zero_duration_and_error_saturation():
    budget = local_history_budget((0,), F(1, 100), speed=0, initial_state_error=F(1, 4),
                                  instrument_errors=(F(1, 3),))
    assert budget['initial_restriction_radius'] == 0
    assert budget['complete_record_error'] == F(7, 12)
    huge = history_error_budget(1, locality_errors=(2,), schedule_errors=(10,),
                               action_errors=(4,), instrument_errors=(1,))
    assert huge['uncapped_error'] == 17 and huge['complete_record_error'] == 1


def test_neutral_reference_counterexample_changes_records_without_any_boundary_error():
    # Low two modes are the field site; upper two are its declared reference.
    # Pair transfer retains total number two and has even CAR parity.
    space = FullFock(4)
    h = (space.number(0)@space.number(1)).toarray()
    number = sum(space.number(j).toarray() for j in range(4))
    state = np.zeros(16, complex)
    state[3] = state[12] = 1/np.sqrt(2)
    plus = np.outer(state, state.conj())
    assert np.linalg.norm(number@state-2*state) == 0
    assert np.linalg.norm(number@plus-plus@number) == 0
    probability = abs(np.vdot(state, expm(-1j*h)@state))**2
    example = coefficient_counterexample()
    difference = 1-probability
    assert float(example['probability_difference_lower']) < difference < float(example['probability_difference_upper'])
    for spectator_energy in (0., .7, 10.):
        # Spectators may change globally; no term couples them to the detector.
        total = np.kron(h, np.eye(2))+np.kron(np.eye(16), np.diag([0., spectator_energy]))
        enlarged = np.kron(state, np.array([1., 1.])/np.sqrt(2))
        final = expm(-1j*total)@enlarged
        assert abs(np.vdot(final, np.kron(plus, np.eye(2))@final)-probability) < 1e-13


def test_exact_certificate_keeps_physical_identifications_open():
    result = certificate()
    assert F(result['example']['complete_record_error']) < F(1, 100)
    assert result['example']['initial_restriction_radius'] == 4
    assert not result['N_to_coefficients_clock_cutoffs_or_state_is_derived']
    assert not result['observed_particle_mass_or_mass_ratio_is_derived']


@pytest.mark.parametrize('bad', [True, .25, float('nan'), complex(1), '1/4'])
def test_inexact_or_invalid_input_does_not_masquerade_as_an_exact_budget(bad):
    with pytest.raises(ValueError):
        _action_bound((1, 1, 1), bad, 0, 1, 1, 0, 1)
    with pytest.raises(ValueError):
        history_error_budget(bad, locality_errors=(0,), schedule_errors=(0,), action_errors=(0,))


def test_wrong_stage_counts_and_out_of_range_errors_are_rejected():
    with pytest.raises(ValueError):
        history_error_budget(0, locality_errors=(0, 0), schedule_errors=(0,), action_errors=(0, 0))
    with pytest.raises(ValueError):
        history_error_budget(0, locality_errors=(3,), schedule_errors=(0,), action_errors=(0,))
    with pytest.raises(ValueError):
        history_error_budget(0, locality_errors=(0,), schedule_errors=(-1,), action_errors=(0,))
    with pytest.raises(ValueError):
        local_history_budget((0,), F(1, 10), simulated_durations=())
    with pytest.raises(ValueError):
        local_history_budget((0,), F(1, 10), simulated_speed=-1)
