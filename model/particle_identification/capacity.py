"""Exact local-history error budgets, without assuming a physical N law.

The action is the existing full-Fock Weyl/contact action at fixed spacing.
Preparation and instrument errors are hypotheses in trace/half-diamond norm;
this module neither estimates them from samples nor chooses a vacuum.
"""
from fractions import Fraction as F
from math import factorial

from ..full_fock_field.bounds import stability_bounds
from ..full_fock_field.local_simulation import (
    local_observer_schedule, source_schedule_counts, split_bound,
)


# Roots consumed by this mathematical adapter; receipt producers should also
# bind their ordinary project import closure, including package initializers.
SOURCE_FILES = (
    'model/particle_identification/capacity.py',
    'model/particle_identification/test_capacity.py',
    'model/particle_identification/CAPACITY.md',
    'model/full_fock_field/bounds.py',
    'model/full_fock_field/local_simulation.py',
    'model/full_fock_field/field.py',
    'model/full_fock_field/LOCAL_SIMULATION.md',
)


def _exact(value, name, *, lower=None, upper=None):
    if type(value) not in (int, F):
        raise ValueError(f'{name} must be an integer or Fraction')
    value = F(value)
    if ((lower is not None and value < lower)
            or (upper is not None and value > upper)):
        raise ValueError(f'{name} lies outside its admitted interval')
    return value


def _sequence(values, name, *, length=None, lower=None, upper=None):
    if type(values) not in (tuple, list) or not values:
        raise ValueError(f'{name} must be a nonempty list or tuple')
    if length is not None and len(values) != length:
        raise ValueError(f'{name} has the wrong number of stages')
    return tuple(_exact(x, name, lower=lower, upper=upper) for x in values)


def duhamel_action_bound(sides, duration, speed, interaction, *,
                         simulated_duration, simulated_speed, simulated_interaction):
    """Bound half diamond distance of two exact finite-box evolutions.

Each spatial edge has two hopping terms of norm c/2, and each site
has a contact projector of norm one. The bound compares tH and sH',
including the clock convention, and holds on the full Fock algebra.
"""
    geometry = stability_bounds(sides, speed, interaction)
    stability_bounds(sides, simulated_speed, simulated_interaction)
    t = _exact(duration, 'target duration')
    s = _exact(simulated_duration, 'simulated duration')
    return (geometry['spatial_edges']*abs(t*speed-s*simulated_speed)
            + geometry['sites']*abs(t*interaction-s*simulated_interaction))


def history_error_budget(initial_state_error, *, locality_errors, schedule_errors,
                         action_errors, instrument_errors=None):
    """Compose uniform errors for all retained outcomes and quantum memory.

Locality entries are operator/cb errors (hence divided by two). Schedule
entries bound unitary operator errors; action entries are Duhamel bounds.
Instrument entries are half-diamond bounds uniform on every incoming
state/reference and every earlier record, not average branch errors.
"""
    initial = _exact(initial_state_error, 'initial trace distance', lower=0, upper=1)
    locality = _sequence(locality_errors, 'locality errors', lower=0, upper=2)
    count = len(locality)
    schedule = _sequence(schedule_errors, 'schedule errors', length=count, lower=0)
    action = _sequence(action_errors, 'action errors', length=count, lower=0)
    instruments = ((F(),)*count if instrument_errors is None else
                   _sequence(instrument_errors, 'instrument errors', length=count,
                             lower=0, upper=1))
    stages = tuple({'locality_operator_error': d, 'schedule_operator_error': e,
                    'duhamel_error': a, 'instrument_half_diamond_error': j,
                    'stage_error': d/2+e+a+j}
                   for d, e, a, j in zip(locality, schedule, action, instruments, strict=True))
    total = initial+sum((row['stage_error'] for row in stages), F())
    return {'initial_joint_trace_distance': initial, 'stages': stages,
            'uncapped_error': total, 'complete_record_error': min(F(1), total)}


def local_history_budget(durations, epsilon, *, speed=F(1, 8), interaction=F(1),
                         simulated_durations=None, simulated_speed=None,
                         simulated_interaction=None, initial_state_error=F(),
                         instrument_errors=None, support_radius=0, support_sites=1,
                         max_radius=4096, refined=False):
    """Reuse the existing nested buffers, then charge actual simulation errors.

The target law chooses the locality buffers and original split counts.
The simulator uses the same counts, with its own coefficients and times;
its split error is recomputed. Consequently epsilon is a planning budget,
not a claim that preparation or mismatched-law errors fit that budget.
"""
    plan = local_observer_schedule(durations, epsilon, speed=speed,
                                   interaction=interaction, support_radius=support_radius,
                                   support_sites=support_sites, max_radius=max_radius,
                                   refined=refined)
    target_times = tuple(row['duration'] for row in plan['stages'])
    simulated_times = (target_times if simulated_durations is None else
                       _sequence(simulated_durations, 'simulated durations',
                                 length=len(target_times)))
    c = speed if simulated_speed is None else simulated_speed
    u = interaction if simulated_interaction is None else simulated_interaction
    stability_bounds((1, 1, 1), c, u)
    locality, schedule, action, details = [], [], [], []
    for row, simulated_time in zip(plan['stages'], simulated_times, strict=True):
        sides = (2*row['input_radius']+1,)*3
        locality.append(row['boundary_operator_bound'])
        schedule.append(split_bound(sides, simulated_time, row['steps'], c, u,
                                    refined=row['coefficient_method'] == 'ordered overlaps'))
        action.append(duhamel_action_bound(
            sides, row['duration'], speed, interaction,
            simulated_duration=simulated_time, simulated_speed=c, simulated_interaction=u))
        details.append({'sides': sides, 'input_radius': row['input_radius'],
                        'output_radius': row['output_radius'], 'steps': row['steps'],
                        'target_duration': row['duration'], 'simulated_duration': simulated_time,
                        'resources': source_schedule_counts(sides, row['steps'], c)})
    budget = history_error_budget(initial_state_error, locality_errors=locality,
                                  schedule_errors=schedule, action_errors=action,
                                  instrument_errors=instrument_errors)
    budget['stages'] = tuple(dict(row, **detail)
                             for row, detail in zip(budget['stages'], details, strict=True))
    budget.update({'initial_restriction_radius': plan['initial_restriction_radius'],
                   'planning_epsilon': F(epsilon), 'target_speed': F(speed),
                   'target_interaction': F(interaction), 'simulated_speed': F(c),
                   'simulated_interaction': F(u), 'ambient_volume_required': False,
                   'initial_joint_state_error_is_supplied': True,
                   'instrument_errors_are_uniform_supplied_bounds': True,
                   'physical_capacity_or_spectral_limit_derived': False})
    return budget


def coefficient_counterexample():
    """A neutral pair-transfer observation separates two onsite laws.

Use a declared two-mode reference and the total-number-two state
    (|pair,empty_ref> + |empty,pair_ref>)/sqrt(2).
The same-state projector commutes with total number. For target c=U=0
and alternate c=0,U=1 at t=1, its probabilities differ by (1-cos(1))/2.
Alternating Taylor bounds give the exact rational interval below.
"""
    cosine_lower = sum((F((-1)**n, factorial(2*n)) for n in range(4)), F())
    cosine_upper = cosine_lower+F(1, factorial(8))
    return {'target_speed': F(), 'target_interaction': F(),
            'alternate_speed': F(), 'alternate_interaction': F(1), 'duration': F(1),
            'probability_difference_lower': (1-cosine_upper)/2,
            'probability_difference_upper': (1-cosine_lower)/2,
            'reference': 'declared compensating pair reference with zero Hamiltonian',
            'total_number': 2, 'all_ambient_volumes_have_same_local_law': True,
            'both_laws_are_asserted_admissible_OPH_universes': False}


def certificate():
    """A compact exact budget; all unproduced hypotheses remain explicit."""
    budget = local_history_budget(
        (F(1, 8),), F(1, 100), refined=True, initial_state_error=F(1, 1000),
        simulated_speed=F(1, 8)+F(1, 10**9),
        simulated_interaction=1-F(1, 10**9), instrument_errors=(F(1, 10000),))
    if not budget['complete_record_error'] < F(1, 100):
        raise ArithmeticError('the finite local example exceeds its declared total budget')

    def serial(value):
        if isinstance(value, F):
            return str(value)
        if isinstance(value, dict):
            return {key: serial(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [serial(item) for item in value]
        return value

    return serial({'scope': 'fixed-spacing full-Fock Weyl/contact local complete histories',
                   'example': budget, 'coefficient_counterexample': coefficient_counterexample(),
                   'N_to_coefficients_clock_cutoffs_or_state_is_derived': False,
                   'observed_particle_mass_or_mass_ratio_is_derived': False})
