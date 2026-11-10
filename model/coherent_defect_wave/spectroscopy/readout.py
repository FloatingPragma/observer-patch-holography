"""Exact four-bin phase instruments for the owned joint unitary.

The clock is an optional coherent apparatus, with its own finite register.
Every output is retained. Its phase labels do not define physical energies.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.model import JointObserver, SparseState, combine
from model.coherent_defect_wave.gauss import gauss_measurement, gamma_amplitudes
from model.spin_shell_readout.bounds import clock_plan


def _real_norm(state):
    return sum((F(value)**2 for _, value in state.amplitudes), F(0))/state.norm_denominator


def _sum(states, coefficients):
    if not states or len(states) != len(coefficients):
        raise ValueError('one coefficient per nonempty source-power list is required')
    denominator = states[0].norm_denominator
    output = {}
    for state, weight in zip(states, coefficients):
        if state.norm_denominator != denominator:
            raise ValueError('source powers must have a common normalization')
        for key, value in state.amplitudes:
            if not isinstance(value, F):
                raise ValueError('this exact execution requires rational real components')
            output[key] = output.get(key, F(0)) + weight*value
    return SparseState(output, denominator)


@dataclass(frozen=True)
class GaussianState:
    """Coherent real plus i times imaginary parts, sharing one normalization."""
    real: SparseState
    imag: SparseState

    def __post_init__(self):
        if self.real.norm_denominator != self.imag.norm_denominator:
            raise ValueError('real and imaginary components must share normalization')
        if any(not isinstance(value, F) for part in (self.real, self.imag)
               for _, value in part.amplitudes):
            raise ValueError('rational components are required')

    def norm_squared(self):
        return _real_norm(self.real)+_real_norm(self.imag)

    def gauss_positive(self):
        return all(gamma_amplitudes(part.as_dict()) == part.as_dict()
                   for part in (self.real, self.imag))

    def as_sparse(self):
        out = self.real.as_dict()
        for key, value in self.imag.amplitudes:
            out[key] = out.get(key, F(0))+s.I*value
        return SparseState(out, self.real.norm_denominator)


def powers(observer, state, maximum=5):
    if type(maximum) is not int or not 0 <= maximum <= 16:
        raise ValueError('bounded exact source power count required')
    observer._validate(state)
    rows = [state]
    for _ in range(maximum):
        rows.append(observer.step(rows[-1]))
    return tuple(rows)


def fourier_branches(source_powers):
    """K_l=(1/4) sum_(t=0)^3 i^(-lt) U^t, all four outcomes."""
    if len(source_powers) != 4:
        raise ValueError('four exact source powers required')
    units = ((F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1)))
    rows = []
    for label in range(4):
        phases = [units[(-label*t) % 4] for t in range(4)]
        real = _sum(source_powers, [a/4 for a, _ in phases])
        imag = _sum(source_powers, [b/4 for _, b in phases])
        rows.append(({'kind': 'joint_phase_clock', 'bin': label,
                      'center_in_quarter_turns': label}, GaussianState(real, imag)))
    return tuple(rows)


def filter_clock_branches(source_powers):
    """Retain both {(I+U²)/2,(I-U²)/2} branches and every clock outcome."""
    if len(source_powers) != 6:
        raise ValueError('source powers zero through five are required')
    rows = []
    for sign, label in ((1, 'even_filter'), (-1, 'nonflat_filter')):
        filtered = tuple(combine(source_powers[t], source_powers[t+2], F(1, 2), F(sign, 2))
                         for t in range(4))
        for clock, state in fourier_branches(filtered):
            rows.append(({'filter': label, **clock}, state))
    return tuple(rows)


def overlap(left, right):
    if left.norm_denominator != right.norm_denominator:
        raise ValueError('common source normalization required')
    first = left.as_dict()
    return sum((first.get(key, F(0))*value for key, value in right.amplitudes), F(0))/left.norm_denominator


@lru_cache(maxsize=1)
def source_experiment():
    observer = JointObserver()
    initial = observer.embedded_port(observer.initial_word, 0, internal=3)
    preparations = gauss_measurement(observer, initial)
    positive = preparations[0][1]
    failure = preparations[1][1]
    source_powers = powers(observer, positive)
    success = _real_norm(positive)
    if success != F(1, 2) or _real_norm(failure) != F(1, 2):
        raise ArithmeticError('complete Gauss preparation changed')
    moments = tuple(overlap(positive, evolved)/success for evolved in source_powers)
    expected = (F(1), F(0), F(55, 72), F(-9571, 54000),
                F(279031, 540000), F(-19868809, 48600000))
    if moments != expected:
        raise ArithmeticError('literal joint clock moments changed')
    if any(_real_norm(state) != success for state in source_powers):
        raise ArithmeticError('joint powers lost probability')
    if any(gamma_amplitudes(state.as_dict()) != state.as_dict() for state in source_powers):
        raise ArithmeticError('joint powers left the positive Gauss sector')
    branches = filter_clock_branches(source_powers)
    if any(not state.gauss_positive() for _, state in branches):
        raise ArithmeticError('a clock branch left the positive constraint')
    probabilities = [state.norm_squared() for _, state in branches]
    if sum(probabilities) + _real_norm(failure) != 1:
        raise ArithmeticError('complete preparation/filter/clock experiment lost an outcome')
    active = sum(p for (label, _), p in zip(branches, probabilities)
                 if label['filter'] == 'nonflat_filter')
    if active != F(17, 288):
        raise ArithmeticError('populated nonflat filter weight changed')
    return observer, initial, preparations, source_powers, branches, moments


def energy_branch_control():
    """Different positive generators and gaps give exactly the same integer ticks."""
    U = s.diag(1, -s.I)
    generators = (s.diag(0, s.pi/2), s.diag(0, 5*s.pi/2))
    if any((-s.I*H).exp() != U for H in generators):
        raise ArithmeticError('energy-branch counterexample changed')
    return {'clock_matrix': [['1', '0'], ['0', '-i']],
            'positive_generator_eigenvalues_in_pi_units': [['0', '1/2'], ['0', '5/2']],
            'same_ground_vector': True,
            'different_positive_gap_ratio': '5',
            'same_integer_tick_unitary_and_all_its_instrument_histories': True,
            'identifies_a_branch_for_the_joint_source': False}


def numerical_free_phase_bound(moments):
    """Nonzero phase support separated from both exact flat eigenvalues."""
    active = (1-moments[2])/2
    threshold = F(1, 16)
    original = max(F(0), (active-threshold)/(1-threshold))
    weighted = max(F(0), active-threshold)
    return {'phase_event': 'abs(sin(theta)) >= 1/4',
            'phase_distance_from_each_of_0_and_pi_lower': 'arcsin(1/4)',
            'spectral_weight_lower_given_positive_Gauss_preparation': str(original),
            'nonflat_filtered_spectral_weight_lower_given_positive_Gauss_preparation': str(weighted),
            'original_trial_weight_lower_after_Gauss_and_nonflat_filter': str(weighted/2),
            'conditional_nonflat_spectral_weight_lower': str(weighted/active),
            'isolated_phase_eigenvalue_or_rest_mass': False}


def resolution_plan():
    plan = clock_plan(F(1, 64), F(1, 16))
    plan['maximum_filter_plus_clock_U_power'] = plan['clock_size']+1
    return {key: str(value) if isinstance(value, F) else value for key, value in plan.items()}
