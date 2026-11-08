"""Complete source-cube localization on a time-invariant response enclosure."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial
from pathlib import Path

import numpy as np

from model.native_ordered_ground import source, certificate
from model.native_joint_spin_spectroscopy import algebra
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds

ROOT = Path(__file__).resolve().parents[2]
BITS = 100
GRID = 1 << BITS


def _down(x):
    return F((x*GRID).numerator//(x*GRID).denominator, GRID)


def _up(x):
    return -_down(-x)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        if type(self.lo) is not F or type(self.hi) is not F or self.lo > self.hi:
            raise ValueError('ordered rational interval endpoints required')

    @classmethod
    def exact(cls, value):
        value = F(value)
        return cls(value, value)

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval.exact(other)
        return Interval(_down(self.lo+other.lo), _up(self.hi+other.hi))

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        other = other if isinstance(other, Interval) else Interval.exact(other)
        return self+-other

    def __rsub__(self, other):
        return Interval.exact(other)+-self

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval.exact(other)
        candidates = (self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi)
        return Interval(_down(min(candidates)), _up(max(candidates)))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError('interval inverse crosses zero')
        return Interval(_down(1/self.hi), _up(1/self.lo))

    def __truediv__(self, other):
        other = other if isinstance(other, Interval) else Interval.exact(other)
        return self*other.reciprocal()

    def square(self):
        lower = F() if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        return Interval(_down(lower), _up(max(self.lo*self.lo, self.hi*self.hi)))

    def widened(self, error):
        if type(error) is not F or error < 0:
            raise ValueError('nonnegative rational error required')
        return Interval(self.lo-error, self.hi+error)

    def decimal(self, places=11):
        return decimal_bounds(self.lo, self.hi, places)


def root(value):
    return Interval(*sqrt_bounds(F(value), BITS))


def inner(left, right):
    return sum(int(a)*int(b) for a, b in zip(left, right, strict=True))


def finite_basis():
    """Exact orthogonal character bases in the full odd-parity source space."""
    actions = algebra.reflections()
    h, scale = source.matrix()
    algebra.audit(h, scale)
    complement = tuple(i ^ 255 for i in range(256))
    if not np.array_equal(h[np.ix_(complement, complement)], h):
        raise ArithmeticError('actual source does not preserve global bit-flip parity')
    representatives = []
    for m in range(256):
        orbit = {action[m] for action in actions} | {action[m] ^ 255 for action in actions}
        if m == min(orbit):
            representatives.append(m)
    bases = []
    for k in range(8):
        columns = []
        for m in representatives:
            v = np.zeros(256, dtype=object)
            v[m], v[m ^ 255] = 1, -1
            v = algebra.character_projection(v, k)
            if any(v):
                columns.append(v)
        if len(columns) != 16:
            raise ArithmeticError('the full odd source character rank is not sixteen')
        norms = tuple(inner(v, v) for v in columns)
        for i, v in enumerate(columns):
            if not np.array_equal(v[list(complement)], -v):
                raise ArithmeticError('a character column left the odd response enclosure')
            for shift, action in enumerate(actions):
                if not np.array_equal(v[list(action)], (-1)**((k & shift).bit_count())*v):
                    raise ArithmeticError('source character basis covariance failed')
            for j, w in enumerate(columns):
                if inner(v, w) != (norms[i] if i == j else 0):
                    raise ArithmeticError('finite source frame is not orthogonal')
        bases.append((tuple(columns), norms))
    # The 8 orthogonal character spaces have total rank128, exactly dim P_odd.
    return h.astype(object), scale, tuple(bases)


def source_frames():
    """Attach the source-probe mode as the first axis of every full character."""
    h, scale, bases = finite_basis()
    with np.load(ROOT/'model/native_ordered_ground/witness.npz', allow_pickle=False) as archive:
        witness = {key: archive[key] for key in archive.files}
    checked = certificate.verify(witness)
    g = witness['vectors'][:, 0].astype(object)
    g = g+g[[i ^ 255 for i in range(256)]]
    norm = inner(g, g)
    probe = g*np.array([1-2*(i & 1) for i in range(256)], dtype=object)
    state_error = 2*F(checked['normalized_ground_vector_error'])
    projected = tuple(algebra.character_projection(probe, k) for k in range(8))
    weights = tuple(F(inner(v, v), 64*norm) for v in projected)
    minimum = min(weights)-2*state_error
    if sum(weights) != 1 or minimum <= F(1, 100):
        raise ArithmeticError('complete source first-axis weights lost their strict lower bound')
    mode_error = 2*state_error/minimum
    v = projected[7]
    hv = h@v
    vnorm = inner(v, v)
    first = F(inner(v, hv), scale*vnorm)
    second = F(inner(hv, hv), scale*scale*vnorm)
    leakage_lower = second-first*first-6*14**2*mode_error
    if leakage_lower <= F(1, 100):
        raise ArithmeticError('source response does not certify enlargement beyond its rank-eight probe code')
    frames, chosen = [], []
    for k, ((columns, norms), projected_vector) in enumerate(zip(bases, projected, strict=True)):
        vnorm = inner(projected_vector, projected_vector)
        dots = tuple(inner(b, projected_vector) for b in columns)
        scores = tuple(F(dot*dot, n*vnorm) for dot, n in zip(dots, norms, strict=True))
        pivot = max(range(16), key=lambda i: scores[i])
        if scores[pivot] < F(1, 16):
            raise ArithmeticError('source-aligned reference axis lost its overlap floor')
        order = (pivot,)+tuple(j for j in range(16) if j != pivot)
        signs = ((1 if dots[pivot] > 0 else -1),)+(1,)*15
        ordered = tuple(sign*columns[j] for sign, j in zip(signs, order, strict=True))
        ordered_norms = tuple(norms[j] for j in order)
        alpha = tuple(Interval.exact(inner(b, projected_vector))/root(n*vnorm)
                      for b, n in zip(ordered, ordered_norms, strict=True))
        if alpha[0].lo-mode_error <= 0:
            raise ArithmeticError('true source Householder chart is not certified away from its singularity')
        u = (alpha[0]+1,)+alpha[1:]
        denominator = 2*(1+alpha[0])
        # F=Householder(e0+alpha)*diag(-1,1,...,1), so F e0=alpha.
        frame = tuple(tuple(((1 if i == j else 0)-2*u[i]*u[j]/denominator)
                            *(-1 if j == 0 else 1)
                            for j in range(16)) for i in range(16))
        frames.append((ordered, ordered_norms, frame))
        chosen.append({'character': k, 'canonical_axis_index': pivot,
                       'squared_axis_overlap_lower': str(scores[pivot]),
                       'true_axis_overlap_strictly_positive': True})
    return {'h': h, 'scale': scale, 'probe': probe, 'norm': norm,
            'frames': tuple(frames), 'weights': weights,
            'state_error': state_error, 'mode_error': mode_error,
            'rank8_code_squared_energy_leakage_lower': leakage_lower,
            'complete_decoder_halfdiamond_error': 6*mode_error,
            'chosen_reference_axes': chosen}


def taylor_state(data, time=F(1, 8), order=28):
    if type(time) is not F or type(order) is not int or not 0 <= order <= 128:
        raise ValueError('rational time and a bounded nonnegative Taylor order required')
    real = [F() for _ in range(256)]
    imag = [F() for _ in range(256)]
    power = data['probe']
    for n in range(order+1):
        coefficient = time**n/F(data['scale']**n*factorial(n))
        sign = (1, -1, -1, 1)[n % 4]
        target = real if n % 2 == 0 else imag
        for j, value in enumerate(power):
            target[j] += sign*coefficient*int(value)
        if n != order:
            power = data['h']@power
    b = 14*abs(time)
    ceiling = -((-b).numerator//(-b).denominator)
    error = F(3**ceiling)*b**(order+1)/factorial(order+1)
    norm = root(data['norm'])
    state = tuple((Interval.exact(a)/norm, Interval.exact(b)/norm)
                  for a, b in zip(real, imag, strict=True))
    return state, error


def decoded_coefficients(data, state):
    all_channels = []
    for columns, norms, frame in data['frames']:
        coefficients = []
        for column, norm in zip(columns, norms, strict=True):
            real = sum((int(c)*v[0] for c, v in zip(column, state, strict=True) if c), Interval.exact(0))/root(norm)
            imag = sum((int(c)*v[1] for c, v in zip(column, state, strict=True) if c), Interval.exact(0))/root(norm)
            coefficients.append((real, imag))
        channels = tuple((sum((frame[j][i]*coefficients[j][0] for j in range(16)), Interval.exact(0)),
                          sum((frame[j][i]*coefficients[j][1] for j in range(16)), Interval.exact(0)))
                         for i in range(16))
        all_channels.append(channels)
    return tuple(all_channels)


def position_probabilities(data, time):
    state, taylor_error = taylor_state(data, time)
    channels = decoded_coefficients(data, state)
    error = (data['state_error']+data['complete_decoder_halfdiamond_error']
             +taylor_error*(2+taylor_error))
    probabilities = []
    for x in range(8):
        probability = Interval.exact(0)
        for i in range(16):
            real = sum(((-1)**((k & x).bit_count())*channels[k][i][0] for k in range(8)), Interval.exact(0))
            imag = sum(((-1)**((k & x).bit_count())*channels[k][i][1] for k in range(8)), Interval.exact(0))
            probability += (real.square()+imag.square())/8
        probability = probability.widened(error)
        probabilities.append([max(F(), probability.lo), min(F(1), probability.hi)])
    if not sum(p[0] for p in probabilities) <= 1 <= sum(p[1] for p in probabilities):
        raise ArithmeticError('complete time-evolved position intervals lost total probability')
    return {'time_in_same_HB_clock': str(time),
            'position_probability_intervals': [decimal_bounds(*p, 11) for p in probabilities],
            'dark_even_parity_probability_exact': '0',
            'taylor_vector_error_upper': str(taylor_error),
            'source_state_and_complete_decoder_error_upper': str(error)}


def locality_control(data):
    # Floating point selects a matrix entry only; rational intervals admit it.
    blocks = []
    for columns, norms, frame in data['frames']:
        basis = np.column_stack([np.asarray(v, dtype=float)/np.sqrt(float(n))
                                 for v, n in zip(columns, norms, strict=True)])
        rotation = np.array([[float((a.lo+a.hi)/2) for a in row] for row in frame])
        blocks.append(basis@rotation)
    eta = sum(blocks)/np.sqrt(8)
    effect = eta@eta.T
    candidates = [(abs(effect[a, b]), a, b) for a in range(256) for b in range(a+1, 256)
                  if (a ^ b).bit_count() >= 2]
    _, a, b = max(candidates)
    entries = []
    for mask in (a, b):
        row = []
        for i in range(16):
            value = Interval.exact(0)
            for columns, norms, frame in data['frames']:
                for j in range(16):
                    if columns[j][mask]:
                        value += int(columns[j][mask])*frame[j][i]/root(norms[j])
            row.append(value)
        entries.append(row)
    coefficient = sum((x*y for x, y in zip(*entries, strict=True)), Interval.exact(0))/8
    # Two decoder unitaries differing by delta give effect error at most2 delta.
    coefficient = coefficient.widened(2*data['complete_decoder_halfdiamond_error'])
    distance = min(abs(coefficient.lo), abs(coefficient.hi)) if coefficient.lo*coefficient.hi > 0 else F()
    if distance <= F(1, 1000):
        raise ArithmeticError('selected multisite effect entry lacks a strict certificate')
    return {'position': 0, 'source_basis_masks': [a, b],
            'differing_sites': (a ^ b).bit_count(),
            'true_effect_matrix_entry_interval': coefficient.decimal(),
            'absolute_entry_strict_lower': '1/1000',
            'effect_is_a_single_site_observable': False,
            'continuum_nonlocality_inferred_from_one_finite_cube': False}


def report():
    data = source_frames()
    initial = position_probabilities(data, F())
    evolved = position_probabilities(data, F(1, 8))
    return {
        'source': 'same open spin cube HB, J=1, B=1/4; unique true ground and local Z probe',
        'full_source_dimension': 256,
        'time_invariant_odd_response_enclosure_dimension': 128,
        'each_graph_reflection_character_dimension': 16,
        'position_outcomes': 8,
        'rank_of_each_position_effect': 16,
        'dark_even_parity_rank': 128,
        'first_probe_mode_in_each_character_preserved': True,
        'complete_position_pvm_and_native_group_covariance_proved': True,
        'dark_projection_commutes_with_actual_HB_exactly': True,
        'actual_minimal_response_dimension_lower': 9,
        'actual_minimal_response_dimension_upper': 128,
        'source_rank8_code_squared_energy_leakage_lower': str(data['rank8_code_squared_energy_leakage_lower']),
        'minimal_response_dimension_exactly_determined': False,
        'all_true_character_weights_strictly_above': '1/100',
        'state_vector_error_upper': str(data['state_error']),
        'complete_decoder_halfdiamond_error_upper': str(data['complete_decoder_halfdiamond_error']),
        'chosen_reference_axes': data['chosen_reference_axes'],
        'same_clock_position_readouts': [initial, evolved],
        'single_site_locality_countercontrol': locality_control(data),
        'native_instrument_tape_executed': False,
        'graph_reflections_identified_as_continuous_momentum': False,
        'field_energy_levels_identified_as_particle_masses': False,
    }
