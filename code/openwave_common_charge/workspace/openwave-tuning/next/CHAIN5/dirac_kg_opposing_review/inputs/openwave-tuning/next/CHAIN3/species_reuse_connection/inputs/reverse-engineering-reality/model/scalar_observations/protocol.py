"""One preparation, two physical field reads, and recorded continuation.

The same seven real Fourier coordinates, smearing, preparation, feedback and
Hamiltonian model times apply on every spatial quadrature.  Occupation is a
separate approximation.  No branch is normalized, discarded or silently reset.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import cosm, expm, sinm

from ..observer import native
from . import source


@dataclass(frozen=True)
class Settings:
    interval: F = F(1, 8)
    mode: int = 1  # sqrt(2/V) cos(x), on the fixed L=2*pi torus
    eta: float = 1.0
    theta: float = math.pi/2
    initial_kick: float = 0.25
    feedback_kick: float = 0.25

    def __post_init__(self):
        if not isinstance(self.interval, (int, F)) or isinstance(self.interval, bool) or self.interval < 0:
            raise ValueError('nonnegative exact Hamiltonian model interval required')
        if type(self.mode) is not int or not 0 <= self.mode < 7:
            raise ValueError('mode must be an integer between zero and six')
        for name in ('eta', 'theta', 'initial_kick', 'feedback_kick'):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be a finite real number')


def evolve(lattice, state, time, steps):
    """Execute symmetric onsite/pair factors; pair generators all commute."""
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer source step count required')
    dt = float(time)/steps
    onsite = [source.onsite(d, i, dt/2, lattice.parameters)[0]
              for i, d in enumerate(lattice.dimensions)]
    pairs = [((i, j), source.pair(lattice.dimensions[i], lattice.dimensions[j],
              lattice.pair_coefficient(i, j), dt, lattice.parameters.omega)[0])
             for i in range(7) for j in range(i+1, 7)]
    out = np.asarray(state, dtype=complex).copy()
    for _ in range(steps):
        for i, matrix in enumerate(onsite):
            out = lattice.apply_local(out, [i], matrix)
        for indices, matrix in pairs:
            out = lattice.apply_local(out, indices, matrix)
        for i, matrix in reversed(list(enumerate(onsite))):
            out = lattice.apply_local(out, [i], matrix)
    return out


def history(lattice, steps=512, sign=1, settings=Settings(), *,
            source_execution=True, feedback=True, reset=False):
    """Return all four joint outcomes on the same continued field state.

    ``reset`` is a separately labelled destructive control: replace each
    disturbed branch by the original prepared state with its branch weight.
    It is never part of the principal history.  The reference uses independent
    sparse Hamiltonian evolution and matrix functions, not source matrices.
    """
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('native preparation sign must be -1 or +1')
    if not isinstance(settings, Settings):
        raise ValueError('Settings instance required')
    m, p = settings.mode, lattice.parameters
    d = lattice.dimensions[m]
    if source_execution:
        prepare = source.feedback(d, sign*settings.initial_kick, p.omega)[0]
        kraus = source.instrument(d, settings.eta, settings.theta, p.omega)[0]
        control = source.feedback(d, settings.feedback_kick, p.omega)[0]
        propagate = lambda v: evolve(lattice, v, settings.interval, steps)
    else:
        q = lattice.operators[m]['q']
        angle = (settings.eta*q+settings.theta*np.eye(d))/2
        prepare = expm(1j*sign*settings.initial_kick*q)
        kraus = cosm(angle), sinm(angle)
        control = expm(1j*settings.feedback_kick*q)
        propagate = lambda v: lattice.evolve(v, float(settings.interval))
    prepared = lattice.apply_local(lattice.vacuum_state(), [m], prepare)
    state = propagate(prepared)
    branches = np.column_stack([lattice.apply_local(state, [m], k) for k in kraus])
    if reset:
        for bit in (0, 1):
            branches[:, bit] = math.sqrt(float(np.vdot(branches[:, bit], branches[:, bit]).real))*prepared
    if feedback:
        branches[:, 1] = lattice.apply_local(branches[:, 1], [m], control)
    branches = propagate(branches)
    return [{'bits': [first, last], 'probability': float(np.vdot(final, final).real)}
            for first in (0, 1) for last, k in enumerate(kraus)
            for final in [lattice.apply_local(branches[:, first], [m], k)]]


def tv(left, right):
    if [row['bits'] for row in left] != [row['bits'] for row in right]:
        raise ValueError('record alphabets differ')
    return math.fsum(abs(a['probability']-b['probability']) for a, b in zip(left, right))/2


def native_records(histories, feedback=True):
    """First native current prepares the kick; last field bit chooses a seam."""
    rows = []
    for initial, loads in enumerate(native.states(binary=True)):
        for coin in (0, 1):
            middle, (_, current) = native.native(loads, native.E, coin)
            phase = int(current > 0)
            for quantum in histories[phase]:
                first, last = quantum['bits']
                seam = native.F if feedback and last else native.E
                for final_coin in (0, 1):
                    final, (_, final_current) = native.native(middle, seam, final_coin)
                    rows.append({'initial': initial,
                        'records': [coin, current, first, last, final_coin, final_current],
                        'final_port_zero': final[0], 'probability': quantum['probability']/264})
    return rows


def source_cost(lattice, steps, settings=Settings()):
    """Quantum compiler costs for both intervals, reads and kicks.

    The two surrounding native repairs and classical routing/control overhead
    are separate.  This counts code load/export services in the emitted words;
    it does not measure host runtime or calibrate a physical clock.
    """
    dt = float(settings.interval)/steps
    p, m = lattice.parameters, settings.mode
    d = lattice.dimensions[m]
    requests = [(4*steps, source.onsite(size, i, dt/2, p)[1])
                for i, size in enumerate(lattice.dimensions)]
    requests += [(2*steps, source.pair(lattice.dimensions[i], lattice.dimensions[j],
                  lattice.pair_coefficient(i, j), dt, p.omega)[1])
                 for i in range(7) for j in range(i+1, 7)]
    requests += [(2, source.instrument(d, settings.eta, settings.theta, p.omega)[1]),
                 (1, source.feedback(d, settings.initial_kick, p.omega)[1]),
                 (1, source.feedback(d, settings.feedback_kick, p.omega)[1])]
    return {key: sum(count*meta['resources'][key] for count, meta in requests)
            for key in source.COUNT_KEYS}


def detector_calibration(lattice):
    """Two independent blank-state read settings, with no fitted parameters.

    The unrestricted Omega Gaussian has <cos(eta Q)>=exp(-eta^2/(4 Omega)).
    The finite-code comparison diagnoses its detector approximation; it is
    neither experimental calibration nor a cutoff certificate for later reads.
    """
    rows = []
    d, omega = lattice.dimensions[1], lattice.parameters.omega
    blank = np.zeros(d, complex)
    blank[0] = 1
    for eta in (1., 2.):
        kraus, metadata = source.instrument(d, eta, 0., omega)
        probability = float(np.linalg.norm(kraus[0]@blank)**2)
        reference = (1+math.exp(-eta**2/(4*omega)))/2
        rows.append({'eta': eta, 'theta': 0, 'source_bit_zero_probability': probability,
                     'unrestricted_gaussian_probability': reference,
                     'absolute_occupation_difference': abs(probability-reference),
                     'source_metadata': metadata})
    return {'settings': rows, 'fitted_parameters': [], 'experimental_data_used': False,
            'scope': 'Supplied preparation and detector consistency; no physical parameter identification'}
