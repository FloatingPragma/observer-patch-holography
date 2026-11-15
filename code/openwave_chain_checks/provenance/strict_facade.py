"""Checked entry points for frozen mathematical helpers; no input correction."""
from __future__ import annotations

from numbers import Real
import math
import numpy as np

STATE_TOL = 2e-12


def scalar(value, name, *, positive=False, upper=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f'{name}: real scalar required')
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f'{name}: finite {"positive" if positive else "nonnegative"} value required')
    if upper is not None and value > upper:
        raise ValueError(f'{name}: exceeds {upper}')


def density(rho, dimension, *, normalized=False, sectors=None):
    """Check the floating representation without projection or normalization.

    Exact mathematical inputs are positive Hermitian subnormalized matrices.
    The checker permits STATE_TOL absolute roundoff in these spectral tests.
    It returns the input object unchanged and never clips an eigenvalue.
    """
    if not isinstance(rho, np.ndarray) or rho.shape != (dimension, dimension):
        raise ValueError('density: exact square ndarray shape required')
    if rho.dtype.kind not in 'fci' or not np.isfinite(rho).all():
        raise ValueError('density: finite numeric entries required')
    if np.max(np.abs(rho - rho.conj().T)) > STATE_TOL:
        raise ValueError('density: Hermitian input required')
    mass = np.trace(rho)
    if abs(mass.imag) > STATE_TOL or mass.real < -STATE_TOL or mass.real > 1 + STATE_TOL:
        raise ValueError('density: subnormalized trace required')
    if normalized and abs(mass - 1) > 1e-12:
        raise ValueError('trajectory: normalized density required')
    if sectors is None:
        blocks = [rho]
    else:
        off = np.array(rho, copy=True)
        blocks = []
        for ids, *_ in sectors.values():
            blocks.append(rho[np.ix_(ids, ids)])
            off[np.ix_(ids, ids)] = 0
        if np.any(off):
            raise ValueError('density: exact number-block-diagonal input required')
    if min(float(np.linalg.eigvalsh(block).min()) for block in blocks) < -STATE_TOL:
        raise ValueError('density: positive input required')
    return rho


class StrictClock:
    """Validate every public state/parameter argument before calling raw code."""
    def __init__(self, raw):
        self._raw = raw

    def evolve(self, rho, duration, h_scale=1.0):
        density(rho, self._raw.n)
        scalar(duration, 'duration'); scalar(h_scale, 'h_scale')
        return self._raw.evolve(rho, duration, h_scale)

    def native(self, rho, detailed=False):
        density(rho, self._raw.n)
        if type(detailed) is not bool:
            raise ValueError('detailed: Boolean required')
        return self._raw.native(rho, detailed=detailed)

    def semigroup(self, rho, duration, rate, *, attempt_mark=1.0, move_mark=1.0,
                  idle_rate=0.0, idle_mark=1.0, h_scale=1.0):
        density(rho, self._raw.n, sectors=self._raw.sectors)
        for name, value in [('duration',duration), ('rate',rate), ('idle_rate',idle_rate), ('h_scale',h_scale)]:
            scalar(value, name)
        for name, value in [('attempt_mark',attempt_mark), ('move_mark',move_mark), ('idle_mark',idle_mark)]:
            scalar(value, name, upper=1)
        return self._raw.semigroup(rho, duration, rate, attempt_mark=attempt_mark,
            move_mark=move_mark, idle_rate=idle_rate, idle_mark=idle_mark, h_scale=h_scale)

    def count_split(self, rho, duration, rate):
        density(rho, self._raw.n, sectors=self._raw.sectors)
        scalar(duration, 'duration'); scalar(rate, 'rate')
        return self._raw.count_split(rho, duration, rate)

    def trajectory(self, rho, duration, rate, uniform_rate, events):
        density(rho, self._raw.n, normalized=True)
        scalar(duration, 'duration'); scalar(rate, 'rate', positive=True)
        scalar(uniform_rate, 'uniform_rate', positive=True)
        if rate > uniform_rate or not isinstance(events, (list, tuple)):
            raise ValueError('ordered events and uniform_rate >= rate required')
        previous = 0.
        for event in events:
            if not isinstance(event, (list, tuple)) or len(event) != 2:
                raise ValueError('event: (time, mark) required')
            at, mark = event
            scalar(at, 'event time', positive=True)
            if not previous < at < duration:
                raise ValueError('strictly ordered interior event times required')
            if mark != 'idle' and (type(mark) is not int or not 0 <= mark < len(self._raw.ops)):
                raise ValueError('invalid native mark')
            if mark == 'idle' and uniform_rate == rate:
                raise ValueError('zero probability idle path')
            previous = at
        return self._raw.trajectory(rho, duration, rate, uniform_rate, events)


class StrictParticle:
    def __init__(self, raw):
        self._raw = raw

    def measure(self, rho):
        density(rho, self._raw.model.n, sectors=self._raw.clock.sectors)
        return self._raw.measure(rho)

    def energy_coefficients(self, rho):
        density(rho, self._raw.model.n, sectors=self._raw.clock.sectors)
        return self._raw.energy_coefficients(rho)


def first_event(raw, state, edges, rates, native, horizon=.001, target_seam=0):
    """Check mathematical domains; source identity requires classical_witness."""
    scalar(horizon, 'horizon')
    if not isinstance(state, (tuple, list)) or len(state) != 2 or any(type(i) is not int or not 0 <= i < 96 for i in state):
        raise ValueError('state: two integer positions in the 96-port source required')
    if not isinstance(edges, (tuple, list)) or len(edges) != 288:
        raise ValueError('edges: 288 seam domain required')
    if any(not isinstance(e, (tuple,list)) or len(e) != 2 or any(not isinstance(i,(int,np.integer)) or not 0 <= i < 96 for i in e) or e[0] == e[1] for e in edges):
        raise ValueError('edges: finite source endpoint pair required')
    if len(set(tuple(e) for e in edges)) != len(edges):
        raise ValueError('edges: duplicate seam')
    if not isinstance(rates, np.ndarray) or rates.shape != (2,288) or rates.dtype.kind not in 'fi' or not np.isfinite(rates).all() or np.any(rates < 0):
        raise ValueError('rates: finite nonnegative array of shape (2,288) required')
    scalar(float(rates.sum()), 'total rate', positive=True)
    if type(target_seam) is not int or not 0 <= target_seam < len(edges) or not callable(native):
        raise ValueError('target seam and callable native endpoint rule required')
    return raw.first_event(state, edges, rates, native, horizon, target_seam)
