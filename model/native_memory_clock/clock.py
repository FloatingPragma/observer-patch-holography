"""Literal finite rail columns and complete reachable-state evolution.

The full carrier is specified by Basis, including dirty buffers and records.
Numerical evolution uses a checked reducing span, never a spectral window.
"""
from dataclasses import dataclass
from itertools import product
from math import comb, cos, sin, sqrt

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from .native import BACKGROUNDS, forward, inverse


def require(condition, message):
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Construction:
    horizon: int
    background: tuple = (0, 0, 0, 0)
    law: str = "native"
    phase: bool = False
    profile: str = "transfer"

    def __post_init__(self):
        require(type(self.horizon) is int and self.horizon > 0, "positive integer horizon")
        require(type(self.background) is tuple and len(self.background) == 4
                and all(type(x) is int for x in self.background)
                and self.background in BACKGROUNDS, "declared integer background")
        require(self.law in ("native", "hold", "uphill"), "unknown event law")
        require(type(self.phase) is bool, "phase must be boolean")
        require(self.profile in ("transfer", "uniform"), "unknown clock profile")

    def coupling(self, j):
        require(type(j) is int and 0 <= j < self.horizon, "edge outside clock rail")
        if self.profile == "uniform":
            return 0.5
        return sqrt((j + 1) * (self.horizon - j)) / self.horizon

    def resources(self):
        t = self.horizon
        data = 2 ** (12 * (t + 1)) * 24 ** t
        return {"source_buffer_qubits": 12 * (t + 1), "old_bit_qubits": t,
                "address_qudits": t, "address_dimension": 12,
                "clock_positions": t + 1, "full_data_dimension": str(data),
                "full_hilbert_dimension": str((t + 1) * data),
                "archive_dimension": str(24 ** t),
                "address_history_dimension_lower_bound": str(12 ** t),
                "clock_carrier": "exactly one particle on an open rail"}


@dataclass(frozen=True, order=True)
class Basis:
    clock: int
    buffers: tuple
    addresses: tuple
    records: tuple


def validate_basis(c, b):
    require(isinstance(b, Basis), "Basis required")
    t = c.horizon
    require(type(b.clock) is int and 0 <= b.clock <= t, "invalid clock position")
    for values, size, upper, name in ((b.buffers, t + 1, 4096, "buffers"),
                                     (b.addresses, t, 12, "addresses"),
                                     (b.records, t, 2, "records")):
        require(type(values) is tuple and len(values) == size
                and all(type(v) is int and 0 <= v < upper for v in values), name)


def step(c, b, direction):
    """One exact gate-labelled hop, also defined on nonblank data.

    Each term touches only neighbouring rail cells. A cell's tetrahedral
    buffer is a supplied bounded-size patch, not emergent physical space.
    """
    validate_basis(c, b)
    require(type(direction) is int and direction in (-1, 1), "hop direction")
    j = b.clock if direction == 1 else b.clock - 1
    require(0 <= j < c.horizon, "hop beyond boundary")
    buffers, records = list(b.buffers), list(b.records)
    if direction == 1:
        old = (buffers[j] >> b.addresses[j]) & 1
        buffers[j], records[j] = forward(buffers[j], b.addresses[j], records[j],
                                         c.background, c.law)
        buffers[j], buffers[j + 1] = buffers[j + 1], buffers[j]
    else:
        buffers[j], buffers[j + 1] = buffers[j + 1], buffers[j]
        old = records[j]
        buffers[j], records[j] = inverse(buffers[j], b.addresses[j], records[j],
                                         c.background, c.law)
    phase = -1 if c.phase and old else 1
    return Basis(b.clock + direction, tuple(buffers), b.addresses, tuple(records)), phase


def column(c, b):
    """Every nonzero column of the complete Hamiltonian, at most two entries."""
    validate_basis(c, b)
    out = {}
    for direction in (-1, 1):
        if 0 <= b.clock + direction <= c.horizon:
            dest, phase = step(c, b, direction)
            j = min(b.clock, dest.clock)
            out[dest] = phase * c.coupling(j)
    return out


def prepared_span(c, source_amplitudes):
    """Enumerate ALL address words and all clock positions of the prepared span.

    No native source-state projection is used: any initial source amplitudes
    can be supplied. Full-carrier statements are proved through column().
    """
    require(isinstance(source_amplitudes, dict) and source_amplitudes, "source amplitudes")
    require(all(type(s) is int and 0 <= s < 4096 for s in source_amplitudes), "source state")
    amps = np.array(list(source_amplitudes.values()), dtype=complex)
    require(np.all(np.isfinite(amps)) and abs(np.vdot(amps, amps).real - 1) < 1e-12,
            "source must have unit norm, with no implicit normalization")
    states, starts = set(), {}
    for source, amplitude in source_amplitudes.items():
        for addresses in product(range(12), repeat=c.horizon):
            b = Basis(0, (source,) + (0,) * c.horizon, addresses, (0,) * c.horizon)
            starts[b] = complex(amplitude) / sqrt(12 ** c.horizon)
            states.add(b)
            for _ in range(c.horizon):
                b, _ = step(c, b, 1)
                states.add(b)
    states = sorted(states)
    index = {b: i for i, b in enumerate(states)}
    require(len(states) == (c.horizon + 1) * 12 ** c.horizon * len(source_amplitudes),
            "histories collide: missing retained information")
    rows, cols, values = [], [], []
    for i, b in enumerate(states):
        for dest, value in column(c, b).items():
            require(dest in index, "nonzero complementary output of proposed reducing span")
            rows.append(index[dest]); cols.append(i); values.append(value)
    h = coo_matrix((values, (rows, cols)), shape=(len(states), len(states))).tocsr()
    require((h - h.T.conjugate()).nnz == 0, "non-Hermitian literal rail")
    psi = np.zeros(len(states), complex)
    for b, amp in starts.items():
        psi[index[b]] = amp
    return states, h, psi


def evolve(h, psi, tau):
    require(np.isscalar(tau) and np.isreal(tau) and np.isfinite(tau), "finite source time")
    out = expm_multiply((-1j * float(tau)) * h, psi, traceA=0)
    require(np.all(np.isfinite(out)), "nonfinite evolution")
    require(abs(np.vdot(out, out).real - np.vdot(psi, psi).real) < 2e-10,
            "unitarity residual exceeds fixed numerical tolerance")
    return out


def clock_matrix(c):
    h = np.zeros((c.horizon + 1, c.horizon + 1))
    for j in range(c.horizon):
        h[j, j + 1] = h[j + 1, j] = c.coupling(j)
    return h


def clock_probabilities(t, tau):
    require(type(t) is int and t > 0 and np.isfinite(tau), "clock arguments")
    return np.array([comb(t, n) * sin(tau / t) ** (2 * n)
                     * cos(tau / t) ** (2 * (t - n)) for n in range(t + 1)])


def readout(states, psi):
    """Complete terminal PVM, coarsened only over unread future registers.

    Outcome: (actual clock position, source buffer there, occurred address/old
    bit prefix). The future tape remains in the state. Zero outcomes outside
    the exact reducing support are specified by this instrument, not deleted
    from the underlying carrier.
    """
    require(len(states) == len(psi), "readout carrier mismatch")
    out = {}
    for b, amp in zip(states, psi):
        key = (b.clock, b.buffers[b.clock], b.addresses[:b.clock], b.records[:b.clock])
        out[key] = out.get(key, 0.0) + float(abs(amp) ** 2)
    return out


def clock_marginal(states, psi, horizon):
    out = np.zeros(horizon + 1)
    for b, amp in zip(states, psi):
        out[b.clock] += abs(amp) ** 2
    return out


def measured_return(c, states, h, initial, read_time, total_time):
    """Nonselective mid-flight clock PVM; retain every branch through return."""
    require(0 <= read_time <= total_time, "ordered readout times")
    mid = evolve(h, initial, read_time)
    outcomes = np.zeros(c.horizon + 1)
    for n in range(c.horizon + 1):
        branch = mid * np.array([b.clock == n for b in states])
        final = evolve(h, branch, total_time - read_time)
        outcomes += clock_marginal(states, final, c.horizon)
    return outcomes
