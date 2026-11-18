"""Uniform local clock, exact packet budgets, and literal joint columns."""
from dataclasses import dataclass
from fractions import Fraction
from math import comb, sqrt

import numpy as np
from scipy.fft import dst
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

from .program import PROGRAM_EDGES, Program, edge, prefixes, preparation, validate_data


def require(condition, message):
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Rail:
    sites: int = 4096
    offset: int = 2048
    order: int = 64

    def __post_init__(self):
        require(all(type(x) is int for x in (self.sites, self.offset, self.order)),
                "integer rail parameters")
        require(self.order >= 1 and self.offset >= 0
                and self.offset + self.order + PROGRAM_EDGES < self.sites,
                "packet and complete program must fit on the rail")

    @property
    def start(self):
        return self.offset + self.order

    @property
    def completed(self):
        return self.start + PROGRAM_EDGES

    def stage(self, position):
        require(type(position) is int and 0 <= position < self.sites, "clock position")
        return min(PROGRAM_EDGES, max(0, position - self.start))


def fraction_record(value):
    value = Fraction(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def moments(order):
    require(type(order) is int and order >= 1, "positive packet order")
    m = order
    return (Fraction(m, 2), Fraction(m * m, 4 * (2 * m - 1)),
            Fraction(m, m + 1), Fraction(2 * m + 1, (m + 1) ** 2 * (m + 2)))


def interval_certificate(rail, t0=256, t1=512):
    """Exact rational Cantelli and finite-boundary certificate, not a time grid."""
    require(type(rail) is Rail and type(t0) is int and type(t1) is int
            and 0 < t0 <= t1, "ordered positive integer interval")
    mean0, variance0, velocity, velocity_variance = moments(rail.order)
    separation = velocity * t0 - (rail.order - mean0 + PROGRAM_EDGES - 1)
    require(separation > 0, "interval begins before the certified packet clears the program")
    variance = variance0 + velocity_variance * t0 * t0
    cantelli = variance / (variance + separation * separation)
    d = min(rail.offset + 1, rail.sites - rail.offset - rail.order)
    ratio = Fraction(3 * t1, d)
    require(ratio < 1 and t1 < d + 1, "finite-boundary tail certificate unavailable")
    boundary = 2 * ratio ** d / (1 - Fraction(t1, d + 1))
    failure = cantelli + boundary
    # A short enclosing dyadic fraction avoids printing the much larger exact tail.
    exponent = max(0, boundary.denominator.bit_length() - boundary.numerator.bit_length() - 1)
    upper = Fraction(1, 2 ** exponent)
    require(boundary <= upper, "boundary enclosure")
    return {
        "window": [t0, t1], "completion_position": rail.completed,
        "position_variance": fraction_record(variance0),
        "mean_velocity": fraction_record(velocity),
        "velocity_variance": fraction_record(velocity_variance),
        "cantelli_failure": fraction_record(cantelli),
        "boundary_distance": d, "boundary_ratio": fraction_record(ratio),
        "boundary_power": d,
        "boundary_prefactor": fraction_record(2 / (1 - Fraction(t1, d + 1))),
        "boundary_upper": fraction_record(upper),
        "failure_upper": fraction_record(cantelli + upper),
        "failure_below_one_in_thousand": failure < Fraction(1, 1000),
        "scope": "one terminal pointer read, arbitrary retained source reference; ideal declared preparation",
    }


def packet(rail, direction=1):
    require(type(rail) is Rail and type(direction) is int and direction in (-1, 1),
            "rail and packet direction")
    out = np.zeros(rail.sites, complex)
    denominator = sqrt(comb(2 * rail.order, rail.order))
    for j in range(rail.order + 1):
        phase = (1, -1j * direction, -1, 1j * direction)[j % 4]
        out[rail.offset + j] = phase * (comb(rail.order, j) / denominator)
    return out


def preparation_error_bound(rail):
    """Exact vector bound for the stored binary amplitudes versus the analytic packet.

    This bounds packet representation only, not SciPy propagation roundoff.
    """
    z = comb(2 * rail.order, rail.order)
    s = Fraction.from_float(sqrt(z))
    denominator_error = abs(s * s - z) / (s * s)
    component_error = Fraction()
    for position, value in enumerate(packet(rail)):
        j = position - rail.offset
        if 0 <= j <= rail.order:
            coefficient = Fraction(comb(rail.order, j), 1) / s
            real, imag = ((1, 0), (0, -1), (-1, 0), (0, 1))[j % 4]
        else:
            coefficient, real, imag = Fraction(), 0, 0
        component_error += abs(Fraction.from_float(float(value.real)) - real * coefficient)
        component_error += abs(Fraction.from_float(float(value.imag)) - imag * coefficient)
    return denominator_error + component_error


def spectrum(rail):
    return np.cos(np.pi * np.arange(1, rail.sites + 1) / (rail.sites + 1))


def evolve_clock(rail, tau, direction=1):
    require(type(tau) in (int, float) and np.isfinite(tau), "finite clock time")
    initial = packet(rail, direction)
    modes = dst(initial, type=1, norm="ortho")
    out = dst(np.exp(-1j * float(tau) * spectrum(rail)) * modes, type=1, norm="ortho")
    require(np.all(np.isfinite(out)), "nonfinite propagated packet")
    return out


def clock_matrix(rail):
    return diags([np.full(rail.sites - 1, 0.5)] * 2, [-1, 1], shape=(rail.sites, rail.sites), format="csr")


def column(rail, config, basis):
    """Every full-carrier nonzero column entry; no prepared-sector assumption."""
    require(type(rail) is Rail and type(config) is Program, "rail and program")
    require(type(basis) is tuple and len(basis) == 2, "joint basis")
    position, data = basis
    rail.stage(position)
    validate_data(data)
    out = {}
    for direction in (-1, 1):
        target = position + direction
        if 0 <= target < rail.sites:
            gate = min(position, target) - rail.start
            changed = edge(config, data, gate, direction) if 0 <= gate < PROGRAM_EDGES else data
            out[(target, changed)] = 0.5
    return out


def prepared_span(rail, config, source_amplitudes):
    """Literal complete reducing span for arbitrary supplied source superpositions."""
    require(type(source_amplitudes) is dict and source_amplitudes, "source amplitudes")
    require(all(type(s) is int and 0 <= s < 4096 for s in source_amplitudes), "source labels")
    require(all(np.isfinite(a) for a in source_amplitudes.values())
            and abs(sum(abs(a) ** 2 for a in source_amplitudes.values()) - 1) < 1e-12,
            "unit source state required, without implicit normalization")
    states, starts = [], {}
    q = packet(rail)
    for source, amplitude in sorted(source_amplitudes.items()):
        for address in range(12):
            history = prefixes(config, preparation(source, address))
            for position in range(rail.sites):
                state = (position, history[rail.stage(position)])
                states.append(state)
                if q[position] != 0:
                    starts[state] = complex(amplitude) * q[position] / sqrt(12)
    index = {state: i for i, state in enumerate(states)}
    require(len(index) == len(states), "lost information in complete prepared histories")
    rows, columns, values = [], [], []
    for i, state in enumerate(states):
        for dest, value in column(rail, config, state).items():
            require(dest in index, "nonzero output outside declared reducing span")
            rows.append(index[dest]); columns.append(i); values.append(value)
    h = coo_matrix((values, (rows, columns)), shape=(len(states), len(states))).tocsr()
    require((h - h.T).nnz == 0, "missing Hermitian column")
    psi = np.zeros(len(states), complex)
    for state, amplitude in starts.items():
        psi[index[state]] = amplitude
    return states, h, psi


def literal_evolve(h, initial, tau):
    require(type(tau) in (int, float) and np.isfinite(tau), "finite clock time")
    return expm_multiply(-1j * tau * h, initial, traceA=0)


def factored_vector(rail, config, source_amplitudes, states, tau):
    """Gauge-theorem prediction, compared against literal column propagation."""
    q = evolve_clock(rail, tau)
    index = {state: i for i, state in enumerate(states)}
    out = np.zeros(len(states), complex)
    for source, amplitude in source_amplitudes.items():
        for address in range(12):
            history = prefixes(config, preparation(source, address))
            for n in range(rail.sites):
                out[index[(n, history[rail.stage(n)])]] = amplitude * q[n] / sqrt(12)
    return out
