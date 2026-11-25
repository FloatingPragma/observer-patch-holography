"""Independent clock and product-detector records after two unconstrained source writes."""

import json
from itertools import combinations, product
from pathlib import Path

import numpy as np
from mpmath import iv

from .....full_algebra import station
from .....full_algebra.enclosure import geometry_jet
from .....full_algebra.reuse import clock_history
from ...geometry import HERE as BALLS, interval_numbers, raw_bounds

HERE = Path(__file__).resolve().parent
ALGEBRA = BALLS.parents[1] / "full_algebra"
TICK_BUDGET = 2**38
SEED = 261126
EDGES = list(combinations(range(4), 2))
# 3 preparations x (12 directed clocks + 2 x 16 Pauli coefficients).
FAMILY = 3 * (12 + 32)
TAIL = 0.01
MEAN_ERROR = np.sqrt(np.log(2 * FAMILY / TAIL) / (2 * TICK_BUDGET))


def log_length_map():
    """Each column corresponds to one unit change in log directed return time."""
    incidence = []
    for axis in np.eye(12):
        raw = station.square(axis)
        gain = (raw - raw.T).sum(axis=1) / 4
        incidence.append(
            [(raw[i, j] + raw[j, i] - gain[i] - gain[j]) / 2 for i, j in EDGES]
        )
    return np.array(incidence).T


def weights():
    cert = json.loads((ALGEBRA / "certificate.json").read_text())
    original = json.loads((ALGEBRA / "results.json").read_text())
    ranges = raw_bounds()[station.ORDER]
    area, volume = geometry_jet(ranges)
    eta = iv.mpf([str(a) for a in cert["gain_interval"]])
    # Center used in the interval multiplier certificate, enclosing its decimal roundoff.
    assert abs(original["price"] - 51.0886428) < 1e-7
    assert cert["multiplier_drift_upper"] < 0.04299
    lam = iv.mpf("51.0886428") + iv.mpf(["-.043", ".043"])
    directed = [
        r * (eta * a - lam * v)
        for r, a, v in zip(ranges, area.g, volume.g, strict=True)
    ]
    by_edge = dict(zip(zip(*station.ORDER), directed, strict=True))
    edge = [by_edge[i, j] + by_edge[j, i] for i, j in EDGES]
    assert all(w.a > 14.71 for w in directed)
    # Euler homogeneity gives a tighter enclosure of the sum than twelve interval additions.
    total = 2 * eta * area.v - 3 * lam * iv.mpf(
        [str(v) for v in cert["volume_interval"]]
    )
    volume_grad = dict(zip(zip(*station.ORDER), ranges * volume.g, strict=True))
    # Both directions scale together on each edge and preserve the enrolled clocks.
    signed_volume = (
        volume_grad[0, 1] + volume_grad[1, 0] - volume_grad[0, 2] - volume_grad[2, 0]
    )
    signed_entropy = edge[1] - edge[0]
    assert signed_volume.a > 0.7494 and signed_entropy.a > 16.95
    return dict(
        shape_control_volume_slope=interval_numbers(signed_volume),
        shape_control_entropy_slope=interval_numbers(signed_entropy),
        directed=[interval_numbers(w) for w in directed],
        edges=[interval_numbers(b) for b in edge],
        total=interval_numbers(total),
    )


def preparations():
    record = json.loads((ALGEBRA / "results.json").read_text())
    # Use a specified positive lottery, separately from the clock and tomography records.
    axes = station.phase.AXES[1:]
    projectors = [(np.eye(2) + sign * a) / 2 for a in axes for sign in (-1, 1)]
    atoms = np.array([np.kron(a, b) for a, b in product(projectors, repeat=2)])
    lottery = np.array([record["preparation"][k] for k in ("kept", "dephased")])
    # Fix an exactly normalized positive dyadic preparation for repeat acquisition.
    numerators = np.rint(lottery * 2**46).astype(np.int64)
    numerators[0, 0] += 2**46 - numerators.sum()
    assert numerators.min() > 0
    kept, dephased = np.einsum("sp,pij->sij", numerators / 2**46, atoms)
    x = np.stack([kept + dephased / 2, dephased / 2])
    total = x.sum(axis=0)
    left = [np.trace(b.reshape(2, 2, 2, 2), axis1=1, axis2=3) for b in x]
    right = np.trace(total.reshape(2, 2, 2, 2), axis1=0, axis2=2)
    return {
        "reference": x,
        "ordinary_write": np.stack([total / 2, total / 2]),
        "independent_right": np.stack([np.kron(m, right) for m in left]),
    }


def clock_effects():
    """Accumulate the charge-pooled, bounded native first-return effects by route and tick."""
    bins = np.zeros((12, 10, 2, 4, 4), dtype=complex)
    for _, parity, route, tick, effect in clock_history.native_events():
        bins[route, tick - 2, parity] += (
            effect[0].astype(float) + 1j * effect[1].astype(float)
        ) / (2 * clock_history.F)
    return bins


def clock_read(state, bins, rng):
    # The twelfth outcome time is also the right-censored outcome. Keep all starts.
    timestamps = np.arange(2, 13) / 4
    probabilities = np.einsum("rtpij,pji->rt", bins, state).real
    rates = np.column_stack([probabilities, 1 - probabilities.sum(axis=1)])
    assert rates.min() >= 0
    counts = np.array([rng.multinomial(TICK_BUDGET, rate) for rate in rates])
    means = counts @ timestamps / TICK_BUDGET
    # The score 0.5 ... 3 has span 2.5. This extra allowance covers dyadic
    # effect approximation and binary probability arithmetic on four-coin blocks.
    assert float(clock_history.timing_error()[1]) < 1e-12
    width = 2.5 * MEAN_ERROR + 1e-10
    ranges = np.array(
        [iv.mpf([str(v - width), str(v + width)]) for v in means], dtype=object
    )
    area, volume = geometry_jet(ranges)
    matrix = log_length_map()
    log_lengths = matrix @ np.array([iv.ln(v) for v in ranges], dtype=object)
    return dict(
        counts=counts.tolist(),
        native_mean_times=means.tolist(),
        length_bands=[interval_numbers(iv.exp(z)) for z in log_lengths],
        area=interval_numbers(area.v),
        volume=interval_numbers(volume.v),
    )


def material_read(state, rng):
    """Two local Pauli outcomes with a parity tag and a separate supply for each setting."""
    axes = [np.eye(2), *station.phase.AXES[1:]]
    coeff = np.zeros((2, 4, 4))
    raw = {}
    signs = np.array(list(product((-1, 1), repeat=2)))
    for i, j in product(range(1, 4), repeat=2):
        effects = np.array(
            [
                np.kron((axes[0] + a * axes[i]) / 2, (axes[0] + b * axes[j]) / 2)
                for a, b in signs
            ]
        )
        rates = np.einsum("kij,pji->pk", effects, state).real
        counts = rng.multinomial(TICK_BUDGET, rates.ravel()).reshape(2, 4)
        raw[f"{i},{j}"] = counts.tolist()
        frequency = counts / TICK_BUDGET
        coeff[:, i, j] = frequency @ np.prod(signs, axis=1)
        if i == 1:
            coeff[:, 0, j] = frequency @ signs[:, 1]
        if j == 1:
            coeff[:, i, 0] = frequency @ signs[:, 0]
        if i == j == 1:
            coeff[:, 0, 0] = frequency.sum(axis=1)
    estimate = sum(
        coeff[:, i, j, None, None] * np.kron(axes[i], axes[j]) / 4
        for i, j in product(range(4), repeat=2)
    )
    assert np.linalg.eigvalsh(estimate).min() > 0
    entropy = station.entropy(estimate)
    # Each coefficient's span is <=2. On this 8-state space, half the
    # trace norm is <=4 times the maximum coefficient error by orthogonality.
    distance = 4 * (2 * MEAN_ERROR + 1e-12)
    entropy_error = distance * np.log(7 / distance) - (1 - distance) * np.log1p(
        -distance
    )
    odd = coeff[1, 0, 0]
    return dict(
        counts=raw,
        entropy=[entropy - entropy_error, entropy + entropy_error],
        odd_energy=[odd - MEAN_ERROR - 1e-12, odd + MEAN_ERROR + 1e-12],
    )


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    bins = clock_effects()
    data = {
        name: {"clocks": clock_read(x, bins, rng), "material": material_read(x, rng)}
        for name, x in preparations().items()
    }
    return dict(
        starts_per_route_or_setting=TICK_BUDGET,
        simultaneous_coverage=1 - TAIL,
        clock_order=list(map(list, zip(*station.ORDER))),
        edge_order=EDGES,
        log_range_weights=weights(),
        acquisitions=data,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
