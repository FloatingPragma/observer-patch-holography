"""Eight-station, inclusive weak flight with unchanged live heads and strain.

The coupled Hamiltonian and four-tick source pump are the registered reciprocal
strain protocol. At each discrete appointment the miss is
exp(-h P_target/2) exp(-ihH), and a click applies sqrt(1-exp(-h)) P_target
after free evolution. The first click switches the target back to the sender.
It retains every internal register. Both possible timeouts pay the cap.

In continuous time the no-click generator is G_b=-iH-P_b/2. The jump from a
to b after time s has density amplitude P_b exp(s G_b) P_a. The reverse
jump after an additional duration u acts on that full post-click amplitude.
The integral over s+u<=T counts every completed trial. The discrete product
formula and finite-dimensional continuity give this monitor as h tends to zero.
The direct generator calculation integrates this density independently.

The first nonzero integer block B=P_b(40H)^mP_a is also the first free
coefficient of the no-click flight since inserting a target projection sooner
requires an earlier arrival. Reverse flight has transposed leading B. Thus
the return amplitude begins with B*B, which is nonzero on a faithful source.
The leading joint probability is proportional to s^(2m) u^(2m), giving its
explicit integral over the deadline triangle in onset(). For all enrolled
station pairs the exact integer calculation finds such an m.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
import json
from math import expm1, factorial, isqrt
from pathlib import Path

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from .. import live_cadence, live_radar, live_survey

HERE = Path(__file__).parent
STATIONS = (0, 1, 38, 39, 5, 7, 40, 42)
DIVISORS = (16, 32)
DEADLINE = 6
KEYS = ("reference", "ordinary", "decorrelated", "hidden", "left", "right")
BINS = ("first_hist", "return_hist")
SCALARS = (
    "capped_range",
    "first_capped_duration",
    "relay",
    "complete",
    "relay_only",
    "no_relay",
)
# Divisions affect appointments or integration, never the native propagation speed.


@lru_cache(maxsize=1)
def states():
    ref, variations, _ = live_survey.retained(True)
    return np.stack([ref] + [variations[k] for k in KEYS[1:]])


@lru_cache(maxsize=2)
def hardware(charge):
    if charge not in (-1, 1):
        raise ValueError("charge must be -1 or +1")
    source = live_radar.source
    board = source.formation.federation(1)
    h = source.system(board, source.formation.region(board), charge, True)[0]
    return h, live_radar.station_rows(board, STATIONS)


@lru_cache(maxsize=2)
def free(charge):
    h, rows = hardware(charge)
    t = np.arange(DEADLINE * max(DIVISORS) + 1) / max(DIVISORS)
    return live_cadence.kernel_grid(h, rows, t)


def receiver_continuum(charge, b):
    """Direct complex-symmetric generator evaluation for one destination bank."""
    h, rows = hardware(charge)
    n = DEADLINE * max(DIVISORS)
    target = rows[b]
    loss = np.zeros(h.shape[0])
    loss[target] = -0.5
    gen = (-1j / 40) * h + diags(loss, format="csr")
    initial = np.zeros((h.shape[0], 8))
    initial[target] = np.eye(8)
    wave = expm_multiply(
        gen, initial, start=0, stop=DEADLINE, num=n + 1, traceA=gen.diagonal().sum()
    )
    # The generator is complex symmetric; T here is a transpose, not adjoint.
    return wave[:, rows, :].swapaxes(-1, -2)


@lru_cache(maxsize=2)
def continuum(charge):
    """Direct no-click jump blocks [t,sender,receiver,new internal,old internal]."""
    return np.stack(
        [receiver_continuum(charge, b) for b in range(len(STATIONS))], axis=2
    )


def discrete_hits(charge, divisor):
    kernels = free(charge)[:: max(DIVISORS) // divisor]
    hits = np.zeros_like(kernels)
    probability = -expm1(-1 / divisor)
    for a, b in permutations(range(8), 2):
        hits[:, a, b] = live_cadence.first_hits(
            kernels[:, b, a], kernels[:, b, b], probability
        )
    return hits


def hermitian_parts(x):
    return np.concatenate([x.real, x.imag], axis=-1).reshape(*x.shape[:-2], -1)


def densities(f, g, rhos, divisor, continuous):
    """Single- and two-click probabilities per appointment or continuous densities."""
    first = np.einsum("tij,zji->zt", live_radar.gram(f), rhos).real
    steered = f[None] @ rhos[:, None] @ f.conj().swapaxes(-1, -2)[None]
    # Directed endpoint trace Tr(E_j sigma_i) and independent time histories.
    # Hermitian real/imag components have their ordinary real inner product.
    nfft = next_fast_len(2 * len(f) - 1, real=True)
    forward = rfft(hermitian_parts(steered), nfft, axis=1)
    backward = rfft(hermitian_parts(live_radar.gram(g)), nfft, axis=0)
    joint = np.einsum("ztd,td->zt", forward, backward)
    done = irfft(joint, nfft, axis=1)[:, : len(f)]
    done[:, :2] = 0  # At t=0 on either leg the off-diagonal endpoint block is zero.
    if continuous:
        # Trapezoid convolution on the two-leg triangle. The zero-duration
        # endpoints vanish exactly for distinct stations.
        done /= divisor
        # Continuous inner integral for t=h is generally nonzero; supply its
        # two-endpoint quadrature value, which is zero on this grid.
    return first, done


def summaries(f, g, rhos, divisor, continuous):
    first, done = densities(f, g, rhos, divisor, continuous)
    t = np.arange(len(f)) / divisor
    trace = np.zeros(len(rhos))
    trace[0] = 1
    step = divisor // 2

    def integrate(x):
        if continuous:
            return np.trapezoid(x, dx=1 / divisor, axis=-1)
        return x.sum(axis=-1)

    def bins(x):
        if continuous:
            return np.stack(
                [integrate(x[:, i : i + step + 1]) for i in range(0, len(t) - 1, step)],
                axis=-1,
            )
        return x[:, 1:].reshape(len(rhos), 12, step).sum(axis=2)

    r, e = integrate(first), integrate(done)
    return dict(
        relay=r,
        complete=e,
        relay_only=r - e,
        no_relay=trace - r,
        first_capped_duration=DEADLINE * trace - integrate((DEADLINE - t) * first),
        capped_range=(DEADLINE * trace - integrate((DEADLINE - t) * done)) / 2,
        first_hist=bins(first),
        return_hist=bins(done),
    )


def measure(divisor, charge, continuous=False):
    if divisor not in DIVISORS:
        raise ValueError("unknown registered divisor")
    amplitudes = (
        continuum(charge)[:: max(DIVISORS) // divisor]
        if continuous
        else discrete_hits(charge, divisor)
    )
    values = {
        k: np.zeros((6, 8, 8, 12) if k in BINS else (6, 8, 8))
        for k in (*SCALARS, *BINS)
    }
    for a, b in permutations(range(8), 2):
        for k, x in summaries(
            amplitudes[:, a, b], amplitudes[:, b, a], states(), divisor, continuous
        ).items():
            values[k][:, a, b] = x
    return dict(
        divisor=divisor,
        charge=charge,
        states={
            name: {k: v[i].tolist() for k, v in values.items()}
            for i, name in enumerate(KEYS)
        },
    )


def pooled(record, divisor, key="capped_range", continuous=False):
    """Fair charge and direction average, [port index,port index,(histogram bin)]."""
    panels = [
        p
        for p in record["continuous_panels" if continuous else "panels"]
        if p["divisor"] == divisor
    ]
    if sorted(p["charge"] for p in panels) != [-1, 1]:
        raise ValueError("expected exactly one panel for each charge")
    result = {}
    for name in KEYS:
        z = sum(np.asarray(p["states"][name][key]) for p in panels)
        result[name] = (z + z.swapaxes(0, 1)) / 4
    return result


def onset(charge):
    """Exact shortest Taylor powers and first continuous cumulative coefficients."""
    h, rows = hardware(charge)
    mtx = h.tocsr()
    row_norm = int(abs(mtx).sum(axis=1).max())
    sites = rows.ravel()
    wave = np.zeros((h.shape[0], len(sites)), dtype=np.int64)
    wave[sites] = np.eye(len(sites), dtype=int)
    order = np.zeros((8, 8), dtype=int)
    fsq = [[0] * 8 for _ in rows]
    values = {
        k: np.zeros((len(KEYS), 8, 8))
        for k in ("first_cumulative", "return_cumulative")
    }
    for degree in range(1, 13):
        if wave.dtype == object or int(np.abs(wave).max()) * row_norm >= 2**62:
            # Integer row reductions for this sparse matrix avoid int64 overflow.
            if np.diff(mtx.indptr).min() <= 0:
                raise ValueError("unexpected zero generator row")
            wave = np.add.reduceat(
                mtx.data[:, None].astype(object) * wave[mtx.indices],
                mtx.indptr[:-1],
                axis=0,
            )
        else:
            wave = mtx @ wave
        for a, b in permutations(range(8), 2):
            if order[a, b]:
                continue
            block = wave[rows[b], 8 * a : 8 * a + 8]
            if not block.any():
                continue
            order[a, b] = degree
            fsq[a][b] = sum(int(z) ** 2 for z in block.ravel())
            lead = np.asarray(block, dtype=float) / (40**degree * factorial(degree))
            intensity = lead.T @ lead
            values["first_cumulative"][:, a, b] = np.einsum(
                "ij,zji->z", intensity, states()
            ).real / (2 * degree + 1)
            triangle = factorial(2 * degree) ** 2 / factorial(4 * degree + 2)
            values["return_cumulative"][:, a, b] = (
                np.einsum("ij,zji->z", intensity @ intensity, states()).real * triangle
            )
        if np.count_nonzero(order) == 56:
            break
    if np.count_nonzero(order) != 56:
        raise ValueError("some port pairs lack a verified leading contact")
    return dict(
        order=order.tolist(),
        integer_block_frobenius_squared=fsq,
        states={
            name: {k: v[i].tolist() for k, v in values.items()}
            for i, name in enumerate(KEYS)
        },
    )


def early_return_floor(order, frobenius_squared):
    """Exact all-time-completion lower bound from one onset and source faithfulness.

    Write the two normalized jump amplitudes as L+E_s and L.T+E_u
    (up to phases). The norm of L is sigma. Each degree-m Taylor remainder
    after division by t^m/m! is at most t*24^(m+1)/(m+1) by contraction.
    On the square [0,a]^2 either error is at most sigma/4. Their product
    differs from L.T L in norm by at most 9 sigma^2/16. Integrating the
    resulting positive density uses only rho >= 0.0014 I, T >= 2a, gamma=1.
    """
    if order < 1 or frobenius_squared < 1:
        raise ValueError("requires an observed nonzero integral contact")
    # An integral 8x8 matrix has norm >= max(1,Frobenius norm/sqrt(8)).
    sigma = Fraction(max(1, isqrt(frobenius_squared // 8)), 40**order)
    a = min(Fraction(DEADLINE, 2), (order + 1) * sigma / (4 * 24 ** (order + 1)))
    normalized = Fraction(7, 16) * sigma**2 / factorial(order) ** 2
    return (
        Fraction(14, 10000)
        * normalized**2
        * (a ** (2 * order + 1) / (2 * order + 1)) ** 2
    )


def run():
    if tuple(live_survey.stations()) != STATIONS:
        raise ValueError("registered station selection changed")
    return dict(
        stations=list(STATIONS),
        calibration_indices=[0, 1, 2, 3, 4],
        held_indices=[5, 6, 7],
        gamma=1,
        deadline=DEADLINE,
        state_order=list(KEYS),
        orientation="first index sends, second relays; diagonal unused and zero",
        derivatives="ordinary-reference, product-reference, normalized XZ+ZX center, left and right marginal-axis rotations",
        histograms="(0,.5],(.5,1],...,(5.5,6]; both timeouts reported separately; first_capped_duration pays T, capped_range pays T/2",
        continuous_method="direct no-click generator; trapezoidal time integrals at the recorded divisor, a numerical quadrature",
        onsets={str(q): onset(q) for q in (-1, 1)},
        panels=[measure(d, q) for d in DIVISORS for q in (-1, 1)],
        continuous_panels=[measure(d, q, True) for d in DIVISORS for q in (-1, 1)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), separators=(",", ":"), allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "flight.json").write_text(data)
    else:
        print(data, end="")
