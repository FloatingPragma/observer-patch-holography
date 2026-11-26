"""Two-observer return clocks, spectral radial scales and held flat-distance tests."""

from math import sqrt
from pathlib import Path
import json
import numpy as np
from scipy.stats import beta
from .. import measure as bearings
from ... import measure as tracer
from . import dynamics

HERE = Path(__file__).resolve().parent
DWELLS = (1 / 512, 1 / 256)
SHOTS = 2**34
STANDARD_SHOTS = 2**44
TAIL = 0.01 / (2 * 1800)


def binomial_band(n, size):
    n = np.asarray(n)
    lo = np.where(n == 0, 0.0, beta.ppf(TAIL, np.maximum(1, n), size - n + 1))
    hi = np.where(n == size, 1.0, beta.ppf(1 - TAIL, n + 1, np.maximum(1, size - n)))
    return lo, hi


def geometry(lo, hi, h, edges):
    # The two independent orientations bracket their symmetrized material transit time.
    lo = (lo + lo.T) / 2
    hi = (hi + hi.T) / 2
    d = (lo + hi) / 2
    error = (hi - lo) / 2
    j = np.eye(12) - np.ones((12, 12)) / 12
    b = -j @ d @ j
    b0 = -j @ h @ j
    values, vectors = np.linalg.eigh(b)
    ideal = np.linalg.eigvalsh(b0)
    gap = ideal[-3] - ideal[-4]
    noise = float(max(error.sum(1)))
    eta = float(np.linalg.norm(b - b0, 2)) + noise
    assert 2 * eta < gap
    _, g = bearings.bearing_frame(edges)
    projection = g / 4
    slow = vectors[:, -3:]
    # Davis-Kahan residual bound relative to the known ideal spectral gap.
    angle = float(np.linalg.norm((np.eye(12) - projection) @ slow, 2))
    flat = -j @ (d * d) @ j / 2
    flat_noise = float(max((error * (2 * abs(d) + error)).sum(1))) / 2
    flats = np.linalg.eigvalsh(flat)
    eig4 = flats[-4]
    distortion = sqrt(1 + 2 * (eig4 - flat_noise) / max((hi * hi).sum(1))) - 1
    rms = float(np.max(abs(d - h)))
    class_bounds = []
    for length in (11, 14, 15):
        mask = h == length
        class_bounds.append([float(lo[mask].min()), float(hi[mask].max())])
    return dict(
        classes=class_bounds,
        max_central_difference_from_limiting_hops=rms,
        centered_transit_operator_error_upper=eta,
        measured_slow_eigenvalues=values[-3:].tolist(),
        slow_spectrum_99_enclosure=[
            float(values[-3] - noise),
            float(values[-1] + noise),
        ],
        fourth_eigenvalue_upper=float(values[-4] + noise),
        nominal_slow_to_bearing_angle_sine=angle,
        population_slow_to_bearing_angle_sine_upper=eta / (gap - eta),
        squared_range_fourth_gram_eigenvalue_lower=float(eig4 - flat_noise),
        squared_range_smallest_centered_eigenvalue_lower=float(flats[1] - flat_noise),
        minimum_relative_distance_error_for_flat_3d=float(distortion),
    )


def count_bank(failure, basis, tag, stock, lam, opts, cal, error, rng):
    gates = []
    L = np.zeros((12, 12))
    U = np.zeros_like(L)
    qs = 5 * np.array(cal)
    for a in range(12):
        source = np.array([location == a and w in stock for location, w in basis])
        admit = sum(source) / len(stock)
        assert admit > 0
        failure_cond = failure[source].mean(0)
        counts = []
        for b, p in enumerate(failure_cond):
            if a == b:
                counts.append(None)
                continue
            # The sender counts all blocked starts. Only a remote detection
            # followed by a home detection before expiry counts as successful.
            raw = rng.multinomial(SHOTS, [admit * p, admit * (1 - p), 1 - admit])
            counts.append(raw.tolist())
            lower, upper = binomial_band(raw[0], sum(raw[:2]))
            L[a, b] = qs[0] * max(0, lower - error) / (2 * lam)
            U[a, b] = qs[1] * min(1, upper + error) / (2 * lam) / (1 - 2 * lam * opts)
        gates.append(counts)
    return gates, L, U


def gravity_free(twice, basis, tag, h, edges, outside):
    stocks = dict(
        uniform=sorted({w for _, w in basis if outside not in w}),
        split=tracer.clock_ball.SPLIT,
        spread=tracer.clock_ball.SPREAD,
    )
    rng = np.random.default_rng(391106)
    reports = []
    for t in DWELLS:
        u, werr = dynamics.wave(twice, t)
        cert = dynamics.certificate(u, werr, tag, h, t)
        rate = cert["jump_rate"]
        bound = 15 / (rate * (1 - cert["conditional_drift_relative_upper"]))
        duration = 1 / (2**14) * t * t
        # Read a separate isolated material load exchange at the same shutter.
        q = abs(tracer._dimer_pulse(t)[1, 0]) ** 2
        assert abs(5 * q - rate) < 1e-13
        clock_hits = int(rng.binomial(STANDARD_SHOTS, float(q)))
        clocks = binomial_band(clock_hits, STANDARD_SHOTS)
        kernel = abs(u) ** 2
        failure, ferr = dynamics.returns(
            kernel, tag, duration, cert["column_probability_error_upper"], bound
        )
        assert 0 <= failure.min() and failure.max() < 0.01 and ferr < 1e-8
        banks = {}
        for name, stock in stocks.items():
            raw, lo, hi = count_bank(
                failure, basis, tag, stock, duration, bound, clocks, ferr, rng
            )
            bank = geometry(lo, hi, h, edges)
            bank["rows_expired_completed_blocked"] = raw
            banks[name] = bank
        reports.append(
            dict(
                dwell=t,
                deadline_probability=duration,
                independently_measured_clock_count=clock_hits,
                clock_transfer_probability_interval=list(map(float, clocks)),
                uniform_first_passage_upper_ticks=bound,
                two_leg_deadline_bias_fraction_upper=2 * duration * bound,
                expired_probability_error_upper=ferr,
                native=cert,
                banks=banks,
            )
        )
    # Uniform Bellman bounds for elapsed oscillator time apply even when a new
    # independent stock is chosen at the other shutter. The range class cancels.
    first, second = [r["native"] for r in reports]
    ratio = DWELLS[0] * second["jump_rate"] / (DWELLS[1] * first["jump_rate"])
    interval = [
        ratio
        * (1 - second["conditional_drift_relative_upper"])
        / (1 + first["conditional_drift_relative_upper"]),
        ratio
        * (1 + second["conditional_drift_relative_upper"])
        / (1 - first["conditional_drift_relative_upper"]),
    ]
    return reports, interval


def run():
    h, basis, tag, times, edges, outside = dynamics.setup()
    shots, clock_ratio = gravity_free(h, basis, tag, times, edges, outside)
    # Change the seam contrast but retain its complete local word alphabet and routing.
    control, order, _, _ = tracer.source(0)
    assert order == basis
    u, e = dynamics.wave(control, DWELLS[-1])
    bound = dynamics.certificate(u, e, tag, times, DWELLS[-1])
    return dict(
        source_words=55,
        shutter_dwells=DWELLS,
        shots_per_ordered_pair=SHOTS,
        standard_starts=STANDARD_SHOTS,
        confidence=0.99,
        reference_jump_counts=[11, 14, 15],
        ideal_slow_eigenvalue=12 / (1 - 1 / sqrt(5)),
        exposures=shots,
        all_resident_states_elapsed_time_ratio=clock_ratio,
        zero_contrast_relative_drift_upper=bound["conditional_drift_relative_upper"],
    )


if __name__ == "__main__":
    HERE.joinpath("results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
