"""Fuel interferometry, a clock read from its ports, and a reader that sees counts only."""

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.stats import binom

from model.candidate_laws.scalar_channel.carrier import require


def return_amplitude(H, times, start=None):
    values, vectors = np.linalg.eigh(H)
    weights = abs(vectors[0])**2 if start is None else abs(vectors.conj().T @ start)**2
    return np.exp(-1j*np.outer(times, values)) @ weights


def clock_amplitude(times_in_dimer_units):
    return return_amplitude(np.array([[0.,1.],[1.,0.]]), times_in_dimer_units,
                            np.array([1.,-1.])/np.sqrt(2))


def fuel_law(amplitude, quarters):
    """Full-fuel probability with the source half occupied, for both initial fuels."""
    phase = np.exp(.5j*np.pi*np.asarray(quarters))
    transfer = abs((phase[None, :] + amplitude[:, None])/2)**2
    return np.stack([.5*(1-transfer), .5*(1+transfer)], axis=-1)


def binomial_counts(rng, shots, probability):
    # One independent open-interval quantile per setting. No parameter-dependent RNG consumption.
    q = (rng.integers(2**52, size=np.shape(probability))+.5) / 2**52
    return np.asarray(binom.ppf(q, shots, probability), dtype=np.int32)


def acquire(rng, probabilities, shots):
    require(np.isfinite(probabilities).all() and probabilities.min() >= -1e-12 and
            probabilities.max() <= 1+1e-12, 'invalid fuel probability')
    # Exact half and endpoint probabilities must take the same RNG branch across BLAS implementations.
    measured = np.rint(np.clip(probabilities, 0, 1)*2**30) / 2**30
    return binomial_counts(rng, shots, measured)


def count_partition(rng, law, shots):
    require(np.isfinite(law).all() and law.min() >= -1e-13 and abs(law.sum()-1) < 1e-11,
            'invalid outcome partition')
    # Quantizing cumulative edges keeps every bin nonnegative and preserves total probability exactly.
    law = np.maximum(law, 0.)
    edges = np.rint(np.cumsum(law/law.sum())*2**30).astype(np.int64)
    parts = np.diff(np.r_[0, edges])
    counts = []
    remaining, total = shots, 2**30
    for part in parts[:-1]:
        count = int(binomial_counts(rng, remaining, part/total)) if total and remaining else 0
        counts.append(count)
        remaining -= count
        total -= part
    return [*counts, remaining]


def quadrature(counts, shots):
    delta = (counts[..., 1].astype(float) - counts[..., 0]) / shots
    return delta[:, 0] - delta[:, 1]


def clock_read(counts, shots):
    require(counts.ndim == 3 and counts.shape[1:] == (4, 2) and len(counts) >= 2
            and counts.dtype.kind in 'iu' and np.all((0 <= counts) & (counts <= shots)),
            'invalid clock counts')
    real = quadrature(counts[:, :2], shots)
    imag = quadrature(counts[:, 2:], shots)
    require(np.max(abs(np.hypot(real, imag)-1)) < .05, 'clock lost its reference contrast')
    phase = np.unwrap(np.angle(real + 1j*imag))
    ticks = phase - phase[0]
    require(np.min(np.diff(ticks)) > 0 and np.max(np.diff(ticks)) < np.pi/2,
            'clock is unresolved or has lost a turn on the declared branch')
    return ticks


def read_lines(times, record, spec):
    """Successive frequencies from the residual, followed by a joint fit. Reserve the tail for prediction."""
    count = spec['fit_samples']
    x, y = times[:count], record[:count]
    dt = (x[-1]-x[0])/(count-1)
    grid = np.linspace(x[0], x[-1], count)
    window = np.hanning(count)
    freq_grid = 2*np.pi*np.fft.rfftfreq(8*count, dt)
    frequencies = np.empty(0)
    weights = np.empty(0)
    residual = y.copy()
    dc = 0.
    for _ in range(spec['maximum_lines']+1):
        centered = residual - residual.mean()
        power = 2*abs(np.fft.rfft(np.interp(grid, x, centered)*window, 8*count))/window.sum()
        peaks, _ = find_peaks(power)
        if not peaks.size or power[peaks].max() < spec['line_floor']:
            dc += residual.mean()
            break
        require(len(frequencies) < spec['maximum_lines'], 'unresolved spectrum exceeds reader capacity')
        candidate = freq_grid[peaks[np.argmax(power[peaks])]]
        n = len(frequencies)+1
        guess = np.r_[frequencies, candidate, weights, power[peaks].max(), dc]

        def error(v):
            return np.cos(np.outer(x, v[:n])) @ v[n:2*n] + v[-1] - y

        def derivative(v):
            angles = np.outer(x, v[:n])
            return np.column_stack([-x[:, None]*np.sin(angles)*v[n:2*n], np.cos(angles), np.ones(len(x))])

        fit = least_squares(error, guess, jac=derivative, max_nfev=120)
        require(fit.success, 'frequency fit failed')
        frequencies, weights, dc = fit.x[:n], fit.x[n:2*n], fit.x[-1]
        residual = -fit.fun
    prediction = np.cos(np.outer(times[count:], frequencies)) @ weights + dc
    rms = float(np.sqrt(np.mean((prediction-record[count:])**2)))
    require(rms < spec['maximum_test_rms'], 'line fit misses the held-out fuel records')
    order = np.argsort(frequencies)
    require(np.all(weights > 0) and np.all(frequencies > 0) and abs(sum(weights)+dc-1) < .015,
            'incomplete or negative return measure')
    return dict(lines=frequencies[order].tolist(), weights=weights[order].tolist(), dc=float(dc), heldout_rms=rms)
