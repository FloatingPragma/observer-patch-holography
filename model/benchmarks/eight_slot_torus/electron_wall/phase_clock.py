"""Finite Fourier clock on the magnetic walk. Both unitary signs enter the folded reader."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

from .cyclotron import HERE, PROTOCOL, WEAK_PERIODS, enrolled, magnetic_block, magnetic_series, oscillator_sources

SPEC = json.loads((HERE / "phase_clock.json").read_text())


def fourier_probabilities(phase, size):
    """Rectangular counter prepared at times 0,...,size-1 then read in its Fourier basis."""
    detuning = (size * phase / (2 * np.pi) - np.arange(size) + size / 2) % size - size / 2
    probabilities = (np.sinc(detuning) / np.sinc(detuning / size)) ** 2
    if abs(probabilities.sum() - 1) > 1e-9:
        raise ArithmeticError("incomplete Fourier clock")
    return probabilities / probabilities.sum()


def check_bins(archive, size, trials, phases, weights):
    observed = dict(archive)
    predicted_bins = size // 4 - np.rint(size * phases / (2 * np.pi)).astype(int)
    bins = np.union1d(list(observed), predicted_bins)
    read = np.array([observed.get(k, 0) for k in bins])
    raw = size // 4 - bins
    contributions = np.zeros((len(raw), len(phases)))
    for orientation in (1, -1):
        d = (size * phases[None, :] / (2 * np.pi) - orientation * raw[:, None] + size / 2) % size - size / 2
        contributions += (np.sinc(d) / np.sinc(d / size)) ** 2
    contributions[(raw == 0) | (raw == size // 2)] /= 2
    p = np.clip(contributions @ weights, 0., 1.)
    tail = np.minimum(binom.cdf(read, trials, p), binom.sf(read - 1, trials, p))
    if not np.isfinite(tail).all() or tail.min() * 2 * (size // 2 + 1) * 24 < 1e-6:
        raise ArithmeticError("folded frequency probabilities reject a retained bin")


def sample(period, length, spec, archives=()):
    rng = np.random.default_rng(spec["seed"])
    records = {}
    size = spec["quanta_per_period"] * period
    for charge in (1, -1):
        for grade in (1, -1):
            matrix = magnetic_block(period, length, charge, grade).toarray()
            roots, vectors = np.linalg.eigh(matrix)
            for n, probe in enumerate(oscillator_sources(period, length, grade)):
                key = f"{charge}/{grade}/{n}"
                weights = abs(vectors.T @ probe) ** 2
                if abs(weights.sum() - 1.) > 1e-12:
                    raise ArithmeticError("incomplete oscillator spectrum")
                # Sampling this spectral mixture gives exactly the Fourier-register marginal.
                # The negative carrier gives the reflected raw bin and the same folded bin.
                trials = spec["trials_per_probe"]
                cdf = np.cumsum(weights / weights.sum())
                cdf[-1] = 1.
                selected = np.searchsorted(cdf, rng.random(trials))
                quantiles = rng.random(trials)
                counts = Counter()
                for idx in np.unique(selected):
                    phase = np.arccos(roots[idx] / 12.)
                    probs = fourier_probabilities(phase, size)
                    clock_cdf = np.cumsum(probs)
                    clock_cdf[-1] = 1.
                    bins = np.searchsorted(clock_cdf, quantiles[selected == idx])
                    folded = size // 4 - np.minimum(bins, size - bins)
                    counts.update(folded.tolist())
                records[key] = sorted([k, v] for k, v in counts.items())
                for archive in (*archives, records):
                    check_bins(archive[key], size, trials, np.arccos(roots / 12), weights)
    return records


def median_rate(counts, size, trials):
    if any(type(k) is not int or type(v) is not int or v <= 0 or abs(k) > size // 4
           for k, v in counts):
        raise ValueError("invalid phase bin")
    bins = [k for k, _ in counts]
    if bins != sorted(set(bins)) or sum(v for _, v in counts) != trials:
        raise ValueError("invalid phase histogram")
    index = np.searchsorted(np.cumsum([v for _, v in counts]), (trials + 1) // 2)
    return float(2 * np.pi * bins[index] / size)


def ratio_interval(rates, errors, charge):
    spin = charge * (rates[1, 0] - rates[0, 0])
    orbital = np.mean(rates[:, 1] - rates[:, 0])
    ds = errors[:, 0].sum()
    do = errors.sum() / 2
    if min(spin - ds, orbital - do) <= 0:
        raise ArithmeticError("the finite clock cannot resolve the two transitions")
    return [2 * (spin - ds) / (orbital + do), 2 * (spin + ds) / (orbital - do)]


def budget(row, spec):
    size = spec["quanta_per_period"] * row["period"]
    if size < 4 or size & (size - 1):
        raise ValueError("the binary phase counter needs a power of two")
    bins = spec["error_bins"]
    radius = spec["raw_root_radius"]
    rates, errors, log_risks = np.zeros((2, 2)), np.zeros((2, 2)), []
    for j, grade in enumerate(row["grades"]):
        for n, line in enumerate(grade["assigned"]):
            rate = grade["rates"][line]
            raw = 12 * np.sin(rate)
            spread = max(abs(np.arcsin((raw + sign * radius) / 12) - rate) for sign in (1, -1))
            rates[j, n] = rate
            errors[j, n] = 2 * np.pi * bins / size + spread
            # The residual bounds the weight of the approximate eigenvector outside the interval.
            leakage_norm = (grade["residual"] + spec["residual_allowance"]) / radius
            weight = grade["weights"][n][line]
            good = max(0., np.sqrt(weight) - leakage_norm) ** 2 * (1 - 1 / (2 * (bins - 1)))
            if good <= .5:
                raise ArithmeticError("intended line lacks a majority after finite readout")
            log_risks.append(-2 * spec["trials_per_probe"] * (good - .5) ** 2)
    return rates, errors, float(logsumexp(log_risks)), ratio_interval(rates, errors, row["charge"])


def verify(records, spec, rows=None):
    dense = spec["dense_period"]
    expected = {f"{q}/{g}/{n}" for q in (-1, 1) for g in (-1, 1) for n in (0, 1)}
    if set(records) != expected:
        raise ValueError("incomplete clock cohort")
    if rows is None:
        rows = magnetic_series(False) + magnetic_series(True)
    if not enrolled(rows, PROTOCOL["periods"] + WEAK_PERIODS):
        raise ValueError("incomplete magnetic clock series")
    summary = []
    risks = []
    for row in rows:
        period, charge = row["period"], row["charge"]
        size = spec["quanta_per_period"] * period
        rates, errors, log_risk, interval = budget(row, spec)
        risks.append(log_risk)
        if period == dense:
            observed = np.array([[median_rate(records[f"{charge}/{g}/{n}"], size, spec["trials_per_probe"])
                                  for n in (0, 1)] for g in (1, -1)])
            if np.any(abs(observed - rates) > errors):
                raise ArithmeticError("a phase median missed its registered line")
            interval = ratio_interval(observed, errors, charge)
            if interval[0] <= 2:
                raise ArithmeticError("coarse field sign unresolved")
        if period == WEAK_PERIODS[-1] and interval[1] >= 2:
            raise ArithmeticError("fine field sign unresolved")
        summary.append(dict(period=period, charge=charge, bits=size.bit_length() - 1,
                            services_per_trial=size - 1, interval=interval, sampled=period == dense))
    log10_risk = float(logsumexp(risks) / np.log(10))
    if not np.isfinite(log10_risk) or log10_risk >= -8:
        raise ArithmeticError("finite-clock confidence gate failed")
    return dict(log10_family_error_bound=log10_risk, series=summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path, help="write the full folded clock histograms to an unused path")
    args = parser.parse_args()
    rows = magnetic_series(False) + magnetic_series(True)
    archives = []
    sources = json.loads((HERE / "clock_sampler.json").read_text())["sources"]
    if set(sources) != {"clock_counts.json", "clock_linux_counts.json"}:
        raise ValueError("incomplete platform count archives")
    for name, digest in sources.items():
        raw = (HERE / name).read_bytes().replace(b"\r\n", b"\n")
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError("source-bound clock counts changed")
        archive = json.loads(raw)
        verify(archive, SPEC, rows)
        archives.append(archive)
    acquired = sample(SPEC["dense_period"], PROTOCOL["private_length"], SPEC, archives)
    if args.acquire:
        with args.acquire.open("x", encoding="utf-8") as stream:
            json.dump(acquired, stream, separators=(",", ":"), sort_keys=True)
            stream.write("\n")
    else:
        stored = json.loads((HERE / "clock_quantiles.json").read_text())
        verify(stored, SPEC, rows)
        if stored != acquired:
            raise ValueError("clock count replay differs")
    print(json.dumps(verify(acquired, SPEC, rows), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
