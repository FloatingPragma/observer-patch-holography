"""A sharp-word spatial dispersion test on five-owner and six-owner polar rings."""

import argparse
import itertools
import json
from math import asin, ceil, floor, log, pi, sin, sqrt
from pathlib import Path

import numpy as np

from model.experiments.material_inertia.pair_spectrum.neutral_guide.polar_ring import fringes, measure, stock


HERE = Path(__file__).resolve().parent
SHOTS = 2**24
WAITS = [0, *(2**i for i in range(16))]


def acquire():
    rng = np.random.default_rng(30192330)
    rows = []
    for case in stock.charts():
        if case["cells"] != 80 or not case["stock"]:
            continue
        n = len(case["stock"]["ring"])
        addr = np.arange(n)
        _, blocks = measure.spectrum(n, 8)
        amplitudes, weight = {}, {}
        for k, (_, energies, _, origin, principal) in enumerate(blocks):
            weights = n * abs(origin)**2
            amplitudes[k] = np.exp(-1j * np.outer(WAITS, energies)) @ weights
            weight[k] = float(weights[principal])
        for harmonic in (1, 2, 3, 4):
            for direction in (1, -1):
                k = direction * harmonic % n
                phase = np.exp(2j * pi * k * addr / n)
                # Translation covariance factorizes the one-upper word amplitude at x.
                one = amplitudes[0][:, None] + phase * amplitudes[k][:, None]
                probs = abs(one)**2 / (2 * n)
                banks = []
                for law in probs:
                    p = np.r_[law, 1 - law.sum()]
                    if min(p) < -1e-10:
                        raise ArithmeticError("Spatial launch exceeds unit probability")
                    p = np.maximum(0, p).round(13)
                    banks.append(rng.multinomial(SHOTS, p / p.sum()).tolist())
                rows.append(dict(n=n, source=direction*harmonic, banks=banks,
                                 weights=[weight[0], weight[k]]))
    return dict(shots=SHOTS, waits=WAITS, rows=rows)


def angle_error(signal, trials):
    # A one-upper word reports a complex unit phasor and a refusal reports zero.
    # This allowance also covers the +/-1 quadrature used in the rest clock.
    cartesian = sqrt(2 * log(2 * 12000 / .01) / trials)
    radius = sqrt(2) * cartesian
    if abs(signal) <= radius:
        raise ValueError("Unresolved complex phase")
    return asin(radius / abs(signal))


def line_error(w):
    if not .5 < w <= 1:
        raise ValueError("No dominant spectral line")
    return asin((1-w)/w)


def unwrap(interval, times, observations):
    lo, hi = interval
    for t, (z, err) in zip(times, observations, strict=True):
        phase = np.angle(z)
        first = ceil((t*lo - err - phase)/(2*pi))
        last = floor((t*hi + err - phase)/(2*pi))
        if first != last:
            raise ValueError("Clock permits zero or multiple phase windings")
        center = (phase + 2*pi*first) / t
        lo, hi = max(lo, center - err/t), min(hi, center + err/t)
    return [lo, hi]


def rest_cost():
    src = Path(measure.__file__).parent
    f = json.loads((src / "fringes.json").read_text())
    old = json.loads((src / "records.json").read_text())
    row = next(r for r in f["panels"] if (r["cells"], r["n"], r["kappa"]) == (80, 10, 8))
    stock_row = next(r for r in old["panels"] if (r["cells"], r["dimer_count"], r["kappa"]) == (80, 10, 8))
    systematic = line_error(row["ground_line_weight"]) + line_error(stock_row["clocks"][0]["one_island_line_weight"])
    observed = []
    for host, excited in zip(row["sharp_lower_return"], row["clocks"][0]["plus_counts"], strict=True):
        h, e = fringes.amplitude(host), fringes.amplitude(excited)
        observed.append((h * e.conjugate(), systematic + angle_error(h, measure.SHOTS) + angle_error(e, measure.SHOTS)))
    return unwrap((0, 50), measure.HOLDS, observed)


def read(raw):
    mass = rest_cost()
    intervals = {}
    if raw["shots"] != SHOTS or raw["waits"] != WAITS:
        raise ValueError("Unrecognized sharp-word schedule")
    for row in raw["rows"]:
        n, source = row["n"], row["source"]
        counts = np.asarray(row["banks"])
        if (counts.shape != (len(WAITS), n+1) or counts.dtype.kind not in "iu"
                or np.any(counts < 0) or np.any(counts.sum(axis=1) != SHOTS)):
            raise ValueError("Incomplete spatial banks")
        harmonics = np.fft.fft(counts[:, :n] / SHOTS, axis=1)
        k = 1 + int(np.argmax(abs(harmonics[0, 1:n//2+1])))
        strengths = sorted(abs(harmonics[0, 1:n//2+1]))
        if strengths[-1] - strengths[-2] < .01 or k != abs(source):
            raise ValueError("Initial word counts do not resolve the launch wavelength")
        signal = harmonics[:, k].conj() if source > 0 else harmonics[:, k]
        systematic = sum(map(line_error, row["weights"]))
        observations = [(z, systematic + angle_error(z, SHOTS)) for z in signal[1:]]
        key = (n, source)
        if key in intervals:
            raise ValueError("Duplicate spatial test")
        # The previously measured narrow island band supplies only the first
        # winding interval. Its size is much wider than the resolved recoils.
        intervals[key] = unwrap((-1, 1), WAITS[1:], observations)
    if set(intervals) != {(n,k) for n in (10,12) for k in (-4,-3,-2,-1,1,2,3,4)}:
        raise ValueError("Incomplete set of spatial tests")
    results = []
    for kind in ("continuum", "lattice", "hopping"):
        misses = []
        for (n,k), (lo,hi) in intervals.items():
            if (n,k) == (10,1):
                continue
            r = (10*k/n)**2 if kind == "continuum" else (sin(pi*k/n) / sin(pi/10))**2
            predictions = []
            for m, d in itertools.product(mass, intervals[10,1]):
                s = d * (2*m+d) * r
                predictions.append(d*r if kind == "hopping" else s/(sqrt(m*m+s)+m))
            misses.append(max(0, min(predictions)-hi, lo-max(predictions)))
        results.append(dict(law=kind, passed=misses.count(0), tests=len(misses), max_excluded_gap=max(misses)))
    return dict(rest=mass, recoil={f"{n}/{k}":v for (n,k),v in intervals.items()}, forecasts=results)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = acquire()
    summary = read(data)
    if args.output:
        with args.output.open("x") as dest:
            json.dump(data, dest, indent=1, allow_nan=False)
            dest.write("\n")
    elif data != json.loads((HERE / "counts.json").read_text()):
        raise ArithmeticError("Polar ring acquisition differs from the archive")
    print(json.dumps(summary, indent=2))
