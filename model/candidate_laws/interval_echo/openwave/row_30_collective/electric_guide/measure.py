"""Counted local returns and measured bulk standing waves in the neutral guide."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

from model.executable_model.carrier import federation
from ...row_20_electric import experiment as electric
from . import release


HERE = Path(__file__).resolve().parent
SAMPLES = 65536
DT = 1/16
SHOTS = 2**20
SPATIAL_SHOTS = 2**23
BENCH = 2.**np.arange(-6, 11)
SCHEDULE = np.sort(np.random.default_rng(203020).choice(SAMPLES, 2048, replace=False))
CASES = ((44, 8, True), (44, 0, True), (44, 8, False), (42, 8, True))


def signs(mean, shots, rng):
    """Number of positive signs out of the specified complete binary population."""
    p = (1 + np.asarray(mean))/2
    if np.min(p) < -1e-10 or np.max(p) > 1+1e-10:
        raise ArithmeticError("Invalid sign expectation")
    return rng.binomial(shots, np.round(np.clip(p, 0, 1), 12)).astype(np.int32)


def quadratures(phases, shots, rng):
    return signs([np.cos(phases), np.sin(phases)], shots, rng)


def modes(n, j, h):
    # In the staggered orientation gauge B_i has a free left boundary at site −1.
    c = 2*(h*np.eye(n) - j*np.eye(n, k=-1))
    _, omega, vh = np.linalg.svd(c)
    return omega, vh.T


def acquire():
    board = federation(2)
    _, path, _ = electric.guide(board)
    longest_even = 2*((len(path)-2)//4)
    if {n for n, _, _ in CASES} != {longest_even, longest_even-2}:
        raise ValueError("Unexpected route from static source")
    rng = np.random.default_rng(203021)
    t = np.arange(SAMPLES)*DT
    raw = dict(protocol=np.frombuffer(hashlib.sha256((HERE/"PROTOCOL.md").read_bytes()).digest(), np.uint8),
               route=np.array(path, dtype=np.int32), spatial_holds=SCHEDULE.astype(np.int32),
               clock=quadratures(t, SHOTS, rng))
    for n, k, moving in CASES:
        j, h = electric.couplings(k)
        if not moving:
            h = 0.
        key = f"{n}_{k}_{int(moving)}"
        raw[key+"_bench"] = np.stack([quadratures(BENCH*freq, SPATIAL_SHOTS, rng)
                                       for freq in (2*j, 2*h)])
        omega, basis = modes(n, j, h)
        cosines = np.cos(np.outer(t, omega))
        raw[key+"_return"] = signs(cosines.mean(axis=1), SHOTS, rng)
        response = (cosines[SCHEDULE] * basis[0]) @ basis.T
        raw[key+"_spatial"] = signs(response, SPATIAL_SHOTS, rng)
    raw["release"] = release.counts(board, path, rng, SHOTS)
    return raw


def read_sign(data, shape, trials):
    values = np.asarray(data)
    if values.shape != shape or values.dtype.kind not in "ui" or np.min(values)<0 or np.max(values)>trials:
        raise ValueError("Incomplete sign counter")
    return 2*values/trials - 1


def bench_frequency(data):
    x, y = read_sign(data, (2, len(BENCH)), SPATIAL_SHOTS)
    angles = np.arctan2(y, x)
    phase = angles[0]
    for a in angles[1:]:
        nxt = a + 2*np.pi*round((2*phase-a)/(2*np.pi))
        if abs(nxt-2*phase) > .02:
            raise ValueError("Unstable bench clock")
        phase = nxt
    return phase/BENCH[-1]


def frequencies(trace, step, n):
    window = np.hanning(SAMPLES)
    a = window*(trace - trace.mean())
    spectrum = 2*abs(np.fft.rfft(a))/window.sum()
    indices = find_peaks(spectrum, height=.15/n)[0]
    times = step*np.arange(SAMPLES)
    width = 2*np.pi/(step*SAMPLES)
    observed = []
    for i in indices:
        # Optimize within a bin in scaled coordinates to retain absolute precision at high frequency.
        shifted = a*np.exp(1j*i*width*times)
        fit = minimize_scalar(lambda u: -abs(shifted @ np.exp(1j*u*width*times)),
                              bounds=(-1., 1.), method="bounded", options={"xatol": 1e-8})
        observed.append(width*(i + fit.x))
    return np.array(observed)


def wavelength(profile):
    x = np.arange(1, len(profile)+1)
    train = x % 3 != 0

    def fit(k):
        f = np.sin(k*x)
        amp = (f[train] @ profile[train])/(f[train] @ f[train])
        return amp, profile-amp*f

    grid = np.linspace(0, np.pi, 513)[1:-1]
    f = np.sin(np.outer(grid, x[train]))
    trained = profile[train]
    scores = (f @ trained)**2/(f*f).sum(axis=1)
    k0 = grid[np.argmax(scores)]
    delta = np.pi/512
    best = minimize_scalar(lambda k: np.linalg.norm(fit(k)[1][train]),
                           bounds=(k0-delta, k0+delta), method="bounded",
                           options={"xatol": 1e-12})
    amp, residual = fit(best.x)
    held = np.sqrt(np.mean(residual[~train]**2))/max(abs(amp), 1e-20)
    return float(best.x), float(held)


def analyze(raw):
    p = hashlib.sha256((HERE/"PROTOCOL.md").read_bytes()).digest()
    if raw["protocol"].tobytes()!=p or not np.array_equal(raw["spatial_holds"], SCHEDULE):
        raise ValueError("Unexpected protocol or spatial appointments")
    _, path, _ = electric.guide(federation(2))
    if not np.array_equal(raw["route"], path):
        raise ValueError("Changed charged-response route")
    x, y = read_sign(raw["clock"], (2, SAMPLES), SHOTS)
    phase = np.unwrap(np.arctan2(y, x))
    if min(np.diff(phase)) <= 0 or max(np.diff(phase)) >= np.pi/2:
        raise ValueError("Unresolved timing increments")
    step, origin = np.polyfit(np.arange(SAMPLES), phase, 1)
    cases = []
    for n, k, moving in CASES:
        key = f"{n}_{k}_{int(moving)}"
        jf, hf = [bench_frequency(panel)/(step/DT) for panel in raw[key+"_bench"]]
        mass, c2 = abs(jf-hf), jf*hf
        trace = read_sign(raw[key+"_return"], (SAMPLES,), SHOTS)
        spatial = read_sign(raw[key+"_spatial"], (len(SCHEDULE), n), SPATIAL_SHOTS)
        omega = frequencies(trace, step, n)
        case = dict(n=n, kappa=k, hopping=moving, twice_j=jf, twice_h=hf, m=mass, c2=c2,
                    dc=float(trace.mean()), lines=omega.tolist(), complete_bulk_count=len(omega) in (n-1, n))
        cases.append(case)
        phases = np.outer(SCHEDULE*step, omega)
        bank = np.column_stack([np.ones(len(phases)), np.cos(phases), np.sin(phases)])
        fits = np.linalg.lstsq(bank, spatial, rcond=None)[0]
        rows = []
        for w, shape in zip(omega, fits[1:len(omega)+1], strict=True):
            q, miss = wavelength(shape)
            s = 2*(1-np.cos(q))
            scores = dict(lattice=np.sqrt(mass**2+c2*s), continuum=np.sqrt(mass**2+c2*q*q),
                          hopping=mass+c2*s/(2*mass), parabola=mass+c2*q*q/(2*mass))
            rows.append(dict(omega=float(w), q=q, held_rms=miss, shape_pass=miss<.05,
                             **{name:float(w-pred) for name,pred in scores.items()}))
        case["rows"] = rows
        if len(rows) >= 3:
            samples = sorted((r for r in rows if r["shape_pass"]), key=lambda r: r["q"])
            if len(samples) >= 3:
                a, b = samples[0], samples[-1]
                slope = (b["omega"]-a["omega"])/(np.cos(b["q"])-np.cos(a["q"]))
                for r in samples[1:-1]:
                    pred = a["omega"]+slope*(np.cos(r["q"])-np.cos(a["q"]))
                    r["refit_hopping"] = float(r["omega"]-pred)
    return dict(clock_step=float(step), clock_max_residual=float(max(abs(phase-origin-step*np.arange(SAMPLES)))),
                cases=cases, release=release.read(raw["release"], SHOTS))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire", type=Path)
    cli.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    if args.acquire and args.counts_only:
        cli.error("Choose acquisition or count replay")
    stored = HERE/"counts.npz"
    if args.counts_only:
        with np.load(stored, allow_pickle=False) as source:
            results = analyze(source)
    else:
        raw = acquire()
        results = analyze(raw)
        if args.acquire:
            with args.acquire.open("xb") as out:
                np.savez_compressed(out, **raw)
        else:
            with np.load(stored, allow_pickle=False) as source:
                if set(source) != set(raw) or any(not np.array_equal(source[k], v) for k,v in raw.items()):
                    raise ArithmeticError("Stored count banks differ")
    print(json.dumps(results, indent=1, allow_nan=False))


if __name__ == "__main__":
    main()
