"""Acquire edge load/current records from the native collision; decode integer records separately."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.sparse.linalg import expm_multiply

from ...row_08_annihilation import experiment, walls
from ...row_08_annihilation.daughter_spectrum.measure import check, reference

HERE = Path(__file__).resolve().parent
N, SHOTS, SEED = 2048, 2**18, 310394
TIME = np.arange(N) / 32
SIGNS = (-1, 1)
IDENTITY = np.eye(8)
CURRENT = np.kron(np.eye(4), [[0, -1j], [1j, 0]])
PROJECTORS = [(IDENTITY + a * CURRENT) / 2 for a in SIGNS]


def sample(probabilities, rng):
    """Sample exhaustive integer outcomes; errors in a complementary outcome are checked."""
    p = np.asarray(probabilities)
    check(np.min(p) > -2e-10 and np.max(abs(np.sum(p, axis=-1)-1)) < 2e-9,
          "incomplete/negative instrument")
    # Fixed-point probabilities make regeneration insensitive to eigensolver roundoff.
    weights = np.rint(np.maximum(p[:, :-1], 0)*2**30).astype(np.int64)
    remaining = np.full(len(p), 2**30, dtype=np.int64)
    left = np.full(len(p), SHOTS, dtype=np.int64)
    counts = np.zeros(p.shape, dtype=np.int32)
    for i, weight in enumerate(weights.T):
        check(np.all(weight <= remaining), "invalid probability complement")
        fraction = np.divide(weight, remaining, out=np.zeros(len(p)), where=remaining > 0)
        counts[:, i] = rng.binomial(left, fraction)
        left -= counts[:, i]
        remaining -= weight
    counts[:, -1] = left
    return counts


def collision():
    boundary, initial, ordered = reference()
    h = sparse.csr_matrix(walls.operators(8, False, 3))
    step = -32j*(h - sparse.eye(len(initial))*np.mean(h.diagonal()))
    return int(boundary[0]), expm_multiply(step, initial).reshape(16, 8, 8, 2)[ordered]


def instruments(states, native, flight):
    """Each native herald acts on the wing and keeps both outcomes."""
    _, v = eigh(native)
    one = v[:, 1:4] @ v[:, 1:4].T
    tags = [IDENTITY-one, one]
    levels, frame = eigh(flight)
    current = frame.T @ CURRENT @ frame
    evolved = np.exp(1j*TIME[:, None, None]*(levels[:, None]-levels[None, :]))*current

    def signs(rho):
        weight = np.trace(rho).real
        mean = np.einsum("tij,ji->t", evolved, frame.T @ rho @ frame).real
        return np.column_stack([(weight + b*mean)/2 for b in SIGNS])

    silent = np.zeros((N, 2, 2))
    started = np.zeros((N, 2, 2, 2))
    disturbance = np.zeros((1, 2, 2, 2))
    rho = states @ states.conj().T
    for tag, q in enumerate(tags):
        r = q @ rho @ q
        silent[:, tag] = signs(r)
        for a, proj in enumerate(PROJECTORS):
            selective = proj @ r @ proj
            started[:, tag, a] = signs(selective)
            disturbance[0, tag, a] = [np.trace(t @ selective).real for t in tags]
    return silent, started, disturbance


def acquire():
    rng = np.random.default_rng(SEED)
    boundary, waves = collision()
    archive = {}
    for src in range(2):
        miss = 1 - np.linalg.norm(waves[:, :, :, src])**2
        for k in (8, 4):
            readings = []
            for phase in range(2):
                pin = boundary if phase == 0 else 2-boundary
                readings.append(instruments(waves[phase, :, :, src], experiment.wing(8, 3, pin),
                                            experiment.wing(k, 3, pin)))
            for i, name in enumerate(("silent", "started", "disturbance")):
                if k == 4 and name == "disturbance":
                    continue
                q = np.stack([phase[i] for phase in readings], axis=1)
                p = np.column_stack([q.reshape(q.shape[0], -1), np.full(q.shape[0], miss)])
                archive[f"{name}_{src}_{k}"] = sample(p, rng)
    return archive


def current_trace(counts, started, tag=None):
    dims = (len(counts), 2, 2, 2, 2) if started else (len(counts), 2, 2, 2)
    observed = counts[:, :-1].reshape(dims).sum(axis=1)
    selected = observed.sum(axis=1) if tag is None else observed[:, tag]
    weights = selected.reshape(len(counts), -1).sum(axis=1)
    check(np.min(weights) > 0, "empty accepted sample")
    signed = (selected[:, 0, 0] + selected[:, 1, 1] - selected[:, 0, 1] - selected[:, 1, 0]
              if started else selected[:, 1] - selected[:, 0])
    return signed/weights, weights


def spectrum(signal):
    """Freeze discovery on even appointments. Odd appointments enter only the prediction test."""
    time, train = TIME[::2], signal[::2]
    spectral = abs(np.fft.rfft((train-train.mean())*np.hanning(len(train))))
    freq = 2*np.pi*np.fft.rfftfreq(len(train), time[1]-time[0])
    peaks = find_peaks(spectral)[0]
    peaks = peaks[(freq[peaks] >= 4) & (freq[peaks] <= 45)]
    check(len(peaks) >= 3, "fewer than three starting peaks")
    seeds = np.sort(freq[peaks[np.argsort(spectral[peaks])[-3:]]])
    def prediction(t, x):
        return np.cos(np.outer(t, x[:3])) @ x[3:]
    result = least_squares(lambda x: prediction(time, x)-train, np.r_[seeds, [.3,.4,.3]],
                           bounds=(np.r_[seeds-.15, [-2]*3], np.r_[seeds+.15, [2]*3]),
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    f = result.x
    held = signal[1::2]
    residual = prediction(TIME[1::2], f)-held
    # The best of the three discovered individual lines has its own fitted amplitude and phase.
    rms = []
    for nu in f[:3]:
        basis = np.c_[np.cos(time*nu), np.sin(time*nu), np.ones_like(time)]
        coeff = np.linalg.lstsq(basis, train, rcond=None)[0]
        odd = TIME[1::2]*nu
        rms.append(np.mean((np.c_[np.cos(odd), np.sin(odd), np.ones_like(odd)] @ coeff-held)**2)**.5)
    return dict(frequencies=f[:3].tolist(), amplitudes=f[3:].tolist(),
                held_rmse=float(np.mean(residual**2)**.5), best_single_discovered_rmse=float(min(rms)))


def decode(archive):
    expected = {f"{name}_{src}_{k}" for src in (0,1) for k in (8,4)
                for name in ("silent","started","disturbance") if k == 8 or name != "disturbance"}
    check(set(archive) == expected, "missing apparatus panel")
    output = {}
    for name, arr in archive.items():
        setting, src, k = name.split("_")
        shape = {"started":(N,17), "silent":(N,9), "disturbance":(1,17)}[setting]
        check(arr.shape == shape and np.issubdtype(arr.dtype,np.integer) and
              np.min(arr) >= 0 and np.all(arr.sum(axis=1)==SHOTS), "bad count archive")
        if setting == "disturbance":
            cells = arr[0,:-1].reshape(2,2,2,2).sum(axis=0)
            eligible = cells[1].sum(axis=0)
            output[name] = dict(one_tag=int(eligible.sum()), repeat_one=int(eligible[1]))
        elif setting == "silent":
            if src == "1":
                signal, n = current_trace(arr, False, 1)
                # Hoeffding conditional on all herald sample sizes, joint 99% across both laws/times.
                bound = np.sqrt(2*np.log(4*N/.01)/min(n))
                output[name] = dict(herald_fraction=float(n.sum()/(N*SHOTS)),
                                    peak_abs_current=float(max(abs(signal))),
                                    simultaneous_99_radius=float(bound))
        else:
            # The launched wall has a substantial one-quantum branch. The sham's *entire*
            # ordered wing is used as the control; its rare one-quantum tags are retained.
            signal, n = current_trace(arr, True, 1 if src == "1" else None)
            fit = spectrum(signal)
            output[name] = dict(**fit, retained_fraction=float(n.sum()/(N*SHOTS)))
    # Counts-only positive sensitivity / zero control, same search as actual records.
    rng = np.random.default_rng(340394)
    nu = [7.1, 14.23, 19.55]
    for name, wave in [("planted",np.cos(np.outer(TIME,nu)) @ [.25,.5,.25]),
                       ("null",np.zeros(N))]:
        output[name] = spectrum(2*rng.binomial(SHOTS,(1+wave)/2)/SHOTS - 1)
    return output


def validate(result):
    registered_path = HERE.parents[1]/"row_08_annihilation"/"daughter_spectrum"/"results.json"
    old = json.loads(registered_path.read_text())["calibrations"][0]["levels"][1:4]
    for k in (8,4):
        sham, arrival = [result[f"started_{s}_{k}"] for s in (0,1)]
        for panel in (sham, arrival):
            check(np.max(abs(np.array(panel["amplitudes"])-[.25,.5,.25])) < .01 and
                  panel["held_rmse"] < .015 and panel["best_single_discovered_rmse"] > .2,
                  "held appointments insufficient to support a three-line correlation")
        check(max(abs(np.array(sham["frequencies"])-arrival["frequencies"]))<.005,
              "source dependence")
        if k == 8:
            check(max(abs(np.asarray(old)-arrival["frequencies"]))<.005, "prior calibration mismatch")
    check(max(result["started_1_4"]["frequencies"]) < min(old)-10, "law intervention unresolved")
    silent = result["silent_1_8"]
    check(silent["peak_abs_current"] < silent["simultaneous_99_radius"], "silent current resolved")
    for case in ("disturbance_1_8",):
        s = result[case]
        check(s["one_tag"] > 100 and abs(s["repeat_one"]/s["one_tag"]-.5)<.025, "disturbance")
    check(np.max(abs(np.array(result["planted"]["frequencies"])-[7.1,14.23,19.55])) < .005 and
          max(abs(np.array(result["null"]["amplitudes"]))) < .005, "spectrometer control")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="acquire fresh counts and write receipt")
    parser.add_argument("--reacquire", action="store_true", help="also repeat deterministic sampling")
    args = parser.parse_args()
    path, summary = HERE/"counts.npz", HERE/"results.json"
    if args.write:
        fresh = acquire()
        result = decode(fresh)
        validate(result)
        np.savez_compressed(path, **fresh)
        summary.write_text(json.dumps(result, indent=2)+"\n")
    else:
        with np.load(path, allow_pickle=False) as stored:
            result = decode(stored)
            validate(result)
            expected = json.loads(summary.read_text())
            check(result.keys() == expected.keys(), "summary panels differ")
            for k, row in result.items():
                check(row.keys() == expected[k].keys(), "summary fields differ")
                for field, val in row.items():
                    check(np.allclose(val, expected[k][field], rtol=3e-7, atol=3e-7), "summary drift")
            if args.reacquire:
                fresh = acquire()
                check(all(np.array_equal(stored[k],v) for k,v in fresh.items()), "reacquisition differs")
    print("Free-end clock records verified.")


if __name__ == "__main__":
    main()
