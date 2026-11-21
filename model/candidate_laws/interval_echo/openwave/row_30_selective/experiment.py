"""Replay physical mode preparation, wavelength analyzers and free propagation."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

from model.candidate_laws.interval_echo.openwave.row_30_arrows.spectrum import check_head_detector
from .analyzer import eigengroups, spatial_basis
from .source import HERE, LAYERS, STEP, Source


N = 262_144
TIMES = np.arange(513) / 32


def portable(value):
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [portable(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return round(float(value), 9)
    if isinstance(value, np.integer):
        return int(value)
    return value


def settings(source, path):
    data = json.loads(path.read_text())
    if (data["ports"], data["step"], data["layers"]) != (source.ports, STEP, LAYERS):
        raise ValueError("Control fixture differs from the registered source")
    return data["pulses"]


def pulse(source, specification):
    angles = np.asarray(specification["angles"], dtype=float)
    if angles.shape != (LAYERS, len(source.ports)) or (angles < 0).any() or (angles >= 2*np.pi).any():
        raise ValueError("Invalid native reference exposure")
    return source.trajectory(angles)[0][-1]


def reverse_read(source, specification, signal):
    phase_rows = np.exp(-1j * (np.asarray(specification["angles"]) @ source.masks))
    for phases in reversed(phase_rows):
        signal = source.step @ (phases * signal)
    return float(abs(signal[0])**2)


def fly(source, state, times):
    coeff = source.modes.T @ state
    return (np.exp(-1j * np.outer(times, source.energies)) * coeff) @ source.modes.T


def fit_frequency(y):
    n = 257
    t = TIMES[:n]
    def fit(w):
        matrix = np.column_stack((np.ones(n), np.cos(w*t), np.sin(w*t)))
        coeff = np.linalg.lstsq(matrix, y[:n], rcond=None)[0]
        return float(np.sum((y[:n] - matrix @ coeff)**2)), coeff
    grid = np.arange(.5, 50.001, .1)
    index = int(np.argmin([fit(w)[0] for w in grid]))
    bounds = (max(.5, grid[index]-.1), min(50., grid[index]+.1))
    w = minimize_scalar(lambda x: fit(x)[0], bounds=bounds, method="bounded").x
    _, c = fit(w)
    yh = c[0] + c[1]*np.cos(w*TIMES[n:]) + c[2]*np.sin(w*TIMES[n:])
    r2 = 1 - np.sum((y[n:]-yh)**2) / np.sum((y[n:]-y[n:].mean())**2)
    amplitude = np.hypot(c[1], c[2])
    # 1/(2 sqrt(N)) bounds every binomial per-exposure standard error.
    return dict(omega=w, held_r2=r2, amplitude=amplitude,
                passes=bool(r2 >= .9 and amplitude > 5/(2*np.sqrt(N))))


def local_records(source, states, seed):
    active = np.unique(source.sites)
    p = np.column_stack([(abs(states[:, source.sites == s])**2).sum(axis=1) for s in active])
    if np.max(abs(p.sum(axis=1) - 1)) > 1e-9:
        raise ArithmeticError("The finite patch records lose probability")
    rng = np.random.default_rng(seed)
    return active, p, np.array([rng.multinomial(N, row/row.sum()) for row in p])


def speed_calibration(source):
    t = np.arange(513)/256
    active, p, counts = local_records(source, fly(source, source.initial, t), 891230)
    path = [0, 1, 82, 83]
    if any(b not in source.car.neighbours[a] for a, b in zip(path, path[1:])):
        raise ValueError("Speed stations are not on consecutive native seams")
    stations = []
    for port in path[1:]:
        index = int(np.where(active == port)[0][0])
        raw = counts[:, index] / N
        smooth = np.convolve(raw, np.ones(9)/9, mode="valid")
        maxima, _ = find_peaks(p[:, index], height=.001)
        estimates, _ = find_peaks(smooth, height=.001)
        estimates += 4
        stations.append(dict(port=port,
                             predicted_peak=float(t[maxima[0]]) if len(maxima) else None,
                             count_peak=float(t[estimates[0]]) if len(estimates) else None))
    speed = {}
    for key in ("predicted_peak", "count_peak"):
        ages = [r[key] for r in stations]
        velocities = []
        if None not in ages and min(np.diff(ages)) > 0:
            velocities = list(1/np.diff(ages))
        speed[key] = dict(segment_speeds=velocities,
                          end_to_end=2/(ages[-1]-ages[0]) if velocities else None,
                          common=bool(velocities and max(velocities)/min(velocities) <= 1.1))
    return dict(stations=stations, speed=speed, ports=active, counts=counts)


def geometry_audit(source, levels, basis):
    def audit(energies, modes, groups, lift):
        rows = []
        for indices in eigengroups(source.energies):
            eigenspace = lift @ source.modes[:, indices]
            weights = [np.linalg.svd(modes[:, g].T @ eigenspace, compute_uv=False)[0]**2 for g in groups]
            best = int(np.argmax(weights))
            rows.append(dict(energy=source.energies[indices[0]], multiplicity=len(indices),
                             best_k2=energies[groups[best][0]], maximum_weight=weights[best],
                             compatible_k2=[energies[g[0]] for g, w in zip(groups, weights, strict=True) if w >= .9],
                             possible_ports=np.unique(source.sites[np.linalg.norm(source.modes[:, indices], axis=1) > 1e-7])))
        return rows
    unique = np.unique(source.sites, return_counts=True)
    lift = np.zeros((source.car.n, len(source.sites)))
    for port, multiplicity in zip(*unique, strict=True):
        lift[port, source.sites == port] = 1/np.sqrt(multiplicity)
    full_levels, full_basis = np.linalg.eigh(source.car.laplacian().toarray())
    return dict(path=audit(levels, basis, eigengroups(levels), np.eye(len(source.sites))),
                screen_equal_history=audit(full_levels, full_basis, eigengroups(full_levels), lift))


def run():
    source = Source()
    check_head_detector()
    controls = settings(source, HERE/"controls.json")
    analyzers = settings(source, HERE/"analyzers.json")
    levels, basis = spatial_basis(source)
    groups = eigengroups(levels)
    states = {key: pulse(source, spec) for key, spec in controls.items()}
    spatial = np.column_stack([pulse(source, row) for row in analyzers])
    fidelity = (abs((basis.conj() * spatial).sum(axis=0)))**2
    rng = np.random.default_rng(891231)
    modes = []
    beats = []
    for i in source.indices:
        name = f"mode_{i}"
        state = states[name]
        target = source.targets[name]
        fid = abs(target @ state)**2
        flight = fly(source, state, TIMES)
        port_range = np.ptp(abs(flight[:, 0])**2)
        # A real-source analyzer can be implemented by transposing the pulse word.
        probs = abs(spatial.T @ state)**2
        check = reverse_read(source, analyzers[0], state)
        if abs(check - probs[0]) > 1e-10:
            raise ArithmeticError("Local root test differs from the reverse-order analyzer")
        clicks = rng.binomial(N, np.clip(probs, 0, 1))
        measured_groups = np.array([clicks[g].sum()/N for g in groups])
        marg = np.sqrt(np.log(2*len(levels)*len(source.indices)/1e-6)/(2*N))
        errors = np.sqrt(np.maximum(0, 1 - fidelity))
        upper = [min(1., val + len(g)*marg + errors[g].sum())
                 for val, g in zip(measured_groups, groups, strict=True)]
        returns = abs(flight @ state)**2
        modes.append(dict(index=i, energy=source.energies[i], source_fidelity=fid,
                          free_root_probability_range=port_range,
                          time_reverse_word_minimum_return=float(returns.min()),
                          analyzer_clicks=clicks, spatial_group_weights=measured_groups,
                          largest_possible_group_weight=max(upper)))
        if i != source.reference:
            beat = states[f"beat_{i}"]
            beat_fid = abs(source.targets[f"beat_{i}"] @ beat)**2
            ports, _, counts = local_records(source, fly(source, beat, TIMES), [891232,i])
            fit = fit_frequency(counts[:, 0]/N)
            outside = counts[:, ports >= 12].sum(axis=1)/N
            a = fly(source, state, TIMES)
            b = fly(source, states[f"mode_{source.reference}"], TIMES)
            choice = rng.binomial(N, .5, len(TIMES))
            mixed_home = (rng.binomial(choice, abs(a[:, 0])**2) +
                          rng.binomial(N-choice, abs(b[:, 0])**2))
            beats.append(dict(index=i, reference=source.reference, source_fidelity=beat_fid,
                              ports=ports, counts=counts, stored_mixture_counts=mixed_home,
                              count_fit=fit,
                              mixture_passes=fit_frequency(mixed_home/N)["passes"],
                              outside_min=outside.min(), outside_max=outside.max(),
                              generator_gap=abs(source.energies[i]-source.energies[source.reference])))
    zero = spatial[:, 0]
    rest_probs = [reverse_read(source, controls[f"mode_{i}"], zero) for i in source.indices]
    rest_clicks = rng.binomial(N, rest_probs)
    rest = dict(zero_mode_fidelity=fidelity[0], energy_indices=source.indices, energy_clicks=rest_clicks,
                complement_estimate=1-rest_clicks.sum()/N,
                minimum_free_return=float(np.min(abs(fly(source, zero, TIMES) @ zero)**2)))
    if min(m["source_fidelity"] for m in modes) < .98 or min(fidelity) < .98:
        raise ArithmeticError("Declared controls failed mode preparation")
    return portable(dict(candidate="interval_echo", scale=dict(length="one eligible native seam",
                         time="native exposure", shots=N), source_ports=source.ports, indices=source.indices,
                         reference=source.reference, spatial_k2=[levels[g[0]] for g in groups],
                         spatial_analyzer_minimum_fidelity=min(fidelity),
                         rest=rest, modes=modes, beats=beats,
                         sharp_packet=speed_calibration(source),
                         spectral_audit=geometry_audit(source, levels, basis)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a fresh evidence file")
    args = parser.parse_args()
    receipt = run()
    if args.output:
        with args.output.open("x") as output:
            json.dump(receipt, output, separators=(",", ":"), allow_nan=False)
            output.write("\n")
    elif receipt != json.loads((HERE/"results.json").read_text()):
        raise SystemExit("Selective free-wave receipt differs from the archive")
    print("Modes", [(m['index'], m['source_fidelity'], m['largest_possible_group_weight'])
                    for m in receipt["modes"]])
    print("Beats", [(b['index'], b['count_fit']['omega'], b['count_fit']['held_r2'],
                     b['outside_max'] - b['outside_min']) for b in receipt["beats"]])
    print("Rest", receipt["rest"])
    print("Speed", receipt["sharp_packet"]["speed"])
