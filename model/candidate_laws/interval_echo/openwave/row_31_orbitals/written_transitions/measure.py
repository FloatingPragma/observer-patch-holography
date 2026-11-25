"""Resolve native written-core transitions against an independent binding reading."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

from model.candidate_laws.interval_echo import candidate
from model.executable_model import EMITTER_SEAM, federation

HERE = Path(__file__).resolve().parent
SHOTS, SEED, BINS, STEP = 2**22, 311953, 4096, 1 / 4
FREQUENCIES = np.linspace(0.44, 0.52, 41)


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def compare(a, b, path="written_transitions"):
    if isinstance(b, dict):
        check(isinstance(a, dict) and a.keys() == b.keys(), path)
        for k in b:
            compare(a[k], b[k], path + "." + k)
    elif isinstance(b, list):
        check(isinstance(a, list) and len(a) == len(b), path)
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            compare(x, y, f"{path}[{i}]")
    elif isinstance(b, float):
        check(np.isfinite(a) and np.isclose(a, b, atol=2e-9, rtol=1e-10), path)
    else:
        check(type(a) is type(b) and a == b, path)


def generator(board, seam, strength):
    h = board.laplacian().toarray() / 2
    a, b = seam
    h[a, a] -= strength + 0.5
    h[b, b] -= strength + 0.5
    h[a, b] += 0.5
    h[b, a] += 0.5
    return h


def kernel(delta):
    phase = STEP * np.asarray(delta)
    return (np.exp(0.5j * (BINS - 1) * phase)
            * np.sinc(BINS * phase / (2 * np.pi)) / np.sinc(phase / (2 * np.pi)))


def record(wave, rng):
    raw = np.r_[abs(wave)**2, 1 - float(np.vdot(wave, wave).real)]
    check(np.all(np.isfinite(raw)) and min(raw) >= 0, "negative or invalid complete record")
    check(abs(sum(raw) - 1) < 2e-11, "record normalization drift")
    order = np.argsort(raw, kind="stable")
    exact = raw[order].copy()
    exact[-1] = 1 - sum(exact[:-1])
    check(exact[-1] >= 0, "negative explicit complement")
    counts = np.empty(len(raw), np.int64)
    counts[order] = rng.multinomial(SHOTS, exact)
    return dict(probabilities_by_port_then_clock_refusal=raw.tolist(),
                counts_by_port_then_clock_refusal=counts.tolist(),
                acceptance=float(sum(raw[:-1])),
                sampling_l1_complement_bound=float(abs(exact[-1] - raw[order[-1]])))


def filtered(e, v, ground, ports, target):
    # This general form also admits the unchanged input under the altered law.
    initial = v.T @ ground
    coeff = kernel(-target) * initial
    for port in ports:
        coeff -= 2 * v[port] * (kernel(e[:, None] - e[None, :] - target)
                                @ (v[port] * initial))
    return v @ coeff


def stationary_filtered(e, v, ground, ports, target):
    pulsed = ground.copy()
    pulsed[list(ports)] *= -1
    return v @ (kernel(e - e[0] - target) * (v.T @ pulsed))


def binding_meter(binding, rng):
    error = float(np.arcsin(2 * np.sqrt(2) * np.sqrt(np.log(32 / 0.01) / (2 * SHOTS))))
    reads, estimate = [], None
    for wait in (0.125, 8):
        z = np.exp(1j * binding * wait)
        probabilities = (1 + np.array([z.real, z.imag])) / 2
        counts = rng.binomial(SHOTS, probabilities)
        angle = float(np.angle(complex(*(2 * counts / SHOTS - 1))))
        if estimate is not None:
            angle += 2 * np.pi * round((wait * estimate - angle) / (2 * np.pi))
        estimate = angle / wait
        reads.append(dict(wait=wait, real_imag_plus=counts.tolist()))
    check(error / 0.125 * 8 + error < np.pi, "binding winding unresolved")
    return dict(readings=reads, measured_binding=float(estimate),
                simultaneous_99_binding_halfwidth=error / 8,
                hydrogen_first_transition=0.75 * float(estimate),
                hydrogen_target_halfwidth=0.75 * error / 8)


def fit_line(records):
    counts = np.array([r["counts_by_port_then_clock_refusal"][-1] for r in records])
    accepted = 1 - counts / SHOTS
    train = np.arange(len(FREQUENCIES)) % 2 == 0

    def model(parameters):
        frequency, amplitude, background = parameters
        return background + amplitude * abs(kernel(frequency - FREQUENCIES))**2

    result = least_squares(lambda p: model(p)[train] - accepted[train],
                           [0.48, 0.08, 1e-6], bounds=([0.47, 0, 0], [0.49, 1, 0.01]),
                           xtol=1e-13, ftol=1e-13, gtol=1e-13)
    residual = float(np.sqrt(np.mean((model(result.x)[~train] - accepted[~train])**2)))
    check(residual < 0.002, "held line profile fails")
    return dict(frequency=float(result.x[0]), amplitude=float(result.x[1]),
                background=float(result.x[2]), held_probability_rms=residual,
                training_indices=np.flatnonzero(train).tolist(),
                held_indices=np.flatnonzero(~train).tolist())


def census(level):
    board = federation(level)
    a, b = EMITTER_SEAM[level]
    # H=L_cut/2 - kappa P. QHQ is nonnegative; PHP is diagonal.
    cut = board.laplacian().astype(float).tolil() / 2
    cut[a, a] -= 0.5
    cut[b, b] -= 0.5
    cut[a, b] += 0.5
    cut[b, a] += 0.5
    cut = cut.tocsr()
    check(max(abs(np.asarray(cut.sum(axis=1)).ravel())) < 1e-13
          and min(cut.diagonal()) >= 0, "cut-graph Laplacian identity fails")
    block = cut[[a, b]][:, [a, b]].toarray() - candidate.KAPPA * np.eye(2)
    check(max(abs(block - np.diag(np.diag(block))).ravel()) < 1e-13, "endpoint compression not diagonal")
    # Cut Laplacian positivity proves H>=-8I, stronger than the rank-two V>=-9I.
    return dict(cells=board.cells, ports=board.n, seam=[int(a), int(b)],
                endpoint_degrees=board.degree[[a, b]].tolist(),
                negative_level_count=2, negative_band=[-8.0, float(max(np.diag(block)))],
                endpoint_compression=block.tolist(), excluded_open_interval=[float(max(np.diag(block))), 0])


def panel(level, rng):
    board = federation(level)
    seam = tuple(EMITTER_SEAM[level])
    h = generator(board, seam, candidate.KAPPA)
    e, v = eigh(h)
    check(sum(e < -1e-9) == 2 and e[2] >= 0, "two-sided level census fails")
    ground = v[:, 0].copy()
    if ground[seam[0]] < 0:
        ground *= -1
    sharp = np.eye(board.n)[:, seam[0]]
    direction = (sharp - ground) / np.linalg.norm(sharp - ground)
    compiled = sharp - 2 * direction * (direction @ sharp)
    check(np.linalg.norm(compiled - ground) < 1e-12, "supplied compiler differs from source")
    meter = binding_meter(-e[0], rng)
    scan = [record(stationary_filtered(e, v, ground, [seam[0]], w), rng) for w in FREQUENCIES]
    fit = fit_line(scan)
    target = meter["hydrogen_first_transition"]
    hydrogen = record(stationary_filtered(e, v, ground, [seam[0]], target), rng)
    common = record(stationary_filtered(e, v, ground, seam, fit["frequency"]), rng)
    native_line = record(stationary_filtered(e, v, ground, [seam[0]], fit["frequency"]), rng)
    planted = e.copy()
    planted[1] = e[0] / 4
    planted_control = record(stationary_filtered(planted, v, ground, [seam[0]], target), rng)
    changed_e, changed_v = eigh(generator(board, seam, 0))
    changed = record(filtered(changed_e, changed_v, ground, [seam[0]], fit["frequency"]), rng)
    target_energy = e[0] + target
    distance = min(-target_energy, target_energy + 5.5) - meter["hydrogen_target_halfwidth"]
    # H=L_cut/2-8P lies in [-8,6]; max cut degree<=6 bounds the upper edge.
    largest_detuning = max(abs(-8 - target_energy), abs(6 - target_energy)) + meter["hydrogen_target_halfwidth"]
    check(distance > 0 and largest_detuning < np.pi / STEP,
          "Fourier sidelobe bound crosses an alias or leaves its monotone sine interval")
    leakage = float((BINS * np.sin(STEP * distance / 2))**-2)
    check(hydrogen["acceptance"] <= leakage, "native hydrogen target exceeds spectral sidelobe bound")
    check(native_line["acceptance"] > 0.07 and common["acceptance"] < 0.001,
          "endpoint/common-phase selection unresolved")
    check(planted_control["acceptance"] > 0.07, "planted first-excited detector sensitivity fails")
    return dict(cells=board.cells, seam=list(map(int, seam)),
                source_compiler=dict(kind="supplied spatial Householder", ideal_acceptance=1,
                                     ideal_refusal=0, local_native_implementation="unestablished"),
                native_bound_energies=e[:2].tolist(), native_next_level=float(e[2]),
                binding_meter=meter, scan=scan, fitted_line=fit, selected_native_line=native_line,
                common_phase_control=common, hydrogen_target=hydrogen,
                hydrogen_target_acceptance_upper_bound=leakage,
                full_native_energy_enclosure=[-8, 6],
                target_minimum_detuning=distance,
                maximum_target_detuning=largest_detuning,
                alias_half_period=float(np.pi / STEP),
                planted_spectrum_detector_control=planted_control,
                zero_contrast_same_input=changed,
                exact_endpoint_transition_weight=float(4 * (v[seam[0], 0] * v[seam[0], 1])**2),
                exact_common_phase_transition_weight=float(4 * sum(v[p, 0] * v[p, 1] for p in seam)**2))


def run():
    rng = np.random.default_rng(SEED)
    return dict(shots=SHOTS, seed=SEED, clock_bins=BINS, clock_step=STEP,
                scan_frequencies=FREQUENCIES.tolist(), censuses=[census(level) for level in range(3)],
                panels=[panel(level, rng) for level in (0, 1)],
                scope="Written-core material transition spectroscopy; no hydrogenic atom or physical nucleus",
                official_status="OW-31 untested at 0df73dc6; full local criterion partial")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--receipt", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    result = run()
    if args.output:
        encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
        with args.output.open("x") as stream:
            stream.write(encoded)
    else:
        compare(result, json.loads(args.receipt.read_text()))
        print("Replayed native written-core line, hydrogenic target refusal and every control.")
