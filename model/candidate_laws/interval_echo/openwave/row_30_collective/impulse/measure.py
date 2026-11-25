"""Reconstruct a real-space wall response from physical, parity-even records."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .. import experiment, formation

HERE = Path(__file__).resolve().parent
TIMES = 0.06 * np.arange(61)
SHOTS = 2**36
HAND = 2**32
SOURCE = (8, 9, 10)
SEED = 312141


def material(n, j, h):
    # i J A_s B_s and i h B_s A_(s+1) are the physical Hamiltonian
    # terms -J tau_z and -h tau_x tau_x. The first-order commutator has no
    # mass, wave number, inferred frequency or square-root dispersion input.
    d = 2 * (n + 1)
    g = np.zeros((d, d))
    for a in range(0, d, 2):
        g[a, a + 1] = 2 * j
        if a + 2 < d:
            g[a + 1, a + 2] = 2 * h
    return g - g.T


def acquire(n):
    pairs = np.array(list(combinations(range(2 * (n + 1)), 2)))
    a, b = pairs.T
    initials = list(combinations(SOURCE, 2))
    rng = np.random.default_rng(183020)
    banks, hands = [], []
    for j in (experiment.K, experiment.K / 2):
        rows, clock = [], []
        r = np.eye(2 * (n + 1))
        step = expm(TIMES[1] * material(n, j, experiment.H))
        for t in TIMES:
            coh = np.stack([r[a, i] * r[b, k] - r[b, i] * r[a, k] for i, k in initials])
            # Two signed product readings of balanced initial dimer orientations.
            # Their lower-order marginals are exactly fair, on every setting.
            plus = rng.binomial(SHOTS, 0.5, size=coh.shape)
            same = np.clip((1 + coh) / 2, 0, 1)
            pp = rng.binomial(plus, same)
            mm = rng.binomial(SHOTS - plus, same)
            rows.append(np.stack([pp, plus - pp, SHOTS - plus - mm, mm], axis=-1))
            clock.append(rng.binomial(HAND, (1 + np.array([np.cos(t), np.sin(t)])) / 2))
            r = step @ r
        banks.append(rows)
        hands.append(clock)
    return np.array(banks, dtype=np.uint64), np.array(hands, dtype=np.uint64)


def read(banks, hand):
    dimension = int((1 + np.sqrt(1 + 8 * banks.shape[-2])) / 2)
    a, b = np.triu_indices(dimension, 1)
    last = np.eye(dimension)[:, SOURCE[0]]
    cols, overlaps, gaps = [], [], []
    for bank in banks:
        observed = (bank.astype(float) @ np.array([1, -1, -1, 1])) / SHOTS
        wedge = np.zeros((3, dimension, dimension))
        wedge[:, a, b] = observed
        wedge[:, b, a] = -observed
        lifted = wedge @ wedge.transpose(0, 2, 1)
        eig, vec = np.linalg.eigh((lifted[0] + lifted[1] - lifted[2]) / 2)
        curr = vec[:, -1]
        overlaps.append(abs(float(last @ curr)))
        gaps.append(float(eig[-1] - eig[-2]))
        if last @ curr < 0:
            curr = -curr
        cols.append(curr)
        last = curr
    cs = 2 * hand.astype(float) / HAND - 1
    times = np.unwrap(np.arctan2(cs[:, 1], cs[:, 0]))
    return np.asarray(cols), times, min(overlaps), min(gaps)


def forecasts(n, j, h, times):
    # Independent second-order initial-value problem: one A component starts
    # at 1, all velocities are zero. End 0 has the wall Robin correction.
    m = n + 1
    mass = 2 * abs(j - h)
    b = 4 * j * h
    wave = np.diag(np.full(m, mass**2 + 2 * b))
    wave += np.diag(np.full(m - 1, -b), 1) + np.diag(np.full(m - 1, -b), -1)
    wave[0, 0] -= 4 * h * h
    squared, modes = np.linalg.eigh(wave)

    def cosine(e):
        return (np.cos(np.outer(times, e)) * modes[SOURCE[0] // 2]) @ modes.T

    kg = cosine(np.sqrt(squared))
    hop = cosine(mass + (squared - mass**2) / (2 * mass))
    stationary = np.zeros_like(kg)
    stationary[:, SOURCE[0] // 2] = np.cos(mass * times)
    return kg, hop, stationary


def ordinary_long_waves(n, j, h, times, field):
    # A geometric spatial basis with ghost values u[-1]=(h/J)u[0] and u[m]=0.
    # Its eigenvalues are 4 sin^2(k/2), independently of the flight readings.
    m = n + 1
    lap = np.diag(np.full(m, 2.0))
    lap += np.diag(-np.ones(m - 1), 1) + np.diag(-np.ones(m - 1), -1)
    lap[0, 0] -= h / j
    lam, basis = np.linalg.eigh(lap)
    k = 2 * np.arcsin(np.sqrt(lam) / 2)
    low = k <= 0.5
    k = k[low]
    basis = basis[:, low]
    rest = 2 * abs(j - h)
    speed_sq = 4 * j * h
    frequencies = np.sqrt(rest**2 + speed_sq * k**2)
    prediction = np.cos(np.outer(times, frequencies)) * basis[SOURCE[0] // 2]
    observed = field @ basis
    error = np.linalg.norm(observed - prediction, axis=1)
    bound = np.max(abs(times)) * speed_sq * max(k) ** 4 / (24 * rest)
    return dict(
        k=k.tolist(),
        max_error=float(max(error)),
        lattice_to_continuum_bound=float(bound),
    )


def report(n, banks, hands, formation_record):
    bench = json.loads((HERE.parent / "coupling_meter.json").read_text())["results"]
    d = 2 * (n + 1)
    # Allow simultaneous 99% coverage for signs and hand means at both contrasts.
    total = len(TIMES) * 2 * (3 * d * (d - 1) // 2 + 2)

    def noise(nshots):
        return np.sqrt(2 * np.log(2 * total / 0.01) / nshots)

    matrix_err = (d - 1) * noise(SHOTS)
    projector_err = 1.5 * (2 * matrix_err + matrix_err**2)
    column_err = np.sqrt(2) * projector_err / (1 - projector_err)
    time_err = np.arcsin(np.sqrt(2) * noise(HAND))
    material_error = 2 * (0.0003 + 0.004)
    rows = []
    for tag, bank, hand in zip(("native", "half_contrast"), banks, hands, strict=True):
        cols, times, overlap, gap = read(bank, hand)
        k = bench["J"]["value"]
        h = bench["h"]["value"]
        norm_budget = 2 * (k + h + 0.0043)
        if np.cos(norm_budget * TIMES[1]) <= 2 * column_err:
            raise ValueError("Continuation could alias a sign")
        kg, hop, stationary = forecasts(n, k, h, times)
        numerical = column_err + norm_budget * time_err
        true_t = np.max(abs(times)) + time_err
        # Duhamel for first-order contact generator. For the hopping foil use
        # independent box bounds on the squared operator and inverse mass.
        cal = material_error * true_t
        mass = 2 * (k - h)
        lower = mass - material_error
        w2err = (2 * norm_budget + material_error) * material_error
        # E_hop=(W^2/m+m)/2. Bound changes in both W and m.
        hop_cal = (
            true_t
            * 0.5
            * (
                w2err / lower
                + norm_budget**2 * material_error / (mass * lower)
                + material_error
            )
        )
        hop_meter_error = column_err + time_err * 0.5 * (
            norm_budget**2 / lower + mass + material_error
        )

        def errors(pred):
            return np.linalg.norm(cols[:, ::2] - pred, axis=1)

        own = forecasts(
            n,
            experiment.K if tag == "native" else experiment.K / 2,
            experiment.H,
            times,
        )[0]
        rows.append(
            dict(
                material=tag,
                min_overlap=overlap,
                min_eigen_gap=gap,
                hand_max_abs_error=float(max(abs(times - TIMES))),
                column_error_bound=float(column_err),
                calibration_box_kg=float(cal),
                kg_count_and_clock=float(numerical),
                hopping_total_allowance=float(hop_cal + hop_meter_error),
                native_kg_max=float(max(errors(kg))),
                hopping_max=float(max(errors(hop))),
                hopping_peak_appointment=int(np.argmax(errors(hop))),
                immobile_max=float(max(errors(stationary))),
                own_law_max=float(max(errors(own))),
                remote_A9_max=float(max(abs(cols[:, 18]))),
                ordinary_long_wave=ordinary_long_waves(n, k, h, times, cols[:, ::2]),
                hand_times=times.tolist(),
                reconstructed_A=cols[:, ::2].tolist(),
                kg_L2_errors=errors(kg).tolist(),
                hopping_L2_errors=errors(hop).tolist(),
            )
        )
    return dict(
        strand=formation_record,
        settings_per_appointment=3 * d * (d - 1) // 2,
        shots_per_pair=SHOTS,
        hand_shots_per_quadrature=HAND,
        time_error_bound=float(time_err),
        confidence=0.99,
        rows=rows,
    )


def main(write):
    experiment.native_checks()
    source = formation.prepare(0, SEED)
    n = source["dimers_enrolled"]
    banks, hand = acquire(n)
    result = report(n, banks, hand, source)
    path, receipt = HERE / "counts.npz", HERE / "results.json"
    if write:
        np.savez_compressed(path, banks=banks, hand=hand)
        receipt.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        with np.load(path, allow_pickle=False) as old:
            for name, arr in [("banks", banks), ("hand", hand)]:
                np.testing.assert_array_equal(old[name], arr)
        old = json.loads(receipt.read_text())
        if old.keys() != result.keys():
            raise ValueError("Summary keys changed")
        for key in result.keys() - {"rows"}:
            if old[key] != result[key]:
                raise ValueError(f"Acquisition metadata changed: {key}")
        for x, y in zip(old["rows"], result["rows"], strict=True):
            if x.keys() != y.keys():
                raise ValueError("Incomplete reconstructed material record")
            for key in x:
                if key == "material":
                    if x[key] != y[key]:
                        raise ValueError("Wrong flight law")
                elif isinstance(x[key], dict):
                    for sub in x[key]:
                        np.testing.assert_allclose(
                            x[key][sub], y[key][sub], atol=1e-10, rtol=1e-10
                        )
                else:
                    np.testing.assert_allclose(x[key], y[key], atol=1e-10, rtol=1e-10)
    for r in result["rows"]:
        print({k: v for k, v in r.items() if not isinstance(v, list)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    main(parser.parse_args().write)
