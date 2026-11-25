"""Coin flips, independent energy bins, and local component reads on a neutral collar."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation

HERE = Path(__file__).parent
SHOTS, TIMER_SHOTS = 2**20, 2**26
TIMES = np.linspace(0.02, 12.02, 193)
ARMS, STEP = 512, 1 / 16


def enroll(rng, p, shots=SHOTS):
    p = np.round(np.asarray(p, dtype=float), 13)
    if min(p) < 0 or abs(sum(p) - 1) > 1e-10:
        raise ValueError("Invalid complete probability law")
    # Place the roundoff on the largest bin, used as the multinomial complement.
    order = np.argsort(p)
    counts = np.empty(p.size, dtype=np.int64)
    counts[order] = rng.multinomial(shots, p[order])
    return counts.tolist()


def fixture():
    echo = law.echo()
    i = law.SRC.index((1, 1))
    ports = [i] + [
        (law.LABELS.index((c, d)) + 1) * law.N_SRC + i for c in (0, 1) for d in (2, -2)
    ]
    basis = np.eye(law.DIM)[:, ports]
    coin = [np.diag([0, 1, 1, 0, 0]), np.diag([0, 0, 0, 1, 1])]
    launch = np.eye(law.DIM)[:, law.SRC.index((2, 0))]
    exchange, contrast = law.local_term_parts()
    split = expm(-0.5j * np.pi * (exchange + candidate.KAPPA * contrast)) @ launch
    ref = np.tile([b * b for _, b in law.SRC], 1 + law.N_LAB)
    phase = {
        "even": np.pi / (8 * candidate.KAPPA),
        "odd": 3 * np.pi / (8 * candidate.KAPPA),
    }
    starts = {}
    for name, v in {
        **{
            s: np.exp(-1j * t * candidate.KAPPA * ref) * split for s, t in phase.items()
        },
        "sharp": launch,
    }.items():
        written = echo @ v
        source = basis.T @ written
        np.testing.assert_allclose(basis @ source, written, atol=1e-13)
        arms = [gate @ source for gate in coin]
        for new in arms:
            np.testing.assert_allclose(np.vdot(new, new), 0.5, atol=1e-13)
        starts[name] = arms
    for name, sign in (("even", 1), ("odd", -1)):
        np.testing.assert_allclose(
            starts[name][0][2], sign * starts[name][0][1], atol=1e-13
        )
    boards = {}
    for level in range(3):
        board = federation(level)
        local = [
            (a, b)
            for edge in board.seams
            for a, b in (edge, edge[::-1])
            if a // 12 == b // 12 == 0 and board.degree[b] < 6
        ]
        a, b = min(local)
        if not np.all(board.degree <= 6):
            raise ValueError("unexpected board")
        boards[str(board.cells)] = [a, b]
    return basis, coin, starts, boards


def energy_weights(h, arms):
    energies, vectors = np.linalg.eigh(h)
    weights = sum(abs(vectors.T.conj() @ w) ** 2 for w in arms)
    return energies, vectors, weights


def qpe(energies, weights):
    phases = np.exp(1j * STEP * np.outer(energies, np.arange(ARMS)))
    return weights @ abs(np.fft.fft(phases, axis=1) / ARMS) ** 2


def energy_read(rng, h, arms, calibration, kappa, basis):
    energies, _, weights = energy_weights(h, arms)
    bins = enroll(rng, qpe(energies, weights))
    gap = 2 * np.arcsin(np.sqrt(calibration / TIMER_SHOTS)) / 16
    energy_axis = 2 * np.pi * np.arange(ARMS) / (ARMS * gap)
    window = np.flatnonzero((energy_axis > 6) & (energy_axis < 46))
    peak = int(window[np.argmax(np.asarray(bins)[window])])
    # Switch on G against idle at pi. At these even integer contrasts the
    # contrast-square phase is unity, so the control measures reversal parity.
    ex, co = law.generator_parts()
    parity = (np.eye(law.DIM) - np.real(expm(-1j * np.pi * (ex + kappa * co)))) / 2
    np.testing.assert_allclose(parity, ex, atol=2e-13)
    transformed = [law.echo() @ basis @ w for w in arms]
    occupied = np.diag(co) == 4
    occup = sum(np.vdot(w[occupied], w[occupied]).real for w in transformed)
    odd = sum(np.vdot(w, parity @ w).real for w in transformed)
    local = enroll(rng, [1 - occup, occup - odd, odd])
    mean = (4 * kappa * sum(local[1:]) + local[2]) / SHOTS
    return {
        "fourier_bins": bins,
        "local_blank_even_odd": local,
        "energy_increment_calibration": calibration,
        "nonzero_window_count": int(sum(np.asarray(bins)[window])),
        "peak_bin": peak,
        "peak": float(energy_axis[peak]),
        "bin_width": float(energy_axis[1]),
        "mean_from_components": mean,
    }


def fit(clock_rows):
    rows = np.array(clock_rows)
    t = 16 * np.arcsin(np.sqrt(rows[:, 0] / TIMER_SHOTS))
    p = rows[:, 2] / SHOTS

    def prediction(w):
        a = np.column_stack((np.ones(len(t)), np.cos(w * t), np.sin(w * t)))
        coeff = np.linalg.lstsq(a[:97], p[:97], rcond=None)[0]
        return a @ coeff, 2 * np.linalg.norm(coeff[1:])

    def loss(w):
        return np.mean((prediction(w)[0][:97] - p[:97]) ** 2)

    grid = np.arange(0.4, 46.01, 0.08)
    best = grid[np.argmin([loss(w) for w in grid])]
    w = float(
        minimize_scalar(
            loss, bounds=(max(0.4, best - 0.08), min(46, best + 0.08)), method="bounded"
        ).x
    )
    pred, amplitude = prediction(w)
    return {
        "frequency": w if amplitude > 0.03 else None,
        "peak_to_trough": float(amplitude),
        "held_rms": float(np.sqrt(np.mean((pred[97:] - p[97:]) ** 2))),
    }


def run():
    basis, coin, starts, boards = fixture()
    rng = np.random.default_rng(509336)
    record = {
        "boards": boards,
        "starts": SHOTS,
        "timer_starts": TIMER_SHOTS,
        "panels": {},
    }
    ex, co = law.local_term_parts()
    clock_cols = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    timer = (ex + candidate.KAPPA * co)[np.ix_(clock_cols, clock_cols)]
    for kappa in (8, 4, 0):
        full = ex + kappa * co
        h = basis.T @ full @ basis
        np.testing.assert_allclose(full @ basis, basis @ h, atol=1e-13)
        eigen, vectors = np.linalg.eigh(h)
        # All neighbouring blank seams have filled loads throughout release.
        filled = np.eye(law.DIM)[:, law.SRC.index((1, 1))]
        np.testing.assert_allclose(full @ filled, 0, atol=1e-13)
        for name, arms in starts.items():
            rows = []
            for t in TIMES:
                phase = vectors @ np.diag(np.exp(-1j * t * eigen)) @ vectors.T.conj()
                flip = sum(
                    np.linalg.norm(coin[1 - c] @ phase @ arms[c]) ** 2 for c in (0, 1)
                )
                blank = sum(abs((phase @ w)[0]) ** 2 for w in arms)
                read = enroll(rng, [1 - flip - blank, flip, blank])
                timer_p = abs(expm(-1j * timer * t / 8)[1, 0]) ** 2
                rows.append([rng.binomial(TIMER_SHOTS, timer_p), *read])
            cal = rng.binomial(
                TIMER_SHOTS, abs(expm(-1j * timer * STEP * 16)[1, 0]) ** 2
            )
            energy = energy_read(rng, h, arms, cal, kappa, basis)
            record["panels"][f"{kappa}_{name}"] = {
                "counts": rows,
                "fit": fit(rows),
                "energy": energy,
            }
    even = starts["even"]
    h = basis.T @ (ex + 8 * co) @ basis
    e, v = np.linalg.eigh(h)
    branches = [
        (v[:, np.isclose(e, val)] @ (v[:, np.isclose(e, val)].T.conj() @ even[c]), c)
        for c in (0, 1)
        for val in (0, 32)
    ]
    control = []
    for t in TIMES:
        u = v @ np.diag(np.exp(-1j * t * e)) @ v.T.conj()
        flip = sum(np.linalg.norm(coin[1 - c] @ u @ w) ** 2 for w, c in branches)
        np.testing.assert_allclose(flip, 0.5, atol=1e-13)
        control.append(enroll(rng, [1 - flip, flip]))
    record["dephased_same_flip"] = control
    sham = basis[:, 0]
    np.testing.assert_allclose((ex + 8 * co) @ sham, 0, atol=1e-13)
    record["uniform_blank_sham_counts"] = [0, 0, SHOTS]
    return record


if __name__ == "__main__":
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
        compare,
    )

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--record", type=Path)
    args = ap.parse_args()
    result = run()
    if args.record:
        with args.record.open("x") as f:
            json.dump(result, f, separators=(",", ":"))
            f.write("\n")
    else:
        compare(result, json.loads((HERE / "records.json").read_text()), "written coin")
    for tag, p in result["panels"].items():
        print(tag, p["fit"], p["energy"]["peak"], p["energy"]["mean_from_components"])
