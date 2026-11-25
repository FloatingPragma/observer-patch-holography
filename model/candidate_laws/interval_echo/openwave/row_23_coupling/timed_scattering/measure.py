"""Two-time, address-resolved scattering from a freely evolving held core."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from model.candidate_laws.interval_echo.openwave.row_21_boost import measure as pair
from model.candidate_laws.interval_echo.openwave.row_23_coupling.mobile_form_factor import (
    measure as static,
)
from model.candidate_laws.interval_echo.openwave.row_30.experiment import native_matrix
from model.executable_model import federation
from model.experiments.count_intervals import clopper_pearson

HERE = Path(__file__).resolve().parent
DELAYS = (0, 1, 4, 16)
TAIL = 0.01 / 300000  # 3 sources x 2 laws x 4 times x 6 x 214 bins x 2 tails
GUARD = 2**-39  # includes complement restoration for <=430 rounded bins
SEED = 231209


def sector(n, entries, winding):
    i, j, shift, amplitude = entries
    phase = np.exp(-2j * np.pi * (winding % n) * shift / n)
    return coo_matrix((amplitude * phase, (i, j)), shape=(2 * n - 1, 2 * n - 1)).tocsr()


def density(n, j):
    idx = np.arange(2 * n - 1)
    a = idx % 2
    q = 2 * np.pi * j / (2 * n)
    return np.exp(-1j * q * a) + np.exp(-1j * q * (a + 1 + idx // 2))


def evolved(state, step):
    states = {0: state}
    for t in range(1, max(DELAYS) + 1):
        state = step(state)
        if t in DELAYS:
            states[t] = state
    return states


def norms(early, late):
    # The preparation is a classical equal mixture of the two starting rails.
    return [
        np.vdot(early, early).real / early.shape[-1],
        np.vdot(late, late).real / early.shape[-1],
        np.vdot(early, late) / early.shape[-1],
    ]


def pair_spectrum(n, winding, state, entries, amp):
    size = 2 * n
    mats = [sector(n, entries, w) for w in range(n)]
    late = evolved(state, mats[winding].dot)
    data = {t: [] for t in DELAYS}
    for j in range(size):
        read = (amp - 1) * density(n, j) / size + (j == 0)
        early = evolved(read[:, None] * state, mats[(winding - j) % n].dot)
        for t in DELAYS:
            data[t].append(norms(early[t], read[:, None] * late[t]))
    return {t: np.array(vals).T for t, vals in data.items()}


def single_step(state, strength):
    for p, theta in enumerate(pair.THETA):
        ev = (np.arange(0, len(state), 2) + p) % len(state)
        od = (ev + 1) % len(state)
        gate = native_matrix(theta, strength)
        moved = state.copy()
        moved[ev] = gate[0, 0] * state[ev] + gate[0, 1] * state[od]
        moved[od] = gate[1, 0] * state[ev] + gate[1, 1] * state[od]
        state = moved
    return state


def single_spectrum(size, strength, amp):
    x = np.arange(size)
    read = (amp - 1) / size * np.exp(-2j * np.pi * x[:, None] * x / size)
    read[:, 0] += 1
    state = np.zeros((size, 1), complex)
    state[0] = 1

    def step(s):
        return single_step(s, strength)

    early = evolved(read * state, step)
    late = evolved(state, step)
    result = {}
    for t in DELAYS:
        early_signal, late_signal = early[t], read * late[t]
        result[t] = np.array(
            [
                np.sum(abs(early_signal) ** 2, axis=0),
                np.sum(abs(late_signal) ** 2, axis=0),
                np.sum(early_signal.conj() * late_signal, axis=0),
            ]
        )
    return result


def acquire(data, admission, rng):
    early_signal, late_signal, overlap = data
    early_signal, late_signal = early_signal.real, late_signal.real
    np.testing.assert_allclose(
        np.sum(early_signal), np.sum(late_signal), atol=2e-12, rtol=0
    )
    spill = max(0.0, admission - early_signal.sum())
    rejected = max(0.0, 1 - admission)

    def observe(values):
        return np.array(static.counts([*values, spill, rejected], rng))

    arms = [observe(early_signal), observe(late_signal)]
    for cross in (overlap.real, -overlap.imag):
        avg = (early_signal + late_signal) / 4
        arms.append(observe(np.r_[avg + cross / 2, avg - cross / 2]))
    return arms


def interval(arr):
    lo, hi = clopper_pearson(arr, static.SHOTS, TAIL)
    return np.maximum(0, lo - GUARD), np.minimum(1, hi + GUARD)


def read_banks(banks, size):
    es, ls, phase0, phase1 = banks
    el, eh = interval(es[:size])
    ll, lh = interval(ls[:size])
    low, high = np.sqrt(el * ll), np.sqrt(eh * lh)
    if np.min(low) <= 0:
        raise ValueError("Unresolved single-arm normalization")
    magnitude = np.sqrt((es[:size] / static.SHOTS) * (ls[:size] / static.SHOTS))
    est, ci = [], []
    for phase, sign in [(phase0, 1), (phase1, -1)]:
        p, m = phase[:size], phase[size : 2 * size]
        pl, ph = interval(p)
        ml, mh = interval(m)
        est.append(sign * (p - m) / static.SHOTS / magnitude)
        ends = np.array([sign * v / d for v in (pl - mh, ph - ml) for d in (low, high)])
        ci.append((ends.min(axis=0), ends.max(axis=0)))
    return np.array(est), np.array(ci)


def compare(standard, switched, radius):
    a, (al, ah) = standard[0], standard[1].transpose(1, 0, 2)
    b, (bl, bh) = switched[0], switched[1].transpose(1, 0, 2)
    lo, hi = al - bh, ah - bl
    nearest = np.maximum(0, np.maximum(lo, -hi))
    farthest = np.maximum(abs(lo), abs(hi))
    below, above = np.linalg.norm(nearest, axis=0), np.linalg.norm(farthest, axis=0)
    estimate = np.linalg.norm(a - b, axis=0)
    size = len(estimate)
    half = size // 2 + 1
    reach = np.flatnonzero(below[1:half] > 0.10) + 1
    first = int(reach[0]) if len(reach) else None
    # A bin before the bracket must be certified below or its earlier crossing
    # cannot be excluded. Start immediately after the certified contiguous prefix.
    free = np.flatnonzero(above[:half] >= 0.10)
    prefix = int(free[0]) - 1 if len(free) else half - 1
    if first and prefix < 0:
        raise ValueError("Forward normalization unresolved")
    last_below = int(np.flatnonzero(above[:first] < 0.10)[-1]) if first else None
    return dict(
        first_above=first,
        last_below_before_first_above=last_below,
        certified_below_prefix_end=prefix,
        crossing_qR=None
        if first is None or radius is None
        else [
            2 * np.pi * prefix / size * radius[0],
            2 * np.pi * first / size * radius[1],
        ],
        detected_bins=int(np.sum(below[1:half] > 0)),
        unresolved_at_threshold=int(
            np.sum((below[1:half] <= 0.10) & (above[1:half] >= 0.10))
        ),
        peak_bin=int(np.argmax(estimate[:half])),
        max_lower=float(max(below[:half])),
        max_upper=float(max(above[:half])),
        difference_all_bins=estimate.tolist(),
        lower_all_bins=below.tolist(),
        upper_all_bins=above.tolist(),
    )


def run():
    prior = json.loads((static.HERE / "records.json").read_text())
    standard = next(s for s in prior["studies"] if s["contrast"] == 8)
    board = federation(1)
    route, _ = static.select_route(board)
    if list(route) != standard["route"] or any(
        neighbour not in set(board.neighbours[p]) - set(route)
        for p, neighbour in zip(route, standard["blank_neighbour"], strict=True)
    ):
        raise ValueError("Missing enrolled off-guide vacancy")
    size = len(route)
    n = size // 2
    amp = native_matrix(static.DWELL, 8)[0, 0]
    flights = {k: pair.transfers(n, k) for k in (0, 8)}
    rng = np.random.default_rng(SEED)
    rows = {}
    counts = {}
    for key in ("single", "wind0", "wind13"):
        if key == "single":
            admission = 1.0
            radius = None
        else:
            w = int(key.removeprefix("wind"))
            _, wave = pair.line(n, flights[8], w)
            state = wave(16)
            admission = float(np.sum(abs(state) ** 2) / 2)
            radius = next(
                p["reading"]["radius_interval"]
                for p in standard["panels"]
                if p["winding"] == w
            )
        reconstructed = {}
        for k in (8, 0):
            spectra = (
                single_spectrum(size, k, amp)
                if key == "single"
                else pair_spectrum(n, w, state, flights[k], amp)
            )
            for t in DELAYS:
                banks = acquire(spectra[t], admission, rng)
                for arm, bank in zip(
                    ("early", "late", "real", "imag"), banks, strict=True
                ):
                    counts[f"{key}_{k}_{t}_{arm}"] = bank
                reconstructed[k, t] = read_banks(banks, size)
        rows[key] = dict(
            radius_interval=radius,
            admission=admission,
            delays={
                t: compare(reconstructed[8, t], reconstructed[0, t], radius)
                for t in DELAYS
            },
        )
    return dict(
        guide_cells=board.cells,
        ports=size,
        shots=static.SHOTS,
        seed=SEED,
        one_sided_tail=TAIL,
        quantization_guard=GUARD,
        rows=rows,
    ), counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results, counts = run()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    np.savez_compressed(args.output / "counts.npz", **counts)
    for key, row in results["rows"].items():
        for t, r in row["delays"].items():
            print(
                key,
                t,
                r["first_above"],
                r["crossing_qR"],
                r["max_lower"],
                r["max_upper"],
                r["detected_bins"],
            )
