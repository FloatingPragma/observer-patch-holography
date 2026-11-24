"""Stroboscopic endpoint separations and center currents of an interacting pair."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.linalg import schur

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30.experiment import (
    Guide,
    native_matrix,
    select_route,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**32
FILTER = 2**14
DETUNING = 0.2
THETA = (np.pi, np.pi - 2 * DETUNING)
# At most 2000 bounded scalar estimates across both laws and all three boards.
LOG_BUDGET = np.log(2 * 2000 / 0.01)


def canonical(n, a, b):
    """Orient the shorter arc; diameters begin at an even port on odd rings."""
    a, b = a % (2 * n), b % (2 * n)
    distance = (b - a) % (2 * n)
    if distance > n or (distance == n and a % 2):
        a, b = b, a
        distance = 2 * n - distance
    return 2 * (distance - 1) + a % 2, a // 2


def transfers(n, strength):
    if n % 2 != 1:
        raise ValueError("The enrolled guides have odd translation order")
    gates = [native_matrix(t, strength) for t in THETA]
    h0, h1 = candidate.native_law.local_term_parts()
    h = h0 + strength * h1
    for load in (0, 2):
        idx = candidate.native_law.SRC.index((load, load))
        np.testing.assert_allclose(h[:, idx], 0, atol=5e-14)
    entries = []
    for source in range(2 * n - 1):
        a = source % 2
        state = {(a, a + 1 + source // 2): 1.0}
        for parity, gate in enumerate(gates):
            nxt = defaultdict(complex)
            for (x, y), amp in state.items():
                xx = (x + (1 if x % 2 == parity else -1)) % (2 * n)
                yy = (y + (1 if y % 2 == parity else -1)) % (2 * n)
                if xx == y:
                    nxt[x, y] += amp
                    continue
                for r, f in ((x, gate[0, 0]), (xx, gate[1, 0])):
                    for s, g in ((y, gate[0, 0]), (yy, gate[1, 0])):
                        nxt[r, s] += amp * f * g
            state = nxt
        for (x, y), amp in state.items():
            dest, shift = canonical(n, x, y)
            entries.append((dest, source, shift, amp))
    i, j, shift, amp = zip(*entries)
    return np.array(i), np.array(j), np.array(shift), np.array(amp)


def matrix(n, entries, winding):
    i, j, shift, amp = entries
    u = np.zeros((2 * n - 1, 2 * n - 1), complex)
    np.add.at(u, (i, j), amp * np.exp(-2j * np.pi * winding * shift / n))
    return u


def line(n, entries, winding):
    u = matrix(n, entries, winding)
    tri, basis = schur(u, output="complex")
    values = np.diag(tri)
    np.testing.assert_allclose(tri, np.diag(values), atol=9e-13)
    angles = np.angle(values)
    nearest = abs(np.angle(values[:, None] / values))
    np.fill_diagonal(nearest, np.inf)
    gaps = nearest.min(axis=1)
    eligible = np.flatnonzero(gaps > 0.02)
    if not len(eligible):
        raise ValueError("No separated line")
    weights = (abs(basis[:2]) ** 2).sum(axis=0)
    peak = eligible[np.argmax(weights[eligible])]
    offset = np.angle(values / values[peak])
    window = (
        np.exp(0.5j * (FILTER - 1) * offset)
        * np.sinc(FILTER * offset / (2 * np.pi))
        / np.sinc(offset / (2 * np.pi))
    )
    # Each column is the accepted state for one adjacent input. The total
    # squared norm retains the actual probability of the accepted timing bin.
    coefficients = window[:, None] * basis[:2, :].conj().T

    def wave(t):
        return basis @ (values[:, None] ** t * coefficients)

    return {
        "phase": float(angles[peak]),
        "nearest_phase_gap": float(gaps[peak]),
        "selected_weight_given_tag": float(
            weights[peak] / np.sum(abs(coefficients) ** 2)
        ),
        "admission_by_orientation": np.sum(abs(coefficients) ** 2, axis=0).tolist(),
        "selected_admission_by_orientation": (abs(basis[:2, peak]) ** 2).tolist(),
    }, wave


def sample(p, rng):
    probs = np.r_[p, max(0.0, 1 - np.sum(p))]
    probs = np.maximum(0.0, np.round(probs, 13))
    return rng.multinomial(SHOTS, probs / probs.sum())


def estimate(counts, values):
    admitted = int(sum(counts[:-1]))
    v = np.asarray(values)
    mean = np.dot(counts[:-1] / admitted, v)
    err = np.ptp(v) * np.sqrt(LOG_BUDGET / (2 * admitted))
    return float(mean), float(err)


def read_panel(waves, n, rng, late):
    d = np.arange(2 * n - 1) // 2 + 1
    start = np.arange(2 * n - 1) % 2
    odd = start + (start + d) % 2
    p0 = np.mean(abs(waves(0)) ** 2, axis=1)
    p1 = np.mean(abs(waves(1)) ** 2, axis=1)
    lengths = {}
    for period in (0, 16) if late else (0,):
        p = np.mean(abs(waves(period)) ** 2, axis=1)
        counts = sample(np.bincount(d - 1, weights=p, minlength=n), rng)
        av, err = estimate(counts, np.arange(1, n + 1))
        lengths[str(period)] = {
            "counts": counts.tolist(),
            "mean_links": av,
            "bound": err,
        }
    rails = []
    for labels, p in ((odd, p0), (2 - odd, p0), (odd, p1)):
        rails.append(sample(np.bincount(labels, weights=p, minlength=3), rng))
    values, errs = zip(*(estimate(row, [0, 1, 2]) for row in rails))
    velocity = (2 * values[1] - values[0] - values[2]) / 4
    bound = (2 * errs[1] + errs[0] + errs[2]) / 4
    return {
        "separation": lengths,
        "odd_occupancy_counts": [v.tolist() for v in rails],
        "velocity_cells_per_period": velocity,
        "velocity_bound": bound,
    }


def fringe(n, winding, cache, rng, orientation):
    times = [0, n // 8]
    d = 1 + np.arange(2 * n - 1) // 2
    start = np.arange(len(d)) % 2
    # The midpoint is recorded to a quarter guide cell. Translating a pair
    # once adds one cell to its midpoint.
    centers = (4 * np.arange(n)[:, None] + 2 * start + d) % (4 * n)
    ks = (winding - 1, winding + 1)
    rows, phases = [], []
    for t in times:
        a, b = (cache[k][1](t)[:, orientation] for k in ks)
        ea, eb = (np.exp(2j * np.pi * k * np.arange(n) / n) for k in ks)
        state = (ea[:, None] * a + eb[:, None] * b) / np.sqrt(2 * n)
        p = np.bincount(
            centers.ravel(), weights=(abs(state) ** 2).ravel(), minlength=4 * n
        )
        counts = sample(p, rng)
        z = (counts[:-1] @ np.exp(4j * np.pi * np.arange(4 * n) / (4 * n))) / sum(
            counts[:-1]
        )
        rows.append(
            {"period": t, "center_counts": counts.tolist(), "visibility": float(abs(z))}
        )
        phases.append(np.angle(z))
    turn = np.angle(np.exp(1j * (phases[1] - phases[0])))
    cal = [cache[k][0] for k in ks]
    admission = sum(s["admission_by_orientation"][orientation] for s in cal)
    purity = (
        sum(s["selected_admission_by_orientation"][orientation] for s in cal)
        / admission
    )
    leakage = 2 * np.sqrt(max(0.0, 1 - purity))
    scale = n / (4 * np.pi * times[1])
    phase_bounds = []
    for row in rows:
        accepted = sum(row["center_counts"][:-1])
        error = 2 * np.sqrt(LOG_BUDGET / accepted) + leakage
        if error >= row["visibility"]:
            phase_bounds = []
            break
        phase_bounds.append(float(np.arcsin(error / row["visibility"])))
    phase_difference = np.angle(np.exp(1j * (cal[1]["phase"] - cal[0]["phase"])))
    return {
        "appointments": rows,
        "speed": float(turn * scale),
        "selected_line_fidelity": float(min(1, purity)),
        "ideal_phase_secant_speed": float(-phase_difference * n / (4 * np.pi)),
        "speed_error_including_filter": float(scale * sum(phase_bounds))
        if phase_bounds
        else None,
    }


def signal(board, route, k, rng):
    guide = Guide(board, route, DETUNING, k)
    records = []
    n = len(route) // 2
    for side in (0, 1):
        state = np.zeros(board.n, complex)
        state[route[side]] = 1
        for t in range(1, 13):
            guide.step(state)
            if t not in (6, 12):
                continue
            direction = 1 - 2 * side
            cell_dist = ((np.arange(n) * direction + n // 2) % n) - n // 2
            p = (abs(state[route]) ** 2).reshape(n, 2).sum(axis=1)
            counts = sample(
                [
                    p[cell_dist == t].sum(),
                    p[cell_dist > t].sum(),
                    p[cell_dist < t].sum(),
                ],
                rng,
            )
            records.append(
                {
                    "direction": direction,
                    "period": t,
                    "front_cells": t,
                    "front_beyond_inside_failed": counts.tolist(),
                }
            )
    return records


def contraction(panels):
    """Compare independent admissible mean-length and endpoint-current intervals."""
    rest = next(p["separation"]["0"] for p in panels if p["winding"] == 0)
    L0, e0 = rest["mean_links"], rest["bound"]
    comparisons = []
    for p in panels:
        if p["winding"] == 0:
            continue
        length = p["separation"]["0"]
        L, e = length["mean_links"], length["bound"]
        v, dv = abs(p["velocity_cells_per_period"]), p["velocity_bound"]
        low, high = (L - e) / (L0 + e0), (L + e) / (L0 - e0)
        vlo, vhi = max(0.0, v - dv), min(1.0, v + dv)
        forecast = [float(np.sqrt(1 - vhi**2)), float(np.sqrt(1 - vlo**2))]
        if not high < 1:
            raise ValueError("The enrolled pair contraction did not resolve")
        comparisons.append(
            {
                "winding": p["winding"],
                "length_ratio": L / L0,
                "length_ratio_interval": [low, high],
                "front_speed_prediction": float(np.sqrt(1 - v * v)),
                "prediction_interval": forecast,
                "compatible_limiting_speed": [
                    float(vlo / np.sqrt(1 - low**2)),
                    float(vhi / np.sqrt(1 - high**2)),
                ],
            }
        )
    return comparisons


def run():
    rng = np.random.default_rng(210517)
    followup = np.random.default_rng(210520)
    studies = []
    for level in range(3):
        board = federation(level)
        route, _ = select_route(board)
        n = len(route) // 2
        ms = sorted({round(n / r) for r in (32, 16, 8)})
        enrolled = sorted({0, *ms, *(-m for m in ms)})
        for strength in (8, 0):
            entries = transfers(n, strength)
            needed = sorted({*enrolled, *(m + j for m in ms for j in (-1, 1))})
            cache = {m: line(n, entries, m) for m in needed}
            panels = []
            for m in enrolled:
                diagnostic, waves = cache[m]
                panel = read_panel(waves, n, rng, m in (0, max(ms)))
                panels.append({"winding": m, **diagnostic, **panel})
            studies.append(
                {
                    "board_cells": board.cells,
                    "guide_cells": n,
                    "route": route,
                    "contrast": strength,
                    "front_reads": signal(board, route, strength, rng),
                    "vacuum_no_double_records": SHOTS,
                    "panels": panels,
                    "contraction": contraction(panels),
                    "first_orientation_center": {
                        str(m): fringe(n, m, cache, rng, 0) for m in ms
                    },
                    "second_orientation_center": {
                        str(m): fringe(n, m, cache, followup, 1) for m in ms
                    },
                }
            )
    return {
        "shots_per_appointment": SHOTS,
        "filter_periods": FILTER,
        "studies": studies,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(result, json.loads((HERE / "results.json").read_text()), "pair")
        print(
            "Reproduced six guided-pair panels, both signal fronts and all endpoint reads."
        )


if __name__ == "__main__":
    main()
