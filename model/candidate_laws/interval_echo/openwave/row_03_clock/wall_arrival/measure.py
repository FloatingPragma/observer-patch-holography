"""Enroll the frozen records, or decode the archived public-word count banks."""

import argparse
import json
from pathlib import Path

import numpy as np

from ...row_30_collective import experiment, packets
from .apparatus import Strand, source, word_law

HERE = Path(__file__).resolve().parent
SHOTS, RECOIL = 2**24, 2**27
DELAYS = np.arange(2048) / 100
EXPOSURES = np.arange(4)
LAWS = [
    (29, "native", 8, 8, 0.5),
    (49, "native", 8, 8, 0.5),
    (49, "half", 4, 4, 0.5),
    (49, "planted", 8, 0, 7.5),
]


def hand(rng, times):
    return rng.binomial(SHOTS, (1 + np.stack((np.cos(times), np.sin(times)))) / 2)


def acquire():
    experiment.native_checks()
    boards = json.loads((experiment.HERE / "results.json").read_text())["boards"]
    rng = np.random.default_rng(31602743)
    raw = {"clock": hand(rng, DELAYS), "recoil_clock": hand(rng, EXPOSURES)}
    for n in (29, 49):
        route = next(
            b["route"] for b in boards if b["seed"] == 312143 and b["dimers"] == n
        )
        raw[f"n{n}_route"] = route
        for i, q in enumerate((0.10, 0.20)):
            prefix = f"n{n}_kick{i}_"
            modes = [source(n, k) for k in (q, -q)]
            raw[prefix + "phase"] = np.array(
                [packets.inspect(f, rng)[1] for f in modes]
            )
            u, e, v = experiment.dual(n, 8)
            p = [[packets.occupations(f, u, e, v, t) for t in EXPOSURES] for f in modes]
            raw[prefix + "density"] = rng.binomial(RECOIL, np.clip(p, 0, 1))
    for n, name, j, after_j, after_h in LAWS:
        tag = f"n{n}_{name}_"
        initial = Strand(n, j)
        near, far = [initial.enroll(source(n, mirror=m)) for m in (False, True)]
        probs = np.array([(1 + g[0, n + 1]) / 2 for g in (near, far)])
        raw[tag + "arrivals"] = rng.binomial(SHOTS, probs)
        bank = word_law(near, Strand(n, after_j, after_h), DELAYS)
        raw[tag + "words"] = rng.multinomial(SHOTS, bank)
    return raw


def read_hand(raw):
    x, y = 2 * raw / SHOTS - 1
    return np.unwrap(np.arctan2(y, x))


def read_phase(raw):
    n = packets.N_SHOTS
    diff = (raw[:, :, 0] - raw[:, :, 1]) / n
    contrast = diff.sum(axis=1) / 2
    variance = ((raw[:, :, :2].sum(axis=2) / n - diff**2) / (4 * n)).sum(axis=1)
    x, y = contrast
    se = np.sqrt(np.array([y * y, x * x]) @ variance) / (x * x + y * y)
    return float(np.arctan2(y, x)), float(se)


def read(raw):
    clock = read_hand(raw["clock"])
    times = read_hand(raw["recoil_clock"])
    times_centered = times - times.mean()
    slope_weights = times_centered / (times_centered @ times_centered)
    c2 = json.loads((experiment.HERE / "coupling_meter.json").read_text())["results"][
        "predicted_c_squared"
    ]
    summary = {
        "provenance": "Seeded multinomial and binomial draws from reduced finite-law Born probabilities. Fit resolution describes the sampler.",
        "word_order": "Lexicographic (trigger, first, seventh), sign order (-1,+1)",
        "shots": SHOTS,
        "recoil_shots": RECOIL,
        "phase_shots": packets.N_SHOTS,
        "material_coefficient_squared": c2,
        "clock_step": float(np.polyfit(np.arange(len(clock)), clock, 1)[0]),
        "recoil": [],
        "arrival_and_clock": [],
    }
    for n in (29, 49):
        # The frozen coordinate has its origin at the endpoint. Wall number can change.
        x = np.arange(n + 1)
        for i, q in enumerate((0.10, 0.20)):
            prefix = f"n{n}_kick{i}_"
            estimates = [read_phase(r) for r in raw[prefix + "phase"]]
            k = (estimates[0][0] - estimates[1][0]) / 2
            se_k = np.hypot(estimates[0][1], estimates[1][1]) / 2
            p = raw[prefix + "density"] / RECOIL
            shifts = (p[0] - p[1]) @ x / 2
            vel = float(slope_weights @ shifts)
            moment_var = (p * (1 - p)).sum(axis=0) @ (x * x) / (4 * RECOIL)
            hand_x, hand_y = 2 * raw["recoil_clock"] / SHOTS - 1
            time_var = (hand_y**2 * (1 - hand_x**2) + hand_x**2 * (1 - hand_y**2)) / (
                SHOTS * (hand_x**2 + hand_y**2) ** 2
            )
            dv_dt = ((shifts - shifts.mean()) - 2 * vel * times_centered) / (
                times_centered @ times_centered
            )
            se_v = float(np.sqrt((slope_weights**2) @ moment_var + dv_dt**2 @ time_var))
            mass = k / vel
            se_mass = mass * np.hypot(se_k / k, se_v / vel)
            summary["recoil"].append(
                dict(
                    dimers=n,
                    nominal_gradient=q,
                    measured_gradient=k,
                    velocity=vel,
                    finite_mass=mass,
                    mass_sampler_se=se_mass,
                    mass_times_material_c2=mass * c2["value"],
                    finite_product_sampler_se=se_mass * c2["value"],
                )
            )
    win = np.hanning(len(clock))
    demod = 2 * win * np.exp(-15j * clock) / win.sum()
    for n, name, *_ in LAWS:
        tag = f"n{n}_{name}_"
        arrivals = raw[tag + "arrivals"] / SHOTS
        words = raw[tag + "words"].reshape(-1, 2, 2, 2)
        counts = words.sum(axis=(2, 3))
        first = (words[:, :, 1, :] - words[:, :, 0, :]).sum(axis=2) / counts
        seventh = (words[:, :, :, 1] - words[:, :, :, 0]).sum(axis=2) / counts
        summary["arrival_and_clock"].append(
            dict(
                dimers=n,
                law=name,
                near_probability=float(arrivals[0]),
                mirror_probability=float(arrivals[1]),
                excess=float(arrivals[0] - arrivals[1]),
                excess_sampler_se=float(
                    np.sqrt((arrivals * (1 - arrivals)).sum() / SHOTS)
                ),
                minimum_count_per_trigger=counts.min(axis=0).tolist(),
                first_amplitude=abs(demod @ first).tolist(),
                seventh_amplitude=abs(demod @ seventh).tolist(),
            )
        )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args()
    if args.acquire:
        data = acquire()
        np.savez_compressed(HERE / "counts.npz", **data)
    with np.load(HERE / "counts.npz") as raw:
        summary = read(raw)
    result = json.dumps(summary, indent=2, allow_nan=False) + "\n"
    if args.acquire:
        (HERE / "results.json").write_text(result)
    print(result)
