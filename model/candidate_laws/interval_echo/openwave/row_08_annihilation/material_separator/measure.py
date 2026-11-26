"""Acquire the frozen trials or decode the archived native and comparison words."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from . import apparatus, neutral_field

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
CLOCK_SHOTS = 2**22
SHAM, WALL = 0, 1


def draw(rng, probabilities):
    """Keep all outcome columns. The final multinomial bin is the largest."""
    if probabilities.min() < -1e-12:
        raise ValueError("negative Born probability")
    np.testing.assert_allclose(probabilities.sum(axis=-1), 1, atol=5e-11, rtol=0)
    out = np.zeros(probabilities.shape, dtype=np.int64)
    for key in np.ndindex(probabilities.shape[:-1]):
        p = np.maximum(0, probabilities[key])
        order = np.argsort(p, kind="stable")
        sorted_p = p[order] / p.sum()
        sorted_p[-1] = 1 - sum(sorted_p[:-1])
        out[key][order] = rng.multinomial(SHOTS, sorted_p)
    return out


def acquire():
    rng = np.random.default_rng(803752)
    h, initial, ordered, omega = apparatus.setup()
    ticks = np.arange(257) / 4
    raw = {
        "hand": rng.binomial(
            CLOCK_SHOTS, (1 + np.stack((np.cos(ticks), np.sin(ticks)))) / 2
        ),
        "routes": apparatus.geometry(),
        "ordered_source_words": ordered,
        "receiver_gap": omega,
    }
    for name, sides in apparatus.ARMS.items():
        laws, energies = apparatus.native(sides, h, initial)
        raw[name + "_probabilities"] = laws
        raw[name + "_counts"] = draw(rng, laws)
        raw[name + "_energy"] = energies
        if len(sides) == 1:
            sparse_state = apparatus.sparse_endpoint(
                apparatus.frozen(h, sides), initial, 64
            )
            error = abs(apparatus.panels(sparse_state[None, :, :])[0] - laws[-1]).max()
            if error > 2e-10:
                raise ValueError(f"independent unilateral propagation error {error}")
            raw[name + "_independent_error"] = error
    local = {}
    for middle, cut in product((0, 1), (False, True)):
        local[middle, cut], energy, hop = neutral_field.calibration(omega, middle, cut)
        raw[f"field_input_energy_{middle}_{int(cut)}"] = energy
    raw["field_hopping"] = hop
    for a, b in product((0, 1), repeat=2):
        for name, cuts in apparatus.ARMS.items():
            left, right = local[a, 0 in cuts], local[b, 1 in cuts]
            law = np.stack(
                (
                    (1 - left) * (1 - right),
                    (1 - left) * right,
                    left * (1 - right),
                    left * right,
                ),
                axis=-1,
            )
            prefix = f"field_{a}{b}_{name}_"
            raw[prefix + "probabilities"] = law
            raw[prefix + "counts"] = draw(rng, law)
    return raw


def observations(measured, levels, modifier=1):
    diff = measured[:, WALL, :] - measured[:, SHAM, :]
    return [
        {
            "time": int(t),
            "left_excess": float(v[2] + v[3]),
            "right_excess": float(v[1] + v[3]),
            "both_excess": float(v[3]),
            "largest_bin_difference": float(max(abs(v))),
            "receiver_energy_excess": float(modifier * (v[1] + v[2] + 2 * v[3])),
        }
        for t, v in zip(levels, diff, strict=True)
    ]


def decode(raw):
    phase = np.unwrap(np.arctan2(*(2 * raw["hand"][::-1] / CLOCK_SHOTS - 1)))
    omega = float(raw["receiver_gap"])
    summary = dict(
        shots=SHOTS,
        clock_shots=CLOCK_SHOTS,
        routes=np.asarray(raw["routes"]).tolist(),
        receiver_gap=omega,
        clock_step=float(np.polyfit(np.arange(len(phase)), phase, 1)[0]),
        clock_reads=phase[4 * apparatus.TIMES].tolist(),
        field_hopping=float(raw["field_hopping"]),
        simultaneous_difference_halfwidth=float(
            np.sqrt(np.log(2 * 1024 / 0.01) / SHOTS)
        ),
        measurement="Seeded draws from calculated finite Born probabilities. Bounds describe simulated sampling.",
        native={},
        field={},
    )
    ordered = np.asarray(raw["ordered_source_words"])
    for name in apparatus.ARMS:

        def reduce(arr):
            shaped = arr.reshape(len(apparatus.TIMES), 2, 16, 4, 4)
            marginal = shaped.sum(axis=(2, 3))
            flagged = shaped[..., 3][:, :, ordered, :].sum(axis=(2, 3))
            recovery = shaped[:, :, ordered].sum(axis=(2, 3, 4))
            return marginal, flagged, recovery

        count, flag, recovery = reduce(raw[name + "_counts"] / SHOTS)
        law, born_flag, born_recovery = reduce(raw[name + "_probabilities"])
        obs = observations(count, apparatus.TIMES, omega)
        for row, observed, pred, source_read, source_law in zip(
            obs, flag, born_flag, recovery, born_recovery, strict=True
        ):
            row["ordered_both_by_source"] = observed.tolist()
            row["ordered_both_born"] = pred.tolist()
            row["ring_recovery_by_source"] = source_read.tolist()
            row["ring_recovery_born"] = source_law.tolist()
        summary["native"][name] = dict(
            observed=obs,
            born=observations(law, apparatus.TIMES, omega),
            initial_source_surplus=float(np.diff(raw[name + "_energy"])[0]),
        )
        if name in ("left", "right"):
            summary["native"][name]["sparse_check_max_abs_error"] = float(
                raw[name + "_independent_error"]
            )
    for a, b in product((0, 1), repeat=2):
        cases = {}
        for name, cuts in apparatus.ARMS.items():
            prefix = f"field_{a}{b}_{name}_"
            launched = [
                np.diff(raw[f"field_input_energy_{pin}_{int(side in cuts)}"])[0]
                for side, pin in enumerate((a, b))
            ]
            cases[name] = dict(
                source_work=list(map(float, launched)),
                observed=observations(
                    raw[prefix + "counts"] / SHOTS, neutral_field.TIMES, omega
                ),
                born=observations(
                    raw[prefix + "probabilities"], neutral_field.TIMES, omega
                ),
            )
        summary["field"][f"{a}{b}"] = cases
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire", action="store_true")
    opts = parser.parse_args()
    if opts.acquire:
        raw = acquire()
        np.savez_compressed(HERE / "counts.npz", **raw)
    with np.load(HERE / "counts.npz") as raw:
        result = decode(raw)
    if opts.acquire:
        (HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    for key, arm in result["native"].items():
        print("native", key, arm["observed"][-1])
    for key, arm in result["field"]["00"].items():
        print(
            "field", key, [row for row in arm["observed"] if row["time"] in (64, 128)]
        )
