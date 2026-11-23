"""Timelike Ricci from the volume of a freely falling radar cohort."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .. import observer_field as classical
from . import measure


def from_radar(tape):
    tick = tape["tick_unit"]
    step = tape["comoving_arm_spacing"]
    trip_pairs = [r["radar_send_return_ticks"] for r in tape["readings"]]
    lengths = np.diff(trip_pairs, axis=-1).squeeze(-1) * tick / (2 * step)
    times = np.array([r["proper_clock_ticks"] * tick for r in tape["readings"]])
    volumes = lengths.prod(axis=1)
    # The rounded appointments need not be exactly equidistant.
    slopes = np.diff(volumes) / np.diff(times)
    radial_acceleration = 2 * np.diff(slopes)[0] / (times[2] - times[0])
    rates = (lengths[2] - lengths[0]) / ((times[2] - times[0]) * lengths[1])
    cross = 2 * (rates[0] * rates[1] + rates[1] * rates[2] + rates[2] * rates[0])
    raw = -radial_acceleration / volumes[1]
    return dict(
        time=float(times[1]),
        volume=float(volumes[1]),
        volume_rate=float((volumes[2] - volumes[0]) / (times[2] - times[0])),
        raw=float(raw),
        expansion_correction=float(cross),
        ricci=float(raw + cross),
    )


def flat_tape(time, step):
    rows = []
    for offset in (-step, 0, step):
        t = time + offset
        distance = step * (1 - t)
        trip = [
            round((t - distance) / classical.TICK),
            round((t + distance) / classical.TICK),
        ]
        rows.append(
            dict(
                proper_clock_ticks=round(t / classical.TICK),
                radar_send_return_ticks=[trip] * 3,
            )
        )
    return dict(tick_unit=classical.TICK, comoving_arm_spacing=step, readings=rows)


def run(receipt):
    rows = []
    for survey in receipt["surveys"]:
        step = survey["spacing"]
        coupling = survey["calibration"]["coupling"]
        vacuum = survey["calibration"]["vacuum"]
        results = {}
        for name in ("y", "z", "even_x", "x", "decoupled_y"):
            comparisons = []
            for tape in survey["records"][name][1 if name == "x" else 0 :]:
                stress = measure.read(tape, receipt["recovered_mass2"])
                response = from_radar(tape["radar"])
                prediction = (
                    0.5 * coupling * (stress["energy"] + sum(stress["pressure"]))
                    - vacuum
                )
                comparisons.append(
                    dict(
                        **response,
                        stress_prediction=float(prediction),
                        error=float(abs(response["ricci"] - prediction)),
                    )
                )
            results[name] = comparisons
        flat = flat_tape(0.35, step)
        rows.append(
            dict(
                spacing=step,
                observations=results,
                flat_tape=flat,
                flat=from_radar(flat),
                max_held_error=max(
                    r["error"]
                    for n, arm in results.items()
                    if n != "decoupled_y"
                    for r in arm
                ),
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "VOLUME_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        source_record_sha256=hashlib.sha256(
            (measure.HERE / "results.json").read_bytes()
        ).hexdigest(),
        refinements=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = run(json.loads((measure.HERE / "results.json").read_text()))
    with args.output.open("x") as out:
        json.dump(results, out, indent=2, allow_nan=False)
        out.write("\n")
    for s in results["refinements"]:
        print(s["spacing"], s["max_held_error"], s["observations"]["y"][0], s["flat"])
