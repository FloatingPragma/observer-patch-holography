"""Test a held one-gain current forecast after both matched sources enter a common survey."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import measure


def propagate(channels, effects, state, steps=2):
    geo = []
    for k in range(steps + 1):
        geo.append(
            measure.station.geometry(measure.station.pair(effects, state)).tolist()
        )
        if k < steps:
            state = np.einsum("bij,bj->bi", channels, state.reshape(2, 16)).reshape(
                2, 4, 4
            )
    return np.array(geo)


def study(channels, effects, receipt):
    rows = []
    uniform = channels.mean(axis=1)
    for mode, item in receipt["comparisons"].items():
        for offset in (0, 2):
            pair = item["endpoints"][offset : offset + 2]
            released = []
            sham = []
            for end in pair:
                weights = np.array(end["routes"])
                state, _ = measure.prepare(channels, weights)
                released.append(propagate(uniform, effects, state))
                sham.append(
                    propagate(
                        np.einsum("d,bdij->bij", weights, channels), effects, state
                    )
                )
            differential = released[1] - released[0]
            # All histories include every launch. Only the known trace drift of
            # the checked integer kernels can leak into a sham's readout.
            sham_drift = float(np.max(abs(np.diff(sham, axis=1))))
            j = (
                pair[1]["central_bridge_energy_current"]
                - pair[0]["central_bridge_energy_current"]
            )
            a, v = (differential[2] - 2 * differential[1] + differential[0]).tolist()
            rows.append(
                dict(
                    matched=mode,
                    radius=pair[1]["route_coordinate"],
                    energy_current_difference=j,
                    two_endpoint_geometry=[group.tolist() for group in released],
                    volume_second_difference_contrast=v,
                    area_second_difference_contrast=a,
                    maximum_stationary_sham_one_cycle_drift=sham_drift,
                )
            )
    gain = (
        rows[0]["volume_second_difference_contrast"]
        / rows[0]["energy_current_difference"]
    )
    for row in rows:
        predicted = gain * row["energy_current_difference"]
        row["forecast_volume_second_difference_contrast"] = predicted
        row["observed_minus_forecast"] = (
            row["volume_second_difference_contrast"] - predicted
        )
    return dict(
        calibration="area_volume at radius 1/64",
        common_evolution="uniform route mixture, fair test charge",
        cycles=[0, 1, 2],
        geometry_order=["area", "volume"],
        gain=gain,
        rows=rows,
    )


def run():
    return study(
        *measure.load_native(), json.loads((measure.HERE / "results.json").read_text())
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(out)
    else:
        print(out, end="")
