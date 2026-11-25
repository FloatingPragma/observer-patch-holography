"""Freeze the ordinary calibration, then compare clock and transport on the collar."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import readout
from .source import measure as sources

HERE = Path(__file__).resolve().parent


def transfer(panel, run):
    if (
        panel["cells"] != run["cells"]
        or run["epsilon"] != readout.EPSILON
        or run["omega"] != readout.OMEGA
    ):
        raise ValueError("source and readout apparatus disagree")
    arms = run["arms"]
    ref = np.r_[1, arms["reference"]["bloch"]]
    ordinary, collar = [
        np.r_[0, np.asarray(arms[arm]["bloch"]) - ref[1:]]
        for arm in ("ordinary", "decorrelated")
    ]
    ratio = (
        arms["decorrelated"]["modular_nats"] - arms["reference"]["modular_nats"]
    ) / (arms["ordinary"]["modular_nats"] - arms["reference"]["modular_nats"])

    def compare(effect):
        observed = float(np.dot(effect, collar))
        calibrated = float(np.dot(effect, ordinary))
        prediction = float(calibrated * ratio)
        return dict(
            ordinary=calibrated,
            collar=observed,
            prediction=prediction,
            residual=observed - prediction,
        )

    # Initial normalized frequency on the isolated bridge. All three initial
    # rates depend on this one scalar, with distinct spatial sensitivities.
    multiplier = np.array([np.cosh(readout.EPSILON), 0, 0, -np.sinh(readout.EPSILON)])
    local = compare(multiplier)
    crossings = []
    for appointment in panel["appointments"]:
        if appointment["time"] <= 1 / 8:
            continue
        row = dict(time=appointment["time"])
        for key in ("band_position", "band_clock"):
            signed = np.array(appointment[key])
            slope = np.dot(appointment[key + "_derivative"], ref)
            # Change of log inverse duration at the *same* reference reading,
            # by local implicit differentiation with respect to source mixing.
            if abs(slope) < 1e-10:
                raise ValueError("a turning point cannot calibrate this crossing")
            row[key] = dict(
                reference_level=float(signed @ ref),
                reference_time_slope=float(slope),
                **compare(signed / (appointment["time"] * slope)),
            )
        crossings.append(row)
    return dict(
        cells=run["cells"],
        forming_charge=run["forming_charge"],
        pump_charge=run["pump_charge"],
        pump_time=run["time"],
        collar_to_ordinary_modular_ratio=ratio,
        initial_bridge_multiplier=local,
        band_crossings=crossings,
    )


def run():
    sources_record = json.loads((sources.HERE / "results.json").read_text())
    panels = [readout.panel(level) for level in (0, 1)]
    _, bond, g, u, f, energy, vals, _, ends = readout.modes(2)
    largest = dict(
        cells=320,
        bridge=list(ends),
        slow_gaps=vals.tolist(),
        rates=readout.initial_reads(bond, g, u, f, energy, ends),
    )
    return dict(
        epsilon=readout.EPSILON,
        omega=readout.OMEGA,
        coefficients=["identity", "X", "Y", "Z"],
        readout_panels=panels,
        independent_320_cell_release=largest,
        source_comparisons=[
            transfer(panels[r["level"]], r) for r in sources_record["runs"]
        ],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
