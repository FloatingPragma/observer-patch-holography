"""Predict two ordinary-calibrated images from the same measured overlap contrast."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from mpmath import iv
import numpy as np

from . import measure as meter
from ..screen_rulers import measure as rulers

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "transfer.json"


def inputs():
    stocks, steps, _ = rulers.preparations()
    ref = stocks["reference"]
    direction, ends, _, coeffs, info = meter.head.select_source(ref)
    steps["held"] = info["step"]
    stocks["held_minus"], stocks["held_plus"] = ends
    archive = json.loads(meter.RECEIPT.read_text())
    gate = archive["gates"]["0.0"]
    np.testing.assert_allclose(
        info["basis_coefficients"], archive["operator_audit"]["source_vector"]
    )
    np.testing.assert_equal(steps["held"], archive["operator_audit"]["selected_step"])
    np.testing.assert_allclose(
        coeffs, gate["pauli_coefficients"]["flagged"], atol=1e-13
    )
    measured = meter.decode(gate["counts"], gate["pauli_coefficients"], 0.0)["flagged"]
    contrast = iv.mpf(measured["estimate"]) + iv.mpf([-1, 1]) * measured["halfwidth"]
    generator = -meter.head.station.log_blocks(ref)
    pair = meter.head.station.pair
    ordinary = pair(generator, stocks["ordinary_plus"] - stocks["ordinary_minus"]) / (
        2 * steps["ordinary"]
    )
    if abs(ordinary) < 1e-12:
        raise ValueError("Unresolved ordinary source calibration")
    np.testing.assert_allclose(
        pair(generator, ends[1] - ends[0]), gate["expected"]["flagged"], atol=1e-14
    )
    np.testing.assert_allclose(meter.head.marginals(ends[1] - ends[0]), 0, atol=1e-13)
    prediction = contrast / (2 * steps["held"] * ordinary)
    return (
        stocks,
        steps,
        prediction,
        {
            "held_contrast_nats": measured,
            "ordinary_modular_slope": float(ordinary),
            "population_source_ratio": float(pair(generator, direction) / ordinary),
            "predicted_ratio_interval": str(prediction),
        },
    )


def read_images(bank, steps, prediction, counted):
    reports = {}
    for probe in rulers.PROBES:
        panels = bank[probe]
        if set(panels) != set(rulers.STOCKS):
            raise ValueError("Incomplete preparation bank")
        areas = {
            key: rulers.panel_read(panel, counted)["image"]
            for key, panel in panels.items()
        }
        ordinary = (areas["ordinary_plus"] - areas["ordinary_minus"]) / (
            2 * steps["ordinary"]
        )
        held = (areas["held_plus"] - areas["held_minus"]) / (2 * steps["held"])
        ratio = None if 0 in ordinary else held / ordinary
        reports[probe] = {
            "ordinary_image_slope": str(ordinary),
            "held_image_slope": str(held),
            "held_over_ordinary": rulers.pack(ratio),
            "held_compatible_with_zero": 0 in held,
            "compatible_with_prediction": None
            if ratio is None
            else 0 in (ratio - prediction),
        }
    return reports


def run():
    stocks, steps, prediction, source = inputs()
    formation = rulers.camera.BEAM.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = rulers.camera.BEAM.screens(board, route)
    banks = {}
    for strength, seed in ((8, 191127), (0, 191128)):
        bank = rulers.collect(board, route, rings, stocks, strength, seed)
        for name, key, counted in (
            ("observed", "counts", True),
            ("calculated", "population", False),
        ):
            bank[name] = read_images(bank[key], steps, prediction, counted)
        banks[str(strength)] = bank
    return {
        "protocol_sha256": sha256((HERE / "TRANSFER.md").read_bytes()).hexdigest(),
        "source_counts_file_sha256": sha256(meter.RECEIPT.read_bytes()).hexdigest(),
        "source": source,
        "half_steps": steps,
        "shots_per_setting": rulers.SHOTS,
        "banks": banks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RECEIPT.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    elif result != json.loads(RECEIPT.read_text()):
        raise ValueError("Source or acquisition differs from archived transfer")
    print(
        json.dumps(
            {
                "source": result["source"],
                "observed": {
                    key: bank["observed"] for key, bank in result["banks"].items()
                },
            },
            indent=2,
        )
    )
