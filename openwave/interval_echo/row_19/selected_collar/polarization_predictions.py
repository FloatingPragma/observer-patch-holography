"""Predict finite flights from observed source positions and ordinary control counts."""

from math import log, sqrt

import numpy as np

from . import full_board_polarization as source

motion = source.motion
local = source.local
SEED = 28122918
STARTS = 2**60


def run():
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = source.steady.symmetric_instrument(board)
    _, _, ports, gates, _, degrees, edges, _ = setup
    occupied = motion.second_instruments(ports, gates, degrees, edges)
    absent = source.empty_table(setup)
    receipt = local.inputs("full_board_polarization")
    # The guarded read has one outside class. Record every position before
    # grouping it. A gateway occupancy is retained as a blocked launch.
    arms = {}
    for label, reads in receipt["outcomes"].items():
        positions = {}
        for tag, arm in reads["arms"].items():
            counts = arm["position_counts"]
            selected = [counts[j] for j in ports]
            selected.append(sum(counts) - sum(selected))
            positions[tag] = np.array(selected) / sum(counts)
        arms[label] = positions["shared"] - positions["fresh"]

    # Each expectation is an independently measured, bounded laboratory
    # score in the model. Allocate 1% family error over the calibration
    # settings and positions in the two source arms, at both read times.
    allotment = 10000
    calibrations, predictions = [], {label: [] for label in arms}
    for it, dt in enumerate(motion.FLIGHTS):
        controls = np.concatenate([occupied[it], absent[it][None]], axis=0)
        settings = []
        scores = []
        for table in controls:
            draws = [
                rng.multinomial(STARTS, np.maximum(0, v) / np.maximum(0, v).sum())
                for v in table
            ]
            settings.append([d.tolist() for d in draws])
            occup = [
                float(d @ np.array([0.0, 1.0, 2.0, 0.0])) / STARTS for d in draws[:4]
            ]
            quadr = [float(int(d[0]) - int(d[1])) / STARTS for d in draws[4:]]
            scores.append(
                [
                    (occup[0] - occup[1]) / dt**2,
                    (occup[2] - occup[3]) / dt**2,
                    (quadr[1] - quadr[3]) / dt,
                ]
            )
        calibrations.append(settings)
        for label, delta in arms.items():
            reads = receipt["outcomes"][label]["arms"]
            position_half = sum(
                sqrt(log(2 * allotment / 0.01) / (2 * sum(arm["position_counts"])))
                for arm in reads.values()
            )
            bands = []
            for k, values in enumerate(np.array(scores).T):
                divisor = dt if k == 2 else dt**2
                calibration_half = (
                    4 * sqrt(log(2 * allotment / 0.01) / (2 * STARTS)) / divisor
                )
                # Apply Hoeffding to the fixed true class score. Its range
                # is enclosed by the independently calibrated range, and
                # replacing it in the empirical contrast costs ||delta||_1 e.
                half = (
                    (np.ptp(values) + 2 * calibration_half) * position_half
                    + np.abs(delta).sum() * calibration_half
                    + 1e-11
                )
                center = float(delta @ values)
                bands.append([center - half, center + half])
            predictions[label].append(bands)
    return {
        "times": list(motion.FLIGHTS),
        "seed": SEED,
        "calibration_starts": STARTS,
        "ports": ports,
        "calibration_counts": calibrations,
        "position_predictions": predictions,
        "inputs": source.exp.save_inputs(
            "full_board_polarization.json", "FULL_BOARD_POLARIZATION.md"
        ),
    }


if __name__ == "__main__":
    result = run()
    source.exp.emit("polarization_predictions", result)
    print(result["position_predictions"])
