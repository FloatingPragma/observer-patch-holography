"""Compare timed routes and freely moving material at identical modular source records."""

from fractions import Fraction as F
from math import fsum

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import parked_energy_charge as source
from . import parked_pair_flight as flight_law


def run():
    local = source.reference.local
    survey = local.survey
    board = survey.bulk.federation(0)
    read = local.inputs("parked_energy_charge")
    flight = local.inputs("parked_pair_flight")
    stations = flight["source_stations"]
    ref, fresh, shared = survey.standards(board)[0][:3]
    baseline = 256 * shared
    denominator = int(baseline.sum())
    # Integer left-null witnesses for the comparison records. Their choice
    # uses the source calibration. This retrospective test reports every
    # station pair and both previously defined free-flight generators.
    u = np.array([-1023, 0, 1088, -1088, 1023, -2080, 2080, 0, 0])
    v = np.eye(9, dtype=int)[3] - np.eye(9, dtype=int)[1]
    change = np.outer(u, v) + np.outer(v, u)
    candidate = baseline + change
    assert min(candidate.min(), baseline.min()) > 0
    assert not np.any(change.sum(axis=1)) and not np.any(change.sum(axis=0))
    _, energy, _ = source.station_records()
    classes = source.categorical(np.array(read["classes_in_calibration_order"]))
    bits = source.categorical(np.array(read["arrival_bit"]))
    station_energy = source.categorical(np.array(read["station_energy"]))
    assert not np.any(u @ energy) and not np.any(u @ classes) and not np.any(u @ bits)
    assert not np.any(v @ station_energy)
    for features in (energy, classes, bits, station_energy):
        assert not np.any(features.T @ change @ features)
    for value in np.unique(ref):
        assert not np.sum(change[ref == value])
    assert np.any(change) and not np.array_equal(shared, fresh)

    code = survey.routes.codebook()[0]
    calibration = survey.routes.calibrate(code)[0]
    active = survey.routes.graph(board, calibration)
    centers = np.arange(0, board.n, 12)
    upper = np.triu_indices(board.cells, 1)
    times = np.array(
        [
            shortest_path(
                active(*stations, i + 1, j + 1)[0], unweighted=True, indices=centers
            )[:, centers][upper]
            for i, j in np.ndindex(9, 9)
        ],
        int,
    )
    changed_histograms = {}
    changed_means = {}
    for k, (left, right) in enumerate(zip(*upper, strict=True)):
        key = f"{left}-{right}"
        timing = times[:, k]
        histogram = {
            str(int(t)): int(np.sum(change.ravel()[timing == t]))
            for t in np.unique(timing)
        }
        if any(histogram.values()):
            changed_histograms[key] = histogram
        shift = F(int(change.ravel() @ timing), denominator)
        if shift:
            changed_means[key] = str(shift)
    material = {}
    for key, row in flight["unit_time_joint_return"].items():
        # The numerical errors were bounded at each of the 81 inputs.
        vals = row["by_source_word_pair"]
        delta = (
            fsum(int(w) * y for w, y in zip(change.ravel(), vals, strict=True))
            / denominator
        )
        eps = (
            flight["single_return_probability_error"]
            * int(abs(change).sum())
            / denominator
            + 2e-15
        )
        material[key] = [delta - eps, delta + eps]
    # Independently test the changing return observable by exact closed
    # walks near the release ports. No finite-time exponential is needed.
    coefficients = [F(0)] * 5
    for i, j in zip(*np.nonzero(change), strict=True):
        marks = flight_law.marks_for(board, stations, (i, j), code)
        sharp = flight_law.return_series(
            board, marks, tuple(flight["release_and_return_ports"]), 4
        )
        for k, c in enumerate(sharp):
            coefficients[k] += F(int(change[i, j]), denominator) * c
    return {
        "source_stations": list(stations),
        "baseline_scale_on_shared": 256,
        "denominator": denominator,
        "u": u.tolist(),
        "v": v.tolist(),
        "minimum_preparation_count": int(candidate.min()),
        "joint_class_counts_on_both": (classes.T @ baseline @ classes).tolist(),
        "common_charge_relative_to_fresh": read["charge_change_nats"]["shared"][
            "exact"
        ],
        "relay_pairs_tested": len(upper[0]),
        "nonzero_relay_histogram_count_differences": changed_histograms,
        "nonzero_mean_relay_changes": changed_means,
        "unit_time_free_pair_return_change": material,
        "native_free_pair_even_taylor_coefficients_through_eighth": list(
            map(str, coefficients)
        ),
        "inputs": source.reference.coherent.steady.exp.save_inputs(
            "parked_energy_charge.json", "parked_pair_flight.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    source.reference.coherent.steady.exp.emit("matched_charge_transport", result)
    print(
        result["nonzero_mean_relay_changes"],
        result["unit_time_free_pair_return_change"],
    )
