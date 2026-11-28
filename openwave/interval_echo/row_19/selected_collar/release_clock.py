"""Can the four stationary radar records determine freely released motion?

Compare explicit continuum completions of the six reconstructed distances.
The clock gradient between the stations is an independent choice in each.
"""

import json
from hashlib import sha256

import numpy as np
from scipy.integrate import solve_ivp

from matched_einstein_limit.source_implication.material_clocks import exterior
from .measure import HERE


def compare(stations, slope):
    stations = np.asarray(stations)
    origin, target = stations[:2]
    separation = np.linalg.norm(target - origin)
    n = (target - origin) / separation
    # An affine function equal to one at the target and zero at the other
    # three stations. Its square keeps the reference station freely falling.
    affine = np.linalg.solve(
        np.column_stack([np.ones(4), stations]), np.array([0, 1, 0, 0])
    )

    def potential(x):
        return slope * (affine[0] + x @ affine[1:]) ** 2 * ((x - target) @ n)

    def gradient(x):
        weight = affine[0] + x @ affine[1:]
        along = (x - target) @ n
        return slope * (weight**2 * n + 2 * weight * along * affine[1:])

    def flight(t, y):
        # Geodesic of exp(2 phi)(-dt^2 + dx^2) in optical coordinates,
        # parameterized by the reference clock's time.
        v = y[3:]
        return np.r_[v, -(1 - v @ v) * gradient(y[:3])]

    times = [0.04, 0.08]
    trajectory = solve_ivp(
        flight,
        (0, times[-1]),
        np.r_[target, np.zeros(3)],
        t_eval=times,
        rtol=1e-12,
        atol=1e-14,
    )
    assert trajectory.success
    distances = np.linalg.norm(trajectory.y[:3].T - origin, axis=1)
    midpoint_acceleration = 2 * (distances - separation) / np.square(times)
    # Two additional supported clocks report their proper ticks per reference
    # tick. Their position along this line is independently tagged by radar.
    h = separation / 8
    rate = np.exp([potential(target + s * h * n) for s in (-1, 1)])
    clock_slope = (np.log(rate[1]) - np.log(rate[0])) / (
        2 * h * (1 + (h / separation) ** 2)
    )
    station_rate = np.exp(potential(stations))
    lengths = np.linalg.norm(stations[:, None] - stations[None, :], axis=2)
    return dict(
        log_rate_gradient=slope,
        stationary_round_trips=(2 * station_rate[:, None] * lengths).tolist(),
        added_clock_rates=rate.tolist(),
        predicted_release_acceleration_from_added_clocks=-clock_slope,
        release_times=times,
        release_accelerations_from_radar=midpoint_acceleration.tolist(),
    )


def run():
    path = HERE / "source_area.json"
    source = json.loads(path.read_bytes())
    counter = HERE / "reference_read.json"
    q = json.loads(counter.read_bytes())["cases"][str(2**38)]["readings"][
        "equilibrium"
    ]["increment"]
    experiments = []
    for tape in source["flat_realizations_of_symmetric_ranges"]:
        stations = np.array(tape["stations"])
        length = np.linalg.norm(stations[1] - stations[0])
        cases = [compare(stations, k * q / length) for k in (-1, 0, 1)]
        for c in cases:
            np.testing.assert_allclose(
                c["stationary_round_trips"],
                cases[1]["stationary_round_trips"],
                rtol=0,
                atol=1e-13,
            )
            np.testing.assert_allclose(
                c["predicted_release_acceleration_from_added_clocks"],
                -c["log_rate_gradient"],
                rtol=0,
                atol=1e-14,
            )
            np.testing.assert_allclose(
                c["release_accelerations_from_radar"],
                -c["log_rate_gradient"],
                rtol=2e-6,
                atol=1e-12,
            )
        experiments.append(cases)
    # Unlike the continuum comparison above these clocks use the native
    # Hamiltonian. Integer columns check their common gap for every input
    # head vector, so both independently prepared source endpoints qualify.
    board = exterior.federation(1)
    route = exterior.native.iso.pairing.formation.region(board)
    clocks, ranges = exterior.clock_cells(board, route)
    for charge in (-1, 1):
        for odd in (False, True):
            hamiltonian = exterior.native.generator(board, route, charge, odd)
            for i, pair in clocks.items():
                exterior.verify_columns(hamiltonian, i, pair)
    return dict(
        area_sha256=sha256(path.read_bytes()).hexdigest(),
        counts_sha256=sha256(counter.read_bytes()).hexdigest(),
        native_exterior_clocks=dict(
            cells=board.cells,
            clock_cells=sorted(clocks),
            relay_hops={str(i): ranges[i] for i in sorted(clocks)},
            common_angular_gap="sqrt(5)",
            signed_endpoint_gap_change=0,
            integer_sector_checks=4 * len(clocks),
        ),
        experiment="kinematic completions of symmetric records, optical speed one",
        resolution="full survey at two endpoints, two additional ideal rate standards",
        signs=[-1, 0, 1],
        endpoints=experiments,
    )


if __name__ == "__main__":
    import argparse

    cmd = argparse.ArgumentParser(description=__doc__)
    cmd.add_argument("--write", action="store_true")
    opts = cmd.parse_args()
    result = run()
    path = HERE / "release_clock.json"
    if opts.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        [
            c["predicted_release_acceleration_from_added_clocks"]
            for c in result["endpoints"][0]
        ]
    )
