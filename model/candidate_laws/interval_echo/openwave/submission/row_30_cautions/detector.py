"""Post-hoc audit of the local OW-30 detector and its archived counts.

Run as a module from the repository root. The default verifies detector.json.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.linalg import eigh

from ...row_30_collective import experiment as ex
from ...row_30_collective import holdout
from ...row_30_collective.rest.physics import propagation
from ...row_30_collective.rest.read import PAIRS, profile


def spatial_check(board, mass, speed_squared):
    """Replay the enrolled settings, checking the unreduced spatial profiles."""
    observation = board["observer"]
    n = board["dimers"]
    wave_checks = []
    probability_ranges = []
    draw_counts = ex.counts

    def inspect_counts(signal, norm, rng, stream, key):
        p = (1 + signal / norm) / 2
        probability_ranges.append([float(np.min(p)), float(np.max(p))])
        return draw_counts(signal, norm, rng, stream, key)

    def inspect_projector(projector, **kwargs):
        values, vectors = eigh(projector)
        wave = vectors[:, -1]
        interior = wave[1:-1]
        neighbors = wave[:-2] + wave[2:]
        multiplier = (interior @ neighbors) / (interior @ interior)
        residual = np.linalg.norm(neighbors - multiplier * interior)
        wave_checks.append(
            [residual / np.linalg.norm(interior), values[-1] - values[-2]]
        )
        return eigh(projector, **kwargs)

    with (
        patch.object(ex, "eigh", inspect_projector),
        patch.object(ex, "counts", inspect_counts),
    ):
        replay = ex.spectroscopy(
            n, np.random.default_rng([board["seed"], board["level"], 47])
        )
    ex.verify(replay, observation)
    nominal = [
        dict(mode, k=np.pi * mode["spatial_launch"] / (n + 2))
        for mode in observation["modes"]
    ]
    score = holdout.score(
        nominal,
        observation,
        board["native_transport"],
        board["changed_transport"],
        mass,
        speed_squared,
    )
    return dict(
        dimers=n,
        measured_lattice_max=board["score"]["lattice_max"],
        launch_label_lattice_max=score["lattice_max"],
        launch_label_shape_ratio=score["shape_rms_lattice"]
        / score["shape_rms_parabola"],
        max_interior_recurrence_residual=max(r[0] for r in wave_checks),
        minimum_projector_eigengap=min(r[1] for r in wave_checks),
        min_click_probability=min(r[0] for r in probability_ranges),
        max_click_probability=max(r[1] for r in probability_ranges),
    )


def packet_check(board, shots):
    raw = np.asarray(board["native"]["density_counts"])
    occupancy = raw.sum(axis=2) / shots
    address = np.arange(board["dimers"] + 1) - board["dimers"] / 2
    moments = raw @ address / shots
    times = board["read_times"]
    slopes = np.polyfit(times, moments, 1)[0]
    normalized = np.polyfit(times, moments / occupancy, 1)[0]
    difference = (slopes[0] - slopes[1]) / 2
    np.testing.assert_allclose(
        difference, board["native"]["velocity"], atol=1e-12, rtol=0
    )
    return dict(
        dimers=board["dimers"],
        wall_occupancy_range=[float(np.min(occupancy)), float(np.max(occupancy))],
        first_moment_slopes=slopes.tolist(),
        normalized_centroid_slopes=normalized.tolist(),
        half_signed_first_moment_slope=float(difference),
        half_signed_normalized_centroid_slope=float(
            (normalized[0] - normalized[1]) / 2
        ),
        forecast=board["envelope_velocity_forecast"],
    )


def source_instrument_counterexample():
    """Build the two physical measurements explicitly on the five-dimer Hilbert space."""
    n = 5
    words = np.arange(1 << n)
    z = 1 - 2 * ((words[:, None] >> np.arange(n)) & 1)
    walls = -z * np.roll(z, 1, axis=1)
    h = np.diag(-ex.K * walls.sum(axis=1))
    for i in range(n):
        h[words, words ^ (1 << i)] = -ex.H
    identity = np.eye(len(words))
    flip = identity[::-1]
    rho = (identity + flip) / len(words)

    def string(a, b):
        op = np.zeros(h.shape, dtype=complex)
        op[words ^ ((1 << b) - (1 << a)), words] = 1j * np.prod(walls[:, a:b], axis=1)
        return op

    source, terminal = string(0, 1), string(2, 3)
    endpoint_z = np.diag(z[:, -1])
    endpoint_y = source @ endpoint_z
    np.testing.assert_allclose(source @ flip, flip @ source)
    np.testing.assert_allclose(endpoint_z @ flip, -flip @ endpoint_z)
    np.testing.assert_allclose(endpoint_y @ flip, -flip @ endpoint_y)

    # Each signed map is the sum of (recorded sign) times the unnormalized
    # postmeasurement state. Terminal factor readings only need their product.
    coarse_signed = np.zeros_like(source)
    local_signed = np.zeros_like(source)
    after_local = np.zeros_like(source)
    for a in (1, -1):
        p = (identity + a * source) / 2
        coarse_signed += a * (p @ rho @ p)
        for b in (1, -1):
            p_local = (identity + a * endpoint_z) @ (identity + b * endpoint_y) / 4
            branch = p_local @ rho @ p_local
            local_signed += a * b * branch
            after_local += branch
    np.testing.assert_allclose(coarse_signed, source @ rho)
    np.testing.assert_allclose(local_signed, source / len(words))
    np.testing.assert_allclose(after_local, identity / len(words))
    values, basis = eigh(h)
    checks = []
    for t in (0.37, 2.15, 5.81):
        flight = (basis * np.exp(-1j * values * t)) @ basis.T
        read = flight.conj().T @ terminal @ flight
        coarse = np.trace(coarse_signed @ read).real
        local = np.trace(local_signed @ read).real
        sector_reads = []
        for f in (1, -1):
            u, energies = propagation(n, ex.K, f)
            kernel = (u * np.cos(energies * t)) @ u.T
            sector_reads.append(np.linalg.det(kernel[np.ix_([0, 1], [2, 3])]))
        np.testing.assert_allclose(
            [coarse, local], [sector_reads[0], np.mean(sector_reads)], atol=1e-12
        )
        checks.append(
            dict(t=t, joint_product=float(coarse), local_factors=float(local))
        )
    return dict(
        dimers=n,
        source=[0, 1],
        terminal=[2, 3],
        j=ex.K,
        flip_before=float(np.trace(flip @ rho)),
        flip_after_local=float(np.trace(flip @ after_local).real),
        responses=checks,
    )


def ring_raw_check(seed, tag, coupling):
    """Compare both instruments with existing binomial counts for every ring pair."""
    with np.load(ex.HERE / "rest" / f"counts_{seed}.npz", allow_pickle=False) as bank:
        n, shots = int(bank["dimers"]), int(bank["shots"])
        waves = profile(n)
        # These are probability audits of the two forward models on the original
        # generator's exposure grid. No observer time is inferred here.
        times = 0.08 * np.arange(len(bank["cos"]))
        means = {pair: [] for pair in PAIRS}
        selected = -((-1) ** n)
        for sector in (selected, -selected):
            u, energies = propagation(n, coupling, sector)
            evolved = np.cos(times[:, None] * energies)
            amplitudes = u.T @ waves
            for a, b in PAIRS:
                fa, fb = amplitudes[:, a], amplitudes[:, b]
                signal = (evolved @ (fa**2)) * (evolved @ (fb**2)) - (
                    evolved @ (fa * fb)
                ) ** 2
                scale = ex.wedge_norm(waves[:, a], waves[:, b])
                means[a, b].append(signal / scale**2)
        sum_sq = np.zeros(2)
        bins = 0
        for (a, b), (sector_mean, opposite_mean) in means.items():
            count = bank[f"{tag}_{a}{b}_pair"].astype(float)
            for i, mean in enumerate((sector_mean, (sector_mean + opposite_mean) / 2)):
                p = (1 + mean) / 2
                sum_sq[i] += np.sum((count - shots * p) ** 2 / (shots * p * (1 - p)))
            bins += len(count)
    rms = np.sqrt(sum_sq / bins)
    return dict(
        seed=seed,
        dimers=n,
        arm=tag,
        pair_time_bins=bins,
        shots_per_bin=shots,
        joint_product_binomial_residual_rms=float(rms[0]),
        local_factors_binomial_residual_rms=float(rms[1]),
    )


def audit():
    held = json.loads((ex.HERE / "holdout.json").read_text())
    packets = json.loads((ex.HERE / "packets.json").read_text())
    return dict(
        scope="Post-hoc local detector audit of archived settings. No official board settings used.",
        open_strands=[
            spatial_check(b, held["m"], held["c_squared"]) for b in held["boards"]
        ],
        packets=[packet_check(b, packets["shots"]) for b in packets["boards"]],
        ring_counterexample=source_instrument_counterexample(),
        ring_raw_counts=[
            ring_raw_check(seed, tag, j)
            for seed in (312142, 927016)
            for tag, j in (("native", ex.K), ("changed", ex.K / 2))
        ],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a new receipt at this path")
    args = cli.parse_args()
    result = audit()
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(result, indent=2) + "\n")
    else:
        ex.verify(result, json.loads(Path(__file__).with_suffix(".json").read_text()))
        print("Reproduced the open-strand, packet, and ring detector diagnostics.")
