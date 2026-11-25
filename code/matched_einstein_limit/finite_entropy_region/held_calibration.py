"""Compare boundary-word capacity to proper-area counts held out of calibration."""

import argparse
from hashlib import sha256
import json
from math import ceil, cos, floor, log, log1p, pi, sin, sqrt, tanh
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, gamma

from ..fixed_volume.measure import spatial_measures
from .measure import COUNTER_BUDGET, TILES

HERE = Path(__file__).resolve().parent
RECORD = HERE / "held_calibration.json"


def rounded(x):
    return round(float(x), 9)


def outward(lo, hi):
    # Seven decimal places, with a guard unit for the numerical quantiles.
    return [(floor(1e7 * lo) - 1) / 1e7, (ceil(1e7 * hi) + 1) / 1e7]


def target_area_bands(sectors, volume):
    # Two tails per event, three gamma stops and three Poisson surface counts.
    tail = 0.01 / (2 * 6)
    n = COUNTER_BUDGET
    lower_volume = volume * gamma.ppf(tail, n) / n
    upper_volume = volume * gamma.isf(tail, n) / n
    result = []
    for sector in sectors:
        k = sector["curvature"]
        y = sector["record"]["regions"]["volume_threshold"]["surface_count"]
        lower_radius, upper_radius = [
            brentq(lambda r: spatial_measures(k, r)[0] - v, 0, 1.5)
            for v in (lower_volume, upper_volume)
        ]
        # dA/dV=2 S'_k/S_k decreases with radius. These charts remain on
        # its positive branch across the full stopping-volume interval.
        if k > 0:
            assert sqrt(k) * upper_radius < pi / 2
            angle = sqrt(k) * lower_radius
            slope = 2 * sqrt(k) * cos(angle) / sin(angle)
        elif k < 0:
            slope = 2 * sqrt(-k) / tanh(sqrt(-k) * lower_radius)
        else:
            slope = 2 / lower_radius
        lower_area = (4 * pi / n) * chi2.ppf(tail, 2 * y) / 2
        upper_area = (4 * pi / n) * chi2.isf(tail, 2 * (y + 1)) / 2
        result.append(
            (
                lower_area - slope * (upper_volume - volume),
                upper_area + slope * (volume - lower_volume),
            )
        )
    return result


def three_sectors(source):
    areas = [
        4
        * pi
        * s["record"]["regions"]["volume_threshold"]["surface_count"]
        / COUNTER_BUDGET
        for s in source["geometry"]
    ]
    # Every allowed word in each of the six sector/bit registers was seen.
    capacities = []
    word_entropy_error = 0.0
    for launch in source["observed"]:
        counts = launch["joint_counts_by_sector_material_bit_edge_word"]
        offset = 0
        launch_capacities = []
        for n in TILES:
            for _ in range(2):
                block = np.array(counts[offset : offset + 2**n])
                assert len(block) == 2**n and min(block) > 0
                offset += len(block)
                frequencies = block / sum(block)
                log_words = log(len(block))
                word_entropy_error = max(
                    word_entropy_error,
                    abs(-frequencies @ np.log(frequencies) - log_words),
                )
                launch_capacities.append(log_words)
        assert offset == len(counts)
        if capacities:
            assert capacities == launch_capacities
        capacities = launch_capacities
    c3, c4, c5 = capacities[::2]
    da_cal, da_held = np.diff(areas)
    calibrated = (c4 - c3) / da_cal
    held = (c5 - c4) / da_held
    areas_ci = target_area_bands(source["geometry"], source["target_volume"])

    def delta(i, j):
        return areas_ci[j][0] - areas_ci[i][1], areas_ci[j][1] - areas_ci[i][0]

    first = delta(0, 1)
    second = delta(1, 2)
    lo, hi = (c4 - c3) / first[1], (c4 - c3) / first[0]
    return {
        "target_proper_volume": rounded(source["target_volume"]),
        "tile_counts": list(TILES),
        "areas": list(map(rounded, areas)),
        "word_capacities_nats": list(map(rounded, capacities[::2])),
        "maximum_sampled_word_entropy_error_nats": rounded(word_entropy_error),
        "eta_fit_3_to_4": rounded(calibrated),
        "eta_held_4_to_5": rounded(held),
        "held_capacity_minus_calibrated_area_nats": rounded(
            c5 - c4 - calibrated * da_held
        ),
        "simultaneous_99_percent": {
            "target_area_bands": [outward(a, b) for a, b in areas_ci],
            "eta_fit_3_to_4": outward(lo, hi),
            "eta_held_4_to_5": outward((c5 - c4) / second[1], (c5 - c4) / second[0]),
            "held_capacity_minus_calibrated_area_nats": outward(
                c5 - c4 - hi * second[1], c5 - c4 - lo * second[0]
            ),
        },
    }


def held_frames(panel, curvatures):
    ell = panel["ell"]
    center = curvatures.index(0)
    low, high = curvatures.index(-1), curvatures.index(1)

    def areas(frame):
        counts = np.array([bank["surface"] for bank in frame["records"]])
        return (4 * pi * ell**2 / panel["geometry_budget"]) * counts

    frames = panel["observers"]
    areas0 = areas(frames[0])
    launches = np.array(panel["terminal_sector_bit_counts"]["reversible"])
    assert frames[0]["frame_id"] == 0 and len(launches) == 2
    tiles = np.array(panel["full_tiles"])
    capacity = (tiles - tiles[center]) * log1p(1 / panel["ordinary_symbols"])
    denom = areas0[high] - areas0[low]
    estimates = {
        "terminal_counts": log(launches[0, high, 0] / launches[0, low, 0]) / denom,
        "constructed_alphabet": (capacity[high] - capacity[low]) / denom,
    }
    result = {}
    for estimator, eta in estimates.items():
        max_capacity, max_stationary = 0.0, 0.0
        for f in frames[1:]:
            af = areas(f)
            area_change = af - af[center]
            max_capacity = max(
                max_capacity, max(abs(capacity - eta * area_change)) / ell**4
            )
            for bank in launches:
                entropy = -np.log(bank / bank[center, 0])
                residual = entropy + eta * area_change[:, None]
                max_stationary = max(max_stationary, np.max(abs(residual)) / ell**4)
        result[estimator] = {
            "eta_fit_rest_frame": rounded(eta),
            "held_max_abs_capacity_minus_eta_area_over_ell4": rounded(max_capacity),
            "held_max_abs_entropy_plus_eta_area_over_ell4": rounded(max_stationary),
        }
    return result


def run():
    src = HERE / "results.json"
    refinement = HERE / "refinement.json"
    data = json.loads(refinement.read_text())
    finest = min(p["ell"] for p in data["panels"])
    return {
        "inputs_sha256": {
            p.name: sha256(p.read_bytes()).hexdigest() for p in (src, refinement)
        },
        "three_sectors": three_sectors(json.loads(src.read_text())),
        "finest_radius": finest,
        "calibration_curvatures": [-1, 1],
        "held_frame_ids": list(range(1, 10)),
        "held_frame_point_estimates": {
            p["branch"]: held_frames(p, data["curvature_sectors"])
            for p in data["panels"]
            if p["ell"] == finest
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    if parser.parse_args().write:
        RECORD.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        assert data == json.loads(RECORD.read_text())
        print("Replayed the three-sector calibration and nine held-frame comparisons.")
