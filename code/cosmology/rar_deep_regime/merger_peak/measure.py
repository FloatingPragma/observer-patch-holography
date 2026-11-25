"""Compare elliptical two-centre convergence with both western and eastern peaks."""

import argparse
import json
import math
from pathlib import Path

from cosmology.rar_deep_regime import bullet_cluster

HERE = Path(__file__).resolve().parent


def difference(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def unit(v):
    length = math.hypot(*v)
    return tuple(x / length for x in v)


def build():
    archival = bullet_cluster.build()
    sky = archival["coordinates_east_north_arcsec"]
    west, east = sky["western_plasma"], sky["eastern_plasma"]
    west_peak = sky["western_lensing_peak"]
    # Section 3 expresses the eastern peak relative to each of two main BCGs.
    # Bracket both readings of the single main BCG position in Table 2.
    bcg = bullet_cluster.sky_offset(35.3, 56.3)
    offsets = {"southern": (-21.2, 17.7), "northern": (2.5, -11.5)}
    east_peaks = {
        label: tuple(round(x + y, 8) for x, y in zip(bcg, offset, strict=True))
        for label, offset in offsets.items()
    }
    # This common direction attains the smaller eastern separation.
    direction = unit(difference(east_peaks["southern"], east))
    west_vector = difference(west, west_peak)
    east_vector = difference(east, west_peak)
    delta = math.atan2(
        abs(west_vector[0] * east_vector[1] - west_vector[1] * east_vector[0]),
        dot(west_vector, east_vector),
    )
    west_gap = dot(west_vector, direction)
    east_gaps = {
        label: dot(difference(peak, east), direction) for label, peak in east_peaks.items()
    }
    margin = min(west_gap, *east_gaps.values())
    axis_limits = []
    for allowance in (5, 10, 20):
        width = delta + sum(
            math.asin(allowance / math.hypot(*v)) for v in (west_vector, east_vector)
        )
        axis_limits.append({
            "combined_peak_and_centre_allowance_arcsec": allowance,
            "upper_bound_on_smaller_axis_ratio": round(math.tan(width / 4), 10),
        })
    return {
        "source": archival["source"],
        "data": "Clowe et al., Section 3, Table 2; 4.413 kpc per arcsec in Section 1.",
        "scope": "Local strict peaks of the total projected convergence. Each component is nonincreasing in one elliptical radius about its centre, with its own fixed positive definite shape. A constant sheet is allowed.",
        "west_peak_east_north_arcsec": west_peak,
        "plasma_centres_east_north_arcsec": [west, east],
        "table_main_bcg_east_north_arcsec": bcg,
        "east_peak_by_table_bcg_east_north_arcsec": east_peaks,
        "common_direction_east_north": [round(x, 12) for x in direction],
        "projected_plasma_separation_arcsec": round(dot(difference(east, west), direction), 8),
        "western_peak_beyond_west_arcsec": round(west_gap, 8),
        "eastern_peak_beyond_east_arcsec": {k: round(v, 8) for k, v in east_gaps.items()},
        "minimum_combined_position_allowance_arcsec": round(margin, 8),
        "minimum_combined_position_allowance_kpc": round(4.413 * margin, 8),
        "west_angle_between_centres_degrees": round(math.degrees(delta), 10),
        "west_only_upper_bound_on_smaller_axis_ratio": round(math.tan(delta / 4), 14),
        "western_peak_position_sensitivity": axis_limits,
        "uncertainty": "Deterministic central-position bound. Uniform peak and uniform centre displacement bounds add. Assign any centroid-to-strict-maximum displacement to the peak bound. The published plasma-offset significance does not give the joint distribution of this four-position statistic.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = json.dumps(build(), indent=2, allow_nan=False) + "\n"
    if args.check:
        if json.loads(text) != json.loads((HERE / "results.json").read_text()):
            raise SystemExit("Published merger coordinates or projections differ from replay.")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
