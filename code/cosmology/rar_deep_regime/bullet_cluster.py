"""Test a two-plasma-centre convergence model against the western Bullet peak.

Clowe et al. (2006), §3 and Table 2, DOI 10.1086/508162. Coordinates
are quoted J2000. The offset of the western peak is relative to the
subcluster BCG; both gas centres were located independently in X-rays.
This is a retrospective comparison of published summaries.
"""

import json
import math
from pathlib import Path


def sky_offset(ra_second, dec_second):
    """East/north arcseconds around 06h58m21.2s, −55°56′30.0″."""
    declination = math.radians(-(55 + 56 / 60 + 30 / 3600))
    return (
        round(15 * (ra_second - 21.2) * math.cos(declination), 8),
        round(30 - dec_second, 8),
    )


def segment_distance(point, start, end):
    direction = [y - x for x, y in zip(start, end)]
    displacement = [p - x for p, x in zip(point, start)]
    fraction = sum(x * y for x, y in zip(displacement, direction)) / sum(
        y * y for y in direction
    )
    clipped = min(1.0, max(0.0, fraction))
    projected = [x + clipped * y for x, y in zip(start, direction)]
    return fraction, math.dist(point, projected)


def build():
    western_gas = sky_offset(21.2, 30.0)
    eastern_gas = sky_offset(30.2, 35.9)
    western_bcg = sky_offset(16.0, 35.1)
    peak = (round(western_bcg[0] + 7.1, 8), round(western_bcg[1] + 6.5, 8))
    fraction, distance = segment_distance(peak, western_gas, eastern_gas)
    return {
        "source": "https://doi.org/10.1086/508162",
        "inputs": "Section 3 western convergence peak relative to the BCG; Table 2 gas and BCG centres; 4.413 kpc per arcsec in Section 1.",
        "domain": "Total convergence is a constant sheet plus two radially nonincreasing profiles centred on the tabulated plasma apertures. Flat sky over 2 arcmin.",
        "coordinates_east_north_arcsec": {
            "western_plasma": western_gas,
            "eastern_plasma": eastern_gas,
            "western_bcg": western_bcg,
            "western_lensing_peak": peak,
        },
        "segment_fraction_before_clipping_west_0_east_1": round(fraction, 6),
        "distance_from_allowed_peak_segment_arcsec": round(distance, 4),
        "distance_from_allowed_peak_segment_kpc": round(distance * 4.413, 3),
        "angular_separation_peak_to_western_bcg_arcsec": round(
            math.dist(peak, western_bcg), 4
        ),
        "central_peak_compatible": distance == 0,
        "uncertainty": "Authors bootstrap the peak and report about eight sigma from its plasma centre. Table 2 supplies no joint centroid covariance to compute a new significance for the segment.",
        "higher_resolution_comparison": "Cha et al. (2025), Section IV.2, report a roughly 150 kpc subcluster mass-peak to gas-core offset in a free-form JWST reconstruction. Their peak and centroid are distinct estimators.",
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.dumps(build(), indent=2, allow_nan=False) + "\n"
    if args.check:
        saved = Path(__file__).with_name("receipts") / "bullet_cluster.json"
        if json.loads(saved.read_text()) != json.loads(result):
            raise ValueError("Published merger comparison differs from the replay")
    else:
        print(result, end="")
