"""Compare NGC 2841's outer stellar inventory with spherical orbital excess mass."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
import math

from . import newton_transport as transport

budget = transport.budget
HERE = budget.HERE
GALAXY = "NGC2841"


def run():
    inputs, _, meta, masses = budget.load_inputs()
    radio_path = budget.local.profiles.maps.HERE / "gas_image_inputs.json"
    gas = json.loads(radio_path.read_bytes())[GALAXY]
    mass = Q(str(masses[GALAXY])) * 10**9
    distance = Q(str(meta[GALAXY]["dist_mpc"]))
    optical_scale = distance / Q("14.1")
    radio_scale = distance / Q(str(gas["radio"]["distance_mpc"]))
    # Zhang et al., section 4.1 (arXiv:1802.02583). The quoted mass is
    # integrated over the entire outer disk. No stellar profile is fitted here.
    stars = (140_000_000 + 2 * 20_000_000) * optical_scale**2
    atomic = Q(str(gas["atomic_He_msun"])) * radio_scale**2
    molecular = Q(str(gas["H2_He_msun"])) * radio_scale**2
    rows = [
        (r, lo * mass, hi * mass)
        for r, lo, hi in transport.mass_intervals(inputs[GALAXY], masses[GALAXY])
        if 30 * optical_scale < r <= 70 * optical_scale
    ]
    assert len(rows) >= 2 and all(a[0] < b[0] for a, b in zip(rows, rows[1:]))

    def witness(pair):
        (a, _, hi), (b, lo, _) = pair
        delta = lo - hi
        return dict(
            inner_kpc=float(a),
            outer_kpc=float(b),
            required_extra_mass_floor_msun=math.floor(delta),
            minimum_extra_per_outer_stellar_mass=math.floor(100 * delta / stars) / 100,
            minimum_extra_per_stars_and_all_measured_gas=math.floor(
                100 * delta / (stars + atomic + molecular)
            )
            / 100,
        )

    largest = max(combinations(rows, 2), key=lambda p: p[1][1] - p[0][2])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "OUTER_DISK.md",
                HERE / "source_profile.json",
                radio_path,
                budget.local.profiles.galaxies.sparc.DATA / "table1.dat",
            )
        },
        optical_source="https://arxiv.org/abs/1802.02583",
        adopted_distance_mpc=float(distance),
        outer_stars_upper_msun=float(stars),
        all_atomic_He_msun=float(atomic),
        all_molecular_He_msun=float(molecular),
        ring_intervals_msun=[
            [float(r), math.floor(lo), math.ceil(hi)] for r, lo, hi in rows
        ],
        strongest=witness(largest),
        first_to_last=witness((rows[0], rows[-1])),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "outer_disk.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    print(
        len(result["ring_intervals_msun"]), result["strongest"], result["first_to_last"]
    )
