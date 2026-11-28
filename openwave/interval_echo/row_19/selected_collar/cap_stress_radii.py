"""Test the fourth-power small-ball premise directly on the measured cap and flight energy."""

import argparse
from hashlib import sha256
import json

from . import distributed_source as dist

HERE, iv = dist.HERE, dist.iv
bounds = dist.lab.six.bounds


def run():
    iv.dps = 35
    cap = dist.read("local_cap")
    flight = dist.read("flight_energy")["result"]["groups"]["mixture"]
    radii = [0.5, 0.25, 0.125]
    rows = []
    for row in cap["balls"]:
        site = str(row["port"])
        local = row["readings"]["overlap"]
        q = iv.mpf(local["modular"])
        area_ratio = iv.mpf(local["area"]) / iv.mpf(row["base_area"])
        energy = iv.mpf(flight["late_minus_dense"]["overlap"][site])
        base = iv.mpf(flight["dense_energy"][site])
        # At fixed proper volume r'/r is the negative relative area rate.
        # On this cone that rate is independent of the enrolled radius.
        result = {}
        for frame, slope in (
            ("fixed_radius", -iv.pi * energy),
            ("fixed_volume", -iv.pi * (energy - base * area_ratio)),
        ):
            assert 0 not in q and 0 not in slope
            safe = min(0.5, bounds(abs(q) / (2 * abs(slope)))[0])
            scores = [q + r * slope for r in radii]
            result[frame] = {
                "coefficient_of_radius": bounds(slope),
                "remainder": [bounds(z) for z in scores],
                "candidate_stress": [
                    bounds(15 * z / (8 * iv.pi**2 * r**4))
                    for z, r in zip(scores, radii, strict=True)
                ],
                "radius_retaining_half_intercept": safe,
                "lower_abs_stress_times_radius_fourth": bounds(
                    15 * abs(q) / (16 * iv.pi**2)
                )[0],
                "zero_crossing_if_within_half": (
                    bounds(-q / slope)
                    if (-q / slope).a > 0 and (-q / slope).b <= 0.5
                    else None
                ),
            }
        rows.append(
            {
                "port": row["port"],
                "triangle": row["triangle"],
                "intercept": bounds(q),
                # For F(r)=r E(r), the signed combination at 1/2,1/4,1/8
                # cancels any central bulk term and is nonnegative for
                # every positive distribution through the rest-ball kernel.
                "three_radius_contrast": bounds(-q / 4),
                "minimum_negative_shell_energy": (
                    bounds(4 * q / (3 * iv.pi))[0] if q.a > 0 else 0.0
                ),
                **result,
            }
        )
    return {
        "radii": radii,
        "balls": rows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "local_cap.json",
                "flight_energy.json",
                "CONICAL_BALL_RADII.md",
                "LOCAL_CAP_PROTOCOL.md",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    path = HERE / "cap_stress_radii.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for row in data["balls"][:6]:
        print(
            row["port"],
            row["fixed_radius"]["zero_crossing_if_within_half"],
            row["fixed_radius"]["candidate_stress"],
        )
