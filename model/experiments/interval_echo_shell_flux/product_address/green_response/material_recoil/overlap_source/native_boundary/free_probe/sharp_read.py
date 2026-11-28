"""Finite occupation readings from independent simulated copies on every port."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import measure as m

SHOTS = 2**32
SEED = 190471


def acquire():
    ideal = m.run()
    _, _, dist, _ = m.make_board()
    rng = np.random.Generator(np.random.PCG64(SEED))
    eps = float(np.sqrt(np.log(2 * 4096 / 0.01) / (2 * SHOTS)))
    output = {}
    for name, data in [*ideal["written"].items(), ("ordinary", ideal["ordinary"])]:
        rows = []
        for probe in data["probes"]:
            if name != "ordinary":
                arms = [
                    probe["with_extra_erased"],
                    data["no_extra_erased"],
                    probe["with_extra_paired"],
                    data["no_extra_paired"],
                ]
            else:
                arms = [
                    probe["source_and_extra"],
                    data["source_alone"],
                    probe["extra_alone"],
                    np.zeros_like(probe["extra_alone"]),
                ]
            prob = np.array(arms)
            assert prob.min() >= -1e-12 and prob.max() <= 1 + 1e-12
            raw = rng.binomial(SHOTS, prob.clip(0, 1))
            contrast = np.einsum("a,atp->tp", [1, -1, -1, 1], raw) / SHOTS
            radius = contrast @ dist
            bound = float(4 * dist.sum() * eps)
            assert max(abs(radius - probe["radial_moment"])) < bound
            rows.append(
                dict(
                    origin=probe["origin"],
                    counts=raw.tolist(),
                    radial=radius.tolist(),
                    radius_99=bound,
                )
            )
        output[name] = rows
    return dict(
        protocol_sha256=sha256((m.HERE / "SHARP_READ.md").read_bytes()).hexdigest(),
        seed=SEED,
        shots_per_probability=SHOTS,
        raw_probability_radius_99=eps,
        enrolled_raw_means=sum(
            np.size(row["counts"]) for rows in output.values() for row in rows
        ),
        results=output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = acquire()
    path = m.HERE / "sharp_read.json"
    if args.write:
        path.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    else:
        m.check(json.loads(path.read_text()), result)
    for name, rows in result["results"].items():
        print(name, [(r["origin"], r["radial"][2], r["radius_99"]) for r in rows])
