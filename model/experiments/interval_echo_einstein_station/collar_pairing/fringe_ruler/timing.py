"""Independent material-clock uncertainty on coherent preparation and flight."""

import argparse
import json
from math import pi

from mpmath import iv
import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.common_receiver import (
    clock_transfer as clock,
)
from . import measure

REPEATS = (1, 8, 1000, 10**6, 10**9, 10**12)


def audit(records=None):
    iv.dps = 55
    rng = np.random.default_rng(823518)
    output = {}
    errors = {}
    for name, duration in (("flight_tick", iv.mpf(1)), ("phase_tick", iv.pi / 96)):
        prior = (0, pi)
        rows = []
        for i, n in enumerate(REPEATS):
            weight = clock.canonical_probability(duration, n)
            hits = (
                int(rng.binomial(clock.SHOTS, weight / clock.GRID))
                if records is None
                else records[name]["counts"][i]
            )
            prior = clock.refine(prior, n, hits, clock.SHOTS)
            rows.append(hits)
        errors[name] = max(max(clock.bounds(abs(iv.mpf(v) - duration))) for v in prior)
        output[name] = dict(counts=rows, single_firing_interval=prior)
    # Both formation and flight have four identical firings of the same
    # stable command.  Two erased marginals can cost twice the formation error.
    error = (
        192 * 4 * errors["flight_tick"] + 2 * (48 / 2 + 24 * 21) * errors["phase_tick"]
    )
    return dict(
        **output,
        hand_starts_each=clock.SHOTS,
        repeats=list(REPEATS),
        per_probability_error=float(np.nextafter(error, np.inf)),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = measure.HERE / "timing.json"
    if args.write:
        with path.open("x") as out:
            out.write(json.dumps(audit(), indent=2) + "\n")
    else:
        stored = json.loads(path.read_text())
        assert audit(stored)["per_probability_error"] < 1e-12


if __name__ == "__main__":
    main()
