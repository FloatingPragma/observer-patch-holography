"""Read the diffusion of one mobile energy excitation around four guarded registers."""

import argparse
import json

import numpy as np

from .measure import BASE, HERE, filter_weights, route_certificate

TRIALS = 2**24
TICKS = 8


def prepare(rng):
    counts = []
    weights = filter_weights(BASE)
    for site in range(4):
        wanted = int(site == 0)
        failures = rng.negative_binomial(TRIALS, weights[wanted])
        alternatives = weights.copy()
        alternatives[wanted] = 0
        bins = rng.multinomial(failures, alternatives / alternatives.sum())
        bins[wanted] = TRIALS
        counts.append(bins.tolist())
    return counts


def run():
    routes = route_certificate()
    required = {(0, 1), (1, 2), (2, 3), (0, 3)}
    if not required <= {tuple(r["pair"]) for r in routes}:
        raise ValueError("Missing exchange route on the diffusion ring")
    rng = np.random.default_rng(398317)
    banks = []
    for parked in (False, True):
        sources = prepare(rng)
        # The key stores just what the observer needs to count signed moves.
        pop = {(0, 0): TRIALS}
        steps = []
        for tick in range(TICKS):
            dest = {}
            offers = np.zeros(8, dtype=np.int64)
            signed = np.zeros((3, 3), dtype=np.int64)
            for (first, displacement), subtotal in sorted(pop.items()):
                owner = displacement % 4
                tags = rng.multinomial(subtotal, np.full(8, 1 / 8))
                offers += tags
                for tag, num in enumerate(tags):
                    receiver = owner
                    if not parked and tag < 4:
                        a, b = tag, (tag + 1) % 4
                        if owner in (a, b):
                            receiver = b if owner == a else a
                    step = (receiver - owner + 1) % 4 - 1
                    first_read = first if tick else step
                    key = (first_read, displacement + step)
                    dest[key] = dest.get(key, 0) + int(num)
                    signed[first_read + 1, step + 1] += num
            pop = dest
            bins = [
                sum(n for (_, x), n in pop.items() if x == i)
                for i in range(-tick - 1, tick + 2)
            ]
            steps.append(
                dict(
                    displacements=bins,
                    first_and_current=signed.tolist(),
                    offers=offers.tolist(),
                )
            )
        banks.append(dict(parked=parked, original_proposals=sources, steps=steps))
    return dict(trials=TRIALS, ticks=TICKS, banks=banks)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    path = HERE / "ring.json"
    s = json.dumps(run(), separators=(",", ":")) + "\n"
    if args.write:
        path.write_text(s)
    elif path.read_text() != s:
        raise SystemExit("Diffusion record mismatch")
