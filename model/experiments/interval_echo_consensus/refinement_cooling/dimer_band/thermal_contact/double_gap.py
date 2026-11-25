"""Compare a two-unit material transition to two independently thermalized dimers."""

import argparse
import json
from itertools import product
from math import pi

import numpy as np
from scipy.linalg import expm
from ..measure import gates, pointer
from .measure import HERE, SHOTS, certify_sign_controls, proposals

TIMES = (0, pi / 6, pi / 4, pi / 2, 2 * pi / 3)
CLOCK_SHOTS = 65536


def experiment():
    pointer()
    certify_sign_controls()
    rng = np.random.default_rng(332735)
    h = gates("native")[0]
    clocks = []
    for t in TIMES:
        settings = []
        for pause in (0, pi / 2):
            p = abs(expm(-1j * (t + pause) * h)[:2, 0]) ** 2
            pair = np.round(np.kron(p, p), 13)
            settings.append(rng.multinomial(CLOCK_SHOTS, pair / pair.sum()).tolist())
        clocks.append(dict(time=t, counts=settings))
    banks = []
    armed = proposals(1, 1)  # both B dimers use the finite-cap repair
    source = proposals(3, 1)  # C uses separate directly filtered cohorts
    for copy, park in ((False, False), (True, False), (True, True)):
        # Replace C's second filter with an addressed copy of its first enrolled
        # sign in the foil. Initialize the second dimer at lower energy.
        locals_ = [armed, armed, source] + ([] if copy else [source])
        refused = []
        surviving = 1.0
        for q in locals_:
            refused.append(surviving * q[2])
            surviving *= q[:2].sum()
        terms, prob = [], []
        for a, b, c in product((0, 1), repeat=3):
            for offer in (0, 1):
                p = armed[a] * armed[b] * source[c] * (1 if copy else source[c]) / 2
                goes = bool(offer and a == b and not park)
                new = (c, c, a) if goes else (a, b, c)
                terms.append([a, b, c, offer, *new])
                prob.append(p)
        refused.append(0 if copy else 2 * armed[:2].sum() ** 2 * source[0] * source[1])
        p = np.array(prob + refused)
        np.testing.assert_allclose(p.sum(), 1, atol=3e-13, rtol=0)
        counts = rng.multinomial(SHOTS, p / p.sum()).tolist()
        banks.append(
            dict(
                copy=copy,
                parked=park,
                events=[
                    [*row, n]
                    for row, n in zip(terms, counts[: len(terms)], strict=True)
                ],
                refusal_counts=counts[len(terms) :],
            )
        )
    return dict(shots=SHOTS, clock_shots=CLOCK_SHOTS, clocks=clocks, contacts=banks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(experiment(), indent=2) + "\n"
    receipt = HERE / "double.json"
    if args.write:
        receipt.write_text(text)
    elif receipt.read_text() != text:
        raise SystemExit("Double-gap receipt has changed")
