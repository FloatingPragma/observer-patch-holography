"""Finite heat-current records under guarded native energy exchanges."""

import argparse
from fractions import Fraction
from functools import cache
import json
from math import acos, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..measure import gates, outcomes, pointer
from ..thermal_contact.measure import certify_sign_controls, route_certificate

BASE = Fraction(1, 2)
DRIVES = [(Fraction(5, 8), Fraction(2, 5)), (Fraction(9, 16), Fraction(4, 9))]
SETTINGS = (BASE, *DRIVES[0], *DRIVES[1])
PATHS, CALIBRATIONS, LAGS, WARMUP = 2**37, 2**26, 64, 1024
HERE = Path(__file__).resolve().parent


@cache
def filter_weights(setting):
    h, decoder, _ = gates("native")
    # Conjugate the AB contact to BC by relabelling the three material ports.
    permutation = [2, 0, 1]
    duration = 2 * acos(sqrt(setting))
    leak = expm(-1j * duration * h[np.ix_(permutation, permutation)])
    raw = sum(outcomes(decoder, leak, rail, 1) for rail in (0, 1)) / 2
    measured = raw[[1, 2, 0]]  # lower, upper, refusal
    exact = np.array([0.5, float(setting / 2), float((1 - setting) / 2)])
    np.testing.assert_allclose(measured, exact, rtol=0, atol=3e-13)
    return exact


def initial(rng, size):
    lower, upper, refusal = filter_weights(BASE)
    census = np.array([PATHS], dtype=np.int64)
    failures = []
    for _ in range(size):
        # Failed probes before the PATHS-th acceptance are independent of
        # accepted labels; this is the joint law of retry-until-accepted.
        failures.append(int(rng.negative_binomial(PATHS, 1 - refusal)))
        high = rng.binomial(census, upper / (lower + upper))
        census = np.r_[census - high, high]
    return census, failures


class Calendar:
    def __init__(self, size, left, right, cut):
        words = np.arange(2**size)
        menu = 2 ** max(2, size.bit_length())
        self.targets = [words]
        self.names = ["idle"]
        weights = [(menu - size - 1) / menu]
        for side, site, setting in (("left", 0, left), ("right", size - 1, right)):
            for sign, weight in enumerate(filter_weights(setting)):
                self.names.append(f"{side}_{('lower', 'upper', 'refused')[sign]}")
                if sign == 2 or (side == "right" and cut):
                    self.targets.append(words)
                else:
                    self.targets.append((words & ~(1 << site)) | (sign << site))
                weights.append(weight / menu)
        for site in range(size - 1):
            different = ((words >> site) ^ (words >> (site + 1))) & 1
            self.targets.append(words ^ ((3 << site) * different))
            self.names.append(f"bond_{site}")
            weights.append(1 / menu)
        self.weights = np.array(weights)
        self.targets = np.array(self.targets).T
        self.left_current = np.zeros_like(self.targets)
        for event in (1, 2):
            self.left_current[:, event] = (self.targets[:, event] & 1) - (words & 1)
        self.nwords = len(words)

    def counts(self, rng, census):
        return rng.multinomial(census, self.weights)

    def propagate(self, counts):
        n = np.zeros(self.nwords, dtype=np.int64)
        np.add.at(n, self.targets, counts)
        return n

    def currents(self, counts):
        return np.bincount(
            self.left_current.ravel() + 1, weights=counts.ravel(), minlength=3
        ).astype(np.int64)


def equilibrium(rng, size, cut):
    cal = Calendar(size, BASE, BASE, cut)
    start, refusals = initial(rng, size)
    counts = cal.counts(rng, start)
    origins = [cal.propagate(counts * (cal.left_current == j)) for j in (-1, 0, 1)]
    offers = counts.sum(axis=0)
    joints = [np.diag(cal.currents(counts)).tolist()]
    for _ in range(LAGS):
        triplet = [cal.counts(rng, n) for n in origins]
        joints.append([cal.currents(c).tolist() for c in triplet])
        origins = [cal.propagate(c) for c in triplet]
        offers += sum(c.sum(axis=0) for c in triplet)
    return dict(
        size=size,
        cut=cut,
        initial_refusals=refusals,
        joints=joints,
        all_offers=dict(zip(cal.names, offers.tolist(), strict=True)),
    )


def driven(rng, size, left, right, cut):
    cal = Calendar(size, left, right, cut)
    pop, refusals = initial(rng, size)
    offers = np.zeros(len(cal.names), dtype=np.int64)
    for _ in range(WARMUP):
        counts = cal.counts(rng, pop)
        offers += counts.sum(axis=0)
        pop = cal.propagate(counts)
    counts = cal.counts(rng, pop)
    offers += counts.sum(axis=0)
    return dict(
        size=size,
        left=str(left),
        right=str(right),
        cut=cut,
        initial_refusals=refusals,
        counts_before_final=pop.tolist(),
        final_offers_by_word=counts.tolist(),
        all_offers=dict(zip(cal.names, offers.tolist(), strict=True)),
        current=cal.currents(counts).tolist(),
    )


def run():
    rng = np.random.default_rng(749103)
    pointer()
    certify_sign_controls()
    routes = route_certificate()
    if len(routes) != 15:
        raise ValueError("Missing route between deployed registers")
    calibration = {
        str(w): rng.multinomial(CALIBRATIONS, filter_weights(w)).tolist()
        for w in SETTINGS
    }
    eq, drive = [], []
    for size in (1, 2, 4):
        for cut in (False, True):
            eq.append(equilibrium(rng, size, cut))
            for hot, cold in DRIVES[:1] if cut else DRIVES:
                for left, right in ((hot, cold), (cold, hot)):
                    drive.append(driven(rng, size, left, right, cut))
    return dict(
        paths=PATHS,
        calibration_proposals=CALIBRATIONS,
        lags=LAGS,
        burn_in=WARMUP,
        thermometer=calibration,
        equilibrium=eq,
        driven=drive,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), separators=(",", ":")) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(output)
    elif path.read_text() != output:
        raise SystemExit("Heat-current receipt mismatch")
