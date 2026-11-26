"""Cool any half-filled input using disjoint one-bond energy intervals."""

import argparse
import json
from math import comb, exp, fsum, log, pi, sqrt
from pathlib import Path
from random import Random

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30.experiment import select_route
from model.executable_model.carrier import federation
from ..measure import ring_modes

HERE = Path(__file__).resolve().parent
SHOTS = 2**30
ALPHA = 0.01
SEEDS = {"inverted": 2901, "mixed": 2902, "alternating": 2903, "reversed": 2904}


def calibration_counts(energies):
    rng = np.random.default_rng(271828)
    n = len(energies)
    return [rng.multinomial(SHOTS, [(1 - e) / n, (1 + e) / n, 1 - 2 / n]).tolist()
            for e in energies]


def read_calibration(counts):
    """Split eigensign uncertainty from the rate of seeing the test load."""
    intervals = []
    for symmetric, antisymmetric, quiet in counts:
        if min(symmetric, antisymmetric, quiet) < 0 or symmetric + antisymmetric + quiet != SHOTS:
            raise ValueError("Every calibration launch needs exactly one terminal bin")
        hits = symmetric + antisymmetric
        if hits == 0:
            raise ValueError("The orbital has no registered bond hits")
        energy = (antisymmetric - symmetric) / hits
        error = sqrt(2 * log(2 * len(counts) / ALPHA) / hits)
        intervals.append((max(-1, energy - error), min(1, energy + error)))
    return intervals


def work_choices(intervals):
    return [[upper < lower for _, upper in intervals] for lower, _ in intervals]


def initial_word(name, n, lower, rng):
    if name == "mixed":
        return set(rng.sample(range(n), n // 2))
    if name == "alternating":
        # A sharp alternating site word has one fermion in each (k,k+π) doublet.
        return {j + (n // 2) * rng.randrange(2) for j in range(n // 2)}
    if name == "inverted":
        return set(range(n)) - lower
    if name == "reversed":
        return lower.copy()
    raise ValueError("Unknown source")


def trajectory(name, energies, choices, lower, intervals):
    rng = Random(SEEDS[name])
    n = len(energies)
    occupied = initial_word(name, n, lower, rng)
    target_group = set(range(n)) - lower if name == "reversed" else lower
    excitations = len(occupied - target_group)
    initial = sorted(occupied)
    checkpoints = []
    accepted = 0
    energy_changes = []
    budget = 16 * comb(n, 2)
    step = 0
    while excitations and step < budget:
        step += 1
        src = rng.randrange(n)
        dst = rng.randrange(n - 1)
        if dst >= src:
            dst += 1
        if dst in occupied:
            src, dst = dst, src
        if src not in occupied or dst in occupied:
            continue
        allowed = choices[dst][src] if name == "reversed" else choices[src][dst]
        if not allowed:
            continue
        occupied.remove(src)
        occupied.add(dst)
        accepted += 1
        energy_changes.append(float(energies[src] - energies[dst]))
        difference = int(dst in target_group) - int(src in target_group)
        excitations -= difference
        if difference and excitations in (64, 32, 16, 8, 4, 2, 1, 0):
            checkpoints.append([step, excitations])
    before = fsum(energies[i] for i in initial)
    after = fsum(energies[i] for i in occupied)
    if abs(before - after - fsum(energy_changes)) > 1e-11:
        raise ArithmeticError("The measured work changes fail to telescope")
    emptied = set(initial) - occupied
    gained = occupied - set(initial)
    read_work = [fsum(intervals[i][s] for i in emptied)
                 - fsum(intervals[i][1 - s] for i in gained) for s in (0, 1)]
    return {"initial_mode_word": initial, "initial_upper_count": len(set(initial) - lower),
            "final_upper_count": len(occupied - lower), "accepted_transfers": accepted,
            "pair_proposals": step, "budget": budget, "unfinished": bool(excitations),
            "initial_energy": before, "final_energy": after,
            "released_hopping_work": fsum(energy_changes), "work_from_count_intervals": read_work,
            "checkpoints": checkpoints}


def build(counts=None):
    car = federation(2)
    route, _ = select_route(car)
    if len(set(route)) != len(route) or any(b not in car.neighbours[a]
            for a, b in zip(route, route[1:] + route[:1], strict=True)):
        raise ValueError("The hopping action needs a cycle of actual neighbor contacts")
    n = len(route)
    momenta, _, filled = ring_modes(n)
    energies = -np.cos(momenta)
    counts = calibration_counts(energies) if counts is None else counts
    intervals = read_calibration(counts)
    if len(intervals) != n or any(not lo <= energy <= hi for (lo, hi), energy
                                  in zip(intervals, energies, strict=True)):
        raise ValueError("An independent calibration interval misses the hopping energy")
    lower = {j for j, (lo, hi) in enumerate(intervals) if hi < 0}
    upper = {j for j, (lo, hi) in enumerate(intervals) if lo > 0}
    if len(lower) != n // 2 or len(upper) != n // 2 or lower != set(np.flatnonzero(filled)):
        raise ValueError("The recorded work intervals fail to separate the half-filled shells")
    choices = work_choices(intervals)
    if not all(choices[j][i] for j in upper for i in lower):
        raise ValueError("A downhill transfer across the gap is unresolved")
    pairs = comb(n, 2)
    return {"stations": n, "load": n // 2, "independent_pair_count": pairs,
            "shots_per_calibration": SHOTS, "counts_symmetric_antisymmetric_quiet": counts,
            "joint_confidence": 1 - ALPHA, "expected_proposals_worst":
                pairs * fsum(1 / j**2 for j in range(1, 1 + n // 2)),
            "finite_failure_bound": (n // 2) * exp(-16),
            "resolved_gap_lower": min(intervals[j][0] for j in upper)
                                  - max(intervals[i][1] for i in lower),
            "exact_gap": 2 * np.sin(pi / n),
            "parked_inverted_upper_count": n // 2,
            "trajectories": {name: trajectory(name, energies, choices, lower, intervals)
                             for name in SEEDS}}


def audit(actual, saved):
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError("The fields of a selection record changed")
        for key in actual:
            audit(actual[key], saved[key])
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise ValueError("The length of a selection record changed")
        for new, old in zip(actual, saved, strict=True):
            audit(new, old)
    elif isinstance(actual, float):
        if not np.isclose(actual, saved, rtol=1e-12, atol=1e-13):
            raise ValueError(f"Selection readout differs: {actual} versus {saved}")
    elif actual != saved:
        raise ValueError(f"Selection record differs: {actual} versus {saved}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(build(), indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        audit(build(saved["counts_symmetric_antisymmetric_quiet"]), saved)
        print("Verified the common work calibration and the full selection trajectories.")


if __name__ == "__main__":
    main()
