"""Infer a downstream-reader bound from counts, then test two held phase settings."""

import argparse
import hashlib
import json
from fractions import Fraction

import numpy as np
from scipy.stats import beta

from . import measure as pair
from ..sequential import compare, grid_counts, require

HERE = pair.HERE
RECEIPT = HERE / "observer_bound.json"
ATTEMPTS = 8192
EXPERIMENTS = ("minus", "quarter", "three_quarters")
GRID_ALLOWANCE = 2 * 3e-12  # A return event is a union of two raw bins.


def lower_probability(hits, attempts, alpha):
    """One-sided exact binomial bound including a zero-success cohort."""
    return 0.0 if hits == 0 else float(beta.ppf(alpha, hits, attempts - hits + 1))


def return_count(bins, start, address, shots):
    """The receiver decodes the second native pointer relative to its first word."""
    require(len({(a, b) for a, b, _ in bins}) == len(bins), "duplicate raw bin")
    require(all(0 <= a < 17 and 0 <= b < 17 and type(n) is int and n > 0
                for a, b, n in bins) and sum(n for _, _, n in bins) == shots,
            "incomplete raw cohort")
    return sum(n for a, b, n in bins if a in (start, start + 1) and b == a ^ address)


def native_joint(head, terminal, start, experiment):
    back = pair.probe_words(start, "port0")[1]
    move = pair.fixed_pointer.inverse(back)
    input_name = "plus" if experiment == "minus" else "plus_i"
    prep = pair.apparatus.preparation_word(input_name)
    angle = {"minus": Fraction(1, 8), "quarter": Fraction(1, 96),
             "three_quarters": Fraction(1, 32)}[experiment]
    phase = [pair.pulses.pulse("reference", 0, 0, angle)]
    if experiment == "minus":
        source = prep + phase + move
        read = back + pair.apparatus.x_analyzer_word()
    else:
        source = prep + move
        read = back + phase + pair.fixed_pointer.inverse(prep)
    vec = pair.pulses.source_unitary(tuple(source), 1)[:, 0]
    analyzer = pair.pulses.source_unitary(tuple(read), 1)
    joint = np.zeros((17, 17))
    for raw, branch in enumerate(head):
        if np.any(branch):
            paths = terminal[raw] @ (analyzer @ branch @ vec)
            joint[raw] = np.sum(abs(paths)**2, axis=1)
    require(abs(joint.sum() - 1) < 3e-12, "missing instrument mass")
    return joint / joint.sum()


def infer(banks, attempts):
    """Only measured public counts enter the bounds."""
    observations = {}
    for r in banks:
        start, name, mode = r["pair_start"], r["probe"], r["instrument"]
        target = int(name in ("port1", "minus"))
        hits = return_count(r["nonzero_joint_counts"], start, target, attempts)
        require((start, mode, name) not in observations, "duplicate probe")
        observations[start, mode, name] = hits
    result = []
    for start in range(0, pair.DIM, 2):
        lower = [
            max(0, lower_probability(observations[start, "pair", name], attempts, 0.005 / 24)
                - GRID_ALLOWANCE)
            for name in ("port0", "port1", "plus", "minus")
        ]
        fidelity = max(0.0, sum(lower) / 2 - 1)
        misses = attempts - observations[start, "address_then_forget", "minus"]
        rival = max(0.0, lower_probability(misses, attempts, 0.0025 / 6) - GRID_ALLOWANCE)
        result.append({"pair_start": start, "entanglement_fidelity_lower": fidelity,
                       "any_event_error_upper": min(1.0, 2 * np.sqrt(1 - fidelity)),
                       "rival_any_event_error_lower": rival})
    return result


def report():
    earlier = json.loads(pair.RECEIPT.read_bytes())
    require(earlier["shots_per_setting"] == ATTEMPTS, "unexpected calibration size")
    calls = pair.address.schedule()
    heads = {"pair": pair.read_branches([c for c in calls if c["bit"] != 0], 0),
             "address_then_forget": pair.read_branches(calls, 0)}
    terminal = {r: pair.read_branches(calls, r) for r in range(pair.DIM)}
    rng = np.random.default_rng(260611)
    banks = []
    for start in range(0, pair.DIM, 2):
        for mode, head in heads.items():
            for name in EXPERIMENTS:
                joint = native_joint(head, terminal, start, name)
                table = grid_counts(rng, ATTEMPTS, joint.ravel()).reshape(17, 17)
                bins = [[int(a), int(b), int(table[a, b])] for a, b in zip(*np.nonzero(table))]
                banks.append(dict(pair_start=start, instrument=mode, probe=name,
                                  nonzero_joint_counts=bins))
    minus_records = [r for r in banks if r["probe"] == "minus"]
    archived = infer(earlier["banks"] + minus_records, ATTEMPTS)
    fresh_rng = np.random.default_rng(260612)
    fresh = []
    for start in range(0, pair.DIM, 2):
        for name in ("port0", "port1", "plus", "minus"):
            if name == "minus":
                joint = native_joint(heads["pair"], terminal, start, name)
            else:
                joint, _ = pair.read_probe(heads["pair"], terminal, start, name)
                joint = joint / joint.sum()
            table = grid_counts(fresh_rng, ATTEMPTS, joint.ravel()).reshape(17, 17)
            bins = [[int(a), int(b), int(table[a, b])] for a, b in zip(*np.nonzero(table))]
            fresh.append(dict(pair_start=start, instrument="pair", probe=name,
                              nonzero_joint_counts=bins))
    rival_minus = [r for r in minus_records if r["instrument"] == "address_then_forget"]
    bounds = infer(fresh + rival_minus, ATTEMPTS)
    radius = float(np.sqrt(np.log(2 * 24 / 0.0025) / (2 * ATTEMPTS))) + GRID_ALLOWANCE
    held = []
    for r in banks:
        if r["probe"] == "minus":
            continue
        start = r["pair_start"]
        ideal = (2 + (1 if r["probe"] == "quarter" else -1) * np.sqrt(2)) / 4
        rate = return_count(r["nonzero_joint_counts"], start, 0, ATTEMPTS) / ATTEMPTS
        limit = bounds[start // 2]["any_event_error_upper"] + radius
        held.append(dict(pair_start=start, instrument=r["instrument"], setting=r["probe"],
                         rate=rate, forecast=ideal, within_coarse_bound=bool(abs(rate - ideal) <= limit)))
        center = ideal if r["instrument"] == "pair" else 0.5
        require(abs(rate - center) <= radius, "held source forecast failed")
    require(all(r["within_coarse_bound"] == (r["instrument"] == "pair") for r in held),
            "held comparison lacks separation")
    sources = [pair.RECEIPT, HERE / "BOUND_PROTOCOL.md", HERE / "measure.py",
               HERE / "observer_bound.py"]
    return {"sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            "seed_pcg64": 260611, "attempts_per_bank": ATTEMPTS, "new_launches": (len(banks) + len(fresh)) * ATTEMPTS,
            "fresh_seed_pcg64": 260612, "fresh_calibration_banks": fresh, "archival_bounds": archived,
            "grid_bits": pair.GRID_BITS, "complete_raw_table_shape": [17, 17],
            "omitted_bins_have_zero_count": True, "combined_confidence": 0.99,
            "held_sampling_allowance": radius, "bounds": bounds, "new_banks": banks, "held": held}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    obs = report()
    if args.write:
        with RECEIPT.open("x") as stream:
            stream.write(json.dumps(obs, indent=2, allow_nan=False) + "\n")
    compare(json.loads(RECEIPT.read_bytes()), obs)
    print("Verified six two-basis certificates and 24 held phase cohorts.")
