"""A count-level execution of each local exchange, preserving original launch totals."""

import argparse
from itertools import combinations
import json
from pathlib import Path
import numpy as np
from energy_release import HERE, certify_gates, filter_labels


def count_run(rng, n, depths, cap, delays, parked):
    failures = []
    starts = n
    for depth in depths:
        bins = rng.multinomial(starts, list(map(float, filter_labels(depth))))
        failures.append(bins[:-2].tolist())
        starts = int(bins[-2:].sum())
    initial = np.arange(64)
    upper = np.array([1 / (1 + 2**d) for d in depths])
    initial_bits = ((initial[:, None] >> np.arange(6)) & 1).astype(bool)
    pword = np.prod(np.where(initial_bits, upper, 1 - upper), axis=1)
    # Filters act independently. Given the total admitted, the six labels
    # have the product conditional law irrespective of the earlier refusals.
    awaiting = rng.multinomial(starts, pword)
    terminal = np.zeros_like(awaiting)
    by_age = []
    interchanges = []
    for i in range(3):
        unlike = ((initial >> (2 * i)) ^ (initial >> (2 * i + 1))) & 1
        interchanges.append(initial ^ (unlike * (3 << (2 * i))))
    for t in range(1, cap + 1):
        # The root read is independent of each initially sharp energy word.
        group = rng.binomial(awaiting, 0.25)
        awaiting -= group
        for step in range(1, t + max(delays) + 1):
            for i, perm in enumerate(interchanges):
                if step <= t + delays[i] and not parked:
                    offered = rng.binomial(group, 0.25)
                    group = group - offered + offered[perm]
        by_age.append(group.tolist())
        terminal += group
    hits = int(terminal.sum())
    assert hits + int(awaiting.sum()) + sum(map(sum, failures)) == n
    singles = [float(terminal[initial_bits[:, 2 * i]].sum() / hits) for i in range(3)]
    joint = [
        float(terminal[initial_bits[:, 2 * i] & initial_bits[:, 2 * j]].sum() / hits)
        for i, j in combinations(range(3), 2)
    ]
    return dict(
        first_filter_leaks=failures,
        timeout=int(awaiting.sum()),
        delivered=hits,
        age_by_terminal_word=by_age,
        probe_upper=singles,
        pair_upper=joint,
        covariance=[
            v - singles[i] * singles[j]
            for v, (i, j) in zip(joint, combinations(range(3), 2))
        ],
    )


def run(out):
    protocol = json.loads((HERE / "precision_protocol.json").read_text())
    source = json.loads((HERE / "protocol.json").read_text())
    certify_gates()
    rng = np.random.default_rng(protocol["seed"])
    banks = {}
    for bank in protocol["banks"]:
        cold = source["bodies"][
            "warm_filter_depth" if bank == "equal_temperature" else "cold_filter_depth"
        ]
        depths = [cold, source["bodies"]["warm_filter_depth"]] * 3
        data = count_run(
            rng,
            protocol["launches_per_bank"],
            depths,
            source["trigger"]["max_reads"],
            source["trigger"]["capture_delay"],
            bank == "parked",
        )
        b = np.sqrt(
            np.log(
                2
                * protocol["checked_moments_per_bank"]
                * len(protocol["banks"])
                / protocol["family_error"]
            )
            / (2 * data["delivered"])
        )
        data["simultaneous_covariance_halfwidth"] = [
            b * (1 + data["probe_upper"][i] + data["probe_upper"][j]) + b * b
            for i, j in combinations(range(3), 2)
        ]
        banks[bank] = data
    out.mkdir(parents=True, exist_ok=True)
    (out / "precision.json").write_text(json.dumps(banks, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE)
    run(p.parse_args().out)
