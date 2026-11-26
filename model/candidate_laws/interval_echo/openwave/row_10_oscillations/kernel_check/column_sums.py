"""Reconstruct the three-label decisions from raw counts and test population balance."""

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

from scipy.stats import beta

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "joint_receivers" / "blind_counts.json"
OWNERS = (12, 17)
WAITS = (256, 512, 2048)
PANELS = (-1, 1)
LABELS = range(3)
TAIL = 0.01 / (2 * 4096)


def read_banks(record):
    expected = set(product(OWNERS, WAITS, LABELS, PANELS, ("training", "validation")))
    banks = {}
    for row in record["banks"]:
        key, counts = tuple(row["key"]), row["counts"]
        if key in banks or key not in expected or len(counts) != 7:
            raise ValueError("Invalid or duplicate original-start bank")
        if (
            any(type(x) is not int or x < 0 for x in counts)
            or sum(counts) != record["starts"]
        ):
            raise ValueError("Inconsistent original-start counts")
        if sum(counts[:3]) == 0:
            raise ValueError("Conditional received-pair probability unavailable")
        banks[key] = counts[:3]
    if set(banks) != expected:
        raise ValueError("Missing original-start bank")
    return banks


def decision(training):
    totals = [sum(r) for r in training]
    # Python's max keeps the first label on an exact tie.
    return [
        max(LABELS, key=lambda a: Fraction(training[a][k], totals[a])) for k in LABELS
    ]


def confidence(successes, total, selected):
    if selected == 0:
        return 0.0, 0.0
    if selected == 3:
        return 1.0, 1.0
    low = beta.ppf(TAIL, successes, total - successes + 1) if successes else 0.0
    high = (
        beta.isf(TAIL, successes + 1, total - successes) if successes < total else 1.0
    )
    return float(low), float(high)


def evaluate(banks, owner, wait):
    # Average the two reference settings with equal weight within each source.
    columns = [{"sum": 0.0, "low": 0.0, "high": 0.0} for _ in LABELS]
    decisions = []
    for panel in PANELS:
        trained = decision([banks[owner, wait, a, panel, "training"] for a in LABELS])
        decisions.append(trained)
        for a in LABELS:
            received = banks[owner, wait, a, panel, "validation"]
            for label in LABELS:
                choices = [k for k in LABELS if trained[k] == label]
                count = sum(received[k] for k in choices)
                total = sum(received)
                low, high = confidence(count, total, len(choices))
                col = columns[label]
                col["sum"] += count / (2 * total)
                col["low"] += low / 2
                col["high"] += high / 2
    # A uniform prior remains uniform under every doubly stochastic kernel.
    # Each band's minimum deviation survives the common confidence event.
    distance = sum(max(0.0, c["low"] - 1, 1 - c["high"]) for c in columns) / 6
    return {
        "owner": owner,
        "wait": wait,
        "decisions": decisions,
        "columns": columns,
        "uniform_output_TV_lower": distance,
    }


def analyze(path=SOURCE):
    raw = path.read_bytes()
    record = json.loads(raw)
    banks = read_banks(record)
    return {
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "original_starts": record["starts"],
        "conditional_single_bin_allocations": 4096,
        "type": "retrospective",
        "results": [
            evaluate(banks, owner, wait) for owner, wait in product(OWNERS, WAITS)
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    destination = parser.parse_args().output
    summary = analyze()
    if destination:
        destination.write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    for result in summary["results"]:
        print(
            result["owner"],
            result["wait"],
            [
                (round(c["sum"], 5), round(c["low"], 5), round(c["high"], 5))
                for c in result["columns"]
            ],
            round(result["uniform_output_TV_lower"], 5),
        )
