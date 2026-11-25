"""Replay the material source followed by eight common-context instruments."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from ..phase_table import measure as analyzer
from ..sequential import compare, require
from . import source

RECEIPT = source.HERE / "results.json"
LAUNCHES = 262_144
SEED = 6190819


def cohort(branches, state, rng):
    """Read yes/no and the tag orientation at the end of each launch."""
    records = [op @ state for op in branches]
    weights = np.asarray([np.sum(abs(out) ** 2, axis=0) for out in records])
    require(abs(float(weights.sum()) - 1) < 2e-11, "instrument lost a branch")
    counts = rng.multinomial(LAUNCHES, (weights / weights.sum()).ravel()).reshape(17, 2)
    return {
        "joint_probability": np.round(weights[[0, 2]], 11).tolist(),
        "joint_counts": counts[[0, 2]].tolist(),
        "all_other_tag_counts": counts[
            [i for i in range(17) if i not in (0, 2)]
        ].tolist(),
        "all_other_probability": round(
            float(weights.sum() - weights[[0, 2]].sum()), 11
        ),
        "yes_probability": round(float(weights[0].sum()), 11),
        "yes_count": int(counts[0].sum()),
    }


def collect():
    rng = np.random.default_rng(SEED)
    inputs = {}
    for label, tag in (("tagged", True), ("copy_omitted", False)):
        inputs[label], duration = source.prepare(tagged=tag)
    # Physical analyzer routes may reuse restored transit loads 72 and 82.
    # Neither the preparation tag nor its reference may be read by a context.
    for mode in ("sharp", "weak"):
        endpoints = {p for _, a, b, _ in analyzer.record_word(mode) for p in (a, b)}
        require(
            endpoints.isdisjoint({76, 77, 79, 239}),
            "analyzer addresses an input spectator",
        )
    outputs = []
    for entry in json.loads(source.TABLE.read_text())["contexts"]:
        target = source.archived_effect(entry)
        weight = float(Fraction(entry["born_weight"]))
        cases = {}
        for mode in ("sharp", "weak"):
            effects = target if mode == "sharp" else np.eye(2) / 4 + target / 2
            cases[mode] = {}
            for j, sign in enumerate(("vacancy", "excess")):
                raw = analyzer.raw_branch_maps(entry["context"], mode, sign=(-1, 1)[j])
                ops = [raw[i] for i in range(17)]
                require(
                    np.linalg.norm(ops[0].conj().T @ ops[0] - effects, 2) < 2e-11,
                    "physical yes effect differs from the named operand",
                )
                require(
                    np.linalg.norm(sum(k.conj().T @ k for k in ops) - np.eye(2), 2)
                    < 2e-11,
                    "incomplete physical branch instrument",
                )
                results = {
                    label: cohort(ops, states[j], rng)
                    for label, states in inputs.items()
                }
                wanted = weight if mode == "sharp" else 0.25 + weight / 2
                require(
                    abs(results["tagged"]["yes_probability"] - wanted) < 2e-10,
                    "common physical preparation failed a table context",
                )
                cases[mode][sign] = results
        outputs.append(
            {
                "context": entry["context"],
                "archived_weight": entry["born_weight"],
                "instruments": cases,
            }
        )
    bound = np.sqrt(np.log(2 * 64 / 0.01) / (2 * LAUNCHES))
    deviations = [
        abs(record["yes_count"] / LAUNCHES - record["yes_probability"])
        for context in outputs
        for mode in context["instruments"].values()
        for sign in mode.values()
        for record in sign.values()
    ]
    return {
        "candidate": "interval_echo",
        "protocol": "Independent multinomial simulation after full source and physical branch propagation",
        "calibration": {
            "source_block_counts": [111, 68],
            "exposure": round(duration, 12),
            "data_ports": [0, 1],
            "tag_ports": [76, 77],
            "host": 79,
            "outbound_swap_path": list(source.ROUTE),
            "charge_reference": 239,
        },
        "terminal_bins": {
            "row": ["yes (raw 0)", "no (raw 2)"],
            "column": ["tag_20", "tag_02"],
            "all_other_raw_rows": [i for i in range(17) if i not in (0, 2)],
        },
        "cohorts": len(deviations),
        "launches_per_cohort": LAUNCHES,
        "seed": SEED,
        "simultaneous_99pct_yes_frequency_halfwidth": round(float(bound), 10),
        "largest_sample_minus_model_absolute": round(max(deviations), 10),
        "inputs_sha256": {
            str(p.relative_to(source.ROOT)): sha256(p.read_bytes()).hexdigest()
            for p in (
                source.ORIGINAL,
                source.TABLE,
                Path(source.__file__),
                Path(source.address.__file__),
                Path(analyzer.__file__),
                Path(__file__),
            )
        },
        "contexts": outputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new receipt at this path")
    args = parser.parse_args()
    result = collect()
    if args.output:
        with args.output.open("x") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare(json.loads(RECEIPT.read_text()), result)
        print(
            "Verified eight native contexts on tagged and coherent inputs across both charge signs."
        )


if __name__ == "__main__":
    main()
