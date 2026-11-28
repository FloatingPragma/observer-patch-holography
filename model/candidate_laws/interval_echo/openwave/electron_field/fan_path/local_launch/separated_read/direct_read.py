"""Uninterrupted signed transfer, with the guide and charge observed apart."""

import argparse
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
from runpy import run_path

import numpy as np

from .measure import (
    HERE,
    ORIGIN,
    PAIR,
    PAULI,
    STEP,
    TARGET,
    TICKS,
    TRIALS,
    WAITS,
    FanPath,
    archived_fixture,
    counter,
    evolution,
    joint_read,
    source_columns,
)

VECTOR = np.array([8, 8 - 5j, 6, -8 + 3j])
TIMES = tuple(sorted({t + w for t in TICKS for w in WAITS}))
COUNTED = (437, 497)
FAMILY = 2 * len(COUNTED)


def run():
    coefficient = np.array(
        [[np.vdot(VECTOR, np.kron(a, b.T) @ VECTOR).real for b in PAULI] for a in PAULI]
    )
    denom = float(4 * np.vdot(VECTOR, VECTOR).real)
    assert set(COUNTED).issubset(TIMES)
    np.testing.assert_array_equal(coefficient, np.rint(coefficient))
    rng = np.random.default_rng(101010451498)
    cases = []
    with archived_fixture() as archive:
        for sign in (1, -1):
            path = FanPath(archive, 4, sign, asymmetric=True)
            mask = sum(path.field.bits.values())
            pair = PAIR if sign == 1 else tuple(mask ^ word for word in PAIR)
            block = path.blocks[2, path.exit[2]]
            assert path.cells[2][block.vertex] == TARGET and path.cells[3][0] == ORIGIN
            readings = []
            for tick, (wave, err) in evolution(
                path.h, source_columns(path), TIMES, sqrt(2), scipy_ticks=COUNTED
            ).items():
                report = joint_read(wave, err, 1.0, path, block, pair)
                records = None
                if tick in COUNTED:
                    records = counter(rng, report, coefficient, denom, family=FAMILY)
                readings.append(
                    {
                        "appointment": float(tick * STEP),
                        "terminal_record": report,
                        "simulated_followup_witness": records,
                    }
                )
            cases.append(
                {
                    "signed_material": sign,
                    "ordered_1045_face_words": list(pair),
                    "uninterrupted_reads": readings,
                }
            )
    return {
        "pre_count_protocol_sha256": sha256(
            (HERE / "DIRECT.md").read_bytes()
        ).hexdigest(),
        "nominal_counted_times": [float(t * STEP) for t in COUNTED],
        "simulated_trials_per_bank": TRIALS,
        "simultaneous_cohort_family_size": FAMILY,
        "witness_numerator_pauli_coefficients": coefficient.astype(int).tolist(),
        "witness_denominator": denom,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parents[2] / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "direct_reads.json").read_text()), "direct_reads"
        )
        print(
            "Replayed uninterrupted signed read and its four simulated witness cohorts."
        )


if __name__ == "__main__":
    main()
