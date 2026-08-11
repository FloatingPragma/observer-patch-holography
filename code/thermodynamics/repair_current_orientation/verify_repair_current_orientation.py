#!/usr/bin/env python3
"""Hermetic independent verifier for the post-hoc repair-current payload.

The verifier intentionally imports neither the simulator producer nor any OPH
library.  It authenticates the exact vendored schema-v3 bytes, pins the source
metadata and count table, and recomputes every designation with standard-library
integer and ``Fraction`` arithmetic.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import permutations
import json
from pathlib import Path
import re
from typing import Any, Sequence


HERE = Path(__file__).resolve().parent
DEFAULT_PAYLOAD = HERE / "repair_current_payload.v3.json"
DEFAULT_LEAN_SOURCE = HERE.parents[2] / "Lean" / "Thermodynamics" / "RepairCurrentOrientation.lean"
EXPECTED_PAYLOAD_SHA256 = (
    "7f8ea7ef9c92a50e23207c2fe85d09ed2bce1c1aa539ae9914a9b9edd0df26d6"
)
EXPECTED_SCHEMA = "oph.sim.repair_current_payload.v3"
EXPECTED_OBSERVER_COUNT = 1028
EXPECTED_TRANSITION_COUNT = 31744
EXPECTED_BUCKETS = tuple(range(8))
EXPECTED_ROW_SUMS = (23552, 43, 172, 1766, 1833, 2332, 1486, 560)

EXPECTED_COUNTS = (
    (23552, 0, 0, 0, 0, 0, 0, 0),
    (0, 2, 32, 0, 1, 5, 0, 3),
    (0, 0, 0, 167, 0, 5, 0, 0),
    (0, 28, 0, 349, 1343, 0, 46, 0),
    (0, 0, 135, 0, 0, 1554, 0, 144),
    (0, 1, 0, 594, 261, 227, 1246, 3),
    (1024, 10, 5, 0, 205, 93, 78, 71),
    (0, 2, 0, 283, 23, 81, 116, 55),
)

EXPECTED_PINNED_INPUT_SHA256 = {
    "conditional_resampling_realization_receipt.json": (
        "d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749"
    ),
    "finite_repair_transition_matrix.npz": (
        "1f8ff2cc413f4bbfeabe9450ae48cb5d0794ec28f284f25913e7c9c61b677163"
    ),
    "finite_repair_transition_matrix_report.json": (
        "e38bb28475fb7b34127b864abd20414efc920d8ef9ddbdf5c78612874ae37f2d"
    ),
    "git_commit.txt": (
        "2a084f768662274335ba7c23afaa18b61bd1184e8d481a4e17ab511c1e9d3293"
    ),
    "observer_views.jsonl": (
        "3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca"
    ),
}

EXPECTED_PROVENANCE = {
    "counting_convention": (
        "ordered pairs of consecutive transition_history_descriptor steps per "
        "observer, unweighted integer multiplicities, under the report state "
        "alphabet; identical to the mixing-chain extraction convention"
    ),
    "designated_cycle_rule": (
        "lexicographically least ordered bucket 3-cycle maximizing "
        "|forward product - backward product|"
    ),
    "designated_pair_rule": (
        "lexicographically least ordered bucket pair maximizing |C(a,b)-C(b,a)|"
    ),
    "falsifiable_checks": (
        "observer, skip, and transition totals equal the pinned report; every "
        "step lies inside the report alphabet with exact integer field values; "
        "the recounted 26-state integer matrix matches the pinned npz weighted "
        "matrix on positivity support and within 1e-6 on weighted recount"
    ),
    "pinned_input_sha256": EXPECTED_PINNED_INPUT_SHA256,
    "projection_field": "repair_load_bucket",
    "recount_identities": (
        "reversing every transition window transposes the count table, and "
        "re-bucketing conserves totals; both hold for any ordered recount by "
        "construction and are identities, not data checks"
    ),
    "run_git_commit": "b39b78faf894894ebe573571e0902ccfaaeac32a",
    "run_id": "b12_prereg_16k_20260806",
}

EXPECTED_TOP_LEVEL_KEYS = {
    "analysis_status",
    "buckets",
    "current_antisymmetric_part",
    "cycle_orientation_nonempty",
    "designated_cycle",
    "designated_pair",
    "designation_rule_preregistered",
    "eligible_as_validation",
    "observer_count",
    "ordered_counts",
    "pair_orientation_nonempty",
    "provenance",
    "schema",
    "statistic_preregistered",
    "transition_count",
}


class VerificationError(ValueError):
    """Raised when the vendored packet fails closed."""


def need(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


@dataclass(frozen=True)
class PairDesignation:
    a: int
    b: int
    forward: Fraction
    backward: Fraction

    @property
    def difference(self) -> Fraction:
        return self.forward - self.backward


@dataclass(frozen=True)
class CycleDesignation:
    a: int
    b: int
    c: int
    forward: Fraction
    backward: Fraction

    @property
    def difference(self) -> Fraction:
        return self.forward - self.backward


EXPECTED_RAW_PAIR = PairDesignation(3, 4, Fraction(1343), Fraction(0))
EXPECTED_NORMALIZED_PAIR = PairDesignation(
    2, 3, Fraction(167, 172), Fraction(0)
)
EXPECTED_RAW_CYCLE = CycleDesignation(
    3, 4, 5, Fraction(1239691068), Fraction(0)
)
EXPECTED_NORMALIZED_CYCLE = CycleDesignation(
    3, 4, 5, Fraction(9391599, 57188378), Fraction(0)
)


def parse_integer_table(
    value: Any, *, size: int = 8, nonnegative: bool = True
) -> tuple[tuple[int, ...], ...]:
    need(isinstance(value, list) and len(value) == size, "table row count mismatch")
    rows: list[tuple[int, ...]] = []
    for row in value:
        need(isinstance(row, list) and len(row) == size, "table column count mismatch")
        parsed: list[int] = []
        for entry in row:
            need(
                isinstance(entry, int) and not isinstance(entry, bool),
                "table contains a non-integer",
            )
            need(not nonnegative or entry >= 0, "count table contains a negative entry")
            parsed.append(entry)
        rows.append(tuple(parsed))
    return tuple(rows)


def parse_lean_count_table(path: Path = DEFAULT_LEAN_SOURCE) -> tuple[tuple[int, ...], ...]:
    """Parse the literal matrix bound to ``repairCounts`` in the Lean mirror."""

    source = path.read_text(encoding="utf-8")
    need(EXPECTED_PAYLOAD_SHA256 in source, "Lean mirror payload hash mismatch")
    need("payload schema v3" in source, "Lean mirror schema marker mismatch")
    match = re.search(
        r"def repairCounts[\s\S]*?:= fun a =>\s*(?P<table>!\[[\s\S]*?\])\s+a\s*\n",
        source,
    )
    need(match is not None, "Lean repairCounts literal not found")
    assert match is not None
    try:
        value = ast.literal_eval(match.group("table").replace("![", "["))
    except (SyntaxError, ValueError) as error:
        raise VerificationError("Lean repairCounts literal is not parseable") from error
    return parse_integer_table(value)


def verify_lean_binding(path: Path = DEFAULT_LEAN_SOURCE) -> None:
    need(
        parse_lean_count_table(path) == EXPECTED_COUNTS,
        "Lean repairCounts literal differs from the vendored payload",
    )


def transpose(table: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    size = len(table)
    need(all(len(row) == size for row in table), "cannot transpose a ragged table")
    return tuple(tuple(table[b][a] for b in range(size)) for a in range(size))


def antisymmetric_part(
    table: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    reversed_table = transpose(table)
    return tuple(
        tuple(table[a][b] - reversed_table[a][b] for b in range(len(table)))
        for a in range(len(table))
    )


def row_normalize(
    table: Sequence[Sequence[int]],
) -> tuple[tuple[Fraction, ...], ...]:
    normalized: list[tuple[Fraction, ...]] = []
    for row in table:
        total = sum(row)
        need(total > 0, "row normalization encountered an empty row")
        normalized.append(tuple(Fraction(entry, total) for entry in row))
    return tuple(normalized)


def select_pair(
    table: Sequence[Sequence[int | Fraction]],
) -> PairDesignation | None:
    best: PairDesignation | None = None
    for a, b in permutations(range(len(table)), 2):
        candidate = PairDesignation(
            a, b, Fraction(table[a][b]), Fraction(table[b][a])
        )
        magnitude = abs(candidate.difference)
        if magnitude == 0:
            continue
        if best is None:
            best = candidate
            continue
        best_magnitude = abs(best.difference)
        if magnitude > best_magnitude or (
            magnitude == best_magnitude and (a, b) < (best.a, best.b)
        ):
            best = candidate
    return best


def select_cycle(
    table: Sequence[Sequence[int | Fraction]],
) -> CycleDesignation | None:
    best: CycleDesignation | None = None
    for a, b, c in permutations(range(len(table)), 3):
        forward = Fraction(table[a][b]) * table[b][c] * table[c][a]
        backward = Fraction(table[b][a]) * table[c][b] * table[a][c]
        candidate = CycleDesignation(a, b, c, forward, backward)
        magnitude = abs(candidate.difference)
        if magnitude == 0:
            continue
        if best is None:
            best = candidate
            continue
        best_magnitude = abs(best.difference)
        if magnitude > best_magnitude or (
            magnitude == best_magnitude
            and (a, b, c) < (best.a, best.b, best.c)
        ):
            best = candidate
    return best


def verify_normalized_designations(table: Sequence[Sequence[int]]) -> None:
    normalized = row_normalize(table)
    need(
        select_pair(normalized) == EXPECTED_NORMALIZED_PAIR,
        "normalized pair maximizer mismatch",
    )
    need(
        select_cycle(normalized) == EXPECTED_NORMALIZED_CYCLE,
        "normalized cycle maximizer mismatch",
    )


def verify_transposition_semantics(
    table: Sequence[Sequence[int]],
    reversed_table: Sequence[Sequence[int]] | None = None,
) -> None:
    expected_reverse = transpose(table)
    actual_reverse = expected_reverse if reversed_table is None else tuple(
        tuple(row) for row in reversed_table
    )
    need(actual_reverse == expected_reverse, "reversal table is not the transpose")
    need(transpose(actual_reverse) == tuple(tuple(row) for row in table), "reversal not involutive")

    for a, b in permutations(range(len(table)), 2):
        original = table[a][b] - table[b][a]
        reversed_value = actual_reverse[a][b] - actual_reverse[b][a]
        need(reversed_value == -original, "pair gap did not flip under transpose")

    for a, b, c in permutations(range(len(table)), 3):
        original = (
            table[a][b] * table[b][c] * table[c][a]
            - table[b][a] * table[c][b] * table[a][c]
        )
        reversed_value = (
            actual_reverse[a][b]
            * actual_reverse[b][c]
            * actual_reverse[c][a]
            - actual_reverse[b][a]
            * actual_reverse[c][b]
            * actual_reverse[a][c]
        )
        need(reversed_value == -original, "cycle gap did not flip under transpose")

    reversed_pair = select_pair(actual_reverse)
    reversed_cycle = select_cycle(actual_reverse)
    need(reversed_pair is not None, "transposed pair designation disappeared")
    need(reversed_cycle is not None, "transposed cycle designation disappeared")
    need(
        (reversed_pair.a, reversed_pair.b) == (3, 4)
        and reversed_pair.difference == -EXPECTED_RAW_PAIR.difference,
        "transposed pair maximizer mismatch",
    )
    need(
        (reversed_cycle.a, reversed_cycle.b, reversed_cycle.c) == (3, 4, 5)
        and reversed_cycle.difference == -EXPECTED_RAW_CYCLE.difference,
        "transposed cycle maximizer mismatch",
    )


def symmetric_control() -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(a + b + 1 for b in range(8)) for a in range(8))


def verify_symmetric_control(
    control: Sequence[Sequence[int]] | None = None,
) -> None:
    table = symmetric_control() if control is None else tuple(
        tuple(row) for row in control
    )
    need(table == transpose(table), "synthetic control is not exactly symmetric")
    need(select_pair(table) is None, "symmetric control has a raw pair orientation")
    need(select_cycle(table) is None, "symmetric control has a raw cycle orientation")

    normalized = row_normalize(table)
    need(
        select_cycle(normalized) is None,
        "symmetric control has a normalized cycle orientation",
    )
    normalized_pair = select_pair(normalized)
    need(
        normalized_pair
        == PairDesignation(0, 7, Fraction(2, 9), Fraction(2, 23)),
        "symmetric-control normalized pair diagnostic mismatch",
    )

    row_sums = tuple(sum(row) for row in table)
    total = sum(row_sums)
    for a, b in permutations(range(8), 2):
        stationary_a = Fraction(row_sums[a], total)
        stationary_b = Fraction(row_sums[b], total)
        need(
            stationary_a * normalized[a][b]
            == stationary_b * normalized[b][a],
            "symmetric control failed exact detailed balance",
        )


def verify_serialized_pair(value: Any, expected: PairDesignation) -> None:
    need(isinstance(value, dict), "serialized pair is not an object")
    need(
        value
        == {
            "bucket_from": expected.a,
            "bucket_to": expected.b,
            "count_forward": int(expected.forward),
            "count_backward": int(expected.backward),
            "difference": int(expected.difference),
        },
        "serialized raw pair designation mismatch",
    )


def verify_serialized_cycle(value: Any, expected: CycleDesignation) -> None:
    need(isinstance(value, dict), "serialized cycle is not an object")
    need(
        value
        == {
            "buckets": [expected.a, expected.b, expected.c],
            "forward_product": int(expected.forward),
            "backward_product": int(expected.backward),
            "difference": int(expected.difference),
        },
        "serialized raw cycle designation mismatch",
    )


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def pair_report(value: PairDesignation) -> dict[str, Any]:
    return {
        "buckets": [value.a, value.b],
        "forward": fraction_text(value.forward),
        "backward": fraction_text(value.backward),
        "difference": fraction_text(value.difference),
    }


def cycle_report(value: CycleDesignation) -> dict[str, Any]:
    return {
        "buckets": [value.a, value.b, value.c],
        "forward": fraction_text(value.forward),
        "backward": fraction_text(value.backward),
        "difference": fraction_text(value.difference),
    }


def verify_payload(payload: dict[str, Any]) -> dict[str, Any]:
    need(set(payload) == EXPECTED_TOP_LEVEL_KEYS, "top-level payload keys mismatch")
    need(payload.get("schema") == EXPECTED_SCHEMA, "schema mismatch")
    need(payload.get("analysis_status") == "post_hoc_diagnostic", "analysis status promoted")
    need(payload.get("statistic_preregistered") is False, "statistic falsely preregistered")
    need(
        payload.get("designation_rule_preregistered") is False,
        "designation rule falsely preregistered",
    )
    need(payload.get("eligible_as_validation") is False, "post-hoc payload promoted to validation")
    need(payload.get("provenance") == EXPECTED_PROVENANCE, "provenance metadata mismatch")
    need(tuple(payload.get("buckets", ())) == EXPECTED_BUCKETS, "bucket alphabet mismatch")
    need(payload.get("observer_count") == EXPECTED_OBSERVER_COUNT, "observer total mismatch")
    need(
        payload.get("transition_count") == EXPECTED_TRANSITION_COUNT,
        "declared transition total mismatch",
    )

    table = parse_integer_table(payload.get("ordered_counts"))
    need(table == EXPECTED_COUNTS, "pinned ordered count table mismatch")
    row_sums = tuple(sum(row) for row in table)
    need(row_sums == EXPECTED_ROW_SUMS, "count-table row sums mismatch")
    need(sum(row_sums) == EXPECTED_TRANSITION_COUNT, "recomputed transition total mismatch")

    serialized_current = parse_integer_table(
        payload.get("current_antisymmetric_part"), nonnegative=False
    )
    need(
        serialized_current == antisymmetric_part(table),
        "serialized antisymmetric current mismatch",
    )

    raw_pair = select_pair(table)
    raw_cycle = select_cycle(table)
    need(raw_pair == EXPECTED_RAW_PAIR, "raw pair maximizer mismatch")
    need(raw_cycle == EXPECTED_RAW_CYCLE, "raw cycle maximizer mismatch")
    need(payload.get("pair_orientation_nonempty") is True, "raw pair flag mismatch")
    need(payload.get("cycle_orientation_nonempty") is True, "raw cycle flag mismatch")
    verify_serialized_pair(payload.get("designated_pair"), EXPECTED_RAW_PAIR)
    verify_serialized_cycle(payload.get("designated_cycle"), EXPECTED_RAW_CYCLE)

    verify_normalized_designations(table)
    verify_transposition_semantics(table)
    verify_symmetric_control()

    return {
        "status": "PASS",
        "payload_sha256": EXPECTED_PAYLOAD_SHA256,
        "analysis_status": payload["analysis_status"],
        "eligible_as_validation": payload["eligible_as_validation"],
        "observer_count": EXPECTED_OBSERVER_COUNT,
        "transition_count": EXPECTED_TRANSITION_COUNT,
        "row_sums": list(EXPECTED_ROW_SUMS),
        "raw_pair_maximizer": pair_report(EXPECTED_RAW_PAIR),
        "normalized_pair_maximizer": pair_report(EXPECTED_NORMALIZED_PAIR),
        "raw_cycle_maximizer": cycle_report(EXPECTED_RAW_CYCLE),
        "normalized_cycle_maximizer": cycle_report(EXPECTED_NORMALIZED_CYCLE),
        "transposition": "raw pair and cycle gaps flip exactly",
        "symmetric_control": (
            "raw pair/cycle and normalized cycle gaps vanish; exact detailed "
            "balance holds"
        ),
    }


def verify_file(path: Path = DEFAULT_PAYLOAD) -> dict[str, Any]:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    need(digest == EXPECTED_PAYLOAD_SHA256, "vendored payload SHA-256 mismatch")
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("vendored payload is not valid JSON") from error
    need(isinstance(payload, dict), "vendored payload root is not an object")
    report = verify_payload(payload)
    verify_lean_binding()
    report["lean_binding"] = "RepairCurrentOrientation.repairCounts matches all 64 payload literals"
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, default=DEFAULT_PAYLOAD)
    args = parser.parse_args()
    try:
        report = verify_file(args.payload)
    except (OSError, VerificationError) as error:
        print(f"FAIL: {error}")
        return 1
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
