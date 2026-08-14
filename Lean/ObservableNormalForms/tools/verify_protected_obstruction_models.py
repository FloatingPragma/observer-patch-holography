#!/usr/bin/env python3
"""Deterministic exact finite explorer for protected-obstruction fixtures.

Public reproduction path:
Lean/ObservableNormalForms/tools/verify_protected_obstruction_models.py

The arithmetic is integer numerator arithmetic over denominator 2. There is
no random seed, floating-point operation, external package, or network input.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Callable, Iterable, Sequence

DENOMINATOR = 2


def compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    if parts == 1:
        return [(total,)]
    return [
        (head,) + tail
        for head in range(total + 1)
        for tail in compositions(total - head, parts - 1)
    ]


def classify(
    kernel: Sequence[Sequence[int]],
    source: int,
    targets: frozenset[int],
    quotient: Callable[[int], int] = lambda x: x,
) -> dict[str, object]:
    states = range(len(kernel))

    def edge(x: int, y: int) -> bool:
        return kernel[x][y] > 0

    if source in targets:
        return {
            "fiber": bool(targets),
            "positive": True,
            "almostSure": True,
            "endpoints": [source],
            "unique": True,
            "tau0": True,
        }

    reachable = {source}
    todo = [source]
    while todo:
        x = todo.pop(0)
        for y in states:
            if y not in targets and edge(x, y) and y not in reachable:
                reachable.add(y)
                todo.append(y)

    endpoints = [
        target for target in targets if any(edge(x, target) for x in reachable)
    ]
    closed_trap = False
    for mask in range(1, 1 << len(kernel)):
        subset = {x for x in states if mask & (1 << x)}
        if (
            all(x in reachable and x not in targets for x in subset)
            and all(not edge(x, y) or y in subset for x in subset for y in states)
        ):
            closed_trap = True
            break

    return {
        "fiber": bool(targets),
        "positive": bool(endpoints),
        "almostSure": bool(endpoints) and not closed_trap,
        "endpoints": endpoints,
        "unique": bool(endpoints)
        and len({quotient(endpoint) for endpoint in endpoints}) == 1,
        "tau0": False,
        "closedTrap": closed_trap,
    }


def stratum(case: dict[str, object]) -> str:
    if not case["fiber"]:
        return "empty-fiber/delta0"
    if not case["positive"]:
        return "unreachable/delta1"
    if not case["almostSure"]:
        return "closed-trap/delta2"
    if not case["unique"]:
        return "ambiguous-endpoint/delta3"
    return "selected-endpoint/L3"


def scheduler(
    kernel: Sequence[Sequence[int]], rewrite: Callable[[int, int], bool]
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    sound = True
    complete = True
    for x, y in itertools.product(range(len(kernel)), repeat=2):
        support = kernel[x][y] > 0
        is_rewrite = bool(rewrite(x, y))
        equality = x == y
        sound_at = not support or is_rewrite or equality
        complete_at = not is_rewrite or support
        sound = sound and sound_at
        complete = complete and complete_at
        rows.append(
            {
                "x": x,
                "y": y,
                "support": support,
                "rewrite": is_rewrite,
                "equality": equality,
                "soundAt": sound_at,
                "completeAt": complete_at,
            }
        )
    return {"sound": sound, "complete": complete, "rows": rows}


def enumerate_kernels(states: int) -> Iterable[tuple[tuple[int, ...], ...]]:
    rows = compositions(DENOMINATOR, states)
    return itertools.product(rows, repeat=states)


def build_output() -> dict[str, object]:
    rows3 = compositions(DENOMINATOR, 3)
    kernels3 = list(enumerate_kernels(3))
    counts: dict[str, int] = {}
    proper_target_cases = 0
    for kernel in kernels3:
        for mask in range(1, 1 << 3):
            if mask & 1:
                continue
            targets = frozenset(i for i in range(3) if mask & (1 << i))
            tag = stratum(classify(kernel, 0, targets))
            counts[tag] = counts.get(tag, 0) + 1
            proper_target_cases += 1

    minima: dict[str, int] = {}
    for state_count in range(1, 4):
        for kernel in enumerate_kernels(state_count):
            for mask in range(1 << state_count):
                if mask & 1:
                    continue
                targets = frozenset(
                    i for i in range(state_count) if mask & (1 << i)
                )
                tag = stratum(classify(kernel, 0, targets))
                if tag != "unreachable/delta1" and tag not in minima:
                    minima[tag] = state_count

    branch3 = ((0, 1, 1), (0, 2, 0), (0, 0, 2))
    deterministic = ((0, 2, 0), (0, 2, 0), (0, 0, 2))
    trap = ((0, 0, 2), (0, 2, 0), (0, 0, 2))
    escape = ((0, 1, 1), (0, 2, 0), (0, 2, 0))

    def rewrite_closed_trap(x: int, y: int) -> bool:
        return (x == 0 and y in (1, 2)) or (x == 2 and y == 2)

    def rewrite_ambiguous_endpoint(x: int, y: int) -> bool:
        return x == 0 and y in (1, 2)

    sched_closed_trap = scheduler(branch3, rewrite_closed_trap)
    sched_ambiguous_endpoint = scheduler(branch3, rewrite_ambiguous_endpoint)
    source00 = sched_closed_trap["rows"][0]
    fixtures = {
        "emptyFiber": classify(branch3, 0, frozenset()),
        "selectedEndpoint": classify(deterministic, 0, frozenset({1})),
        "closedTrap": classify(branch3, 0, frozenset({1})),
        "ambiguousEndpoint": classify(branch3, 0, frozenset({1, 2})),
        "timeZeroHit": classify(branch3, 0, frozenset({0, 1, 2})),
        "quotientControl": classify(
            branch3, 0, frozenset({1, 2}), lambda x: 1 if x else 0
        ),
    }

    checks = [
        (stratum(fixtures["emptyFiber"]) == "empty-fiber/delta0", "empty fiber"),
        (
            stratum(fixtures["selectedEndpoint"]) == "selected-endpoint/L3",
            "selected endpoint",
        ),
        (stratum(fixtures["closedTrap"]) == "closed-trap/delta2", "closed trap"),
        (
            stratum(fixtures["ambiguousEndpoint"]) == "ambiguous-endpoint/delta3",
            "ambiguous endpoint",
        ),
        (fixtures["timeZeroHit"]["tau0"], "time-zero hit control"),
        (fixtures["quotientControl"]["unique"], "quotient control"),
        (
            stratum(classify(deterministic, 0, frozenset({1})))
            == "selected-endpoint/L3",
            "remove trap",
        ),
        (
            stratum(classify(trap, 0, frozenset({1})))
            == "unreachable/delta1",
            "remove hit",
        ),
        (classify(escape, 0, frozenset({1}))["almostSure"], "break closure"),
        (
            sched_closed_trap["sound"] and sched_closed_trap["complete"],
            "closed-trap scheduler",
        ),
        (
            sched_ambiguous_endpoint["sound"]
            and sched_ambiguous_endpoint["complete"],
            "ambiguous-endpoint scheduler",
        ),
        (
            source00["equality"]
            and not source00["support"]
            and not source00["rewrite"],
            "source00",
        ),
        (
            source00["support"]
            != (source00["rewrite"] or source00["equality"]),
            "reject old iff",
        ),
        (
            not scheduler(
                ((0, 1, 1), (0, 2, 0), (1, 1, 0)), rewrite_closed_trap
            )["sound"],
            "sound mutant",
        ),
        (
            not scheduler(deterministic, rewrite_closed_trap)["complete"],
            "complete mutant",
        ),
        (
            minima
            == {
                "empty-fiber/delta0": 1,
                "selected-endpoint/L3": 2,
                "closed-trap/delta2": 3,
                "ambiguous-endpoint/delta3": 3,
            },
            "minima",
        ),
    ]
    failures = [name for passed, name in checks if not passed]
    return {
        "schema": "protected-obstruction-model-explorer-v1",
        "arithmetic": "integer numerators / 2",
        "deterministic": True,
        "seed": None,
        "enumeration": {
            "rows": len(rows3),
            "kernels": len(kernels3),
            "properTargetCases": proper_target_cases,
            "counts": counts,
            "exhaustiveMinimumStates": minima,
        },
        "schedulerRegression": {
            "source00": source00,
            "allPairs": {
                "closedTrap": {
                    "sound": sched_closed_trap["sound"],
                    "complete": sched_closed_trap["complete"],
                },
                "ambiguousEndpoint": {
                    "sound": sched_ambiguous_endpoint["sound"],
                    "complete": sched_ambiguous_endpoint["complete"],
                },
            },
            "oldFalseIffAtSource00": source00["support"]
            == (source00["rewrite"] or source00["equality"]),
        },
        "assertions": {
            "total": len(checks),
            "passed": len(checks) - len(failures),
            "failures": failures,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = build_output()
    rendered = json.dumps(output, indent=2, sort_keys=False) + "\n"
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")
    return 1 if output["assertions"]["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
