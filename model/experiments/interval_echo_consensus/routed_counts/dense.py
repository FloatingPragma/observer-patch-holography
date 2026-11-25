"""Complete routed scans from an independently prepared extensive load texture."""

import argparse
from itertools import combinations
import json

import numpy as np

from model.executable_model import federation
from .material import Routing, attempt
from .measure import HERE


def run():
    board = federation(0)
    net = Routing(board)
    rng = np.random.default_rng(22131011)
    initial = rng.integers(0, 3, board.n).tolist()
    dirty = sorted(rng.choice(len(net.chords), 11, replace=False).tolist())
    available = [e for i, e in enumerate(net.chords) if i not in dirty]
    pairs = list(combinations(range(board.n), 2))
    banks = []
    for rev in (False, True):
        chords = available[::-1] if rev else available
        for disabled in (False, True):
            loads = list(initial)
            passed = []
            pulses = 0
            for index in reversed(range(len(pairs))) if rev else range(len(pairs)):
                chord = chords[len(passed)]
                marked, n = attempt(loads, net, pairs[index], chord, disabled)
                pulses += n
                if marked:
                    passed.append([index, list(chord)])
            banks.append(
                dict(
                    reverse=rev,
                    exchange_disabled=disabled,
                    accepted=passed,
                    quiet=len(pairs) - len(passed),
                    routing_pulses=pulses,
                    vacancies=loads.count(0),
                    doubles=loads.count(2),
                    cost=sum(x * x for x in loads),
                    marks=len(dirty) + len(passed),
                )
            )
    return dict(
        initial_word=initial,
        initially_occupied_chords=dirty,
        occupied_memory_labels=[
            dict(chord=list(net.chords[idx]), coin=i % 2, contrast=(i % 5) - 2)
            for i, idx in enumerate(dirty)
        ],
        total_load=sum(initial),
        initial_cost=sum(x * x for x in initial),
        total_offers=len(pairs),
        banks=banks,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    rec = run()
    if args.write:
        (HERE / "dense.json").write_text(json.dumps(rec, indent=2) + "\n")
    elif rec != json.loads((HERE / "dense.json").read_bytes()):
        raise ArithmeticError("Dense complete-scan receipt changed")
    for b in rec["banks"]:
        print(
            b["reverse"],
            b["exchange_disabled"],
            len(b["accepted"]),
            b["marks"],
            b["cost"],
        )
