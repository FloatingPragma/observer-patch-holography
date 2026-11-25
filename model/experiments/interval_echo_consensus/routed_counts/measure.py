"""Replay ternary-census repair records, retaining quiet attempts and terminal loads."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from .flight import cohort
from .material import Routing, attempt, certify

HERE = Path(__file__).parent
SOURCE_WORDS = [
    "222220000011",
    "220000000111",
    "222222201111",
    "222221111111",
    "000001111111",
]


def sharp_bank(level):
    board = federation(level)
    routing = Routing(board)
    rng = np.random.default_rng(880621 + level)
    sites = list(map(int, rng.choice(board.n, 12, replace=False)))
    menu = list(combinations(sites, 2))
    banks = []
    for word in SOURCE_WORDS:
        for reverse in (False, True):
            offers = menu[::-1] if reverse else menu
            chords = routing.chords[::-1] if reverse else routing.chords
            for disabled in (False, True):
                loads = [1] * board.n
                for s, n in zip(sites, word, strict=True):
                    loads[s] = int(n)
                committed = []
                pulses = 0
                history = []
                partial = None
                for k, pair in enumerate(offers):
                    chord = chords[len(committed)]
                    success, work = attempt(loads, routing, pair, chord, disabled)
                    pulses += work
                    if success:
                        committed.append(list(chord))
                        history.append(k)
                    if k == 9:
                        partial = dict(
                            offers=k + 1,
                            marks=len(committed),
                            vacancy=loads.count(0),
                            double=loads.count(2),
                        )
                v, d = loads.count(0), loads.count(2)
                banks.append(
                    dict(
                        stopped_after_ten=partial,
                        input_word=word,
                        reverse=reverse,
                        exchange_disabled=disabled,
                        offers=len(offers),
                        positive_at=history,
                        marks=committed,
                        vacancy=v,
                        double=d,
                        pi_contacts=pulses,
                        quiet_remaining_pairs_when_normal=(
                            (board.n * (board.n - 1) // 2 - len(offers))
                            if min(v, d) == 0
                            else None
                        ),
                    )
                )
    return dict(
        cells=board.cells,
        ports=board.n,
        tree_edges=board.n - 1,
        chords=len(routing.chords),
        source_ports=sites,
        banks=banks,
    )


def run():
    native = certify()
    geometry = [sharp_bank(level) for level in range(3)]
    results, counts = cohort(federation(0), np.random.default_rng(67093211))
    return dict(native=native, sharp_boards=geometry, coherent_flight=results), counts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    summary, counts = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
        np.savez_compressed(HERE / "counts.npz", **counts)
    else:
        if summary != json.loads((HERE / "results.json").read_bytes()):
            raise ArithmeticError("Published summary changed")
        with np.load(HERE / "counts.npz") as saved:
            if set(saved.files) != set(counts):
                raise ArithmeticError("Changed bank list")
            for k, v in counts.items():
                np.testing.assert_array_equal(v, saved[k])
    print(
        "Checked three federations and", len(counts), "complete coherent-flight banks."
    )


if __name__ == "__main__":
    main()
