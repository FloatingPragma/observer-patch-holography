"""Sample complete survivor histories conditional on a successful common-source comparison."""

import argparse
from bisect import bisect_right
from fractions import Fraction
from itertools import accumulate
import json
from pathlib import Path
from random import Random

from .modular_standard import standard_count
from .offer_service import finite_service

HERE = Path(__file__).resolve().parent
OFFERS = 17
SEED = 307267
DRAWS = 2**21
STEPS = {"0": 0, "1": 0, "-": -1, "+": 1}
MODES = ((0, "even"), (1, "even"), (1, "odd"), (2, "even"), (2, "odd"))


def survivor_table():
    counts = [(1, 1, 1)]  # One equally weighted initial preparation per load.
    for _ in range(OFFERS):
        previous = counts[-1]
        counts.append(
            tuple(
                sum(
                    previous[end - step]
                    for step in STEPS.values()
                    if 0 <= end - step <= 2
                )
                for end in range(3)
            )
        )
    return counts


def load_history(table, end, rank):
    """Unrank one equally weighted surviving chain, including the two distinct idle symbols."""
    final = end
    reversed_symbols = []
    for counts in reversed(table[:-1]):
        for symbol, step in STEPS.items():
            predecessor = end - step
            size = counts[predecessor] if 0 <= predecessor <= 2 else 0
            if rank < size:
                reversed_symbols.append(symbol)
                end = predecessor
                break
            rank -= size
        else:
            raise ArithmeticError("History rank exceeds its endpoint block")
    symbols = "".join(reversed(reversed_symbols))
    load = end
    for symbol in symbols:
        load += STEPS[symbol]
        if not 0 <= load <= 2:
            raise ArithmeticError("Forward replay encountered a rejected offer")
    if rank or load != final:
        raise ArithmeticError("Decoded history has the wrong endpoint")
    return [end, symbols, final]


def pair_history(table, contrast, rank):
    counts = table[-1]
    for left in range(3):
        right = left - contrast
        if 0 <= right <= 2:
            size = counts[left] * counts[right]
            if rank < size:
                x, y = divmod(rank, counts[right])
                return [load_history(table, left, x), load_history(table, right, y)]
            rank -= size
    raise ArithmeticError("Pair rank exceeds its signed-word block")


def material_example(table, sizes, m, signed_word, parity, rank):
    reporter_word = int(parity == "odd")
    rank, reporter_id = divmod(rank, sizes[reporter_word])
    reporter = pair_history(table, reporter_word, reporter_id)
    witnesses = []
    for _ in range(m):
        rank, pair_id = divmod(rank, sizes[abs(signed_word)])
        witnesses.append(pair_history(table, signed_word, pair_id))
    if rank or [l[-1] - r[-1] for l, r in witnesses] != [signed_word] * m:
        raise ArithmeticError("Replayed witness words do not agree")
    if reporter[0][-1] - reporter[1][-1] != reporter_word:
        raise ArithmeticError("Replayed parity reporter rejected the selected outcome")
    return dict(
        signed_word=signed_word,
        material_parity=parity,
        witnesses=witnesses,
        reporter=reporter,
    )


def run():
    m = standard_count()
    service = finite_service(OFFERS, m, cap=2048)
    table = survivor_table()
    a, b, c = table[-1]
    sizes = [a * a + b * b + c * c, a * b + b * c, a * c]
    blocks = [sizes[d] ** m * sizes[int(parity == "odd")] for d, parity in MODES]
    boundaries = list(accumulate(blocks))
    if tuple(Fraction(w, boundaries[-1]) for w in blocks) != service["modes"]:
        raise ArithmeticError(
            "Independent integer blocks differ from the native mode law"
        )
    rng = Random(SEED)
    totals = [0] * len(MODES)
    signs = [{} for _ in MODES]
    rare = {}
    for draw in range(DRAWS):
        index = rng.randrange(boundaries[-1])
        mode = bisect_right(boundaries, index)
        magnitude, parity = MODES[mode]
        # A sharp nonzero word projects to either exchange parity with
        # probability 1/2. The two equally likely original signs cancel this
        # factor in each mode block; restore the orientation for material replay.
        d = -magnitude if magnitude and rng.getrandbits(1) else magnitude
        totals[mode] += 1
        signs[mode][str(d)] = signs[mode].get(str(d), 0) + 1
        if magnitude == 2 and parity not in rare:
            within = index - boundaries[mode - 1]
            rare[parity] = dict(
                draw=draw + 1,
                global_index=str(index),
                **material_example(table, sizes, m, d, parity, within),
            )
    if set(rare) != {"even", "odd"}:
        raise ArithmeticError("The enrolled conditional draw lacked a held parity")
    return dict(
        seed=SEED,
        draws=DRAWS,
        witnesses=m,
        offers=OFFERS,
        offer_symbols=STEPS,
        survivor_path_counts=table[-1],
        signed_pair_counts=sizes,
        mode_order=MODES,
        counts=totals,
        signed_counts=signs,
        conditional_probabilities=list(map(float, service["modes"])),
        first_held_material_examples=rare,
        source_pair_survival_probability=float(service["pair_source_success"]),
        first_saturated_offer_probabilities=list(
            map(float, service["first_failure_by_offer"])
        ),
        unconditioned_epoch={
            k: float(format(service[k], ".12g"))
            for k in ("overflow", "disagreement", "parity_reject", "pointer_flip")
        },
        mean_raw_pairs_per_held_odd=float(
            format(service["mean_pairs_per_held_odd"], ".12g")
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--record", action="store_true", help="create the conditional sample receipt"
    )
    args = cli.parse_args()
    report = run()
    receipt = HERE / "offer_sampler.json"
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.record:
        with receipt.open("x") as stream:
            stream.write(encoded)
    elif json.loads(receipt.read_text()) != json.loads(encoded):
        raise ArithmeticError("The conditional sample differs from its receipt")
    print("Draws:", report["draws"], "mode counts:", report["counts"])
