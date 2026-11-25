"""Fresh complete word counts for a ternary inventory and the cubic valuation witness."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from . import circuit as c

HERE = Path(__file__).resolve().parent
N = {
    "zero": 111 * 2**14,
    "one": 68 * 2**14,
    "plus_x": 2**20,
    "plus_y": 2**20,
    "minus_y": 2**20,
}
LABELS = ("fine", "union", "single")


def committed_coin():
    """Certify fair blank write and read the two definite native material outcomes."""
    law = candidate.native_law
    v = law.echo()[:, law.SRC.index((1, 0))]
    p = (abs(v.reshape(-1, 9)) ** 2).sum(axis=0)
    active = [law.SRC.index(x) for x in ((0, 1), (1, 0))]
    np.testing.assert_allclose(p[active], [0.5, 0.5], rtol=0, atol=1e-14)
    np.testing.assert_allclose(p.sum(), 1, atol=1e-14)
    board = federation(0)
    seams = [(12 * owner, 12 * owner + 1) for owner in (1, 2, 3)]
    if any(b not in board.neighbours[a] for a, b in seams):
        raise ValueError("Missing isolated coin contact")
    return p[active], seams


def classical_labels(y, word):
    if y >= 2:
        return 3
    quota = 3 if y == 0 else 1
    return 0 if word < quota else 1 if word < 2 * quota else 2


def counts(rng, shots, mass):
    mass = np.round(mass, 13)
    return rng.multinomial(shots, mass / mass.sum())


def score(record):
    if set(record["banks"]) != {*LABELS, "split"} or set(record["source_reads"]) != set(
        N
    ):
        raise ValueError("Missing instrument or source preparation")
    for name, banks in [*record["banks"].items(), ("source", record["source_reads"])]:
        if set(banks) != set(N):
            raise ValueError("Incomplete preparation enrollment")
        for src, bank in banks.items():
            if (
                len(bank) != (56 if name == "split" else 7)
                or any(type(n) is not int or n < 0 for n in bank)
                or sum(bank) != N[src]
            ):
                raise ValueError("Incomplete or invalid raw count bank")

    def grouped(menu, bank):
        if menu == "split":
            out = np.zeros(4, dtype=np.int64)
            for y, word in product(range(7), range(8)):
                out[classical_labels(y, word)] += bank[8 * y + word]
            return out
        return np.r_[np.array(bank[:6]).reshape(3, 2).sum(axis=1), bank[-1]]

    parsed = {
        name: {s: grouped(name, v) for s, v in bank.items()}
        for name, bank in record["banks"].items()
    }
    for bank in parsed.values():
        bank["fixture"] = bank["zero"] + bank["one"]

    def rate(menu, src, start, end):
        x = parsed[menu][src]
        n = int(x.sum())
        return sum(map(int, x[start:end])) / n, np.sqrt(
            np.log(2 * 512 / 0.01) / (2 * n)
        )

    def diff(x, y):
        return {
            "delta": round(x[0] - y[0], 9),
            "joint_radius": round(float(x[1] + y[1]), 9),
            "distinguished": bool(abs(x[0] - y[0]) > x[1] + y[1]),
        }

    comparisons = {}
    for src in (*N, "fixture"):
        one = rate("single", src, 0, 1)
        total = rate("union", src, 0, 1)
        comparisons[src] = {
            "first_half": diff(rate("fine", src, 0, 1), one),
            "second_half": diff(rate("fine", src, 1, 2), one),
            "union": diff(rate("fine", src, 0, 2), total),
            "independent_coin_union": diff(rate("split", src, 0, 2), total),
            "independent_coin_half": diff(rate("split", src, 0, 1), one),
        }
    # Round up after including the nine-decimal rounding in both displayed inputs.
    common = [
        comparisons["fixture"][name] for name in ("first_half", "second_half", "union")
    ]
    upper = max(abs(x["delta"]) + x["joint_radius"] for x in common) + 1e-9
    return {
        "maximum_cross_context_disagreement_upper_99": float(
            np.ceil(upper * 1e9) / 1e9
        ),
        "cubic_minimum_disagreement": 1 / 256,
        "grouped_fixture": {m: v["fixture"].tolist() for m, v in parsed.items()},
        "historical_single_against_cubic": diff(
            rate("single", "fixture", 0, 1), (143 / 512, 0)
        ),
        "historical_union_against_cubic": diff(
            rate("union", "fixture", 0, 1), (35 / 64, 0)
        ),
        "comparisons": comparisons,
    }


def run():
    rng = np.random.default_rng(113031207)
    coin, seams = committed_coin()
    partition = np.array([np.prod(coin[list(w)]) for w in product(range(2), repeat=3)])
    program = {menu: c.compile_reader(c.target(menu)) for menu in LABELS}
    program["split"] = c.word_rotation(0, np.pi / 4, -np.pi / 2)
    output = {}
    for k in (c.K, 0):
        initial = {src: c.prepare(src, k) for src in N}
        record = {
            "source_reads": {
                s: counts(rng, N[s], c.read(v)).tolist() for s, v in initial.items()
            },
            "banks": {},
        }
        joined = np.column_stack(list(initial.values()))
        for menu, pulses in program.items():
            final = c.play(pulses, joined, k)
            bank = {}
            for i, src in enumerate(N):
                p = c.read(final[:, i])
                if menu == "split":
                    p = np.kron(p, partition)
                bank[src] = counts(rng, N[src], p).tolist()
            record["banks"][menu] = bank
        record["readback"] = score(record)
        output[str(k)] = record
    return {
        "ports": list(c.ports()),
        "data_words": [list(w) for w in c.WORDS],
        "readout": "word 0..5 plus all other outputs; split uses word-major coin 000..111",
        "coin_seams": seams,
        "pulses": {m: len(v) for m, v in program.items()},
        "shots": N,
        "seed": 113031207,
        "law": output,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_bytes()):
        raise ValueError("Incomplete or changed experiment receipt")
    for k, row in result["law"].items():
        print(
            k,
            row["readback"]["grouped_fixture"],
            row["readback"]["comparisons"]["plus_y"],
        )


if __name__ == "__main__":
    main()
