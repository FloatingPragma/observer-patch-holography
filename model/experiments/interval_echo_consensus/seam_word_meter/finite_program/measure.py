"""Acquire the exact sharp records of the frozen finite observer programs."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

from .machine import Bench, HEADS, REFS, WRITER, decode_word
from . import programs as p

HERE = Path(__file__).resolve().parent


def check(result):
    for key in ("adaptive", "real_menu"):
        for row in result[key]:
            flag = {None: 0, "low": 1, "high": 2}[row["chosen"]]
            if (
                any(decode_word(raw) != flag for raw in row["public"])
                or row["final"][0] != row["final"][1]
            ):
                raise ArithmeticError(
                    "The material publication did not reach both readers"
                )
    if Counter(x["chosen"] for x in result["adaptive"]) != {
        None: 3,
        "low": 3,
        "high": 3,
    }:
        raise ArithmeticError("Adaptive predictions failed")
    if Counter(x["chosen"] for x in result["real_menu"]) != {
        None: 1,
        "low": 5,
        "high": 3,
    }:
        raise ArithmeticError("Incommensurate calibration failed")
    for row in result["real_menu"]:
        if set(map(tuple, row["inputs"])) == {p.STOCK[0], p.STOCK[2]} and row["final"][
            0
        ][4:6] != [4, 12]:
            raise ArithmeticError("Mixed standard routed incorrectly")
    bus = [
        (v - w) // 2
        for v, w in zip(result["bus"]["final"][0][::2], result["bus"]["final"][0][1::2])
    ]
    if bus != [24, -8, 8, -8, 0, 0]:
        raise ArithmeticError("Wrong reusable bus record")
    if result["accumulator"]["final"][0] != [66, 62, 64, 64, 60, 68]:
        raise ArithmeticError("Wrong scaled accumulator record")
    for row in result["hop"]:
        if row["stages"][-1] != [[row["source"], 0, 0, row["source"]]] * 2:
            raise ArithmeticError("Routed archive or retirement error")
    refusals = result["controls"]["refusals"]
    if [r["outcome"] for r in refusals] != [
        "read_disagreement",
        "invalid_standard",
        "off_grid",
    ]:
        raise ArithmeticError("A refused attempt reached a different exit")
    if not all(r["preserved"] for r in refusals):
        raise ArithmeticError("Refusal changed a material port")


def report():
    b = Bench()
    values = (48, 64, 80)
    result = dict(
        provenance="Deterministic simulated native write, load-read, inverse and swap effects. No sampled launches.",
        adaptive=[p.adaptive(b, x, y) for x in values for y in values],
        real_menu=[p.real_menu(b, x, y) for x in p.STOCK for y in p.STOCK],
        bus=p.bus(b),
        accumulator=p.accumulator(b),
        hop=[p.hop(b, x) for x in range(0, 9, 2)],
        controls=dict(
            fold=[p.adaptive(b, x, y, "fold") for x in values for y in values],
            force_low=p.adaptive(b, 80, 80, "force_low"),
            equal_weights=p.real_menu(b, p.STOCK[0], p.STOCK[2], "equal_weights"),
            no_reset=p.hop(b, 6, True),
            refusals=[
                p.refusal(b, 48, 80, "parked_write"),
                p.refusal(b, 48, 80, "parked_sensor"),
                p.refusal(b, 64, 65, "off_grid"),
            ],
        ),
    )
    sources = [
        HERE / name
        for name in ("machine.py", "programs.py", "measure.py", "PROTOCOL.md")
    ]
    from .. import apparatus, chart

    sources.extend(
        Path(module.__file__).resolve() for module in (apparatus, chart, apparatus.LAW)
    )
    root = HERE.parents[4]
    result["sources"] = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sources
    }
    result["layout"] = dict(
        words=b.words,
        flags=b.flags,
        writing_seam=WRITER,
        heads=HEADS,
        references=REFS,
        paths=[
            [u, v, list(b.route(u, v))]
            for u in sorted(
                set(x for word in [*b.words, b.flags] for pair in word for x in pair)
                | set(WRITER)
            )
            for v in (*HEADS, *WRITER)
            if u != v and not (u in WRITER and v in WRITER)
        ],
    )
    check(result)
    return result


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", type=Path)
    opts = cli.parse_args()
    payload = (
        json.dumps(report(), separators=(",", ":"), allow_nan=False) + "\n"
    ).encode()
    if opts.output is None:
        if gzip.decompress((HERE / "records.json.gz").read_bytes()) != payload:
            raise ValueError("The native program records differ")
    else:
        with opts.output.open("xb") as dest:
            dest.write(gzip.compress(payload, mtime=0))


if __name__ == "__main__":
    main()
