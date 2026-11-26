"""Interpret a finite read/repair tree on independently prepared sharp interventions."""

import argparse
import gzip
import hashlib
import json

from .. import apparatus
from .machine import Bench
from .measure import HERE
from .programs import accepted, read_family


def execute(bench, tree):
    """A node is a bounded mean, a finite observation menu, or a material publication."""
    public = []
    node = tree
    while node["op"] != "publish":
        if node["op"] == "mean":
            outcome = bench.mean_axes(node["axes"])
            if outcome != "accepted":
                return dict(prefix=public, refusal=outcome, flags=bench.public(None))
            node = node["next"]
        elif node["op"] == "inspect":
            station = node["station"]
            read = tuple(bench.read(r, station) for r in node["registers"])
            public.append(read)
            menu = node["cases"]
            if read not in menu:
                return dict(
                    prefix=public, refusal="outside_read_menu", flags=bench.public(None)
                )
            node = menu[read]
        else:
            raise ValueError("Unknown finite program step")
    choices = node["possible"]
    choice = accepted(choices, tuple(public))
    return dict(prefix=public, classification=choice, flags=bench.public(choice))


def sample(next_node, possible):
    return dict(
        op="inspect",
        station=0,
        registers=(0,),
        cases={(v,): next_node for v in possible},
    )


def move(pair, next_node):
    return dict(op="mean", axes=[pair], next=next_node)


def program(first):
    # The full possible-observation relation is fixed without a source label on input.
    menu = {
        ((64,), (48,)): {"low"},
        ((64,), (80,)): {"high"},
        ((64,), (64,), (64,), (64,)): {"low", "high"},
    }
    finish = dict(op="publish", possible=menu)
    if first == (0, 1):
        after = sample(finish, [48, 80])
    else:
        after = sample(
            move((0, 1), sample(move((0, 2), sample(finish, [64])), [64])), [64]
        )
    return sample(move(first, after), [64])


def report():
    board = Bench()
    # Complementary outputs of the native half exchange have exact mass 1/2.
    # Check its registered matrix and archive the exact dyadic law on every platform.
    checked = apparatus.calendar_bit(0.5)
    if any(abs(value - 0.5) > apparatus.TOL for value in checked):
        raise ArithmeticError(
            "The physical calendar disagrees with its exact half-exchange law"
        )
    calendar = [0.5, 0.5]
    runs = []
    for first in ((0, 1), (1, 2)):
        for src in ((64, 32, 96), (64, 96, 32)):
            board.start([*src] + [0] * (len(board.words) - 3))
            returned = execute(board, program(first))
            runs.append(
                dict(
                    first=first,
                    source=src,
                    returned=returned,
                    final=read_family(board, range(3)),
                    tape=board.tape,
                )
            )
    sources = {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest()
        for f in (
            HERE / "CHECKPOINT_PROTOCOL.md",
            HERE / "checkpoint.py",
            HERE / "records.json.gz",
        )
    }
    return dict(
        provenance="Exhaustive conditional sharp branch effects. Calendar probabilities from a separate "
        "native half-exchange and two-reference read. No finite-shot sampling.",
        sources=sources,
        calendar_stay_cross=calendar,
        runs=runs,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", type=str)
    options = cli.parse_args()
    payload = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if options.output:
        with open(options.output, "xb") as handle:
            handle.write(gzip.compress(payload, mtime=0))
    elif payload != gzip.decompress((HERE / "checkpoint.json.gz").read_bytes()):
        raise ValueError("Checkpoint records differ")
