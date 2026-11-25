import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from .material import Machine

HERE = Path(__file__).resolve().parent
SMALL = ("nand", ("var", 0), ("var", 1))
LARGE = (
    "nand",
    ("nand", ("var", 0), ("const", True)),
    ("nand", ("var", 1), ("nand", ("var", 0), ("var", 2))),
)


def compile_formula(form, inputs):
    nodes = []

    def visit(f):
        op, x, *rest = f
        if op == "nand":
            args = (visit(x), visit(rest[0]))
        elif op == "var":
            op, args = "copy", (x,)
        else:
            op, args = ("true" if x else "false"), ()
        addr = inputs + len(nodes)
        nodes.append((op, args))
        return addr

    visit(form)
    return nodes


def reference(form, inputs):
    if form[0] == "var":
        return inputs[form[1]]
    if form[0] == "const":
        return form[1]
    return not (reference(form[1], inputs) and reference(form[2], inputs))


def offer_at(machine, state, inputs, nodes, index, reverse=False, mode="native"):
    """One offer of the indefinitely callable address calendar."""
    pos = index % len(nodes)
    i = len(nodes) - 1 - pos if reverse else pos
    return machine.attempt(state, inputs + i, nodes[i], mode)


def run():
    grouped = []
    for level, form, n in ((0, SMALL, 2), (1, LARGE, 3)):
        nodes = compile_formula(form, n)
        machine = Machine(level, n + len(nodes))
        textures = (
            list(product(range(2), repeat=len(nodes)))
            if n == 2
            else [
                tuple((i >> j) & 1 for j in range(len(nodes)))
                for i in (0, 341, (1 << len(nodes)) - 1)
            ]
        )
        for mode in ("native", "parked", "swap_query", "reversed_output"):
            for reverse in (False, True):
                summary = []
                hist = np.zeros((3, 3), int)
                machine.rail_swaps = 0
                for inputs in product(range(2), repeat=n):
                    for stale in textures:
                        state = machine.fresh(inputs + stale)
                        attempts = changes = 0
                        tickets = []
                        for turn in range(len(nodes) + 1):
                            before = changes
                            for _ in nodes:
                                q, r = offer_at(
                                    machine, state, n, nodes, attempts, reverse, mode
                                )
                                hist[q, r] += 1
                                tickets.append(3 * q + r)
                                attempts += 1
                                changes += int(r == 1)
                            if before == changes:
                                break
                        bits = [int(state[row[0][0]] // 2) for row in machine.rows]
                        wanted = list(inputs)
                        for op, args in nodes:
                            if op == "nand":
                                wanted.append(int(not all(wanted[j] for j in args)))
                            elif op == "copy":
                                wanted.append(wanted[args[0]])
                            else:
                                wanted.append(int(op == "true"))
                        if bool(wanted[-1]) != bool(reference(form, inputs)):
                            raise ArithmeticError("Compiler disagreement")
                        summary.append(
                            dict(
                                inputs=inputs,
                                initial=stale,
                                final=bits[n:],
                                tickets=tickets,
                                settled=before == changes,
                                attempts=attempts,
                                writes=changes,
                                all_consistent=bits == wanted,
                                root_correct=bits[-1] == wanted[-1],
                            )
                        )
                grouped.append(
                    dict(
                        owners=machine.car.cells,
                        mode=mode,
                        reverse=reverse,
                        reads=hist.tolist(),
                        routing_swaps=machine.rail_swaps,
                        runs=summary,
                    )
                )
    return grouped


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), separators=(",", ":")) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(out)
    elif path.read_text() != out:
        raise ValueError("Boolean node records changed")
    for bank in json.loads(out):
        print(
            bank["owners"],
            bank["mode"],
            bank["reverse"],
            sum(r["all_consistent"] for r in bank["runs"]),
            len(bank["runs"]),
        )
