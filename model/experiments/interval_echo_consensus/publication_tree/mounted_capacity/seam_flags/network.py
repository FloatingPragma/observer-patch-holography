"""Count independently controllable public seam flags using native sharp interventions."""

import argparse
import gzip
import json
from functools import cache
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np

from .apparatus import exchange_word, federation, local_calibration, route

HERE = Path(__file__).resolve().parent


@cache
def station(ports):
    """Ports use the local cell chart. Audit every conditional sharp source word."""
    graph = nx.Graph((a, b) for a, b in federation(0).seams if a < 12 and b < 12)
    other = sorted(set(graph) - set(ports))
    spares = other[: len(ports)]
    standard = other[len(ports)]
    steps = [
        exchange_word(route(graph, a, b)) for a, b in zip(ports, spares, strict=True)
    ]
    for p in ports:
        route(graph, standard, min(graph[p]), avoid=(p,))
    table = []
    for order in product(range(2), repeat=len(ports)):
        loads = [0] * 12
        loads[standard] = 1
        for p in spares:
            loads[p] = 2
        for bit, swaps in zip(order, steps, strict=True):
            if bit:
                for a, b in swaps:
                    loads[a], loads[b] = loads[b], loads[a]
        table.append(tuple(loads[p] // 2 for p in ports))
        for bit, swaps in zip(order, steps, strict=True):
            if bit:
                for a, b in swaps:
                    loads[a], loads[b] = loads[b], loads[a]
        if (
            any(loads[p] for p in ports)
            or any(loads[p] != 2 for p in spares)
            or loads[standard] != 1
        ):
            raise ValueError("Stock not restored after a second drive")
    return np.array(table, dtype=np.uint8)


@cache
def layout(level):
    b = federation(level)
    edges = np.array(sorted(tuple(sorted(e)) for e in b.external))
    owners = edges // 12
    n = b.cells
    local = [
        sorted(int(k % 12) for k in edges.ravel() if k // 12 == a) for a in range(n)
    ]
    if any(len(set(row)) != 3 for row in local):
        raise ValueError("Need three physical contacts per owner")
    graph0 = {tuple(sorted((p, q))) for p, q in b.seams if p < 12 and q < 12}
    expected = {
        (p % 12, q % 12) if p % 12 < q % 12 else (q % 12, p % 12)
        for p, q in federation(0).seams
        if p < 12 and q < 12
    }
    physical = {
        (p // 12, *sorted((p % 12, q % 12))) for p, q in b.seams if p // 12 == q // 12
    }
    if graph0 != expected or physical != {
        (a, p, q) for a in range(n) for p, q in expected
    }:
        raise ValueError("Station routes do not match the installed local contacts")
    table = np.stack([station(tuple(row)) for row in local])
    slots = np.array([[local[p // 12].index(p % 12) for p in pair] for pair in edges])
    return owners, slots, table


def encode(bits):
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def apply(level, menu, drive, writer="repair", twice=False):
    owners, slots, table = layout(level)
    n = table.shape[0]
    bits = np.zeros((n, 3), dtype=np.uint8)
    if menu == "owner":
        bits[:] = np.array([(drive >> i) & 1 for i in range(n)], dtype=np.uint8)[
            :, None
        ]
    else:
        bits[owners[:, 0], slots[:, 0]] = [(drive >> j) & 1 for j in range(len(owners))]
    # Every stock drive is an audited endpoint transposition (including two offers).
    selector = bits @ np.array([4, 2, 1], dtype=np.uint8)
    if twice:
        selector = np.zeros_like(
            selector
        )  # exact second-drive cancellation in station()
    addressed = table[np.arange(n), selector]
    ends = addressed[owners, slots]
    flags = local_calibration()[1][writer][ends[:, 0], ends[:, 1]]
    return [hex(encode(col)) for col in flags.T]


def acquire():
    rng = np.random.Generator(np.random.PCG64(1282612))
    rows = []
    for level in range(3):
        owners, _, tbl = layout(level)
        for menu in ("owner", "seam"):
            size = len(tbl) if menu == "owner" else len(owners)
            held = [
                encode(rng.integers(2, size=size, dtype=np.uint8)) for _ in range(11)
            ]
            held.append((1 << size) - 1)
            rows.append(
                dict(
                    level=level,
                    menu=menu,
                    width=len(owners),
                    empty=apply(level, menu, 0),
                    singles=[apply(level, menu, 1 << i) for i in range(size)],
                    held=[[hex(x), apply(level, menu, x)] for x in held],
                    twice=apply(level, menu, held[-1], twice=True),
                    changed=apply(level, menu, held[-1], writer="unchanged_loads"),
                )
            )
    return rows


def codeword(raw, width):
    if not isinstance(raw, list) or len(raw) != 2:
        raise ValueError("Missing endpoint")
    a, b = (int(v, 16) for v in raw)
    if a != b or not 0 <= a < 1 << width:
        raise ValueError("Inconsistent endpoint words")
    return a


def reduce_word(word, basis):
    for pivot in sorted(basis, reverse=True):
        if word & (1 << pivot):
            word ^= basis[pivot]
    return word


def score(rows):
    result = {}
    for row in rows:
        level = row["level"]
        menu = row["menu"]
        key = level, menu
        n = 20 * 4**level
        width = 3 * n // 2
        if (
            level not in range(3)
            or menu not in ("owner", "seam")
            or key in result
            or row["width"] != width
        ):
            raise ValueError("Unexpected panel")

        def read(raw):
            return codeword(raw, width)

        singles = [read(x) for x in row["singles"]]
        size = n if menu == "owner" else width
        if (
            len(singles) != size
            or len(row["held"]) != 12
            or any(read(row[k]) for k in ("empty", "changed", "twice"))
        ):
            raise ValueError("Incomplete calibration or failed control")
        if int(row["held"][-1][0], 16) != (1 << size) - 1:
            raise ValueError("Missing full drive")
        for source, output in row["held"]:
            target = int(source, 16)
            if not 0 <= target < 1 << size:
                raise ValueError("Out-of-range intervention")
            forecast = 0
            for i, word in enumerate(singles):
                if target >> i & 1:
                    forecast ^= word
            if forecast != read(output):
                raise ValueError("Held word disagrees with calibration")
        basis = {}
        for word in singles:
            remaining = reduce_word(word, basis)
            if remaining:
                basis[remaining.bit_length() - 1] = remaining
        # In the owner menu the source identity carried by each row recovers the measured graph.
        if menu == "owner":
            ends = [
                [i for i, w in enumerate(singles) if w >> j & 1] for j in range(width)
            ]
            if any(len(pair) != 2 for pair in ends) or any(
                w.bit_count() != 3 for w in singles
            ):
                raise ValueError(
                    "Flag does not have two owners and three ends per station"
                )
            if not nx.is_connected(nx.Graph(ends)):
                raise ValueError("Reconstructed graph is disconnected")
        unit_count = sum(reduce_word(1 << j, basis) == 0 for j in range(width))
        result[key] = dict(
            rank=len(basis),
            recoverable_single_flags=unit_count,
            held=len(row["held"]),
            all_ports_flip=read(row["held"][-1][1]).bit_count(),
        )
    if len(result) != 6:
        raise ValueError("Missing network")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = acquire()
    result = score(data)
    path = HERE / "network.json.gz"
    raw = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    if args.write:
        with path.open("xb") as f:
            f.write(gzip.compress(raw, mtime=0))
    elif gzip.decompress(path.read_bytes()) != raw:
        raise ValueError("Network archive mismatch")
    for k, v in result.items():
        print(k, v)


if __name__ == "__main__":
    main()
