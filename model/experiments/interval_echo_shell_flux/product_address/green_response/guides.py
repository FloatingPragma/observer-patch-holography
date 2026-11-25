"""Construct three vertex-disjoint address guides on the deposited spherical boards."""

import argparse
from collections import deque
from functools import cache
import hashlib
from itertools import product
import json
from pathlib import Path

from model.executable_model import CELL_SEAMS, federation
from model.executable_model.refinement import DEPOSIT_SHA256, load_refinement
from .. import measure as rails

RECEIPT = Path(__file__).with_name("guides.json")
INTERNAL = set(CELL_SEAMS)


@cache
def visit_cell(entrance, exit_port):
    """Return one spanning route with prescribed endpoints in the fixed 12-port cell."""
    if entrance == exit_port:
        raise ValueError("Entrance and exit must be different")

    def extend(path):
        if len(path) == 12:
            return path
        for v in range(12):
            if v in path or tuple(sorted((path[-1], v))) not in INTERNAL:
                continue
            if (v == exit_port) != (len(path) == 11):
                continue
            answer = extend((*path, v))
            if answer:
                return answer
        return ()

    result = extend((entrance,))
    if not result:
        raise ValueError("No spanning path for these cell ports")
    return result


def corner_strip(adjacent, face, per_face):
    """Find the shortest owner path between two corner triangles inside this face."""
    owners = range(face * per_face, (face + 1) * per_face)
    inside = {c: sorted(v for v in adjacent[c] if v in owners) for c in owners}
    corners = [c for c in owners if len(inside[c]) == 1]
    if len(corners) != 3:
        raise ValueError("Deposited face lacks the three corner owners")
    start, end = corners[:2]
    queue = deque([start])
    parent = {start: None}
    while queue:
        current = queue.popleft()
        if current == end:
            break
        for following in inside[current]:
            if following not in parent:
                parent[following] = current
                queue.append(following)
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return path[::-1]


def verify(level, external, guides, flag):
    """Check every physical contact, owner, address and directed native pulse."""
    per_face = 4**level
    owners_per_guide = 2 * 2**level - 1
    addresses = 12 * owners_per_guide
    if len(guides) != 3 or any(len(path) != addresses for path in guides):
        raise ValueError("Guide length differs from the full triangular strip")
    if len(flag) != 2 or any(port // (12 * per_face) != 3 for port in flag):
        raise ValueError("Source flag must be in the fourth coarse face")
    seen = set()
    directed = 0
    for face, path in enumerate([*guides, flag]):
        if any(p in seen for p in path) or len(set(path)) != len(path):
            raise ValueError("Reused guide or flag address")
        seen.update(path)
        if face < 3:
            owners = [p // 12 for p in path]
            entry_order = [
                v for i, v in enumerate(owners) if i == 0 or owners[i - 1] != v
            ]
            if (
                len(set(entry_order)) != owners_per_guide
                or any(v // per_face != face for v in owners)
                or owners != [v for v in entry_order for _ in range(12)]
            ):
                raise ValueError(
                    "Guide leaves its face or fails to span each owner once"
                )
        for u, v in zip(path, path[1:]):
            edge = tuple(sorted((u, v)))
            if u // 12 == v // 12:
                native = (edge[0] % 12, edge[1] % 12) in INTERNAL
            else:
                native = edge in external
            if not native:
                raise ValueError("Missing deposited seam")
            for before, after in ((u, v), (v, u)):
                if rails.pulse((u, v), before) != after:
                    raise ValueError("The sharp native pulse failed on a guide or flag")
                if rails.pulse((u, v), before, parked=True) != before:
                    raise ValueError("The parked control changed an address")
                directed += 1
    return directed


def construct(level, external):
    per_face = 4**level
    cells = 20 * per_face
    gateways = {}
    adjacent = [set() for _ in range(cells)]
    seams = {tuple(sorted(map(int, e))) for e in external}
    for u, v in seams:
        a, b = u // 12, v // 12
        if a == b or (a, b) in gateways:
            raise ValueError("Unexpected external cell incidence")
        adjacent[a].add(b)
        adjacent[b].add(a)
        gateways[a, b], gateways[b, a] = u % 12, v % 12
    paths, strips = [], []
    for face in range(3):
        strip = corner_strip(adjacent, face, per_face)
        path = []
        for j, cell in enumerate(strip):
            entry = gateways[cell, strip[j - 1]] if j else None
            leave = gateways[cell, strip[j + 1]] if j + 1 < len(strip) else None
            if entry is None:
                entry = int(leave == 0)
            if leave is None:
                leave = int(entry == 0)
            path.extend(12 * cell + v for v in visit_cell(entry, leave))
        paths.append(path)
        strips.append(strip)
    flag = [36 * per_face, 36 * per_face + 1]
    directed = verify(level, seams, paths, flag)
    middle = len(paths[0]) // 2
    central = [path[middle - 1 : middle + 1] for path in paths]
    released = []
    for tag in product((0, 1), repeat=3):
        if rails.pulse(flag, flag[0]) != flag[1]:
            raise ValueError("Source release flag failed")
        outcome = [
            rails.pulse(pair, pair[0], parked=not bit)
            for pair, bit in zip(central, tag, strict=True)
        ]
        if outcome != [pair[bit] for pair, bit in zip(central, tag, strict=True)]:
            raise ValueError("Corner release failed on a surveyed address")
        released.append(outcome)
    return dict(
        released_port_triples=released,
        level=level,
        cells=cells,
        openwave_source_sha256=DEPOSIT_SHA256[level],
        external_label_sha256=hashlib.sha256(
            json.dumps(sorted(seams)).encode()
        ).hexdigest(),
        owners=strips,
        ports=paths,
        source_flag=flag,
        addresses_per_guide=len(paths[0]),
        directed_sharp_transfers=directed,
        directed_parked_controls=directed,
    )


def check_cell_paths():
    for a in range(12):
        for b in range(12):
            if a == b:
                continue
            p = visit_cell(a, b)
            if (
                set(p) != set(range(12))
                or p[0] != a
                or p[-1] != b
                or any(tuple(sorted((u, v))) not in INTERNAL for u, v in zip(p, p[1:]))
            ):
                raise ValueError("Invalid endpoint-prescribed cell traversal")
    return 12 * 11


def run(deposits=None):
    # The first two results can be reproduced from the bundled labelled seams.
    records = {}
    for level in (1, 2) if deposits is None else (1, 2, 3):
        if deposits is None:
            external = federation(level).external
        else:
            ref = load_refinement(
                deposits / f"FEDERATION_LEVEL{level}.npz", level=level
            )
            external = ref.external
        records[str(level)] = construct(level, external)
    return dict(ordered_cell_endpoint_pairs=check_cell_paths(), levels=records)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--deposits", type=Path, help="Directory with the pinned L1, L2, L3 archives"
    )
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.write and args.deposits is None:
        p.error("--write requires the original three deposits")
    records = run(args.deposits)
    if args.write:
        RECEIPT.write_text(json.dumps(records, indent=2) + "\n")
    else:
        expected = json.loads(RECEIPT.read_text())
        for level, observed in records["levels"].items():
            if observed != expected["levels"][level]:
                raise ValueError("Actual guide differs from recorded hardware")
        if (
            expected["ordered_cell_endpoint_pairs"]
            != records["ordered_cell_endpoint_pairs"]
        ):
            raise ValueError("Cell traversal count differs")
    print(
        [
            (r["cells"], r["addresses_per_guide"], r["directed_sharp_transfers"])
            for r in records["levels"].values()
        ]
    )


if __name__ == "__main__":
    main()
