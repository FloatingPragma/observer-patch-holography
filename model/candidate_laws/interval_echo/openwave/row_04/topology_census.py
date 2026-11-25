"""Exact connected-pair certificates over every cell of the three native boards.

This is a retrospective extension of the prepared eleven-record experiment.
It enumerates operator columns, not native event calls or observed trajectories.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from .corridor_pair import certify_native_collar, seal_path
from .open_pair_release import source_bound


HERE = Path(__file__).parent
ARCHIVE = HERE / "topology_census.json"


def edge(a, b):
    return tuple(sorted((int(a), int(b))))


def layouts(board, home):
    """Enumerate every eligible motif before applying a certificate or ranking."""
    local = range(12 * home, 12 * (home + 1))
    for center in local:
        if board.degree[center] != 5:
            continue
        tips = [p for p in board.neighbours[center]
                if p // 12 == home and board.degree[p] == 6]
        for left, right in combinations(tips, 2):
            if right in board.neighbours[left]:
                continue
            receivers = (set(board.neighbours[left]) & set(board.neighbours[right])) - {center}
            for receiver in sorted(p for p in receivers if p // 12 == home):
                yield (left, center, right), receiver


def neutral_edges(board, path, receiver):
    ports = set(path)
    openings = {edge(path[0], receiver), edge(path[2], receiver)}
    return tuple(e for e in board.seams if len(ports.intersection(e)) == 1 and e not in openings)


def connected(board, written):
    reached, pending = {0}, [0]
    while pending:
        port = pending.pop()
        for neighbour in board.neighbours[port]:
            if neighbour not in reached and edge(port, neighbour) not in written:
                reached.add(neighbour)
                pending.append(neighbour)
    return len(reached) == board.n


def native_blocks():
    """Check every blank binary seam column, including complementary outputs."""
    law = candidate.native_law
    kinetic, contrast = law.local_term_parts()
    blank = np.eye(law.DIM)[:, :law.N_SRC]
    binary = ((0, 0), (0, 1), (1, 0), (1, 1))
    for word in binary:
        index = law.SRC.index(word)
        sharp = blank[:, index]
        opposite = blank[:, law.SRC.index(word[::-1])]
        if (not np.allclose(kinetic @ sharp, (sharp - opposite) / 2, rtol=0, atol=5e-14)
                or not np.allclose(contrast @ sharp, (word[0] - word[1]) ** 2 * sharp,
                                   rtol=0, atol=5e-14)):
            raise ArithmeticError("Native blank binary column differs from the certified law")
    certify_native_collar()


def action(board, written, state, *, preparation=False, path=None):
    """Integer action (2K, C) with every possible two-load destination retained."""
    kinetic, contrast = Counter(), Counter()
    addressed = {edge(path[0], path[1]), edge(path[1], path[2])} if preparation else None
    for occupied, amplitude in state.items():
        occupied = set(occupied)
        incident = {edge(p, q) for p in occupied for q in board.neighbours[p]}
        for contact in incident:
            if contact in written or (addressed is not None and contact not in addressed):
                continue
            if len(occupied.intersection(contact)) == 1:
                kinetic[tuple(sorted(occupied))] += amplitude
                kinetic[tuple(sorted(occupied.symmetric_difference(contact)))] -= amplitude
                contrast[tuple(sorted(occupied))] += amplitude
        reference = [path[1]] if preparation else occupied
        contrast[tuple(sorted(occupied))] += amplitude * sum(
            6 - int(board.degree[p]) for p in reference if p in occupied)
    return ({p: v for p, v in kinetic.items() if v},
            {p: v for p, v in contrast.items() if v})


def certify_layout(board, path, receiver):
    left, center, right = path
    if (len(set((*path, receiver))) != 4 or len({p // 12 for p in (*path, receiver)}) != 1
            or tuple(board.degree[list(path)]) != (6, 5, 6)
            or any(edge(p, center) not in board.seams for p in (left, right))
            or edge(left, right) in board.seams
            or any(edge(p, receiver) not in board.seams for p in (left, right))):
        raise ValueError("Ineligible connected-pair motif")
    written = set(neutral_edges(board, path, receiver))
    if len(written) != 11 or not connected(board, written):
        raise ArithmeticError("Eleven-record preparation disconnects the native board")
    installation = seal_path(path, tuple(sorted(written)))
    if len(installation["live_swaps"]) != 4 or set(installation["neutral_edges"]) != written:
        raise ArithmeticError("The same two loads cannot install the declared records")

    x, y, z = edge(left, center), edge(left, right), edge(center, right)
    basis = (x, y, z)
    expected = ((1, -1, 0), (-1, 2, -1), (0, -1, 1))
    for column, pair in enumerate(basis):
        kinetic, contrast = action(board, written, {pair: 1}, preparation=True, path=path)
        if (kinetic != {word: expected[row][column] for row, word in enumerate(basis)
                        if expected[row][column]}
                or contrast != {pair: 2}):
            raise ArithmeticError("Local preparation differs from the three-state source")

    dark = {x: 1, z: -1}
    kinetic, contrast = action(board, written, dark)
    if kinetic != {x: 2, z: -2} or contrast != {x: 3, z: -3}:
        raise ArithmeticError("Full-board dark vector is not an exact common eigenvector")
    home = left // 12
    if any(p // 12 != home or q // 12 != home or edge(p, q) not in board.seams for p, q in dark):
        raise ArithmeticError("Dark support misses the original home-seam detector")

    # These are certificate-sensitivity controls. They do not award emergence
    # credit to a property deliberately engineered by the eleven writes.
    if action(board, written, {x: 1, z: 1})[0] == {x: 2, z: 2}:
        raise ArithmeticError("The certificate cannot distinguish the bright relative sign")
    for contact in written:
        if action(board, written - {contact}, dark) == (kinetic, contrast):
            raise ArithmeticError("The certificate misses one released neutral contact")
    return {"path": list(path), "receiver": receiver,
            "neutral_edges": [list(e) for e in sorted(written)]}


def census():
    native_blocks()
    boards = []
    for level in range(3):
        board = federation(level)
        cells = []
        for home in range(board.cells):
            eligible = list(layouts(board, home))
            certificates = [certify_layout(board, path, receiver) for path, receiver in eligible]
            cells.append({"home": home, "eligible_layouts": len(eligible),
                          "certified_layouts": len(certificates),
                          "first_layout": certificates[0] if certificates else None})
        boards.append({"level": level, "cells": board.cells, "ports": board.n,
                       "board_sha256": sha256((CARRIER_HOME / "carriers" /
                                               f"FEDERATION_LEVEL{level}.npz").read_bytes()).hexdigest(),
                       "eligible_layouts": sum(c["eligible_layouts"] for c in cells),
                       "certified_layouts": sum(c["certified_layouts"] for c in cells),
                       "eligible_cells": sum(bool(c["eligible_layouts"]) for c in cells),
                       "ineligible_cells": [c["home"] for c in cells if not c["eligible_layouts"]],
                       "cell_results": cells})
    floor, exhaustion = source_bound()
    sources = {str(Path(module.__file__).resolve().relative_to(HERE.parents[4])):
               sha256(Path(module.__file__).read_bytes()).hexdigest()
               for module in (candidate, candidate.native_law)}
    sources[str(Path(__file__).resolve().relative_to(HERE.parents[4]))] = sha256(
        Path(__file__).read_bytes()).hexdigest()
    return {"schema": "oph.interval_echo.connected_pair_topology.v1",
            "source_sha256": sources,
            "scope": "Retrospective operator certificate for prepared eleven-record material; no flight acquisition, ordinary formation, translating particle or OpenWave award.",
            "source_rounds": 3, "source_attempts": 8,
            "joint_home_floor": str(floor), "source_exhaustion": str(exhaustion),
            "flight_twice_exchange_eigenvalue": 2, "flight_contrast_eigenvalue": 3,
            "native_column_check_tolerance": 5e-14,
            "certificate_arithmetic": "Integer operator action, with every destination from each occupied port retained. The native 99-state column cross-check uses the stated floating-point tolerance.",
            "sensitivity_controls": {"bright_sign": sum(b["eligible_layouts"] for b in boards),
                                     "single_contact_unseal": 11 * sum(b["eligible_layouts"] for b in boards)},
            "all_real_couplings": True, "every_free_exposure": True, "boards": boards}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path, help="explicit first record; refuses an existing path")
    args = cli.parse_args()
    result = census()
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            json.dump(result, destination, indent=2, allow_nan=False)
            destination.write("\n")
    elif json.loads(ARCHIVE.read_text(encoding="utf-8")) != result:
        raise ArithmeticError("Exact topology census differs from its archived certificate")
    for board in result["boards"]:
        print(f"{board['cells']} cells: {board['eligible_cells']}/{board['cells']} eligible; "
              f"{board['certified_layouts']}/{board['eligible_layouts']} layouts certified")


if __name__ == "__main__":
    main()
