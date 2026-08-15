#!/usr/bin/env python3
"""Exact discrete Coulomb Green receipt for the committed seam table.

The producer re-derives the thirty-seam incidence from the pinned committed
Lean sources, computes the exact rational Green matrix of the twelve-port
seam Laplacian with ``Fraction`` arithmetic, checks every identity that the
kernel-checked Lean receipts state (Green identity, symmetry, zero row sums,
kernel rank, hop-distance structure, dipole values, Thomson decomposition,
and the ``I - L/60`` uniform-repair bridge), and writes one canonical
receipt with a design-only handoff interface for the instrument lane.

The receipt is finite exact mathematics on the committed table.  It does
not identify a port load with a physical charge, a seam value with a
physical field or flux, or the table with physical space; those
attachments are the open premise-register rows PR-53 and PR-54.  The word
Coulomb names the canonical discrete Gauss solution only.  No public
measurement is read and no decision rule is armed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
CARRIER_LEAN_PATH = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentCarrierQuotient.lean"
GREEN_LEAN_PATH = REPO_ROOT / "Lean" / "Screen" / "DiscreteCoulombGreen.lean"
RECEIPT_PATH = RUNTIME / "discrete_coulomb_green_receipt.json"

SCHEMA = "oph.discrete_coulomb_green.v1"
STATUS = (
    "EXACT_DISCRETE_COULOMB_GREEN_THOMSON_AND_REPAIR_RECEIPT__"
    "PHYSICAL_ATTACHMENTS_PR53_PR54_OPEN"
)

SOURCE_CONTRACT = (
    {
        "path": CARRIER_LEAN_PATH,
        "role": "committed thirty-seam incidence table (seamLeft/seamRight)",
        "bytes": 19162,
        "sha256": "e5f712cbccc5a1462945f0cb511a47064fe10818748498021d5c70ae169a933f",
    },
    {
        "path": GREEN_LEAN_PATH,
        "role": "kernel-checked discrete Coulomb Green receipts (issue 733)",
        "bytes": 50759,
        "sha256": "fe337b4fae484c3da24cff4f0e9fabe56018b64b31b6d2be254c42bdb70df0bf",
    },
)

PORTS = 12
SEAMS = 30
TREE_SEAMS = (0, 1, 2, 3, 4, 7, 8, 11, 14, 17, 20)

GREEN_LEAN_THEOREMS = (
    "laplacian_mul_green",
    "greenMatrix_symm",
    "greenMatrix_row_sum",
    "laplacian_eq_zero_iff",
    "coulombField_gauss",
    "thomson_energy_decomposition",
    "thomson_unique_orthogonal",
    "thomson_minimum",
    "dipole_receipt_values",
    "green_single_source",
    "green_distance_strictly_decreasing",
    "uniformSeamRepair_eq_id_sub_realLaplacian",
    "uniformSeamRepair_fixed_iff_constant",
    "constant_iff_orthogonal_to_boundary_range",
)


class CoulombGreenError(ValueError):
    """The discrete Coulomb Green receipt failed closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CoulombGreenError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def sha256_hex(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def tagged_sha256(raw: bytes) -> str:
    return "sha256:" + sha256_hex(raw)


def load_pinned_source(contract: dict[str, Any]) -> str:
    path = contract["path"]
    raw = path.read_bytes()
    require(len(raw) == contract["bytes"], f"pin byte-count drift: {path.name}")
    require(sha256_hex(raw) == contract["sha256"], f"pin digest drift: {path.name}")
    text = raw.decode("utf-8")
    require("sorry" not in text and "admit" not in text, f"placeholder in {path.name}")
    return text


def parse_fin_vector(text: str, name: str, count: int) -> list[int]:
    match = re.search(
        rf"def {name} : Fin {count} → Fin 12 :=\s*!\[([^\]]*)\]",
        text,
        re.DOTALL,
    )
    require(match is not None, f"missing committed table {name}")
    entries = [int(tok) for tok in match.group(1).replace("\n", " ").split(",")]
    require(len(entries) == count, f"table {name} arity drift")
    return entries


def parse_int_matrix(text: str, name: str, stop: str) -> list[list[int]]:
    start = text.find(f"def {name}")
    require(start >= 0, f"missing integer matrix {name}")
    end = text.find(stop, start)
    require(end > start, f"unterminated integer matrix {name}")
    segment = text[start:end]
    rows = re.findall(r"!\[([\-0-9,\s]+)\]", segment)
    require(len(rows) == PORTS, f"matrix {name} row-count drift")
    parsed = [[int(tok) for tok in row.replace("\n", " ").split(",")] for row in rows]
    require(all(len(row) == PORTS for row in parsed), f"matrix {name} column drift")
    return parsed


def seam_table(carrier_text: str) -> tuple[list[int], list[int]]:
    left = parse_fin_vector(carrier_text, "seamLeft", SEAMS)
    right = parse_fin_vector(carrier_text, "seamRight", SEAMS)
    require(
        all(0 <= l < r < PORTS for l, r in zip(left, right, strict=True)),
        "seam orientation drift",
    )
    require(
        len({(l, r) for l, r in zip(left, right, strict=True)}) == SEAMS,
        "duplicate seam rows",
    )
    return left, right


def laplacian_from_incidence(
    left: list[int], right: list[int]
) -> list[list[Fraction]]:
    lap = [[Fraction(0)] * PORTS for _ in range(PORTS)]
    for l, r in zip(left, right, strict=True):
        lap[l][l] += 1
        lap[r][r] += 1
        lap[l][r] -= 1
        lap[r][l] -= 1
    return lap


def invert(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(matrix)
    work = [
        row[:] + [Fraction(1) if i == j else Fraction(0) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for col in range(n):
        pivot_row = next(
            (r for r in range(col, n) if work[r][col] != 0), None
        )
        require(pivot_row is not None, "regularized Laplacian is singular")
        work[col], work[pivot_row] = work[pivot_row], work[col]
        pivot = work[col][col]
        work[col] = [x / pivot for x in work[col]]
        for r in range(n):
            if r != col and work[r][col] != 0:
                factor = work[r][col]
                work[r] = [x - factor * y for x, y in zip(work[r], work[col])]
    return [row[n:] for row in work]


def green_matrix(lap: list[list[Fraction]]) -> list[list[Fraction]]:
    """The unique symmetric solution of ``L G = I - J/12`` with ``G 1 = 0``,
    computed as ``(L + J)^{-1} - J/144``."""

    regularized = [
        [lap[i][j] + Fraction(1) for j in range(PORTS)] for i in range(PORTS)
    ]
    inverse = invert(regularized)
    return [
        [inverse[i][j] - Fraction(1, 144) for j in range(PORTS)]
        for i in range(PORTS)
    ]


def verify_green_identities(
    lap: list[list[Fraction]], green: list[list[Fraction]]
) -> None:
    for i in range(PORTS):
        require(sum(green[i]) == 0, "green row sum drift")
        for j in range(PORTS):
            require(green[i][j] == green[j][i], "green symmetry drift")
            product = sum(lap[i][k] * green[k][j] for k in range(PORTS))
            expected = (Fraction(1) if i == j else Fraction(0)) - Fraction(1, 12)
            require(product == expected, "green identity drift")


def matrix_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rank = 0
    n_rows, n_cols = len(work), len(work[0])
    for col in range(n_cols):
        pivot_row = next(
            (r for r in range(rank, n_rows) if work[r][col] != 0), None
        )
        if pivot_row is None:
            continue
        work[rank], work[pivot_row] = work[pivot_row], work[rank]
        pivot = work[rank][col]
        work[rank] = [x / pivot for x in work[rank]]
        for r in range(n_rows):
            if r != rank and work[r][col] != 0:
                factor = work[r][col]
                work[r] = [x - factor * y for x, y in zip(work[r], work[rank])]
        rank += 1
    return rank


def bfs_distances(left: list[int], right: list[int]) -> list[list[int]]:
    adjacency: list[set[int]] = [set() for _ in range(PORTS)]
    for l, r in zip(left, right, strict=True):
        adjacency[l].add(r)
        adjacency[r].add(l)
    require(all(len(nbrs) == 5 for nbrs in adjacency), "degree drift")
    table = []
    for source in range(PORTS):
        dist = {source: 0}
        frontier = [source]
        while frontier:
            nxt = []
            for u in frontier:
                for w in adjacency[u]:
                    if w not in dist:
                        dist[w] = dist[u] + 1
                        nxt.append(w)
            frontier = nxt
        require(len(dist) == PORTS, "seam graph is not connected")
        table.append([dist[p] for p in range(PORTS)])
    return table


def hop_classification(left: list[int], right: list[int]) -> list[list[int]]:
    adjacency: list[set[int]] = [set() for _ in range(PORTS)]
    for l, r in zip(left, right, strict=True):
        adjacency[l].add(r)
        adjacency[r].add(l)
    table = []
    for u in range(PORTS):
        row = []
        for v in range(PORTS):
            if u == v:
                row.append(0)
            elif v in adjacency[u]:
                row.append(1)
            elif any(w in adjacency[v] for w in adjacency[u]):
                row.append(2)
            else:
                row.append(3)
        table.append(row)
    return table


def boundary(
    left: list[int], right: list[int], current: list[Fraction]
) -> list[Fraction]:
    load = [Fraction(0)] * PORTS
    for e in range(SEAMS):
        load[right[e]] += current[e]
        load[left[e]] -= current[e]
    return load


def coboundary(
    left: list[int], right: list[int], potential: list[Fraction]
) -> list[Fraction]:
    return [potential[right[e]] - potential[left[e]] for e in range(SEAMS)]


def spanning_tree_check(left: list[int], right: list[int]) -> None:
    require(len(TREE_SEAMS) == PORTS - 1, "tree seam count drift")
    reached = {0}
    changed = True
    while changed:
        changed = False
        for e in TREE_SEAMS:
            l, r = left[e], right[e]
            if l in reached and r not in reached:
                reached.add(r)
                changed = True
            elif r in reached and l not in reached:
                reached.add(l)
                changed = True
    require(len(reached) == PORTS, "declared tree does not span the ports")


def tree_solution(
    left: list[int], right: list[int], load: list[Fraction]
) -> list[Fraction]:
    """The unique Gauss solution supported on the committed spanning tree."""

    require(sum(load) == 0, "tree solve requires a neutral load")
    adjacency: dict[int, list[tuple[int, int]]] = {p: [] for p in range(PORTS)}
    for e in TREE_SEAMS:
        adjacency[left[e]].append((right[e], e))
        adjacency[right[e]].append((left[e], e))
    parent: dict[int, tuple[int, int] | None] = {0: None}
    order = [0]
    queue = [0]
    while queue:
        node = queue.pop(0)
        for neighbour, seam in adjacency[node]:
            if neighbour not in parent:
                parent[neighbour] = (node, seam)
                order.append(neighbour)
                queue.append(neighbour)
    subtree = {p: Fraction(load[p]) for p in range(PORTS)}
    for node in reversed(order):
        entry = parent[node]
        if entry is not None:
            subtree[entry[0]] += subtree[node]
    current = [Fraction(0)] * SEAMS
    for node in order:
        entry = parent[node]
        if entry is None:
            continue
        _, seam = entry
        if right[seam] == node:
            current[seam] = subtree[node]
        else:
            current[seam] = -subtree[node]
    require(boundary(left, right, current) == load, "tree solution drift")
    return current


def fundamental_cycles(
    left: list[int], right: list[int]
) -> list[list[Fraction]]:
    adjacency: dict[int, list[tuple[int, int, int]]] = {
        p: [] for p in range(PORTS)
    }
    for e in TREE_SEAMS:
        adjacency[left[e]].append((right[e], e, 1))
        adjacency[right[e]].append((left[e], e, -1))
    parent: dict[int, tuple[int, int, int] | None] = {0: None}
    queue = [0]
    while queue:
        node = queue.pop(0)
        for neighbour, seam, sign in adjacency[node]:
            if neighbour not in parent:
                parent[neighbour] = (node, seam, sign)
                queue.append(neighbour)

    def path_from_root(node: int) -> list[tuple[int, int]]:
        steps: list[tuple[int, int]] = []
        while parent[node] is not None:
            previous, seam, sign = parent[node]  # type: ignore[misc]
            steps.append((seam, sign))
            node = previous
        return steps

    cycles = []
    for chord in range(SEAMS):
        if chord in TREE_SEAMS:
            continue
        cycle = [Fraction(0)] * SEAMS
        cycle[chord] = Fraction(1)
        for seam, sign in path_from_root(right[chord]):
            cycle[seam] -= sign
        for seam, sign in path_from_root(left[chord]):
            cycle[seam] += sign
        require(
            all(x == 0 for x in boundary(left, right, cycle)),
            "fundamental cycle drift",
        )
        cycles.append(cycle)
    require(len(cycles) == 19, "chord count drift")
    return cycles


def seam_inner(a: list[Fraction], b: list[Fraction]) -> Fraction:
    return sum((x * y for x, y in zip(a, b, strict=True)), start=Fraction(0))


def seam_energy(current: list[Fraction]) -> Fraction:
    return seam_inner(current, current)


def thomson_check(
    left: list[int],
    right: list[int],
    green: list[list[Fraction]],
    cycles: list[list[Fraction]],
    load: list[Fraction],
) -> dict[str, Any]:
    require(sum(load) == 0, "thomson check requires a neutral load")
    potential = [
        sum(green[p][q] * load[q] for q in range(PORTS)) for p in range(PORTS)
    ]
    coulomb = coboundary(left, right, potential)
    require(boundary(left, right, coulomb) == load, "coulomb gauss drift")
    require(
        all(seam_inner(coulomb, cycle) == 0 for cycle in cycles),
        "coulomb cycle-orthogonality drift",
    )
    tree = tree_solution(left, right, load)
    difference = [t - c for t, c in zip(tree, coulomb, strict=True)]
    coulomb_energy = seam_energy(coulomb)
    tree_energy = seam_energy(tree)
    difference_energy = seam_energy(difference)
    require(
        tree_energy == coulomb_energy + difference_energy,
        "thomson decomposition drift",
    )
    require(coulomb_energy <= tree_energy, "thomson minimality drift")
    return {
        "load": {str(p): str(load[p]) for p in range(PORTS) if load[p] != 0},
        "coulomb_energy": str(coulomb_energy),
        "tree_solution_energy": str(tree_energy),
        "cycle_part_energy": str(difference_energy),
        "gauss_boundary_check_passed": True,
        "cycle_orthogonality_check_passed": True,
        "energy_decomposition_exact": True,
    }


def pair_average_identity(
    left: list[int], right: list[int], lap: list[list[Fraction]],
    weight: Fraction,
) -> None:
    """The uniform average of the thirty pair-averaging updates equals
    ``I - weight * L`` exactly; only ``weight = 1/60`` passes."""

    average = [[Fraction(0)] * PORTS for _ in range(PORTS)]
    for e in range(SEAMS):
        for i in range(PORTS):
            for j in range(PORTS):
                if i in (left[e], right[e]):
                    entry = (
                        Fraction(1, 2) if j in (left[e], right[e]) else Fraction(0)
                    )
                else:
                    entry = Fraction(1) if i == j else Fraction(0)
                average[i][j] += Fraction(1, SEAMS) * entry
    for i in range(PORTS):
        for j in range(PORTS):
            expected = (
                Fraction(1) if i == j else Fraction(0)
            ) - weight * lap[i][j]
            require(
                average[i][j] == expected,
                "uniform repair weight drift",
            )


def build_receipt() -> dict[str, Any]:
    carrier_text = load_pinned_source(dict(SOURCE_CONTRACT[0]))
    green_text = load_pinned_source(dict(SOURCE_CONTRACT[1]))
    left, right = seam_table(carrier_text)
    lap = laplacian_from_incidence(left, right)
    require(
        all(lap[p][p] == 5 for p in range(PORTS)),
        "degree-five regularity drift",
    )

    lean_lap = parse_int_matrix(green_text, "laplacianEntriesZ", "/-- One hundred")
    require(
        all(
            lap[i][j] == lean_lap[i][j]
            for i in range(PORTS)
            for j in range(PORTS)
        ),
        "Lean laplacianEntriesZ drift against derived incidence",
    )

    green = green_matrix(lap)
    verify_green_identities(lap, green)
    green_times_180 = [[green[i][j] * 180 for j in range(PORTS)] for i in range(PORTS)]
    require(
        all(x.denominator == 1 for row in green_times_180 for x in row),
        "green denominator drift",
    )
    lean_green = parse_int_matrix(green_text, "greenNumZ", "set_option")
    require(
        all(
            int(green_times_180[i][j]) == lean_green[i][j]
            for i in range(PORTS)
            for j in range(PORTS)
        ),
        "Lean greenNumZ drift against computed Green matrix",
    )

    for fragment in GREEN_LEAN_THEOREMS:
        require(
            f"theorem {fragment}" in green_text,
            f"claim-bearing Lean theorem missing: {fragment}",
        )

    require(matrix_rank(lap) == PORTS - 1, "laplacian rank drift")
    kernel_image = [
        sum(lap[i][j] for j in range(PORTS)) for i in range(PORTS)
    ]
    require(all(x == 0 for x in kernel_image), "constant kernel drift")

    distances = bfs_distances(left, right)
    hops = hop_classification(left, right)
    require(distances == hops, "hop classification / BFS distance mismatch")
    single_source = {0: Fraction(7, 36), 1: Fraction(1, 90),
                     2: Fraction(-7, 180), 3: Fraction(-1, 18)}
    for p in range(PORTS):
        for q in range(PORTS):
            require(
                green[p][q] == single_source[distances[p][q]],
                "single-source distance value drift",
            )
    require(
        Fraction(7, 36) > Fraction(1, 90) > Fraction(-7, 180) > Fraction(-1, 18),
        "single-source ordering drift",
    )

    spanning_tree_check(left, right)
    cycles = fundamental_cycles(left, right)

    dipole = [Fraction(0)] * PORTS
    dipole[0], dipole[1] = Fraction(1), Fraction(-1)
    dipole_potential = [
        sum(green[p][q] * dipole[q] for q in range(PORTS)) for p in range(PORTS)
    ]
    expected_class_values = {
        (0, 1): Fraction(11, 60),
        (1, 0): Fraction(-11, 60),
        (1, 1): Fraction(0),
        (1, 2): Fraction(1, 20),
        (2, 1): Fraction(-1, 20),
        (2, 2): Fraction(0),
        (2, 3): Fraction(1, 60),
        (3, 2): Fraction(-1, 60),
    }
    seen_classes = set()
    for p in range(PORTS):
        key = (distances[0][p], distances[1][p])
        require(key in expected_class_values, "unexpected dipole class")
        require(
            dipole_potential[p] == expected_class_values[key],
            "dipole class value drift",
        )
        seen_classes.add(key)
    require(
        seen_classes == set(expected_class_values),
        "dipole class census drift",
    )

    dipole_thomson = thomson_check(left, right, green, cycles, dipole)
    require(
        dipole_thomson["coulomb_energy"] == "11/30"
        and dipole_thomson["tree_solution_energy"] == "1"
        and dipole_thomson["cycle_part_energy"] == "19/30",
        "dipole energy drift",
    )

    antipodal = [Fraction(0)] * PORTS
    antipodal[3], antipodal[8] = Fraction(1), Fraction(-1)
    require(distances[3][8] == 3, "antipodal fixture drift")
    three_port = [Fraction(0)] * PORTS
    three_port[0], three_port[5], three_port[9] = (
        Fraction(2), Fraction(-1), Fraction(-1),
    )
    spot_checks = [
        dipole_thomson,
        thomson_check(left, right, green, cycles, antipodal),
        thomson_check(left, right, green, cycles, three_port),
    ]

    pair_average_identity(left, right, lap, Fraction(1, 60))

    chords = [e for e in range(SEAMS) if e not in TREE_SEAMS]

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 733,
        "scope": {
            "type": "exact discrete Green/Thomson/repair receipt",
            "statement": (
                "exact rational Green matrix, canonical Coulomb solution, "
                "Thomson decomposition, hop-distance receipt, and the "
                "I - L/60 repair bridge on the committed thirty-seam table"
            ),
            "physical_prediction": False,
            "comparison_permitted": False,
        },
        "source_contracts": [
            {
                "path": contract["path"].relative_to(REPO_ROOT).as_posix(),
                "role": contract["role"],
                "bytes": contract["bytes"],
                "sha256": f"sha256:{contract['sha256']}",
            }
            for contract in SOURCE_CONTRACT
        ],
        "seam_table": {"seam_left": left, "seam_right": right},
        "degree_sequence": [int(lap[p][p]) for p in range(PORTS)],
        "laplacian_matrix": [
            [int(lap[i][j]) for j in range(PORTS)] for i in range(PORTS)
        ],
        "green_matrix_times_180": [
            [int(green_times_180[i][j]) for j in range(PORTS)]
            for i in range(PORTS)
        ],
        "green_identities": {
            "defining_equations": (
                "L G = 1 - J/12, G symmetric, G 1 = 0, over exact rationals"
            ),
            "lg_identity_exact": True,
            "symmetric": True,
            "row_sums_zero": True,
            "laplacian_rank": 11,
            "kernel": "constants",
            "uniqueness": (
                "forced: any two symmetric zero-row-sum solutions differ by "
                "a kernel matrix, and the kernel is the constants"
            ),
        },
        "distance_receipt": {
            "hop_classification_matches_bfs_distance": True,
            "single_source_values": {
                "0": "7/36", "1": "1/90", "2": "-7/180", "3": "-1/18"
            },
            "strictly_decreasing": True,
        },
        "dipole_receipt": {
            "poles": [0, 1],
            "seam": 0,
            "class_values": {
                f"{a},{b}": str(v)
                for (a, b), v in sorted(expected_class_values.items())
            },
            "class_census_complete": True,
            "thomson": dipole_thomson,
        },
        "thomson_spot_checks": spot_checks,
        "spanning_tree": {
            "tree_seams": list(TREE_SEAMS),
            "is_spanning_tree": True,
            "chords": chords,
            "chord_count": len(chords),
        },
        "repair_bridge": {
            "identity": "uniform pair averaging over 30 seams = I - L/60",
            "pair_average_identity_exact": True,
            "fixed_space": "constants = kernel of L",
            "kernel_dimension": 1,
            "wrong_weight_guard": "exercised by the test suite, not here",
        },
        "lean_receipts": {
            "file": GREEN_LEAN_PATH.relative_to(REPO_ROOT).as_posix(),
            "theorems": list(GREEN_LEAN_THEOREMS),
            "sorry_free": True,
        },
        "physical_boundary": {
            "port_load_identified_with_physical_charge": False,
            "seam_value_identified_with_physical_field_or_flux": False,
            "green_matrix_identified_with_laboratory_potential": False,
            "seam_table_identified_with_physical_space": False,
            "continuum_limit_proved": False,
            "laboratory_readout_proved": False,
            "open_premise_rows": ["PR-53", "PR-54"],
            "statement": (
                "the word Coulomb names the canonical discrete Gauss "
                "solution of the committed finite table only"
            ),
        },
        "exposure_boundary": {
            "public_measurement_read": False,
            "comparison_data_read": False,
            "comparison_inputs": [],
            "comparison_permitted": False,
            "score_emitted": False,
            "verdict_emitted": False,
        },
        "handoff_interface": {
            "schema": "oph.discrete_coulomb_green.handoff.v1",
            "consumer": "instrument lane (issue 737)",
            "design_only": True,
            "frozen": False,
            "statement": (
                "interface design for a later preregistration; no decision "
                "rule is armed and no comparison budget exists here"
            ),
            "observables": [
                {
                    "name": "port_potential_difference",
                    "definition": (
                        "phi = greenMatrix.mulVec rho for a neutral rational "
                        "load rho; observable(p, q) = phi p - phi q"
                    ),
                    "exact_type": "rational",
                    "lean_anchor": "OPH.DiscreteCoulombGreen.dipole_receipt_values",
                },
                {
                    "name": "seam_flux",
                    "definition": (
                        "coulombField rho e on the thirty committed seams"
                    ),
                    "exact_type": "rational",
                    "lean_anchor": "OPH.DiscreteCoulombGreen.coulombField_gauss",
                },
                {
                    "name": "chord_field_strength_component",
                    "definition": (
                        "fieldStrength A on the nineteen chord seams of the "
                        "committed spanning tree"
                    ),
                    "chords": chords,
                    "exact_type": "real; exact rational on rational inputs",
                    "lean_anchor": (
                        "OPH.PositionSpaceMaxwellAction."
                        "fieldStrength_determines_up_to_gauge"
                    ),
                },
            ],
            "replay": {
                "producer": (
                    ".venv/bin/python code/electromagnetism/"
                    "discrete_coulomb_green.py --verify"
                ),
                "independent": (
                    ".venv/bin/python code/electromagnetism/"
                    "verify_discrete_coulomb_green_independent.py"
                ),
                "pytest": (
                    ".venv/bin/python -m pytest code/electromagnetism/"
                    "test_discrete_coulomb_green.py -q"
                ),
            },
            "decision_rule_template": {
                "template_only": True,
                "not_a_freeze": True,
                "verdict_labels": ["REPLICATED", "FAILED", "INCONCLUSIVE"],
                "comparison": (
                    "exact rational equality of the declared observables "
                    "against an independently recomputed value"
                ),
                "REPLICATED": "every declared observable matches exactly",
                "FAILED": (
                    "at least one declared observable mismatches under a "
                    "valid replay"
                ),
                "INCONCLUSIVE": (
                    "replay could not be completed (pin drift, missing "
                    "file, or tool failure)"
                ),
                "arming": (
                    "a future instrument-lane freeze must copy this "
                    "template into a registered preregistration; this "
                    "receipt does not arm it"
                ),
            },
        },
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise CoulombGreenError("invalid committed receipt JSON") from error
    require(isinstance(committed, dict), "receipt root is not an object")
    require(raw == canonical_json_bytes(committed), "noncanonical receipt")
    require(committed.get("schema") == SCHEMA, "receipt schema drift")
    require(committed.get("status") == STATUS, "receipt status drift")
    body = {k: v for k, v in committed.items() if k != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == tagged_sha256(canonical_json_bytes(body)),
        "receipt self-digest drift",
    )
    require(
        raw == canonical_json_bytes(build_receipt()),
        "receipt ancestry or result drift",
    )
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if not args.write and not args.verify:
        parser.error("choose --write or --verify")
    if args.write:
        RUNTIME.mkdir(parents=True, exist_ok=True)
        RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    if args.verify:
        verify_committed_receipt()
        print("DISCRETE_COULOMB_GREEN_RECEIPT_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
