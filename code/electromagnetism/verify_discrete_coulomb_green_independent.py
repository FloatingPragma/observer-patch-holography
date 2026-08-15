#!/usr/bin/env python3
"""Separately replay the discrete Coulomb Green receipt.

The verifier does not import the producer.  It re-parses the pinned
committed Lean sources with its own scanner, recomputes the exact Green
matrix by a different route (per-column solves of ``L x = e_j - 1/12``
under the zero-sum constraint, instead of the producer's regularized
inverse), recomputes distances by Floyd-Warshall instead of breadth-first
search, rebuilds the tree Gauss solution from the committed
``rationalBoundarySection`` formula instead of subtree sums, reconstructs
the complete expected receipt, and checks the committed receipt byte for
byte: canonical serialization, self-digest, and full semantic equality.

The replay is separate code, not an independent scientific implementation.
It attaches no physical charge, field, space, or readout; PR-53 and PR-54
stay open.  No public measurement is read and no decision rule is armed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_CARRIER = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentCarrierQuotient.lean"
DEFAULT_GREEN_LEAN = REPO_ROOT / "Lean" / "Screen" / "DiscreteCoulombGreen.lean"
DEFAULT_RECEIPT = HERE / "runtime" / "discrete_coulomb_green_receipt.json"

SCHEMA = "oph.discrete_coulomb_green.v1"
STATUS = (
    "EXACT_DISCRETE_COULOMB_GREEN_THOMSON_AND_REPAIR_RECEIPT__"
    "PHYSICAL_ATTACHMENTS_PR53_PR54_OPEN"
)

CARRIER_BYTES = 19162
CARRIER_SHA256 = "e5f712cbccc5a1462945f0cb511a47064fe10818748498021d5c70ae169a933f"
GREEN_LEAN_BYTES = 50759
GREEN_LEAN_SHA256 = (
    "fe337b4fae484c3da24cff4f0e9fabe56018b64b31b6d2be254c42bdb70df0bf"
)

PORTS = 12
SEAMS = 30
TREE_SEAMS = (0, 1, 2, 3, 4, 7, 8, 11, 14, 17, 20)

LEAN_THEOREMS = (
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


class VerificationError(ValueError):
    """The discrete Coulomb Green replay failed closed."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


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


def self_digest(value: dict[str, Any], field: str) -> str:
    body = {key: item for key, item in value.items() if key != field}
    return "sha256:" + hashlib.sha256(canonical_json_bytes(body)).hexdigest()


def load_pinned(path: Path, expected_bytes: int, expected_sha: str) -> str:
    raw = path.read_bytes()
    check(len(raw) == expected_bytes, f"pin byte-count drift: {path.name}")
    check(
        hashlib.sha256(raw).hexdigest() == expected_sha,
        f"pin digest drift: {path.name}",
    )
    text = raw.decode("utf-8")
    check("sorry" not in text and "admit" not in text, f"placeholder: {path.name}")
    return text


def bracket_ints(text: str, opener: str) -> list[int]:
    """Scan integers inside the first ``![ ... ]`` block after ``opener``."""

    start = text.index(opener)
    start = text.index("![", start) + 2
    depth = 1
    out: list[str] = []
    token = ""
    index = start
    while depth > 0:
        char = text[index]
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                break
        if char in "-0123456789":
            token += char
        else:
            if token:
                out.append(token)
                token = ""
        index += 1
    if token:
        out.append(token)
    return [int(tok) for tok in out]


def nested_int_matrix(text: str, opener: str) -> list[list[int]]:
    """Scan a ``![![...], ...]`` matrix after ``opener`` row by row."""

    start = text.index(opener)
    start = text.index("![![", start) + 2
    rows: list[list[int]] = []
    index = start
    while True:
        row_start = text.index("![", index) + 2
        row_end = text.index("]", row_start)
        rows.append(
            [int(tok) for tok in text[row_start:row_end].replace("\n", " ").split(",")]
        )
        after = text[row_end + 1 : row_end + 3]
        if after.startswith("]"):
            break
        index = row_end + 1
        if len(rows) > PORTS:
            raise VerificationError(f"matrix overflow after {opener}")
    check(len(rows) == PORTS, f"matrix row drift after {opener}")
    check(all(len(row) == PORTS for row in rows), f"matrix column drift after {opener}")
    return rows


def solve(matrix: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction]:
    """Gaussian elimination on a (possibly overdetermined) exact system."""

    n_rows = len(matrix)
    n_cols = len(matrix[0])
    work = [matrix[i][:] + [rhs[i]] for i in range(n_rows)]
    pivots: list[int] = []
    rank = 0
    for col in range(n_cols):
        pivot_row = next((r for r in range(rank, n_rows) if work[r][col] != 0), None)
        if pivot_row is None:
            continue
        work[rank], work[pivot_row] = work[pivot_row], work[rank]
        pivot = work[rank][col]
        work[rank] = [x / pivot for x in work[rank]]
        for r in range(n_rows):
            if r != rank and work[r][col] != 0:
                factor = work[r][col]
                work[r] = [x - factor * y for x, y in zip(work[r], work[rank])]
        pivots.append(col)
        rank += 1
    check(rank == n_cols, "underdetermined solve")
    for r in range(rank, n_rows):
        check(work[r][n_cols] == 0, "inconsistent solve")
    solution = [Fraction(0)] * n_cols
    for row, col in enumerate(pivots):
        solution[col] = work[row][n_cols]
    return solution


def independent_green(lap: list[list[Fraction]]) -> list[list[Fraction]]:
    """Column-by-column solve of ``L g_j = e_j - 1/12`` with ``sum g_j = 0``."""

    columns = []
    ones_row = [Fraction(1)] * PORTS
    for j in range(PORTS):
        rhs = [
            (Fraction(1) if i == j else Fraction(0)) - Fraction(1, 12)
            for i in range(PORTS)
        ]
        system = [row[:] for row in lap] + [ones_row]
        column = solve(system, rhs + [Fraction(0)])
        columns.append(column)
    return [[columns[j][i] for j in range(PORTS)] for i in range(PORTS)]


def floyd_warshall(left: list[int], right: list[int]) -> list[list[int]]:
    big = 10**6
    dist = [[0 if i == j else big for j in range(PORTS)] for i in range(PORTS)]
    for l, r in zip(left, right, strict=True):
        dist[l][r] = 1
        dist[r][l] = 1
    for k in range(PORTS):
        for i in range(PORTS):
            for j in range(PORTS):
                through = dist[i][k] + dist[k][j]
                if through < dist[i][j]:
                    dist[i][j] = through
    check(all(dist[i][j] < big for i in range(PORTS) for j in range(PORTS)),
          "graph is not connected")
    return dist


def boundary(left: list[int], right: list[int],
             current: list[Fraction]) -> list[Fraction]:
    load = [Fraction(0)] * PORTS
    for e in range(SEAMS):
        load[right[e]] += current[e]
        load[left[e]] -= current[e]
    return load


def committed_section(load: list[Fraction]) -> list[Fraction]:
    """The committed ``rationalBoundarySection`` formula, transcribed."""

    x = load
    current = [Fraction(0)] * SEAMS
    current[0] = x[1] + x[5] + x[7] + x[11]
    current[1] = x[2] + x[8]
    current[2] = x[3] + x[9]
    current[3] = x[4] + x[10]
    current[4] = x[6]
    current[7] = x[5] + x[11]
    current[8] = x[7]
    current[11] = x[8]
    current[14] = x[9]
    current[17] = x[10]
    current[20] = x[11]
    return current


def inner(a: list[Fraction], b: list[Fraction]) -> Fraction:
    return sum((x * y for x, y in zip(a, b, strict=True)), start=Fraction(0))


def chord_cycles(left: list[int], right: list[int],
                 green: list[list[Fraction]]) -> list[list[Fraction]]:
    """Fundamental chord cycles built from tree paths found by search."""

    tree_adj: dict[int, list[tuple[int, int, int]]] = {p: [] for p in range(PORTS)}
    for e in TREE_SEAMS:
        tree_adj[left[e]].append((right[e], e, 1))
        tree_adj[right[e]].append((left[e], e, -1))
    parent: dict[int, tuple[int, int, int] | None] = {0: None}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbour, seam, sign in tree_adj[node]:
            if neighbour not in parent:
                parent[neighbour] = (node, seam, sign)
                stack.append(neighbour)
    check(len(parent) == PORTS, "declared tree does not span")

    def root_path(node: int) -> list[tuple[int, int]]:
        steps = []
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
        for seam, sign in root_path(right[chord]):
            cycle[seam] -= sign
        for seam, sign in root_path(left[chord]):
            cycle[seam] += sign
        check(all(x == 0 for x in boundary(left, right, cycle)),
              "fundamental cycle drift")
        cycles.append(cycle)
    check(len(cycles) == 19, "chord count drift")
    del green
    return cycles


def thomson_block(left: list[int], right: list[int],
                  green: list[list[Fraction]],
                  cycles: list[list[Fraction]],
                  load: list[Fraction]) -> dict[str, Any]:
    check(sum(load) == 0, "spot load is not neutral")
    potential = [
        sum(green[p][q] * load[q] for q in range(PORTS)) for p in range(PORTS)
    ]
    coulomb = [potential[right[e]] - potential[left[e]] for e in range(SEAMS)]
    check(boundary(left, right, coulomb) == load, "coulomb boundary drift")
    check(all(inner(coulomb, cycle) == 0 for cycle in cycles),
          "coulomb orthogonality drift")
    tree = committed_section(load)
    check(boundary(left, right, tree) == load, "committed section drift")
    difference = [t - c for t, c in zip(tree, coulomb, strict=True)]
    coulomb_energy = inner(coulomb, coulomb)
    tree_energy = inner(tree, tree)
    difference_energy = inner(difference, difference)
    check(tree_energy == coulomb_energy + difference_energy,
          "thomson decomposition drift")
    return {
        "load": {str(p): str(load[p]) for p in range(PORTS) if load[p] != 0},
        "coulomb_energy": str(coulomb_energy),
        "tree_solution_energy": str(tree_energy),
        "cycle_part_energy": str(difference_energy),
        "gauss_boundary_check_passed": True,
        "cycle_orthogonality_check_passed": True,
        "energy_decomposition_exact": True,
    }


def repair_identity(left: list[int], right: list[int],
                    lap: list[list[Fraction]], weight: Fraction) -> bool:
    average = [[Fraction(0)] * PORTS for _ in range(PORTS)]
    for e in range(SEAMS):
        endpoints = (left[e], right[e])
        for i in range(PORTS):
            for j in range(PORTS):
                if i in endpoints:
                    entry = Fraction(1, 2) if j in endpoints else Fraction(0)
                else:
                    entry = Fraction(1) if i == j else Fraction(0)
                average[i][j] += Fraction(1, SEAMS) * entry
    return all(
        average[i][j]
        == (Fraction(1) if i == j else Fraction(0)) - weight * lap[i][j]
        for i in range(PORTS)
        for j in range(PORTS)
    )


def expected_receipt(carrier_text: str, green_text: str) -> dict[str, Any]:
    left = bracket_ints(carrier_text, "def seamLeft")
    right = bracket_ints(carrier_text, "def seamRight")
    check(len(left) == SEAMS and len(right) == SEAMS, "seam table arity drift")
    check(all(0 <= l < r < PORTS for l, r in zip(left, right, strict=True)),
          "seam orientation drift")
    check(len({pair for pair in zip(left, right, strict=True)}) == SEAMS,
          "duplicate seam rows")

    lap = [[Fraction(0)] * PORTS for _ in range(PORTS)]
    for l, r in zip(left, right, strict=True):
        lap[l][l] += 1
        lap[r][r] += 1
        lap[l][r] -= 1
        lap[r][l] -= 1
    check(all(lap[p][p] == 5 for p in range(PORTS)), "degree drift")

    lean_lap = nested_int_matrix(green_text, "def laplacianEntriesZ")
    check(all(lap[i][j] == lean_lap[i][j] for i in range(PORTS)
              for j in range(PORTS)), "Lean laplacianEntriesZ drift")

    green = independent_green(lap)
    for i in range(PORTS):
        check(sum(green[i]) == 0, "green row sum drift")
        for j in range(PORTS):
            check(green[i][j] == green[j][i], "green symmetry drift")
            product = sum(lap[i][k] * green[k][j] for k in range(PORTS))
            check(product == (Fraction(1) if i == j else Fraction(0))
                  - Fraction(1, 12), "green identity drift")
    green_180 = [[green[i][j] * 180 for j in range(PORTS)] for i in range(PORTS)]
    check(all(x.denominator == 1 for row in green_180 for x in row),
          "green denominator drift")
    lean_green = nested_int_matrix(green_text, "def greenNumZ")
    check(all(int(green_180[i][j]) == lean_green[i][j] for i in range(PORTS)
              for j in range(PORTS)), "Lean greenNumZ drift")

    for fragment in LEAN_THEOREMS:
        check(f"theorem {fragment}" in green_text,
              f"claim-bearing Lean theorem missing: {fragment}")

    distances = floyd_warshall(left, right)
    single_source = {0: Fraction(7, 36), 1: Fraction(1, 90),
                     2: Fraction(-7, 180), 3: Fraction(-1, 18)}
    for p in range(PORTS):
        for q in range(PORTS):
            check(green[p][q] == single_source[distances[p][q]],
                  "single-source value drift")

    cycles = chord_cycles(left, right, green)
    chords = [e for e in range(SEAMS) if e not in TREE_SEAMS]

    dipole = [Fraction(0)] * PORTS
    dipole[0], dipole[1] = Fraction(1), Fraction(-1)
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
    dipole_potential = [
        sum(green[p][q] * dipole[q] for q in range(PORTS)) for p in range(PORTS)
    ]
    seen = set()
    for p in range(PORTS):
        key = (distances[0][p], distances[1][p])
        check(key in expected_class_values, "unexpected dipole class")
        check(dipole_potential[p] == expected_class_values[key],
              "dipole class value drift")
        seen.add(key)
    check(seen == set(expected_class_values), "dipole class census drift")

    dipole_thomson = thomson_block(left, right, green, cycles, dipole)
    check(dipole_thomson["coulomb_energy"] == "11/30", "dipole energy drift")

    antipodal = [Fraction(0)] * PORTS
    antipodal[3], antipodal[8] = Fraction(1), Fraction(-1)
    check(distances[3][8] == 3, "antipodal fixture drift")
    three_port = [Fraction(0)] * PORTS
    three_port[0], three_port[5], three_port[9] = (
        Fraction(2), Fraction(-1), Fraction(-1),
    )
    spot_checks = [
        dipole_thomson,
        thomson_block(left, right, green, cycles, antipodal),
        thomson_block(left, right, green, cycles, three_port),
    ]

    check(repair_identity(left, right, lap, Fraction(1, 60)),
          "uniform repair identity drift")

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
                "path": "Lean/Screen/SeamCurrentCarrierQuotient.lean",
                "role": (
                    "committed thirty-seam incidence table "
                    "(seamLeft/seamRight)"
                ),
                "bytes": CARRIER_BYTES,
                "sha256": f"sha256:{CARRIER_SHA256}",
            },
            {
                "path": "Lean/Screen/DiscreteCoulombGreen.lean",
                "role": (
                    "kernel-checked discrete Coulomb Green receipts "
                    "(issue 733)"
                ),
                "bytes": GREEN_LEAN_BYTES,
                "sha256": f"sha256:{GREEN_LEAN_SHA256}",
            },
        ],
        "seam_table": {"seam_left": left, "seam_right": right},
        "degree_sequence": [int(lap[p][p]) for p in range(PORTS)],
        "laplacian_matrix": [
            [int(lap[i][j]) for j in range(PORTS)] for i in range(PORTS)
        ],
        "green_matrix_times_180": [
            [int(green_180[i][j]) for j in range(PORTS)] for i in range(PORTS)
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
            "file": "Lean/Screen/DiscreteCoulombGreen.lean",
            "theorems": list(LEAN_THEOREMS),
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
                    "lean_anchor": (
                        "OPH.DiscreteCoulombGreen.dipole_receipt_values"
                    ),
                },
                {
                    "name": "seam_flux",
                    "definition": (
                        "coulombField rho e on the thirty committed seams"
                    ),
                    "exact_type": "rational",
                    "lean_anchor": (
                        "OPH.DiscreteCoulombGreen.coulombField_gauss"
                    ),
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
    receipt["receipt_sha256"] = self_digest(receipt, "receipt_sha256")
    return receipt


def verify(carrier_path: Path, green_lean_path: Path,
           receipt_path: Path) -> dict[str, Any]:
    carrier_text = load_pinned(carrier_path, CARRIER_BYTES, CARRIER_SHA256)
    green_text = load_pinned(green_lean_path, GREEN_LEAN_BYTES,
                             GREEN_LEAN_SHA256)
    raw = receipt_path.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid receipt JSON") from error
    check(isinstance(committed, dict), "receipt root is not an object")
    check(raw == canonical_json_bytes(committed), "noncanonical receipt")
    check(committed.get("schema") == SCHEMA, "receipt schema drift")
    check(committed.get("status") == STATUS, "receipt status drift")
    check(
        committed.get("receipt_sha256") == self_digest(committed, "receipt_sha256"),
        "receipt self-digest drift",
    )
    check(committed == expected_receipt(carrier_text, green_text),
          "receipt semantic drift")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carrier-lean", type=Path, default=DEFAULT_CARRIER)
    parser.add_argument("--green-lean", type=Path, default=DEFAULT_GREEN_LEAN)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        verify(args.carrier_lean, args.green_lean, args.receipt)
    except (OSError, VerificationError, ValueError) as error:
        print(f"DISCRETE_COULOMB_GREEN_REPLAY_FAIL: {error}", file=sys.stderr)
        return 1
    print("DISCRETE_COULOMB_GREEN_SEPARATE_REPLAY_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
