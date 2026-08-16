#!/usr/bin/env python3
"""Independent exact replay of the local face-Maxwell receipt.

This file does not import the producer.  It uses a separate source scanner,
fraction-free integer rank elimination, nullspace certificates for the exact
spectra, and an augmented zero-sum solve for the scalar Green potential.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DEFAULT_RECEIPT = HERE / "runtime/local_face_maxwell_action_receipt.json"
SOURCE_ROWS = (
    ("Lean/Screen/SeamCurrentCarrierQuotient.lean", 19162,
     "e5f712cbccc5a1462945f0cb511a47064fe10818748498021d5c70ae169a933f",
     "committed 30-seam table"),
    ("Lean/ObserverPatchHolography/CoreAxioms.lean", 8196,
     "addcea08fe55c0ccc3edbeb90ebd9df4b640c96bd20e9755417bcdb84f093fe4",
     "committed 20 oriented faces"),
    ("Lean/Screen/DiscreteCoulombGreen.lean", 50954,
     "a58830d60688b7df505d2c8cbc1ff2290b057f9be055def25db7006fae7bd078",
     "exact scalar Green operator"),
    ("Lean/Screen/LocalFaceMaxwellAction.lean", 32766,
     "c54da8884d3a3d19ed73baab2218378795fa5f8f2681344d545f285a647cf763",
     "kernel-checked local action and typed join"),
)
THEOREMS = (
    "faceVertices_matches_committed", "face_port_incidence_product_zero",
    "unsigned_face_incidence_breaks_boundary", "localKineticZ_five_per_row",
    "localKineticZ_total_support", "localKineticZ_diagonal_two",
    "ker_faceCurvature_eq_gradient", "range_localMaxwellOperator",
    "localSourcedAction_gauge_invariant_iff", "localStationary_solvable_iff",
    "localSourced_global_minimum", "greenPotential_stationary",
    "scalarCoulombAction_constant_gauge_iff", "scalarCoulomb_global_minimum",
    "staticAction_gauge_invariant", "staticAction_global_minimum",
)
PORTS, SEAMS, FACES = 12, 30, 20


class VerificationError(ValueError):
    pass


def check(ok: bool, msg: str) -> None:
    if not ok:
        raise VerificationError(msg)


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def pinned_text(relative: str, size: int, sha: str) -> str:
    raw = (ROOT / relative).read_bytes()
    check(len(raw) == size, f"byte pin drift: {relative}")
    check(hashlib.sha256(raw).hexdigest() == sha, f"digest pin drift: {relative}")
    text = raw.decode()
    check("sorry" not in text and "admit" not in text, f"placeholder: {relative}")
    return text


def scan_vector(text: str, declaration: str, expected: int) -> list[int]:
    pos = text.index(declaration)
    pos = text.index("![", pos) + 2
    end = text.index("]", pos)
    values = []
    token = ""
    for char in text[pos:end] + " ":
        if char.isdigit() or (char == "-" and not token):
            token += char
        elif token:
            values.append(int(token)); token = ""
    check(len(values) == expected, f"vector arity drift: {declaration}")
    return values


def scan_faces(text: str) -> list[tuple[int, int, int]]:
    pos = text.index("def orientedFaces")
    pos = text.index("[", pos)
    end = text.index("]", pos)
    numbers = []
    token = ""
    for char in text[pos:end] + " ":
        if char.isdigit(): token += char
        elif token: numbers.append(int(token)); token = ""
    check(len(numbers) == 60, "face table arity drift")
    return [tuple(numbers[i:i + 3]) for i in range(0, 60, 3)]


def transpose(a):
    return [list(x) for x in zip(*a, strict=True)]


def mm(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(r, c, strict=True)) for c in bt] for r in a]


def integer_rank(a: list[list[int]]) -> int:
    """Fraction-free elimination with row gcd reduction."""
    w = [row[:] for row in a]
    r = 0
    for c in range(len(w[0])):
        p = next((i for i in range(r, len(w)) if w[i][c]), None)
        if p is None: continue
        w[r], w[p] = w[p], w[r]
        pivot = w[r][c]
        for i in range(r + 1, len(w)):
            factor = w[i][c]
            if not factor: continue
            w[i] = [pivot * x - factor * y for x, y in zip(w[i], w[r], strict=True)]
            g = math.gcd(*map(abs, w[i]))
            if g > 1: w[i] = [x // g for x in w[i]]
        r += 1
    return r


def shifted(a, q: int):
    return [[a[i][j] - q * int(i == j) for j in range(len(a))]
            for i in range(len(a))]


def qmatrix(a, linear: int, constant: int):
    square = mm(a, a)
    return [[square[i][j] + linear * a[i][j] + constant * int(i == j)
             for j in range(len(a))] for i in range(len(a))]


def solve_augmented(matrix: list[list[int]], rhs: list[Fraction]) -> list[Fraction]:
    w = [[Fraction(x) for x in row] + [rhs[i]] for i, row in enumerate(matrix)]
    rows, cols, r = len(w), len(matrix[0]), 0
    pivots = []
    for c in range(cols):
        p = next((i for i in range(r, rows) if w[i][c]), None)
        if p is None: continue
        w[r], w[p] = w[p], w[r]
        q = w[r][c]; w[r] = [x / q for x in w[r]]
        for i in range(rows):
            if i != r and w[i][c]:
                q = w[i][c]
                w[i] = [x - q * y for x, y in zip(w[i], w[r], strict=True)]
        pivots.append(c); r += 1
    check(r == cols, "augmented Green solve underdetermined")
    out = [Fraction(0)] * cols
    for i, c in enumerate(pivots): out[c] = w[i][-1]
    return out


def build_matrices(left, right, faces):
    b = [[0] * SEAMS for _ in range(PORTS)]
    lookup = {}
    for e, (l, r) in enumerate(zip(left, right, strict=True)):
        check(l < r and (l, r) not in lookup, "seam table invalid")
        lookup[l, r] = e; b[l][e] = -1; b[r][e] = 1
    c = [[0] * SEAMS for _ in range(FACES)]
    for f, face in enumerate(faces):
        for u, v in zip(face, face[1:] + face[:1], strict=True):
            key = min(u, v), max(u, v)
            check(key in lookup, "unmatched face edge")
            c[f][lookup[key]] = 1 if u < v else -1
    return b, c


def verify(receipt_path: Path = DEFAULT_RECEIPT) -> dict[str, Any]:
    texts = {Path(rel).name: pinned_text(rel, size, sha)
             for rel, size, sha, _ in SOURCE_ROWS}
    left = scan_vector(texts["SeamCurrentCarrierQuotient.lean"],
                       "def seamLeft", SEAMS)
    right = scan_vector(texts["SeamCurrentCarrierQuotient.lean"],
                        "def seamRight", SEAMS)
    faces = scan_faces(texts["CoreAxioms.lean"])
    local = texts["LocalFaceMaxwellAction.lean"]
    check(all(("theorem " + name) in local for name in THEOREMS),
          "Lean theorem inventory drift")
    b, c = build_matrices(left, right, faces)
    bt, ct = transpose(b), transpose(c)
    h, dual, vertex = mm(ct, c), mm(c, ct), mm(b, bt)
    check(mm(c, bt) == [[0] * PORTS for _ in range(FACES)], "C B^T != 0")
    check(integer_rank(b) == 11 and integer_rank(c) == 19, "rank drift")
    check(sum(x != 0 for row in c for x in row) == 60, "C support drift")
    hrows = [sum(x != 0 for x in row) for row in h]
    check(hrows == [5] * SEAMS and sum(hrows) == 150, "H support drift")
    check(all(h[i][i] == 2 for i in range(SEAMS)), "H diagonal drift")
    unsigned_boundary = mm([[abs(x) for x in row] for row in c], bt)
    check(unsigned_boundary[0][0] == -2 and
          any(x != 0 for row in unsigned_boundary for x in row),
          "orientation-erasure control drift")

    face_dims = {q: FACES - integer_rank(shifted(dual, q))
                 for q in (0, 2, 3, 5)}
    check(face_dims == {0: 1, 2: 5, 3: 4, 5: 4}, "face spectrum drift")
    face_quadratic_dim = FACES - integer_rank(qmatrix(dual, -6, 4))
    check(face_quadratic_dim == 6, "face quadratic spectrum drift")
    vertex_dims = {q: PORTS - integer_rank(shifted(vertex, q)) for q in (0, 6)}
    check(vertex_dims == {0: 1, 6: 5}, "vertex spectrum drift")
    vertex_quadratic_dim = PORTS - integer_rank(qmatrix(vertex, -10, 20))
    check(vertex_quadratic_dim == 6, "vertex quadratic spectrum drift")

    # Independent column solve of L phi=rho with sum(phi)=0.
    rho = [Fraction(0)] * PORTS; rho[0] = 1; rho[1] = -1
    system = vertex + [[1] * PORTS]
    phi = solve_augmented(system, rho + [Fraction(0)])
    field = [sum(Fraction(b[p][e]) * phi[p] for p in range(PORTS))
             for e in range(SEAMS)]
    recovered = [sum(Fraction(b[p][e]) * field[e] for e in range(SEAMS))
                 for p in range(PORTS)]
    check(recovered == rho, "scalar Green boundary drift")

    receipt = json.loads(receipt_path.read_bytes())
    body = {k: v for k, v in receipt.items() if k != "receipt_sha256"}
    check(receipt.get("receipt_sha256") ==
          "sha256:" + hashlib.sha256(canonical(body)).hexdigest(), "self digest drift")
    check(receipt_path.read_bytes() == canonical(receipt), "noncanonical JSON")
    check(receipt.get("schema") == "oph.local_face_maxwell_action.v1", "schema drift")
    check(receipt.get("issue") == 733, "issue drift")
    check(receipt.get("status") ==
          "EXACT_LOCAL_FACE_CURVATURE_AND_SCALAR_COULOMB_ACTION__PHYSICAL_ATTACHMENTS_PR53_PR54_OPEN",
          "status drift")
    expected_contracts = [{"path": rel, "bytes": size, "sha256": "sha256:" + sha,
                           "role": role} for rel, size, sha, role in SOURCE_ROWS]
    check(receipt.get("source_contracts") == expected_contracts, "contract drift")
    check(receipt.get("lean_receipts") == {
        "file": "Lean/Screen/LocalFaceMaxwellAction.lean",
        "theorems": list(THEOREMS), "sorry_free": True}, "Lean block drift")
    check(receipt.get("incidence") == {
        "ports": 12, "seams": 30, "faces": 20, "rank_B": 11, "rank_C": 19,
        "kernel_C_dimension": 11, "kernel_C_equals_image_B_transpose": True,
        "C_times_B_transpose_zero": True,
        "C_support": {"nonzero": 60, "total": 600}}, "incidence receipt drift")
    check(receipt.get("local_hessian") == {
        "definition": "H = C^T C", "rank": 19,
        "support": {"nonzero": 150, "total": 900, "nonzero_per_row": [5] * 30},
        "diagonal": [2] * 30}, "H receipt drift")
    check(receipt.get("exact_spectra") == {
        "C_C_transpose": {
            "rational_eigenspace_dimensions": {
                str(q): face_dims[q] for q in (0, 2, 3, 5)
            },
            "quadratic_polynomial": "x^2-6x+4",
            "quadratic_sector_dimension": face_quadratic_dim,
            "eigenvalues": ["0^1", "2^5", "3^4", "5^4",
                            "(3-sqrt(5))^3", "(3+sqrt(5))^3"],
        },
        "B_B_transpose": {
            "rational_eigenspace_dimensions": {
                str(q): vertex_dims[q] for q in (0, 6)
            },
            "quadratic_polynomial": "x^2-10x+20",
            "quadratic_sector_dimension": vertex_quadratic_dim,
            "eigenvalues": ["0^1", "6^5", "(5-sqrt(5))^3",
                            "(5+sqrt(5))^3"],
        },
    }, "spectral receipt drift")
    check(receipt.get("scalar_green_join") == {
        "charge_type": "Fin 12 -> rational",
        "field_type": "Fin 30 -> rational",
        "test_charge": {"0": "1", "1": "-1"},
        "green_potential": [str(x) for x in phi],
        "boundary_d_green_rho_equals_rho": True,
        "typed_separately_from_seam_current_source": True,
    }, "scalar Green join drift")
    check(receipt.get("adversarial_controls") == {
        "orientation_erasure_rejected": True,
        "unsigned_face_zero_port_boundary_value": unsigned_boundary[0][0],
        "dense_hessian_rewrite_rejected_by_exact_support": sum(hrows) < SEAMS ** 2,
    }, "adversarial-control receipt drift")
    check(receipt.get("handoff_interface") == {
        "schema": "oph.local_face_maxwell_action.handoff.v1",
        "consumer": "instrument lane (issue 737)",
        "design_only": True,
        "frozen": False,
        "statement": (
            "Expose separately typed scalar charge/potential and conserved "
            "seam-current/vector-potential observables; no physical "
            "identification or comparison is armed."
        ),
        "source_types": {
            "charge": "Fin 12 -> real",
            "seam_current": "Fin 30 -> real",
        },
        "comparison_permitted": False,
    }, "handoff drift")
    check(receipt.get("physical_boundary") == {
        "open_premise_rows": ["PR-53", "PR-54"],
        "position_attachment_proved": False,
        "physical_source_identification_proved": False,
        "laboratory_readout_proved": False,
        "continuum_limit_proved": False,
    }, "physical boundary drift")
    check(set(receipt) == {"schema", "status", "issue", "source_contracts",
          "lean_receipts", "incidence", "local_hessian", "exact_spectra",
          "scalar_green_join", "adversarial_controls", "handoff_interface",
          "physical_boundary", "receipt_sha256"}, "top-level field drift")
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    verify(args.receipt)
    print("LOCAL_FACE_MAXWELL_INDEPENDENT_VERIFICATION_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
