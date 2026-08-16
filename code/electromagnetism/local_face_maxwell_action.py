#!/usr/bin/env python3
"""Exact receipt for the committed local face-curvature Maxwell packet.

This producer parses the committed seam and oriented-face tables, constructs
the port--seam incidence ``B`` and face--seam incidence ``C``, and performs
exact integer/rational checks.  It does not use floating-point eigensolvers.
The physical attachment remains open: "Maxwell" and "Coulomb" name the two
finite action sectors formalized in Lean, not a laboratory identification.
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
ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT = RUNTIME / "local_face_maxwell_action_receipt.json"

SOURCES = (
    (ROOT / "Lean/Screen/SeamCurrentCarrierQuotient.lean", 19162,
     "e5f712cbccc5a1462945f0cb511a47064fe10818748498021d5c70ae169a933f",
     "committed 30-seam table"),
    (ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean", 8196,
     "addcea08fe55c0ccc3edbeb90ebd9df4b640c96bd20e9755417bcdb84f093fe4",
     "committed 20 oriented faces"),
    (ROOT / "Lean/Screen/DiscreteCoulombGreen.lean", 50954,
     "a58830d60688b7df505d2c8cbc1ff2290b057f9be055def25db7006fae7bd078",
     "exact scalar Green operator"),
    (ROOT / "Lean/Screen/LocalFaceMaxwellAction.lean", 32766,
     "c54da8884d3a3d19ed73baab2218378795fa5f8f2681344d545f285a647cf763",
     "kernel-checked local action and typed join"),
)

SCHEMA = "oph.local_face_maxwell_action.v1"
STATUS = (
    "EXACT_LOCAL_FACE_CURVATURE_AND_SCALAR_COULOMB_ACTION__"
    "PHYSICAL_ATTACHMENTS_PR53_PR54_OPEN"
)
PORTS, SEAMS, FACES = 12, 30, 20

THEOREMS = (
    "faceVertices_matches_committed",
    "face_port_incidence_product_zero",
    "unsigned_face_incidence_breaks_boundary",
    "localKineticZ_five_per_row",
    "localKineticZ_total_support",
    "localKineticZ_diagonal_two",
    "ker_faceCurvature_eq_gradient",
    "range_localMaxwellOperator",
    "localSourcedAction_gauge_invariant_iff",
    "localStationary_solvable_iff",
    "localSourced_global_minimum",
    "greenPotential_stationary",
    "scalarCoulombAction_constant_gauge_iff",
    "scalarCoulomb_global_minimum",
    "staticAction_gauge_invariant",
    "staticAction_global_minimum",
)


class ReceiptError(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ReceiptError(message)


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_sources() -> dict[str, str]:
    out = {}
    for path, size, sha, _ in SOURCES:
        raw = path.read_bytes()
        require(len(raw) == size, f"source byte pin drift: {path.name}")
        require(digest(raw) == sha, f"source digest pin drift: {path.name}")
        text = raw.decode()
        require("sorry" not in text and "admit" not in text,
                f"placeholder in {path.name}")
        out[path.name] = text
    return out


def parse_vector(text: str, name: str, count: int) -> list[int]:
    m = re.search(rf"def {name} : Fin {count} → Fin 12 :=\s*!\[([^\]]*)\]",
                  text, re.S)
    require(m is not None, f"missing {name}")
    values = [int(x) for x in m.group(1).replace("\n", " ").split(",")]
    require(len(values) == count, f"{name} arity drift")
    return values


def parse_faces(text: str) -> list[tuple[int, int, int]]:
    start = text.index("def orientedFaces")
    end = text.index("def faceEdges", start)
    faces = [tuple(map(int, x)) for x in
             re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", text[start:end])]
    require(len(faces) == FACES, "oriented face count drift")
    return faces


def transpose(a: list[list[int | Fraction]]) -> list[list[int | Fraction]]:
    return [list(row) for row in zip(*a, strict=True)]


def mul(a: list[list[int | Fraction]], b: list[list[int | Fraction]]):
    bt = transpose(b)
    return [[sum((x * y for x, y in zip(row, col, strict=True)), start=0)
             for col in bt] for row in a]


def identity(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def sub_scalar_identity(a: list[list[int]], scalar: int) -> list[list[int]]:
    return [[a[i][j] - scalar * int(i == j) for j in range(len(a))]
            for i in range(len(a))]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def scale(a, k: int):
    return [[k * x for x in row] for row in a]


def rank(a: list[list[int | Fraction]]) -> int:
    work = [[Fraction(x) for x in row] for row in a]
    r = 0
    for c in range(len(work[0])):
        pivot = next((i for i in range(r, len(work)) if work[i][c]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        q = work[r][c]
        work[r] = [x / q for x in work[r]]
        for i in range(len(work)):
            if i != r and work[i][c]:
                q = work[i][c]
                work[i] = [x - q * y for x, y in zip(work[i], work[r], strict=True)]
        r += 1
    return r


def nullity(a) -> int:
    return len(a[0]) - rank(a)


def inverse(a: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(a)
    w = [row[:] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for c in range(n):
        p = next((i for i in range(c, n) if w[i][c]), None)
        require(p is not None, "singular regularized Laplacian")
        w[c], w[p] = w[p], w[c]
        q = w[c][c]
        w[c] = [x / q for x in w[c]]
        for i in range(n):
            if i != c and w[i][c]:
                q = w[i][c]
                w[i] = [x - q * y for x, y in zip(w[i], w[c], strict=True)]
    return [row[n:] for row in w]


def build_incidence(left: list[int], right: list[int], faces):
    b = [[0] * SEAMS for _ in range(PORTS)]
    edge_index = {}
    for e, (l, r) in enumerate(zip(left, right, strict=True)):
        require(l < r, "seam orientation drift")
        require((l, r) not in edge_index, "duplicate seam")
        edge_index[l, r] = e
        b[l][e], b[r][e] = -1, 1
    c = [[0] * SEAMS for _ in range(FACES)]
    for f, (x, y, z) in enumerate(faces):
        for u, v in ((x, y), (y, z), (z, x)):
            key = (min(u, v), max(u, v))
            require(key in edge_index, "face edge absent from seam table")
            c[f][edge_index[key]] = 1 if u < v else -1
    return b, c


def spectrum_receipt(matrix: list[list[int]], rational_parts,
                     quadratic: tuple[int, int], quadratic_dim: int):
    checks = {}
    for eigenvalue, multiplicity in rational_parts:
        got = nullity(sub_scalar_identity(matrix, eigenvalue))
        require(got == multiplicity, f"eigenspace drift at {eigenvalue}")
        checks[str(eigenvalue)] = got
    # q(M)=M^2+aM+bI for q(x)=x^2+a x+b.
    a, b = quadratic
    qmat = add(add(mul(matrix, matrix), scale(matrix, a)), scale(identity(len(matrix)), b))
    require(nullity(qmat) == quadratic_dim, "quadratic spectral sector drift")
    return {"rational_eigenspace_dimensions": checks,
            "quadratic_polynomial": f"x^2{a:+d}x{b:+d}",
            "quadratic_sector_dimension": quadratic_dim}


def scalar_green_join(b: list[list[int]]) -> dict[str, Any]:
    lap = mul(b, transpose(b))
    regular = [[Fraction(lap[i][j] + 1) for j in range(PORTS)]
               for i in range(PORTS)]
    inv = inverse(regular)
    green = [[inv[i][j] - Fraction(1, 144) for j in range(PORTS)]
             for i in range(PORTS)]
    rho = [Fraction(0)] * PORTS
    rho[0], rho[1] = 1, -1
    phi = [sum((green[p][q] * rho[q] for q in range(PORTS)), Fraction(0))
           for p in range(PORTS)]
    field = [sum((Fraction(b[p][e]) * phi[p] for p in range(PORTS)), Fraction(0))
             for e in range(SEAMS)]
    recovered = [sum((Fraction(b[p][e]) * field[e] for e in range(SEAMS)), Fraction(0))
                 for p in range(PORTS)]
    require(recovered == rho, "scalar Green join boundary drift")
    return {"charge_type": "Fin 12 -> rational", "field_type": "Fin 30 -> rational",
            "test_charge": {"0": "1", "1": "-1"},
            "green_potential": [str(x) for x in phi],
            "boundary_d_green_rho_equals_rho": True,
            "typed_separately_from_seam_current_source": True}


def build_receipt() -> dict[str, Any]:
    texts = load_sources()
    left = parse_vector(texts["SeamCurrentCarrierQuotient.lean"], "seamLeft", SEAMS)
    right = parse_vector(texts["SeamCurrentCarrierQuotient.lean"], "seamRight", SEAMS)
    faces = parse_faces(texts["CoreAxioms.lean"])
    local_text = texts["LocalFaceMaxwellAction.lean"]
    require(all(f"theorem {name}" in local_text for name in THEOREMS),
            "Lean theorem inventory drift")
    b, c = build_incidence(left, right, faces)
    h = mul(transpose(c), c)
    dual = mul(c, transpose(c))
    vertex = mul(b, transpose(b))
    require(mul(c, transpose(b)) == [[0] * PORTS for _ in range(FACES)],
            "C B^T is nonzero")
    require(rank(b) == 11 and rank(c) == 19, "incidence rank drift")
    require(nullity(c) == 11, "curvature kernel dimension drift")
    c_support = sum(x != 0 for row in c for x in row)
    h_rows = [sum(x != 0 for x in row) for row in h]
    require(c_support == 60, "C support drift")
    require(h_rows == [5] * SEAMS, "H row support drift")
    require(sum(h_rows) == 150, "H support drift")
    require([h[e][e] for e in range(SEAMS)] == [2] * SEAMS,
            "H diagonal drift")
    face_spectrum = spectrum_receipt(
        dual, ((0, 1), (2, 5), (3, 4), (5, 4)), (-6, 4), 6)
    vertex_spectrum = spectrum_receipt(
        vertex, ((0, 1), (6, 5)), (-10, 20), 6)
    value = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 733,
        "source_contracts": [
            {"path": str(path.relative_to(ROOT)), "bytes": size,
             "sha256": "sha256:" + sha, "role": role}
            for path, size, sha, role in SOURCES],
        "lean_receipts": {"file": "Lean/Screen/LocalFaceMaxwellAction.lean",
                          "theorems": list(THEOREMS), "sorry_free": True},
        "incidence": {
            "ports": PORTS, "seams": SEAMS, "faces": FACES,
            "rank_B": 11, "rank_C": 19,
            "kernel_C_dimension": 11,
            "kernel_C_equals_image_B_transpose": True,
            "C_times_B_transpose_zero": True,
            "C_support": {"nonzero": 60, "total": 600},
        },
        "local_hessian": {
            "definition": "H = C^T C", "rank": 19,
            "support": {"nonzero": 150, "total": 900,
                        "nonzero_per_row": h_rows},
            "diagonal": [2] * SEAMS,
        },
        "exact_spectra": {
            "C_C_transpose": {**face_spectrum,
                "eigenvalues": ["0^1", "2^5", "3^4", "5^4",
                                "(3-sqrt(5))^3", "(3+sqrt(5))^3"]},
            "B_B_transpose": {**vertex_spectrum,
                "eigenvalues": ["0^1", "6^5", "(5-sqrt(5))^3",
                                "(5+sqrt(5))^3"]},
        },
        "scalar_green_join": scalar_green_join(b),
        "adversarial_controls": {
            "orientation_erasure_rejected": True,
            "unsigned_face_zero_port_boundary_value": -2,
            "dense_hessian_rewrite_rejected_by_exact_support": True,
        },
        "handoff_interface": {
            "schema": "oph.local_face_maxwell_action.handoff.v1",
            "consumer": "instrument lane (issue 737)",
            "design_only": True, "frozen": False,
            "statement": "Expose separately typed scalar charge/potential and conserved seam-current/vector-potential observables; no physical identification or comparison is armed.",
            "source_types": {"charge": "Fin 12 -> real",
                             "seam_current": "Fin 30 -> real"},
            "comparison_permitted": False,
        },
        "physical_boundary": {
            "open_premise_rows": ["PR-53", "PR-54"],
            "position_attachment_proved": False,
            "physical_source_identification_proved": False,
            "laboratory_readout_proved": False,
            "continuum_limit_proved": False,
        },
    }
    value["receipt_sha256"] = "sha256:" + digest(canonical(value))
    return value


def verify_committed() -> dict[str, Any]:
    committed = json.loads(RECEIPT.read_bytes())
    expected = build_receipt()
    require(committed == expected, "committed receipt semantic drift")
    require(RECEIPT.read_bytes() == canonical(expected), "noncanonical receipt bytes")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        RUNTIME.mkdir(parents=True, exist_ok=True)
        RECEIPT.write_bytes(canonical(build_receipt()))
    if args.verify:
        verify_committed()
    if not args.write and not args.verify:
        print(canonical(build_receipt()).decode(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
