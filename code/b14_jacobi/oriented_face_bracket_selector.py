#!/usr/bin/env python3
"""Exact oriented-face bracket discriminator for issue B14 / GitHub #705.

The twenty coherently oriented faces in the pinned echosahedral source packet
define the most direct equal-weight alternating bracket: every face contributes
its three positive cyclic structure constants.  This program independently
locates that tensor in the pinned fourteen-dimensional Reynolds basis, computes
all Jacobi residual components, and projects it orthogonally onto the already
classified compact P/F/G loci.

BOUNDARY: the face tensor itself is source-derived.  Hilbert--Schmidt distance,
orthogonal projection, and any proposed Jacobi repair/minimization rule are not
source-derived by this calculation.  The nearest-locus result is therefore an
exact discriminator conditional on that explicitly added metric rule, not a
selection or closure theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json"
COMPACT = HERE / "b14_compact_locus.certificate.json"
OUTPUT = HERE / "oriented_face_bracket_selector.certificate.json"
VERIFIER = HERE / "verify_oriented_face_bracket_selector.py"

N = 12
M = 14
CHANNELS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m",
    "t_ff_to_p", "t_ff_to_m", "t_pm_to_f", "t_pf_to_p",
    "t_pf_to_m", "t_pf_to_f", "t_mf_to_p", "t_mf_to_m",
    "t_mf_to_f",
]
CH = {name: i for i, name in enumerate(CHANNELS)}


class CertificateError(RuntimeError):
    pass


def require(value: bool, message: str) -> None:
    if not value:
        raise CertificateError(message)


class Q5:
    """Exact element a + b sqrt(5)."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self + (-other)

    def __mul__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a * other.a + 5 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def reciprocal(self):
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "division by zero in Q(sqrt(5))")
        return Q5(self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self * other.reciprocal()

    def __eq__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self.a == other.a and self.b == other.b

    def __bool__(self):
        return bool(self.a or self.b)


ZERO = Q5()
ONE = Q5(1)
ROOT5 = Q5(0, 1)


def enc(x: Q5) -> list[int]:
    return [x.a.numerator, x.a.denominator, x.b.numerator, x.b.denominator]


def sign(x: Q5) -> int:
    """Exact real-embedding sign, using sqrt(5)>0 and rational squares."""
    if not x:
        return 0
    if x.a == 0:
        return 1 if x.b > 0 else -1
    if x.b == 0:
        return 1 if x.a > 0 else -1
    if (x.a > 0) == (x.b > 0):
        return 1 if x.a > 0 else -1
    comparison = x.a * x.a - 5 * x.b * x.b
    require(comparison != 0, "unexpected zero in Q(sqrt(5)) sign test")
    if x.a > 0:  # positive rational part, negative radical part
        return 1 if comparison > 0 else -1
    return 1 if comparison < 0 else -1


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def upper_index() -> list[tuple[int, int, int]]:
    return [(o, i, j) for o in range(N) for i in range(N) for j in range(i + 1, N)]


COORDS = upper_index()


def load_inputs():
    return tuple(json.loads(path.read_text()) for path in (MANIFEST, BASIS, RECEIPT, SYSTEM, COMPACT))


def canonical_faces(manifest) -> list[tuple[int, int, int]]:
    ports = manifest["carrier"]["ports"]
    require(ports == [f"p{i:02d}" for i in range(N)], "unexpected source port order")
    index = {port: i for i, port in enumerate(ports)}
    faces = [tuple(index[port] for port in face) for face in manifest["carrier"]["oriented_faces"]]
    require(len(faces) == 20 and len(set(faces)) == 20, "expected twenty distinct oriented faces")
    return faces


def face_tensor(faces):
    tensor = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a, b, c in faces:
        for out, left, right in ((c, a, b), (a, b, c), (b, c, a)):
            tensor[out][left][right] += 1
            tensor[out][right][left] -= 1
    return tensor


def reynolds_tensors(raw_basis):
    result = [[[[Fraction(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in range(M)]
    require([row["basis_id"] for row in raw_basis["basis"]] == [f"R{i:02d}" for i in range(M)],
            "unexpected Reynolds basis order")
    for a, row in enumerate(raw_basis["basis"]):
        for out, left, right, numerator, denominator in row["entries"]:
            require(left < right, "noncanonical Reynolds entry")
            value = Fraction(numerator, denominator)
            result[a][out][left][right] = value
            result[a][out][right][left] = -value
    return result


def q5_matrix(rows_key, system):
    rows = system["fixed_line_reduction"]["channel_decomposition"][rows_key]
    matrix = [[ZERO for _ in range(M)] for _ in range(M)]
    for i, row in enumerate(rows):
        for j, an, ad, bn, bd in row:
            matrix[i][j] = Q5(Fraction(an, ad), Fraction(bn, bd))
    return matrix


def vector_from_x(reynolds, x):
    return [sum((x[a] * reynolds[a][o][i][j] for a in range(M)), ZERO)
            for o, i, j in COORDS]


def x_from_channels(inverse, values: dict[str, Q5]):
    y = [ZERO for _ in range(M)]
    for name, value in values.items():
        y[CH[name]] = value
    return [sum((inverse[i][j] * y[j] for j in range(M)), ZERO) for i in range(M)]


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), ZERO)


def solve(matrix, rhs):
    n = len(rhs)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col]), None)
        require(pivot is not None, "singular projection Gram matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col].reciprocal()
        aug[col] = [scale * x for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col]:
                factor = aug[r][col]
                aug[r] = [x - factor * y for x, y in zip(aug[r], aug[col])]
    return [aug[i][-1] for i in range(n)]


def family_bases(inverse, reynolds):
    specs = {
        "P": [
            {"t_pp_to_p": ONE},
            {"t_mm_to_m": ONE},
        ],
        "F": [
            {"t_pp_to_p": ONE, "t_pf_to_f": -ROOT5},
            {"t_ff_to_p": ONE},
            {"t_mm_to_m": ONE},
        ],
        "G": [
            {"t_mm_to_m": ONE, "t_mf_to_f": ROOT5},
            {"t_ff_to_m": ONE},
            {"t_pp_to_p": ONE},
        ],
    }
    return {
        family: [vector_from_x(reynolds, x_from_channels(inverse, channels)) for channels in rows]
        for family, rows in specs.items()
    }


def project(target, basis):
    gram = [[dot(x, y) for y in basis] for x in basis]
    rhs = [dot(x, target) for x in basis]
    coefficients = solve(gram, rhs)
    residual = [target[i] - sum((coefficients[j] * basis[j][i] for j in range(len(basis))), ZERO)
                for i in range(len(target))]
    return gram, rhs, coefficients, dot(residual, residual)


def jacobi_census(tensor):
    entries = []
    witnesses = []
    for out in range(N):
        for i in range(N):
            for j in range(i + 1, N):
                for k in range(j + 1, N):
                    value = sum(
                        tensor[m][i][j] * tensor[out][m][k]
                        + tensor[m][j][k] * tensor[out][m][i]
                        + tensor[m][k][i] * tensor[out][m][j]
                        for m in range(N)
                    )
                    if value:
                        entries.append(value)
                        if len(witnesses) < 12:
                            witnesses.append([out, i, j, k, value])
    return {
        "coordinate_domain": "output o and strictly increasing input triple i<j<k",
        "domain_size": N * (N * (N - 1) * (N - 2) // 6),
        "nonzero_count": len(entries),
        "positive_count": sum(value > 0 for value in entries),
        "negative_count": sum(value < 0 for value in entries),
        "value_histogram": {str(value): entries.count(value) for value in sorted(set(entries))},
        "squared_norm": sum(value * value for value in entries),
        "first_witnesses": witnesses,
    }


def build_certificate():
    manifest, raw_basis, receipt, system, compact = load_inputs()
    require(compact["schema"] == "oph.b14_compact_locus.certificate.v1", "unexpected compact-locus schema")
    faces = canonical_faces(manifest)
    face = face_tensor(faces)
    reynolds = reynolds_tensors(raw_basis)
    target = [Q5(face[o][i][j]) for o, i, j in COORDS]

    exact_r13 = all(Q5(face[o][i][j]) == Q5(60 * reynolds[13][o][i][j])
                    for o in range(N) for i in range(N) for j in range(N))
    require(exact_r13, "canonical oriented-face bracket is not exactly 60 R13")
    require(dot(target, target) == Q5(60), "unexpected face-bracket Hilbert--Schmidt norm")

    transform = q5_matrix("transform_rows", system)
    inverse = q5_matrix("inverse_transform_rows", system)
    x = [ZERO for _ in range(M)]
    x[13] = Q5(60)
    channels = [sum((transform[i][j] * x[j] for j in range(M)), ZERO) for i in range(M)]
    expected_channels = [
        ZERO, ZERO, ZERO, ROOT5 / 20, -ROOT5 / 20,
        ROOT5 / 20, -ROOT5 / 20, ZERO, ZERO, ZERO,
        ROOT5 / 20, ZERO, ZERO, -ROOT5 / 20,
    ]
    require(channels == expected_channels, "unexpected oriented-face channel coordinates")

    jacobi = jacobi_census(face)
    require(jacobi["nonzero_count"] == 240, "unexpected Jacobi failure count")
    require(jacobi["positive_count"] == jacobi["negative_count"] == 120, "unexpected Jacobi sign census")
    require(jacobi["value_histogram"] == {"-1": 120, "1": 120}, "unexpected Jacobi values")
    require(jacobi["squared_norm"] == 240, "unexpected Jacobi residual norm")

    bases = family_bases(inverse, reynolds)
    names = {"P": ["p", "m"], "F": ["a", "b", "e"], "G": ["a", "b", "e"]}
    projections = {}
    for family in ("P", "F", "G"):
        gram, rhs, coefficients, residual = project(target, bases[family])
        projections[family] = {
            "parameter_order": names[family],
            "gram": [[enc(value) for value in row] for row in gram],
            "rhs": [enc(value) for value in rhs],
            "orthogonal_projection_parameters": [enc(value) for value in coefficients],
            "squared_distance": enc(residual),
        }
        if family in ("F", "G"):
            product = coefficients[0] * coefficients[1]
            projections[family]["compact_sign_product_a_times_b"] = enc(product)
            projections[family]["compact_sign"] = "negative" if sign(product) < 0 else "positive"
            projections[family]["projection_lies_in_compact_open_stratum"] = (
                sign(product) < 0 if family == "F" else sign(product) > 0
            )

    expected_parameters = {
        "P": [ROOT5 / 20, -ROOT5 / 20],
        "F": [Q5(Fraction(-1, 22), Fraction(1, 220)), ROOT5 / 20, -ROOT5 / 20],
        "G": [Q5(Fraction(-1, 22), Fraction(-1, 220)), -ROOT5 / 20, ROOT5 / 20],
    }
    expected_distances = {
        "P": Q5(45),
        "F": Q5(Fraction(615, 22), Fraction(123, 22)),
        "G": Q5(Fraction(615, 22), Fraction(-123, 22)),
    }
    for family in ("P", "F", "G"):
        require(projections[family]["orthogonal_projection_parameters"] == [enc(x) for x in expected_parameters[family]],
                f"unexpected {family} projection")
        require(projections[family]["squared_distance"] == enc(expected_distances[family]),
                f"unexpected {family} distance")
    require(projections["F"]["projection_lies_in_compact_open_stratum"], "F projection not compact")
    require(projections["G"]["projection_lies_in_compact_open_stratum"], "G projection not compact")
    require(sign(expected_distances["P"] - expected_distances["G"]) > 0, "G not closer than P")
    require(sign(expected_distances["F"] - expected_distances["G"]) > 0, "G not closer than F")

    group = receipt["proper_port_action"]["permutation_rows"]
    require(len(group) == 60 and len({tuple(g) for g in group}) == 60, "proper action order mismatch")

    body = {
        "schema": "oph.b14.oriented_face_bracket_selector.v1",
        "issue": 705,
        "claim_boundary": (
            "The equal-weight cyclic bracket is derived from the pinned oriented-face incidence and is exactly 60*R13. "
            "It fails Jacobi. The squared-distance comparison is conditional on adding the displayed coefficient-space "
            "Hilbert--Schmidt metric. Neither minimum-Hilbert--Schmidt repair nor any Jacobi-repair law is source-derived; "
            "therefore the unique nearest G result is a discriminator, not an OPH source-selection or B14 closure theorem."
        ),
        "field": "Q(sqrt(5)); [a_num,a_den,b_num,b_den] encodes a+b*sqrt(5)",
        "source_face_bracket": {
            "definition": "for every oriented face (a,b,c), add [e_a,e_b]=e_c, [e_b,e_c]=e_a, [e_c,e_a]=e_b and antisymmetry",
            "oriented_face_count": len(faces),
            "nonzero_upper_structure_constants": sum(bool(x) for x in target),
            "upper_coordinate_squared_norm": enc(dot(target, target)),
            "reynolds_coordinates": [[13, *enc(Q5(60))]],
            "identity": "B_face = 60 * R13 exactly in all 12*12*12 tensor coordinates",
            "proper_action_order": len(group),
        },
        "channel_coordinates": {name: enc(value) for name, value in zip(CHANNELS, channels)},
        "jacobi_failure": jacobi,
        "orthogonal_compact_locus_discriminator": {
            "metric_premise": "Euclidean/Hilbert--Schmidt inner product on upper-triangular structure-constant coordinates in the pinned port basis",
            "metric_is_source_derived": False,
            "minimum_hs_or_jacobi_repair_is_source_derived": False,
            "families": projections,
            "unique_nearest_family": "G",
            "exact_gaps": {
                "distance_F_minus_distance_G": enc(expected_distances["F"] - expected_distances["G"]),
                "distance_P_minus_distance_G": enc(expected_distances["P"] - expected_distances["G"]),
            },
            "conclusion_status": "TARGET_CLEAN_CONDITIONAL_DISCRIMINATOR__NOT_SOURCE_SELECTION",
        },
        "upstream": {
            "source_manifest_path": str(MANIFEST.relative_to(REPO)),
            "source_manifest_sha256": file_hash(MANIFEST),
            "stage1_basis_path": str(BASIS.relative_to(REPO)),
            "stage1_basis_sha256": file_hash(BASIS),
            "stage1_receipt_path": str(RECEIPT.relative_to(REPO)),
            "stage1_receipt_sha256": file_hash(RECEIPT),
            "stage2_system_path": str(SYSTEM.relative_to(REPO)),
            "stage2_system_sha256": file_hash(SYSTEM),
            "compact_locus_path": str(COMPACT.relative_to(REPO)),
            "compact_locus_file_sha256": file_hash(COMPACT),
            "compact_locus_self_sha256": compact["certificate_sha256"],
        },
        "implementation_pins": {
            "producer_sha256": file_hash(Path(__file__).resolve()),
            "verifier_sha256": file_hash(VERIFIER),
        },
    }
    body["certificate_sha256"] = object_hash(body)
    return body


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    if args.check:
        require(args.output.exists(), f"missing certificate {args.output}")
        require(json.loads(args.output.read_text()) == result, "committed certificate is stale")
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "certificate": str(args.output),
        "identity": result["source_face_bracket"]["identity"],
        "jacobi_nonzero": result["jacobi_failure"]["nonzero_count"],
        "unique_nearest": result["orthogonal_compact_locus_discriminator"]["unique_nearest_family"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
