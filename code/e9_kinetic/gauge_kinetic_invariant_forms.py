#!/usr/bin/env python3
"""Exact E9 item 1a certificate: carrier-invariant versus ad-invariant forms.

The computation starts from the pinned fourteen-dimensional Reynolds basis,
the pinned channel transform, and the certified compact samples.  It does not
read the existing E9 kinetic-form output.  For each sample it reconstructs the
bracket over Q(sqrt(5)), the A5 spectral projectors, and the linear
ad-invariance equations on the carrier-projector form space.

The result is deliberately representation-level.  In particular, requiring
one form to be invariant for both mirror brackets is an *additional
simultaneous-invariance premise*, not a source-selection theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASIS_PATH = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT_PATH = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM_PATH = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json"
COMPACT_PATH = REPO / "code/b14_jacobi/b14_compact_locus.certificate.json"
OUTPUT_PATH = HERE / "gauge_kinetic_invariant_forms.certificate.json"
VERIFIER_PATH = HERE / "verify_gauge_kinetic_invariant_forms.py"

PORTS = 12
PARAMS = 14
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
CHANNELS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m",
    "t_ff_to_p", "t_ff_to_m", "t_pm_to_f", "t_pf_to_p",
    "t_pf_to_m", "t_pf_to_f", "t_mf_to_p", "t_mf_to_m",
    "t_mf_to_f",
]
CIDX = {name: i for i, name in enumerate(CHANNELS)}


class CertificateError(RuntimeError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise CertificateError(message)


class Q5:
    """Exact a + b sqrt(5)."""

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
        return self + (-other if isinstance(other, Q5) else -Q5(other))

    def __mul__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a * other.a + 5 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def inv(self):
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "division by zero in Q(sqrt(5))")
        return Q5(self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self * other.inv()

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self.a == other.a and self.b == other.b


ZERO = Q5()
ONE = Q5(1)
ROOT5 = Q5(0, 1)


def enc(x: Q5) -> list[int]:
    return [x.a.numerator, x.a.denominator, x.b.numerator, x.b.denominator]


def dec(x: list[int]) -> Q5:
    return Q5(Fraction(x[0], x[1]), Fraction(x[2], x[3]))


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def rref(rows: list[list[Q5]], width: int) -> tuple[list[list[Q5]], list[int]]:
    work = [list(row) for row in rows if any(row)]
    pivots: list[int] = []
    pivot_row = 0
    for col in range(width):
        hit = next((i for i in range(pivot_row, len(work)) if work[i][col]), None)
        if hit is None:
            continue
        work[pivot_row], work[hit] = work[hit], work[pivot_row]
        scale = work[pivot_row][col].inv()
        work[pivot_row] = [scale * x for x in work[pivot_row]]
        for i in range(len(work)):
            if i != pivot_row and work[i][col]:
                factor = work[i][col]
                work[i] = [x - factor * y for x, y in zip(work[i], work[pivot_row])]
        pivots.append(col)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return work[:pivot_row], pivots


def load_inputs():
    basis = json.loads(BASIS_PATH.read_text())
    receipt = json.loads(RECEIPT_PATH.read_text())
    system = json.loads(SYSTEM_PATH.read_text())
    compact = json.loads(COMPACT_PATH.read_text())
    require(compact["upstream"]["stage1_basis_file_sha256"] == file_hash(BASIS_PATH), "compact certificate basis pin differs")
    require(compact["upstream"]["stage1_receipt_file_sha256"] == file_hash(RECEIPT_PATH), "compact certificate receipt pin differs")
    require(compact["upstream"]["stage2_system_file_sha256"] == file_hash(SYSTEM_PATH), "compact certificate system pin differs")
    return basis, receipt, system, compact


def reynolds_tensors(basis):
    tensors = [[[[Fraction(0) for _ in range(PORTS)] for _ in range(PORTS)]
                for _ in range(PORTS)] for _ in range(PARAMS)]
    require(len(basis["basis"]) == PARAMS, "expected fourteen Reynolds rows")
    for a, row in enumerate(basis["basis"]):
        for out, left, right, num, den in row["entries"]:
            require(left < right, "Reynolds entry is not upper triangular")
            value = Fraction(num, den)
            tensors[a][out][left][right] = value
            tensors[a][out][right][left] = -value
    return tensors


def inverse_transform(system):
    rows = system["fixed_line_reduction"]["channel_decomposition"]["inverse_transform_rows"]
    inverse = [[ZERO for _ in range(PARAMS)] for _ in range(PARAMS)]
    for i, row in enumerate(rows):
        for j, an, ad, bn, bd in row:
            inverse[i][j] = Q5(Fraction(an, ad), Fraction(bn, bd))
    return inverse


def bracket(reynolds, inverse, channels: dict[int, Q5]):
    x = [sum((inverse[i][j] * value for j, value in channels.items()), ZERO)
         for i in range(PARAMS)]
    result = [[[ZERO for _ in range(PORTS)] for _ in range(PORTS)] for _ in range(PORTS)]
    for a in range(PARAMS):
        if not x[a]:
            continue
        for out in range(PORTS):
            for left in range(PORTS):
                for right in range(PORTS):
                    if reynolds[a][out][left][right]:
                        result[out][left][right] += x[a] * reynolds[a][out][left][right]
    return result


def spectral_projectors(group):
    remaining = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        remaining -= orbit
        orbits.append(orbit)
    sixty = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))
    require(len(sixty) == 2, "unexpected ordered-pair orbit structure")
    adjacency = [[ZERO for _ in range(PORTS)] for _ in range(PORTS)]
    for i, j in sixty[0]:
        adjacency[i][j] = ONE
    eigen = {"fixed": Q5(5), "three_plus": ROOT5,
             "three_minus": -ROOT5, "five": Q5(-1)}
    projectors = {}
    for name, value in eigen.items():
        matrix = [[ONE if i == j else ZERO for j in range(PORTS)] for i in range(PORTS)]
        for other, other_value in eigen.items():
            if other == name:
                continue
            factor = [[adjacency[i][j] - (other_value if i == j else ZERO)
                       for j in range(PORTS)] for i in range(PORTS)]
            matrix = [[sum((matrix[i][k] * factor[k][j] for k in range(PORTS)), ZERO)
                       for j in range(PORTS)] for i in range(PORTS)]
            divisor = value - other_value
            matrix = [[entry / divisor for entry in row] for row in matrix]
        projectors[name] = matrix
        trace = sum((matrix[i][i] for i in range(PORTS)), ZERO)
        require(trace == Q5(DIMS[name]), f"wrong projector rank for {name}")
    return projectors


def sector_bases(projectors):
    result = {}
    for name in SECTORS:
        rows, _ = rref([list(row) for row in projectors[name]], PORTS)
        require(len(rows) == DIMS[name], f"wrong basis rank for {name}")
        result[name] = rows
    return result


def lie(B, x, y):
    return [sum((B[o][i][j] * x[i] * y[j]
                 for i in range(PORTS) for j in range(PORTS)
                 if B[o][i][j] and x[i] and y[j]), ZERO)
            for o in range(PORTS)]


def form(projector, x, y):
    py = [sum((projector[i][j] * y[j] for j in range(PORTS) if y[j]), ZERO)
          for i in range(PORTS)]
    return sum((x[i] * py[i] for i in range(PORTS) if x[i]), ZERO)


def constraint_space(B, names, projectors, bases):
    vectors = [vector for name in names for vector in bases[name]]
    equations = []
    for x in vectors:
        for u in vectors:
            xu = lie(B, x, u)
            for v in vectors:
                xv = lie(B, x, v)
                row = [form(projectors[name], xu, v) + form(projectors[name], u, xv)
                       for name in names]
                if any(row):
                    equations.append(row)
    reduced, pivots = rref(equations, len(names))
    return reduced, pivots, len(equations)


SAMPLES = {
    "P": ("PLANE", ["1", "1"]),
    "F": ("SU3PF", ["1", "-1", "1"]),
    "G": ("SU3MF", ["1", "1", "1"]),
}


def build_certificate() -> dict:
    basis, receipt, system, compact = load_inputs()
    reynolds = reynolds_tensors(basis)
    inverse = inverse_transform(system)
    group = [tuple(row) for row in receipt["proper_port_action"]["permutation_rows"]]
    require(len(set(group)) == 60, "expected sixty proper carrier maps")
    projectors = spectral_projectors(group)
    bases = sector_bases(projectors)
    family_records = {}
    for name, (family, parameters) in SAMPLES.items():
        sample = next(s for s in compact["samples"]
                      if s["family"] == family and s["parameters"] == parameters)
        values = {CIDX[channel]: dec(value) for channel, value in sample["channel_values"].items()}
        B = bracket(reynolds, inverse, values)
        sectors = list(sample["derived_sectors"])
        reduced, pivots, raw_count = constraint_space(B, sectors, projectors, bases)
        family_records[name] = {
            "certified_family": family,
            "certified_parameters": parameters,
            "weight_order": sectors,
            "carrier_invariant_dimension": len(sectors),
            "ad_invariance_equation_count_before_reduction": raw_count,
            "ad_invariance_constraint_rank": len(reduced),
            "ad_invariant_dimension": len(sectors) - len(reduced),
            "rref_constraints": [[enc(value) for value in row] for row in reduced],
            "pivot_columns": pivots,
        }

    require(family_records["F"]["rref_constraints"] == [[[1, 1, 0, 1], [0, 1, 0, 1], [0, 1, -1, 1]]],
            "F relation is not w_three_plus - sqrt(5) w_five = 0")
    require(family_records["G"]["rref_constraints"] == [[[0, 1, 0, 1], [1, 1, 0, 1], [0, 1, -1, 1]]],
            "G relation is not w_three_minus - sqrt(5) w_five = 0")
    require(family_records["P"]["rref_constraints"] == [], "P should impose no carrier-weight relation")

    body = {
        "schema": "oph.e9.gauge_kinetic_invariant_forms.v1",
        "issue": 716,
        "claim_boundary": (
            "On the supplied certified generic brackets, this exact finite record compares the A5-carrier-invariant "
            "projector form space with the ad-invariant subspace. It selects no bracket, coupling scale, relative "
            "factor coefficient, source action, continuum field, or laboratory current. Simultaneous F/G invariance "
            "is an extra mirror-common premise and is not inferred from the carrier or source."
        ),
        "field": "Q(sqrt(5)); [a_num,a_den,b_num,b_den] encodes a+b*sqrt(5)",
        "families": family_records,
        "derived_relations": {
            "F": "w_three_plus = sqrt(5) * w_five; w_three_minus is free",
            "G": "w_three_minus = sqrt(5) * w_five; w_three_plus is free",
            "P": "two carrier weights remain two ad-invariant weights",
        },
        "mirror_common_control": {
            "extra_premise": "one and the same carrier-projector form is ad-invariant for both supplied mirror brackets F and G",
            "weight_order": ["three_plus", "three_minus", "five"],
            "combined_constraint_rank": 2,
            "solution_dimension": 1,
            "generator": [enc(ROOT5), enc(ROOT5), enc(ONE)],
            "relation": "(w_three_plus,w_three_minus,w_five)=t*(sqrt(5),sqrt(5),1)",
            "not_source_selected": True,
        },
        "upstream": {
            "stage1_basis_path": str(BASIS_PATH.relative_to(REPO)),
            "stage1_basis_sha256": file_hash(BASIS_PATH),
            "stage1_receipt_path": str(RECEIPT_PATH.relative_to(REPO)),
            "stage1_receipt_sha256": file_hash(RECEIPT_PATH),
            "stage2_system_path": str(SYSTEM_PATH.relative_to(REPO)),
            "stage2_system_sha256": file_hash(SYSTEM_PATH),
            "compact_locus_path": str(COMPACT_PATH.relative_to(REPO)),
            "compact_locus_file_sha256": file_hash(COMPACT_PATH),
            "compact_locus_self_sha256": compact["certificate_sha256"],
        },
        "implementation_pins": {
            "producer_sha256": file_hash(Path(__file__).resolve()),
            "verifier_sha256": file_hash(VERIFIER_PATH),
        },
    }
    body["certificate_sha256"] = object_hash(body)
    return body


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    if args.check:
        require(args.output.exists(), f"missing certificate {args.output}")
        require(json.loads(args.output.read_text()) == result, "checked certificate is stale")
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"certificate": str(args.output), "families": {
        name: [record["carrier_invariant_dimension"], record["ad_invariant_dimension"]]
        for name, record in result["families"].items()}, "mirror_common_dimension": 1}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
