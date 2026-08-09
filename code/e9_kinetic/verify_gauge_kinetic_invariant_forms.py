#!/usr/bin/env python3
"""Independent replay for gauge_kinetic_invariant_forms.certificate.json.

This verifier intentionally does not import the producer.  It rebuilds the
Reynolds tensors, channel transform, A5 spectral projectors, compact sample
brackets, and carrier-weight ad-invariance equations with separate exact
Q(sqrt(5)) arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_CERT = HERE / "gauge_kinetic_invariant_forms.certificate.json"
PRODUCER = HERE / "gauge_kinetic_invariant_forms.py"
BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json"
COMPACT = REPO / "code/b14_jacobi/b14_compact_locus.certificate.json"

N = 12
M = 14
CHANNELS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m",
    "t_ff_to_p", "t_ff_to_m", "t_pm_to_f", "t_pf_to_p",
    "t_pf_to_m", "t_pf_to_f", "t_mf_to_p", "t_mf_to_m",
    "t_mf_to_f",
]
CH = {name: i for i, name in enumerate(CHANNELS)}
DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}


class VerificationError(RuntimeError):
    pass


def check(value: bool, message: str) -> None:
    if not value:
        raise VerificationError(message)


class F5:
    __slots__ = ("r", "s")

    def __init__(self, r=0, s=0):
        self.r = r if isinstance(r, Fraction) else Fraction(r)
        self.s = s if isinstance(s, Fraction) else Fraction(s)

    def __add__(self, x):
        x = x if isinstance(x, F5) else F5(x)
        return F5(self.r + x.r, self.s + x.s)

    __radd__ = __add__

    def __neg__(self):
        return F5(-self.r, -self.s)

    def __sub__(self, x):
        x = x if isinstance(x, F5) else F5(x)
        return self + (-x)

    def __mul__(self, x):
        x = x if isinstance(x, F5) else F5(x)
        return F5(self.r * x.r + 5 * self.s * x.s,
                  self.r * x.s + self.s * x.r)

    __rmul__ = __mul__

    def reciprocal(self):
        n = self.r * self.r - 5 * self.s * self.s
        check(n != 0, "zero divisor")
        return F5(self.r / n, -self.s / n)

    def __truediv__(self, x):
        x = x if isinstance(x, F5) else F5(x)
        return self * x.reciprocal()

    def __bool__(self):
        return bool(self.r or self.s)

    def __eq__(self, x):
        x = x if isinstance(x, F5) else F5(x)
        return self.r == x.r and self.s == x.s


Z = F5()
U = F5(1)
S5 = F5(0, 1)


def decode(x):
    return F5(Fraction(x[0], x[1]), Fraction(x[2], x[3]))


def encode(x):
    return [x.r.numerator, x.r.denominator, x.s.numerator, x.s.denominator]


def file_sha(path):
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def obj_sha(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def echelon(rows, width):
    rows = [list(row) for row in rows if any(row)]
    pivots = []
    at = 0
    for col in range(width):
        hit = next((i for i in range(at, len(rows)) if rows[i][col]), None)
        if hit is None:
            continue
        rows[at], rows[hit] = rows[hit], rows[at]
        rows[at] = [x * rows[at][col].reciprocal() for x in rows[at]]
        for i, row in enumerate(rows):
            if i != at and row[col]:
                c = row[col]
                rows[i] = [x - c * y for x, y in zip(row, rows[at])]
        pivots.append(col)
        at += 1
        if at == len(rows):
            break
    return rows[:at], pivots


def replay():
    raw_basis = json.loads(BASIS.read_text())
    raw_receipt = json.loads(RECEIPT.read_text())
    raw_system = json.loads(SYSTEM.read_text())
    raw_compact = json.loads(COMPACT.read_text())

    reynolds = [[[[Fraction(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in range(M)]
    for a, row in enumerate(raw_basis["basis"]):
        for o, i, j, p, q in row["entries"]:
            check(i < j, "noncanonical Reynolds entry")
            reynolds[a][o][i][j] = Fraction(p, q)
            reynolds[a][o][j][i] = Fraction(-p, q)

    inv = [[Z for _ in range(M)] for _ in range(M)]
    inv_rows = raw_system["fixed_line_reduction"]["channel_decomposition"]["inverse_transform_rows"]
    for i, row in enumerate(inv_rows):
        for j, an, ad, bn, bd in row:
            inv[i][j] = F5(Fraction(an, ad), Fraction(bn, bd))

    group = [tuple(row) for row in raw_receipt["proper_port_action"]["permutation_rows"]]
    check(len(set(group)) == 60, "proper group does not have order sixty")
    unvisited = {(i, j) for i in range(N) for j in range(N)}
    orbits = []
    while unvisited:
        seed = min(unvisited)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unvisited.difference_update(orbit)
        orbits.append(orbit)
    orbital = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))[0]
    adjacency = [[Z for _ in range(N)] for _ in range(N)]
    for i, j in orbital:
        adjacency[i][j] = U

    eigen = {"fixed": F5(5), "three_plus": S5, "three_minus": -S5, "five": F5(-1)}
    projectors = {}
    sector_bases = {}
    for label, ev in eigen.items():
        P = [[U if i == j else Z for j in range(N)] for i in range(N)]
        for other, ev2 in eigen.items():
            if label == other:
                continue
            factor = [[adjacency[i][j] - (ev2 if i == j else Z) for j in range(N)] for i in range(N)]
            P = [[sum((P[i][k] * factor[k][j] for k in range(N)), Z) for j in range(N)] for i in range(N)]
            P = [[x / (ev - ev2) for x in row] for row in P]
        projectors[label] = P
        sector_bases[label], _ = echelon(P, N)
        check(len(sector_bases[label]) == DIMS[label], f"bad {label} projector")

    def make_bracket(sample):
        y = {CH[name]: decode(value) for name, value in sample["channel_values"].items()}
        x = [sum((inv[i][j] * value for j, value in y.items()), Z) for i in range(M)]
        B = [[[Z for _ in range(N)] for _ in range(N)] for _ in range(N)]
        for a in range(M):
            for o in range(N):
                for i in range(N):
                    for j in range(N):
                        if x[a] and reynolds[a][o][i][j]:
                            B[o][i][j] += x[a] * reynolds[a][o][i][j]
        return B

    def bracket_value(B, x, y):
        return [sum((B[o][i][j] * x[i] * y[j] for i in range(N) for j in range(N)
                    if B[o][i][j] and x[i] and y[j]), Z) for o in range(N)]

    def qform(P, x, y):
        return sum((x[i] * P[i][j] * y[j] for i in range(N) for j in range(N)
                    if x[i] and P[i][j] and y[j]), Z)

    wanted = {
        "P": ("PLANE", ["1", "1"]),
        "F": ("SU3PF", ["1", "-1", "1"]),
        "G": ("SU3MF", ["1", "1", "1"]),
    }
    records = {}
    for short, (family, params) in wanted.items():
        sample = next(s for s in raw_compact["samples"] if s["family"] == family and s["parameters"] == params)
        names = sample["derived_sectors"]
        vectors = [v for name in names for v in sector_bases[name]]
        B = make_bracket(sample)
        equations = []
        for x in vectors:
            for u in vectors:
                xu = bracket_value(B, x, u)
                for v in vectors:
                    xv = bracket_value(B, x, v)
                    row = [qform(projectors[name], xu, v) + qform(projectors[name], u, xv) for name in names]
                    if any(row):
                        equations.append(row)
        rr, piv = echelon(equations, len(names))
        records[short] = {
            "certified_family": family,
            "certified_parameters": params,
            "weight_order": names,
            "carrier_invariant_dimension": len(names),
            "ad_invariance_equation_count_before_reduction": len(equations),
            "ad_invariance_constraint_rank": len(rr),
            "ad_invariant_dimension": len(names) - len(rr),
            "rref_constraints": [[encode(x) for x in row] for row in rr],
            "pivot_columns": piv,
        }
    return records, raw_compact


def verify_certificate(path: Path = DEFAULT_CERT) -> dict:
    cert = json.loads(path.read_text())
    check(cert.get("schema") == "oph.e9.gauge_kinetic_invariant_forms.v1", "schema mismatch")
    body = dict(cert)
    stored = body.pop("certificate_sha256", None)
    check(stored == obj_sha(body), "certificate self-hash mismatch")
    pins = cert["implementation_pins"]
    check(pins["producer_sha256"] == file_sha(PRODUCER), "producer pin mismatch")
    check(pins["verifier_sha256"] == file_sha(Path(__file__).resolve()), "verifier pin mismatch")
    upstream = cert["upstream"]
    check(upstream["stage1_basis_sha256"] == file_sha(BASIS), "basis pin mismatch")
    check(upstream["stage1_receipt_sha256"] == file_sha(RECEIPT), "receipt pin mismatch")
    check(upstream["stage2_system_sha256"] == file_sha(SYSTEM), "system pin mismatch")
    check(upstream["compact_locus_file_sha256"] == file_sha(COMPACT), "compact pin mismatch")
    replayed, compact = replay()
    check(cert["families"] == replayed, "family replay differs")
    check(upstream["compact_locus_self_sha256"] == compact["certificate_sha256"], "compact self pin mismatch")
    check(replayed["F"]["rref_constraints"] == [[[1, 1, 0, 1], [0, 1, 0, 1], [0, 1, -1, 1]]], "F relation mismatch")
    check(replayed["G"]["rref_constraints"] == [[[0, 1, 0, 1], [1, 1, 0, 1], [0, 1, -1, 1]]], "G relation mismatch")
    check(replayed["P"]["rref_constraints"] == [], "P control mismatch")
    common = cert["mirror_common_control"]
    check(common["combined_constraint_rank"] == 2 and common["solution_dimension"] == 1, "mirror-common dimension mismatch")
    check(common["generator"] == [encode(S5), encode(S5), encode(U)], "mirror-common generator mismatch")
    check(common["not_source_selected"] is True, "mirror-common boundary weakened")
    return {"verified": str(path), "dimensions": {k: replayed[k]["ad_invariant_dimension"] for k in replayed}}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERT)
    args = parser.parse_args()
    print(json.dumps(verify_certificate(args.certificate), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
