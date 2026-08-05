#!/usr/bin/env python3
"""Independent replay of the issue #705 compact-real-locus certificate.

This verifier reads only the pinned stage-one basis, the pinned stage-two
reduction, and the emitted certificate. It recomputes, with its own exact
arithmetic over Q(sqrt(5)),

* the 2640 coordinate Jacobi equations, their 44 signed-orbit representatives,
  the coefficient-row rank 38, and the exact 11+27 split;
* the spectral identities behind the compact-implies-derivation-free exclusion;
* the restricted quadric spans for every analyzed channel support, the claimed
  span bases, and the coexistence facts;
* every family Killing monomial matrix, sector block, and signature;
* every sample point: Jacobi vanishing evaluated directly on the bracket
  tensor, the Killing matrix, its signature, the derived algebra, the center,
  and the compact-type criterion.

The verifier fails closed: any mismatch raises VerifyError.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STAGE1_BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
STAGE1_RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
STAGE2_SYSTEM = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json"
STAGE2_RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json"
CERTIFICATE_PATH = HERE / "b14_compact_locus.certificate.json"
PRODUCER_PATH = HERE / "classify.py"

PORTS = 12
PARAMS = 14
MONOMIALS = [(a, b) for a in range(PARAMS) for b in range(a, PARAMS)]
MIDX = {m: i for i, m in enumerate(MONOMIALS)}
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
CHANNEL_IDS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m", "t_ff_to_p",
    "t_ff_to_m", "t_pm_to_f", "t_pf_to_p", "t_pf_to_m", "t_pf_to_f",
    "t_mf_to_p", "t_mf_to_m", "t_mf_to_f",
]
CHANNEL_INDEX = {cid: i for i, cid in enumerate(CHANNEL_IDS)}


class VerifyError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def check(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise VerifyError(code, message)


class F5:
    """Exact a + b*sqrt(5); replay arithmetic written for this verifier."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return F5(-self.a, -self.b)

    def __sub__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - 5 * self.b * self.b
        check(n != 0, "F5_DIV", "division by zero")
        return F5(self.a / n, -self.b / n)

    def __truediv__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return self * o.inv()

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def sign(self):
        if not self:
            return 0
        if self.b == 0:
            return 1 if self.a > 0 else -1
        if self.a == 0:
            return 1 if self.b > 0 else -1
        if self.a * self.a > 5 * self.b * self.b:
            return 1 if self.a > 0 else -1
        return 1 if self.b > 0 else -1


ZERO = F5()
ONE = F5(1)
RT5 = F5(0, 1)


def decode_f5(entry):
    return F5(Fraction(entry[0], entry[1]), Fraction(entry[2], entry[3]))


def sha_file(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def sha_object(value) -> str:
    return "sha256:" + hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    ).hexdigest()


def echelon(rows):
    work = [list(r) for r in rows if any(r)]
    if not work:
        return [], []
    pivots = []
    pr = 0
    for col in range(len(work[0])):
        p = next((r for r in range(pr, len(work)) if work[r][col]), None)
        if p is None:
            continue
        work[pr], work[p] = work[p], work[pr]
        inv = work[pr][col].inv() if isinstance(work[pr][col], F5) else 1 / work[pr][col]
        work[pr] = [v * inv for v in work[pr]]
        for r in range(len(work)):
            if r != pr and work[r][col]:
                f = work[r][col]
                work[r] = [x - f * y for x, y in zip(work[r], work[pr])]
        pivots.append(col)
        pr += 1
        if pr == len(work):
            break
    return work[:pr], pivots


def member(basis, pivots, row) -> bool:
    residual = list(row)
    for base, col in zip(basis, pivots):
        if residual[col]:
            f = residual[col]
            residual = [x - f * y for x, y in zip(residual, base)]
    return not any(residual)


def inertia(M):
    n = len(M)
    W = [row[:] for row in M]
    pos = neg = 0
    size = n
    while size:
        d = next((i for i in range(size) if W[i][i]), None)
        if d is None:
            hit = None
            for i in range(size):
                for j in range(i + 1, size):
                    if W[i][j]:
                        hit = (i, j)
                        break
                if hit:
                    break
            if hit is None:
                break
            i, j = hit
            for k in range(size):
                W[i][k] = W[i][k] + W[j][k]
            for k in range(size):
                W[k][i] = W[k][i] + W[k][j]
            d = i
        if d:
            W[0], W[d] = W[d], W[0]
            for row in W:
                row[0], row[d] = row[d], row[0]
        piv = W[0][0]
        if piv.sign() > 0:
            pos += 1
        else:
            neg += 1
        W = [[W[i][j] - (W[i][0] / piv) * W[0][j] for j in range(1, size)] for i in range(1, size)]
        size -= 1
    return [pos, neg, n - pos - neg]


def parity(vals):
    return sum(vals[i] > vals[j] for i in range(len(vals)) for j in range(i + 1, len(vals))) % 2


def main() -> dict:
    basis_raw = json.loads(STAGE1_BASIS.read_text(encoding="utf-8"))
    stage1_receipt = json.loads(STAGE1_RECEIPT.read_text(encoding="utf-8"))
    system = json.loads(STAGE2_SYSTEM.read_text(encoding="utf-8"))
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))

    # pins and self-hash
    up = certificate["upstream"]
    check(up["stage1_basis_file_sha256"] == sha_file(STAGE1_BASIS), "PIN1", "stage-one basis hash differs")
    check(up["stage1_receipt_file_sha256"] == sha_file(STAGE1_RECEIPT), "PIN2", "stage-one receipt hash differs")
    check(up["stage2_system_file_sha256"] == sha_file(STAGE2_SYSTEM), "PIN3", "stage-two system hash differs")
    check(up["stage2_receipt_file_sha256"] == sha_file(STAGE2_RECEIPT), "PIN4", "stage-two receipt hash differs")
    body = dict(certificate)
    stored_hash = body.pop("certificate_sha256")
    check(stored_hash == sha_object(body), "SELF_HASH", "certificate self-hash differs")
    pins = certificate["implementation_pins"]
    check(pins["producer_sha256"] == sha_file(PRODUCER_PATH), "PIN_PRODUCER", "producer hash differs")
    check(pins["verifier_sha256"] == sha_file(Path(__file__).resolve()), "PIN_VERIFIER", "verifier hash differs")

    # structure constants
    C = [[[[Fraction(0)] * PORTS for _ in range(PORTS)] for _ in range(PORTS)] for _ in range(PARAMS)]
    check(len(basis_raw["basis"]) == PARAMS, "BASIS14", "expected fourteen basis rows")
    for a, row in enumerate(basis_raw["basis"]):
        for o, l, r, n, d in row["entries"]:
            C[a][o][l][r] = Fraction(n, d)
            C[a][o][r][l] = -Fraction(n, d)
    group = [tuple(r) for r in stage1_receipt["proper_port_action"]["permutation_rows"]]
    check(len(set(group)) == 60, "GROUP60", "expected sixty group rows")

    # 2640 Jacobi rows
    pair_support: dict = {}
    col_index: dict = {}
    for a in range(PARAMS):
        for o in range(PORTS):
            for l in range(PORTS):
                for r in range(l + 1, PORTS):
                    if C[a][o][l][r]:
                        pair_support.setdefault((l, r), []).append((a, o, C[a][o][l][r]))
    for b in range(PARAMS):
        for o in range(PORTS):
            for m in range(PORTS):
                for k in range(PORTS):
                    if C[b][o][m][k]:
                        col_index.setdefault((m, k), []).append((b, o, C[b][o][m][k]))

    def pair_val(i, j):
        if i < j:
            return pair_support.get((i, j), [])
        return [(a, m, -v) for (a, m, v) in pair_support.get((j, i), [])]

    raw = {}
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        rows = [[Fraction(0)] * len(MONOMIALS) for _ in range(PORTS)]
        for (pi, pj, pk) in ((i, j, k), (j, k, i), (k, i, j)):
            for (a, m, v) in pair_val(pi, pj):
                for (b, o, w) in col_index.get((m, pk), []):
                    rows[o][MIDX[tuple(sorted((a, b)))]] += v * w
        for o in range(PORTS):
            check(any(rows[o]), "RAW_ZERO", "identically zero coordinate equation")
            raw[(o, i, j, k)] = tuple(rows[o])
    check(len(raw) == 2640, "RAW_COUNT", "wrong raw equation count")

    # orbits, representatives, rank 38
    unassigned = set(raw)
    orbits = []
    while unassigned:
        seed = min(unassigned)
        signs = {}
        for g in group:
            o, i, j, k = seed
            images = [g[i], g[j], g[k]]
            image = (g[o], *sorted(images))
            signs[image] = -1 if parity(images) else 1
        unassigned -= set(signs)
        orbits.append((seed, signs))
    check(len(orbits) == 44, "ORBITS44", "expected forty-four orbits")
    for seed, signs in orbits:
        for coord, s in signs.items():
            check(raw[coord] == tuple(s * v for v in raw[seed]), "ORBIT_CONST", "orbit sign constancy fails")
    reps = [raw[seed] for seed, _ in orbits]
    rep_ech, rep_piv = echelon([[F5(v) for v in row] for row in reps])
    check(len(rep_ech) == 38, "RANK38", "representative rank is not 38")
    jt = system["jacobi_tensor"]
    check(jt["representative_coefficient_rank_over_Q"] == 38, "PINNED38", "pinned rank differs")
    pinned_by_coord = {}
    for p in jt["representative_equations"]:
        row = [Fraction(0)] * len(MONOMIALS)
        for a, b, n, d in p["equation"]:
            row[MIDX[(a, b)]] = Fraction(n, d)
        pinned_by_coord[tuple(p["target_coordinate"])] = tuple(row)
    check(sum(1 for seed, _ in orbits if pinned_by_coord.get(seed) == raw[seed]) == 44, "PINNED_REPS", "pinned representatives differ")

    # contraction rank 11 and split 11+27
    contraction = []
    for o in range(PORTS):
        for i, j in itertools.combinations(range(PORTS), 2):
            total = [Fraction(0)] * len(MONOMIALS)
            for s in range(PORTS):
                if len({s, i, j}) < 3:
                    continue
                sign = -1 if parity((s, i, j)) else 1
                for idx, v in enumerate(raw[(o, *sorted((s, i, j)))]):
                    if v:
                        total[idx] += sign * v
            contraction.append(total)
    check(len(contraction) == 792, "CON792", "wrong contracted row count")
    con_ech, con_piv = echelon([[F5(v) for v in row] for row in contraction])
    check(len(con_ech) == 11, "RANK11", "contracted rank is not 11")

    fl = system["fixed_line_reduction"]
    transform = [[ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(fl["channel_decomposition"]["transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            transform[i][j] = F5(Fraction(n1, d1), Fraction(n2, d2))
    inverse = [[ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(fl["channel_decomposition"]["inverse_transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            inverse[i][j] = F5(Fraction(n1, d1), Fraction(n2, d2))
    for i in range(PARAMS):
        for j in range(PARAMS):
            v = sum((transform[i][k] * inverse[k][j] for k in range(PARAMS) if transform[i][k] and inverse[k][j]), ZERO)
            check(v == (ONE if i == j else ZERO), "TINV", "pinned transform times inverse is not identity")

    weights = [(1, 0, 0), (0, 1, 0), (-1, 0, 2), (0, -1, 2), (1, 1, -1), (0, 0, 1), (1, -1, 1), (1, 0, 0), (-1, 1, 1), (0, 0, 1), (0, 1, 0)]
    product_rows = []
    for number, weight in enumerate(weights):
        wform = [sum((F5(w) * transform[t][col] for t, w in enumerate(weight) if w), ZERO) for col in range(PARAMS)]
        cform = transform[3 + number]
        row = [ZERO] * len(MONOMIALS)
        for i, va in enumerate(wform):
            if not va:
                continue
            for j, vb in enumerate(cform):
                if vb:
                    row[MIDX[tuple(sorted((i, j)))]] = row[MIDX[tuple(sorted((i, j)))]] + va * vb
        product_rows.append(row)
    prod_ech, prod_piv = echelon(product_rows)
    check(len(prod_ech) == 11, "PROD11", "product rank is not 11")
    for row in con_ech:
        check(member(prod_ech, prod_piv, row), "SPLIT_A", "contracted row leaves the product rowspace")
    for row in prod_ech:
        check(member(con_ech, con_piv, row), "SPLIT_B", "product row leaves the contracted rowspace")
    residual = []
    for entries in fl["residual_integer_equations_in_x"]:
        row = [ZERO] * len(MONOMIALS)
        for a, b, v in entries:
            row[MIDX[(a, b)]] = F5(v)
        residual.append(row)
    check(len(residual) == 27 and len(echelon(residual)[0]) == 27, "RES27", "residual rank is not 27")
    for row in residual:
        check(member(rep_ech, rep_piv, row), "RES_IN_FULL", "residual row leaves the Jacobi rowspace")
    check(len(echelon(prod_ech + residual)[0]) == 38, "SPLIT38", "product plus residual rank is not 38")
    check(fl["full_rank_split"] == [11, 27, 38], "PINNED_SPLIT", "pinned split differs")
    rep = certificate["reproduction"]
    check(
        rep["representative_rank_over_Q"] == 38
        and rep["contracted_rank_over_Q"] == 11
        and rep["residual_rank_over_Qsqrt5"] == 27
        and rep["split_matches_pinned"] == [11, 27, 38],
        "CERT_REPRO",
        "certificate reproduction block differs",
    )

    # spectral identities
    orbital = fl["channel_decomposition"]["spectral_decomposition"]["canonical_valency_five_orbital"]
    A = [[ZERO] * PORTS for _ in range(PORTS)]
    for i, j in orbital:
        A[i][j] = ONE
    eigenvalues = {"fixed": F5(5), "three_plus": RT5, "three_minus": -RT5, "five": F5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        M = [[ONE if i == j else ZERO for j in range(PORTS)] for i in range(PORTS)]
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else ZERO) for j in range(PORTS)] for i in range(PORTS)]
            scale = (ev - ev2).inv()
            M = [
                [scale * sum((M[i][k] * factor[k][j] for k in range(PORTS) if M[i][k] and factor[k][j]), ZERO) for j in range(PORTS)]
                for i in range(PORTS)
            ]
        projectors[label] = M
    for label, P in projectors.items():
        sq = [[sum((P[i][k] * P[k][j] for k in range(PORTS) if P[i][k] and P[k][j]), ZERO) for j in range(PORTS)] for i in range(PORTS)]
        check(sq == P, "PROJ_IDEM", f"{label} projector is not idempotent")
        check(sum((P[i][i] for i in range(PORTS)), ZERO) == F5(SECTOR_DIMS[label]), "PROJ_TRACE", f"{label} trace differs")
    for i in range(PORTS):
        for j in range(PORTS):
            total = sum((projectors[l][i][j] for l in SECTORS), ZERO)
            check(total == (ONE if i == j else ZERO), "PROJ_SUM", "projectors do not resolve the identity")

    derivations = []
    for a in range(PARAMS):
        derivations.append([[F5(sum((C[a][o][s][c] for s in range(PORTS)), Fraction(0))) for c in range(PORTS)] for o in range(PORTS)])
    for a in range(PARAMS):
        claimed = [
            [
                transform[0][a] * projectors["three_plus"][i][j]
                + transform[1][a] * projectors["three_minus"][i][j]
                + transform[2][a] * projectors["five"][i][j]
                for j in range(PORTS)
            ]
            for i in range(PORTS)
        ]
        check(derivations[a] == claimed, "MP4", "ad_u spectral identity fails")
    for a in range(PARAMS):
        for b in range(a, PARAMS):
            gram = sum(
                (derivations[a][i][j] * derivations[b][j][i] for i in range(PORTS) for j in range(PORTS) if derivations[a][i][j] and derivations[b][j][i]),
                ZERO,
            )
            claimed = F5(3) * transform[0][a] * transform[0][b] + F5(3) * transform[1][a] * transform[1][b] + F5(5) * transform[2][a] * transform[2][b]
            check(gram == claimed, "MP3", "K(u,u) Gram identity fails")
    for a in range(PARAMS):
        for i in range(PORTS):
            for j in range(i + 1, PORTS):
                check(sum((C[a][o][i][j] for o in range(PORTS)), Fraction(0)) == 0, "MP2", "fixed-sector output appears")
    pair_orbit_seen = set()
    orbit_count = 0
    unassigned_pairs = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    while unassigned_pairs:
        seed = min(unassigned_pairs)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned_pairs -= orbit
        orbit_count += 1
    check(orbit_count == 4, "MP1", "ordered-pair orbit count is not four")

    # restricted spans
    reps_f5 = [[F5(v) for v in row] for row in reps]

    def span_of(support):
        pairs = [(support[i], support[j]) for i in range(len(support)) for j in range(i, len(support))]
        pidx = {p: i for i, p in enumerate(pairs)}
        rows = []
        for row in reps_f5:
            out = [ZERO] * len(pairs)
            for (alpha, beta), coefficient in zip(MONOMIALS, row):
                if not coefficient:
                    continue
                for ci in support:
                    va = inverse[alpha][ci]
                    if not va:
                        continue
                    for cj in support:
                        vb = inverse[beta][cj]
                        if vb:
                            slot = pidx[tuple(sorted((ci, cj)))]
                            out[slot] = out[slot] + coefficient * va * vb
            rows.append(out)
        ech, piv = echelon(rows)
        return pairs, ech, piv

    for case_id, support in [("A", [3]), ("B", [4]), ("C", [3, 4])]:
        _, ech, _ = span_of(support)
        check(len(ech) == 0, f"CASE_{case_id}", "restricted system does not vanish")
        check(certificate["case_analysis"][case_id]["span_dimension"] == 0, f"CERT_CASE_{case_id}", "certificate case record differs")

    def decode_vec(entries, length):
        row = [ZERO] * length
        for entry in entries:
            row[entry[0]] = decode_f5(entry[1:])
        return row

    for case_id, support, forced_zero_channel, quad_channel, partner, branch_scale, free_channel in [
        ("D", [3, 5, 8, 10], 8, 10, 3, -RT5, 5),
        ("E", [4, 6, 12, 13], 12, 13, 4, RT5, 6),
    ]:
        pairs, ech, piv = span_of(support)
        check(len(ech) == 5, f"CASE_{case_id}_DIM", "span dimension is not five")
        record = certificate["case_analysis"][case_id]
        check(record["support"] == support and record["span_dimension"] == 5, f"CERT_CASE_{case_id}", "certificate case record differs")
        claimed = [decode_vec(entries, len(pairs)) for entries in record["claimed_basis"]]
        for row in claimed:
            check(member(ech, piv, row), f"CASE_{case_id}_MEMBER", "claimed basis row leaves the span")
        check(len(echelon(claimed)[0]) == 5, f"CASE_{case_id}_RANK", "claimed basis is dependent")
        # structural shape of the claimed basis rows, replayed independently:
        # rows 0..2 are the pure product monomials with the forced-zero channel.
        expected_products = {tuple(sorted((forced_zero_channel, other))) for other in support if other != forced_zero_channel}
        seen_products = set()
        for row in claimed[:3]:
            nz = [(pairs[i], v) for i, v in enumerate(row) if v]
            check(len(nz) == 1 and nz[0][1] == ONE, f"CASE_{case_id}_PRODUCT", "product row is not a single monomial")
            seen_products.add(nz[0][0])
        check(seen_products == expected_products, f"CASE_{case_id}_PRODUCTSET", "product monomials differ")
        # row 3 is quad_channel*(quad_channel - branch_scale*partner).
        row3 = claimed[3]
        for i, v in enumerate(row3):
            pair = pairs[i]
            if pair == (quad_channel, quad_channel):
                check(v == ONE, f"CASE_{case_id}_Q1", "quadratic row diagonal differs")
            elif pair == tuple(sorted((partner, quad_channel))):
                check(v == -branch_scale, f"CASE_{case_id}_Q2", "quadratic row cross term differs")
            else:
                check(not v, f"CASE_{case_id}_Q3", "quadratic row has an extra term")
        # row 4: on the branch quad_channel = branch_scale*partner with the forced channel zero,
        # the row vanishes identically; at quad_channel=0 it reduces to a nonzero multiple of
        # partner*free_channel; its forced-channel-squared coefficient is nonzero.
        row4 = claimed[4]
        coeff = {pairs[i]: v for i, v in enumerate(row4) if v}
        check(set(coeff) <= {tuple(sorted((partner, free_channel))), tuple(sorted((quad_channel, free_channel))), (forced_zero_channel, forced_zero_channel)}, f"CASE_{case_id}_R4", "fifth row support differs")
        c_pf = coeff.get(tuple(sorted((partner, free_channel))), ZERO)
        c_qf = coeff.get(tuple(sorted((quad_channel, free_channel))), ZERO)
        c_ff = coeff.get((forced_zero_channel, forced_zero_channel), ZERO)
        check(bool(c_ff), f"CASE_{case_id}_R4A", "forced-channel square coefficient vanishes")
        check(bool(c_pf), f"CASE_{case_id}_R4B", "partner-free coefficient vanishes")
        check(c_pf + branch_scale * c_qf == ZERO, f"CASE_{case_id}_R4C", "fifth row does not vanish on the branch")

    for case_id, support, absent in [("F", [3, 4, 5, 8, 10], 4), ("G", [3, 4, 6, 12, 13], 3)]:
        pairs, ech, _ = span_of(support)
        check(len(ech) == 5, f"CASE_{case_id}_DIM", "span dimension is not five")
        for row in ech:
            for (i, j), v in zip(pairs, row):
                if v:
                    check(absent not in (i, j), f"CASE_{case_id}_CROSS", "span contains a cross monomial")

    # families and samples
    def channel_values_for(family_id, params):
        a, b, e = (list(params) + [0, 0, 0])[:3]
        a, b, e = F5(a), F5(b), F5(e)
        if family_id == "PLANE":
            return {3: a, 4: b}
        if family_id == "SU3PF":
            return {3: a, 5: b, 10: F5(0, -1) * a, 4: e}
        if family_id == "SU3MF":
            return {4: a, 6: b, 13: F5(0, 1) * a, 3: e}
        if family_id == "DERIV":
            return {0: a, 1: b, 2: e}
        raise VerifyError("FAMILY", f"unknown family {family_id}")

    def bracket_of(channel_values):
        x = [sum((inverse[i][j] * v for j, v in channel_values.items()), ZERO) for i in range(PARAMS)]
        B = [[[ZERO] * PORTS for _ in range(PORTS)] for _ in range(PORTS)]
        for a in range(PARAMS):
            if not x[a]:
                continue
            for o in range(PORTS):
                for i in range(PORTS):
                    row = C[a][o][i]
                    for j in range(PORTS):
                        if row[j]:
                            B[o][i][j] = B[o][i][j] + x[a] * row[j]
        return x, B

    def killing_of(B):
        K = [[ZERO] * PORTS for _ in range(PORTS)]
        for i in range(PORTS):
            for j in range(i, PORTS):
                v = sum(
                    (B[m][i][o] * B[o][j][m] for m in range(PORTS) for o in range(PORTS) if B[m][i][o] and B[o][j][m]),
                    ZERO,
                )
                K[i][j] = v
                K[j][i] = v
        return K

    sector_bases = {label: echelon([[projectors[label][i][j] for j in range(PORTS)] for i in range(PORTS)])[0] for label in SECTORS}

    def restrict(K, rows):
        return [
            [
                sum((rows[r][i] * K[i][j] * rows[s][j] for i in range(PORTS) for j in range(PORTS) if rows[r][i] and K[i][j] and rows[s][j]), ZERO)
                for s in range(len(rows))
            ]
            for r in range(len(rows))
        ]

    param_counts = {"PLANE": 2, "SU3PF": 3, "SU3MF": 3, "DERIV": 3}
    expected_blocks = {
        "PLANE": {(0, 0): ("three_plus", [0, 3, 0]), (1, 1): ("three_minus", [0, 3, 0])},
        "SU3PF": {(0, 0): ("three_plus", [0, 3, 0]), (0, 1): ("five", [5, 0, 0]), (2, 2): ("three_minus", [0, 3, 0])},
        "SU3MF": {(0, 0): ("three_minus", [0, 3, 0]), (0, 1): ("five", [0, 5, 0]), (2, 2): ("three_plus", [0, 3, 0])},
        "DERIV": {(0, 0): ("fixed", [1, 0, 0]), (1, 1): ("fixed", [1, 0, 0]), (2, 2): ("fixed", [1, 0, 0])},
    }
    for family_id, n in param_counts.items():
        points = []
        for i in range(n):
            e = [0] * n
            e[i] = 1
            points.append(tuple(e))
        for i in range(n):
            for j in range(i + 1, n):
                e = [0] * n
                e[i] = e[j] = 1
                points.append(tuple(e))
        K_at = {}
        for point in points:
            x, B = bracket_of(channel_values_for(family_id, point))
            for row in reps_f5:
                v = sum(
                    (coefficient * x[alpha] * x[beta] for (alpha, beta), coefficient in zip(MONOMIALS, row) if coefficient and x[alpha] and x[beta]),
                    ZERO,
                )
                check(not v, "FAM_JACOBI", f"{family_id} violates Jacobi at {point}")
            K_at[point] = killing_of(B)
        monomials = {}
        for i in range(n):
            e = [0] * n
            e[i] = 1
            monomials[(i, i)] = K_at[tuple(e)]
        for i in range(n):
            for j in range(i + 1, n):
                e = [0] * n
                e[i] = e[j] = 1
                monomials[(i, j)] = [
                    [K_at[tuple(e)][r][c] - monomials[(i, i)][r][c] - monomials[(j, j)][r][c] for c in range(PORTS)]
                    for r in range(PORTS)
                ]
        nonzero = {m: M for m, M in monomials.items() if any(any(row) for row in M)}
        check(set(nonzero) == set(expected_blocks[family_id]), "FAM_MONO", f"{family_id} monomial support differs")
        cert_monos = certificate["families"][family_id]["killing_monomials"]
        check(len(cert_monos) == len(nonzero), "CERT_MONO_COUNT", f"{family_id} certificate monomial count differs")
        for mono, M in nonzero.items():
            key = "*".join(f"p{t}" for t in mono)
            record = cert_monos[key]
            stored = [[ZERO] * PORTS for _ in range(PORTS)]
            for i, j, n1, d1, n2, d2 in record["matrix"]:
                stored[i][j] = F5(Fraction(n1, d1), Fraction(n2, d2))
                stored[j][i] = stored[i][j]
            check(stored == M, "CERT_MONO_MATRIX", f"{family_id} {key} stored matrix differs")
            label, sig = expected_blocks[family_id][mono]
            check(record["sector"] == label and record["block_signature"] == sig, "CERT_MONO_META", f"{family_id} {key} metadata differs")
            block = restrict(M, sector_bases[label])
            check(inertia(block) == sig, "FAM_SIG", f"{family_id} {key} block signature differs")
            for ia, la in enumerate(SECTORS):
                for lb in SECTORS[ia + 1:]:
                    cross = [
                        [sum((va[i] * M[i][j] * vb[j] for i in range(PORTS) for j in range(PORTS) if va[i] and M[i][j] and vb[j]), ZERO) for vb in sector_bases[lb]]
                        for va in sector_bases[la]
                    ]
                    check(not any(any(row) for row in cross), "FAM_CROSS", f"{family_id} has a cross-sector block")
                for other in SECTORS:
                    if other == label:
                        continue
                    block_other = restrict(M, sector_bases[other])
                    check(not any(any(row) for row in block_other), "FAM_SUPPORT", f"{family_id} {key} leaks outside {label}")

    # samples
    for record in certificate["samples"]:
        family_id = record["family"]
        params = tuple(Fraction(p) for p in record["parameters"])
        channel_values = channel_values_for(family_id, params)
        stored_channels = {CHANNEL_INDEX[cid]: decode_f5(entry) for cid, entry in record["channel_values"].items()}
        check(stored_channels == {k: v for k, v in channel_values.items()}, "SAMPLE_CH", "stored channel values differ")
        x, B = bracket_of(channel_values)
        stored_x = [ZERO] * PARAMS
        for entry in record["x_vector"]:
            stored_x[entry[0]] = decode_f5(entry[1:])
        check(stored_x == x, "SAMPLE_X", "stored x vector differs")
        for (i, j, k) in itertools.combinations(range(PORTS), 3):
            for o in range(PORTS):
                v = sum(
                    (B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j] for m in range(PORTS)),
                    ZERO,
                )
                check(not v, "SAMPLE_JACOBI", f"{family_id}{params} violates Jacobi directly")
        K = killing_of(B)
        check(inertia(K) == record["killing_signature"], "SAMPLE_SIG", "stored Killing signature differs")
        image_rows = []
        for i in range(PORTS):
            for j in range(i + 1, PORTS):
                col = [B[o][i][j] for o in range(PORTS)]
                if any(col):
                    image_rows.append(col)
        derived, _ = echelon(image_rows) if image_rows else ([], [])
        check(len(derived) == record["derived_dimension"], "SAMPLE_DERIVED", "derived dimension differs")
        flat_rows = []
        for v in range(PORTS):
            flat_rows.append([B[o][v][j] for o in range(PORTS) for j in range(PORTS)])
        tr = [[flat_rows[i][j] for i in range(PORTS)] for j in range(len(flat_rows[0]))]
        work, piv = echelon(tr)
        pivset = set(piv)
        kernel = []
        for fc in (c for c in range(PORTS) if c not in pivset):
            vec = [ZERO] * PORTS
            vec[fc] = ONE
            for row, pc in zip(work, piv):
                vec[pc] = -row[fc]
            kernel.append(vec)
        for vec in kernel:
            for o in range(PORTS):
                for j in range(PORTS):
                    s = sum((vec[v] * B[o][v][j] for v in range(PORTS) if vec[v]), ZERO)
                    check(not s, "SAMPLE_CENTER", "center vector fails centrality")
        check(len(kernel) == record["center_dimension"], "SAMPLE_CENTER_DIM", "center dimension differs")
        direct = False
        if derived or kernel:
            union_rank = len(echelon(derived + kernel)[0])
            direct = union_rank == len(derived) + len(kernel) == PORTS
        else:
            direct = PORTS == 0
        check(direct == record["center_plus_derived_is_direct_and_full"], "SAMPLE_DIRECT", "direct-sum flag differs")
        negdef = True
        if derived:
            negdef = inertia(restrict(K, derived)) == [0, len(derived), 0]
        check(negdef == record["killing_negative_definite_on_derived"], "SAMPLE_NEGDEF", "negative-definiteness flag differs")
        check(record["compact"] == (direct and negdef), "SAMPLE_COMPACT", "compact verdict differs")

    # classification block consistency
    families = {f["id"]: f for f in certificate["classification"]["compact_families"]}
    check(set(families) == {"P", "F", "G"}, "CLASS_IDS", "compact family identifiers differ")
    check(families["F"]["sign_condition"] == "a*b < 0", "CLASS_F_SIGN", "family F sign condition differs")
    check(families["G"]["sign_condition"] == "a*b > 0", "CLASS_G_SIGN", "family G sign condition differs")
    check(families["P"]["dimension"] == 2 and families["F"]["dimension"] == 3 and families["G"]["dimension"] == 3, "CLASS_DIMS", "family dimensions differ")

    compact_samples = [r for r in certificate["samples"] if r["compact"]]
    check(len(compact_samples) == 7, "SAMPLE_COUNT", "expected seven compact sample points")

    return {
        "verified": True,
        "representative_rank": 38,
        "split": [11, 27, 38],
        "compact_families": sorted(families),
        "samples_checked": len(certificate["samples"]),
    }


if __name__ == "__main__":
    result = main()
    print(json.dumps(result, indent=1, sort_keys=True))
