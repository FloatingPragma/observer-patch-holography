#!/usr/bin/env python3
"""Independent exact replay of the B14 oriented-face discriminator packet.

No code is imported from the producer.  This verifier separately rebuilds the
face tensor, Reynolds/channel coordinates, Jacobi census, P/F/G least-square
normal equations, and the exact L1/Linfinity primal-dual certificates over
Q(sqrt(5)).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT = HERE / "oriented_face_bracket_selector.certificate.json"
PRODUCER = HERE / "oriented_face_bracket_selector.py"
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM = REPO / "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json"
COMPACT = HERE / "b14_compact_locus.certificate.json"

N = 12
M = 14
CHANNELS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m",
    "t_ff_to_p", "t_ff_to_m", "t_pm_to_f", "t_pf_to_p",
    "t_pf_to_m", "t_pf_to_f", "t_mf_to_p", "t_mf_to_m",
    "t_mf_to_f",
]
CH = {name: i for i, name in enumerate(CHANNELS)}
COORDS = [(o, i, j) for o in range(N) for i in range(N) for j in range(i + 1, N)]


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

    def __add__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return F5(self.r + other.r, self.s + other.s)

    __radd__ = __add__

    def __neg__(self):
        return F5(-self.r, -self.s)

    def __sub__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self + (-other)

    def __mul__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return F5(self.r * other.r + 5 * self.s * other.s,
                  self.r * other.s + self.s * other.r)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.r * self.r - 5 * self.s * self.s
        check(norm != 0, "division by zero")
        return F5(self.r / norm, -self.s / norm)

    def __truediv__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self * other.inverse()

    def __eq__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self.r == other.r and self.s == other.s

    def __bool__(self):
        return bool(self.r or self.s)


Z = F5()
U = F5(1)
S = F5(0, 1)


def encode(x):
    return [x.r.numerator, x.r.denominator, x.s.numerator, x.s.denominator]


def file_sha(path):
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def obj_sha(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def exact_sign(x):
    if not x:
        return 0
    if x.r == 0:
        return 1 if x.s > 0 else -1
    if x.s == 0:
        return 1 if x.r > 0 else -1
    if (x.r > 0) == (x.s > 0):
        return 1 if x.r > 0 else -1
    delta = x.r * x.r - 5 * x.s * x.s
    check(delta != 0, "unexpected zero sign boundary")
    if x.r > 0:
        return 1 if delta > 0 else -1
    return 1 if delta < 0 else -1


def f_abs(value):
    return value if exact_sign(value) >= 0 else -value


def sum_abs(values):
    return sum((f_abs(value) for value in values), Z)


def max_abs(values):
    result = Z
    for value in values:
        absolute = f_abs(value)
        if exact_sign(absolute - result) > 0:
            result = absolute
    return result


def pairing(left, right):
    return sum((a * b for a, b in zip(left, right)), Z)


def family_point(vectors, parameters):
    check(len(vectors) == len(parameters), "family parameter length mismatch")
    return [sum((parameters[j] * vectors[j][i] for j in range(len(vectors))), Z)
            for i in range(len(COORDS))]


def coordinate_image(coordinate, permutation):
    out, left, right = coordinate
    image_out = permutation[out]
    image_left = permutation[left]
    image_right = permutation[right]
    orientation = 1
    if image_right < image_left:
        image_left, image_right = image_right, image_left
        orientation = -1
    return (image_out, image_left, image_right), orientation


def orbit_covector(seeds, group):
    coordinate_index = {coordinate: i for i, coordinate in enumerate(COORDS)}
    values = [Z for _ in COORDS]
    assigned = [False for _ in COORDS]
    for coordinate, seed_value in seeds:
        check(coordinate in coordinate_index, "dual seed outside coordinate domain")
        for permutation in group:
            image, orientation = coordinate_image(coordinate, permutation)
            index = coordinate_index[image]
            value = orientation * seed_value
            if assigned[index]:
                check(values[index] == value, "inconsistent signed dual orbit")
            else:
                values[index] = value
                assigned[index] = True
    return values


def sparse_covector(entries):
    coordinate_index = {coordinate: i for i, coordinate in enumerate(COORDS)}
    values = [Z for _ in COORDS]
    for coordinate, value in entries:
        index = coordinate_index[coordinate]
        check(not values[index], "duplicate sparse dual coordinate")
        values[index] = value
    return values


def lp_record(target, vectors, parameters, distance, dual, norm_kind):
    point = family_point(vectors, parameters)
    residual = [target[i] - point[i] for i in range(len(target))]
    if norm_kind == "l1":
        primal = sum_abs(residual)
        dual_norm = max_abs(dual)
    elif norm_kind == "linfinity":
        primal = max_abs(residual)
        dual_norm = sum_abs(dual)
    else:
        raise VerificationError(f"unsupported norm {norm_kind}")
    check(primal == distance, f"{norm_kind} primal mismatch")
    check(exact_sign(F5(1) - dual_norm) >= 0, f"{norm_kind} dual infeasible")
    annihilations = [pairing(dual, vector) for vector in vectors]
    check(not any(annihilations), f"{norm_kind} dual misses annihilator")
    objective = pairing(dual, target)
    check(objective == distance, f"{norm_kind} duality gap")
    return {
        "parameters": [encode(value) for value in parameters],
        "distance": encode(distance),
        "primal_residual_support": sum(bool(value) for value in residual),
        "primal_active_max_coordinates": (
            sum(f_abs(value) == distance for value in residual)
            if norm_kind == "linfinity" else None
        ),
        "dual_support": sum(bool(value) for value in dual),
        "dual_norm": encode(dual_norm),
        "dual_objective": encode(objective),
        "dual_annihilates_family": [encode(value) for value in annihilations],
    }


def matrix_from_sparse(rows):
    matrix = [[Z for _ in range(M)] for _ in range(M)]
    for i, row in enumerate(rows):
        for j, an, ad, bn, bd in row:
            matrix[i][j] = F5(Fraction(an, ad), Fraction(bn, bd))
    return matrix


def solve_normal(gram, rhs):
    n = len(rhs)
    a = [list(row) + [rhs[i]] for i, row in enumerate(gram)]
    for c in range(n):
        pivot = next((r for r in range(c, n) if a[r][c]), None)
        check(pivot is not None, "singular normal equations")
        a[c], a[pivot] = a[pivot], a[c]
        inv = a[c][c].inverse()
        a[c] = [inv * value for value in a[c]]
        for r in range(n):
            if r != c and a[r][c]:
                multiplier = a[r][c]
                a[r] = [x - multiplier * y for x, y in zip(a[r], a[c])]
    return [a[i][-1] for i in range(n)]


def replay():
    manifest = json.loads(MANIFEST.read_text())
    raw_basis = json.loads(BASIS.read_text())
    receipt = json.loads(RECEIPT.read_text())
    system = json.loads(SYSTEM.read_text())

    ports = manifest["carrier"]["ports"]
    ix = {port: i for i, port in enumerate(ports)}
    faces = [tuple(ix[p] for p in face) for face in manifest["carrier"]["oriented_faces"]]
    check(len(faces) == 20 and len(set(faces)) == 20, "source faces differ")
    B = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for i, j, k in faces:
        for o, u, v in ((k, i, j), (i, j, k), (j, k, i)):
            B[o][u][v] += 1
            B[o][v][u] -= 1

    reynolds = [[[[Fraction(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in range(M)]
    for a, row in enumerate(raw_basis["basis"]):
        for o, i, j, num, den in row["entries"]:
            value = Fraction(num, den)
            reynolds[a][o][i][j] = value
            reynolds[a][o][j][i] = -value
    check(all(Fraction(B[o][i][j]) == 60 * reynolds[13][o][i][j]
              for o in range(N) for i in range(N) for j in range(N)), "60*R13 identity fails")
    target = [F5(B[o][i][j]) for o, i, j in COORDS]
    norm = sum((x * x for x in target), Z)

    decomposition = system["fixed_line_reduction"]["channel_decomposition"]
    transform = matrix_from_sparse(decomposition["transform_rows"])
    inverse = matrix_from_sparse(decomposition["inverse_transform_rows"])
    x13 = [Z for _ in range(M)]
    x13[13] = F5(60)
    channel_values = [sum((transform[i][j] * x13[j] for j in range(M)), Z) for i in range(M)]

    entries = []
    witnesses = []
    for o in range(N):
        for i in range(N):
            for j in range(i + 1, N):
                for k in range(j + 1, N):
                    value = sum(
                        B[m][i][j] * B[o][m][k]
                        + B[m][j][k] * B[o][m][i]
                        + B[m][k][i] * B[o][m][j]
                        for m in range(N)
                    )
                    if value:
                        entries.append(value)
                        if len(witnesses) < 12:
                            witnesses.append([o, i, j, k, value])
    jacobi = {
        "coordinate_domain": "output o and strictly increasing input triple i<j<k",
        "domain_size": N * (N * (N - 1) * (N - 2) // 6),
        "nonzero_count": len(entries),
        "positive_count": sum(v > 0 for v in entries),
        "negative_count": sum(v < 0 for v in entries),
        "value_histogram": {str(v): entries.count(v) for v in sorted(set(entries))},
        "squared_norm": sum(v * v for v in entries),
        "first_witnesses": witnesses,
    }

    def vector_for_channels(spec):
        y = [Z for _ in range(M)]
        for name, value in spec.items():
            y[CH[name]] = value
        x = [sum((inverse[i][j] * y[j] for j in range(M)), Z) for i in range(M)]
        return [sum((x[a] * reynolds[a][o][i][j] for a in range(M)), Z) for o, i, j in COORDS]

    specifications = {
        "P": (["p", "m"], [{"t_pp_to_p": U}, {"t_mm_to_m": U}]),
        "F": (["a", "b", "e"], [
            {"t_pp_to_p": U, "t_pf_to_f": -S}, {"t_ff_to_p": U}, {"t_mm_to_m": U},
        ]),
        "G": (["a", "b", "e"], [
            {"t_mm_to_m": U, "t_mf_to_f": S}, {"t_ff_to_m": U}, {"t_pp_to_p": U},
        ]),
    }
    records = {}
    raw_distances = {}
    family_vectors = {}
    for family, (names, specs) in specifications.items():
        vectors = [vector_for_channels(spec) for spec in specs]
        family_vectors[family] = vectors
        gram = [[sum((a * b for a, b in zip(left, right)), Z) for right in vectors] for left in vectors]
        rhs = [sum((a * b for a, b in zip(vector, target)), Z) for vector in vectors]
        coefficients = solve_normal(gram, rhs)
        residual = [target[i] - sum((coefficients[j] * vectors[j][i] for j in range(len(vectors))), Z)
                    for i in range(len(target))]
        distance = sum((value * value for value in residual), Z)
        raw_distances[family] = distance
        record = {
            "parameter_order": names,
            "gram": [[encode(value) for value in row] for row in gram],
            "rhs": [encode(value) for value in rhs],
            "orthogonal_projection_parameters": [encode(value) for value in coefficients],
            "squared_distance": encode(distance),
        }
        if family in ("F", "G"):
            product = coefficients[0] * coefficients[1]
            record.update({
                "compact_sign_product_a_times_b": encode(product),
                "compact_sign": "negative" if exact_sign(product) < 0 else "positive",
                "projection_lies_in_compact_open_stratum": exact_sign(product) < 0 if family == "F" else exact_sign(product) > 0,
            })
        records[family] = record

    group = receipt["proper_port_action"]["permutation_rows"]
    check(len({tuple(row) for row in group}) == 60,
          "proper group pin differs")

    l1_parameters = {
        "P": [Z, Z],
        "F": [Z, Z, Z],
        "G": [
            F5(Fraction(-1, 8), Fraction(1, 40)),
            F5(Fraction(1, 8), Fraction(-1, 8)),
            S / 10,
        ],
    }
    l1_distances = {"P": F5(60), "F": F5(60), "G": F5(-30, 30)}
    l1_seeds = {
        "P": [
            ((0, 1, 8), F5(Fraction(3, 10))),
            ((0, 1, 9), F5(Fraction(4, 15))),
            ((0, 1, 10), F5(Fraction(-13, 30))),
            ((0, 2, 8), U),
        ],
        "F": [
            ((0, 1, 5), F5(Fraction(-1, 2), Fraction(1, 2))),
            ((0, 1, 8), F5(Fraction(-2, 5), Fraction(1, 4))),
            ((0, 1, 9), F5(Fraction(-1, 4), Fraction(1, 4))),
            ((0, 1, 10), F5(Fraction(-3, 20))),
            ((0, 2, 8), U),
        ],
        "G": [
            ((0, 1, 4), F5(Fraction(-17, 40), Fraction(1, 4))),
            ((0, 1, 5), F5(-1)),
            ((0, 1, 6), F5(Fraction(-11, 60))),
            ((0, 1, 9), F5(Fraction(-1, 60))),
            ((0, 2, 4), F5(Fraction(1, 4), Fraction(-29, 120))),
            ((0, 2, 8), F5(Fraction(-1, 2), Fraction(1, 2))),
        ],
    }
    l1 = {}
    for family in ("P", "F", "G"):
        dual = orbit_covector(l1_seeds[family], group)
        l1[family] = lp_record(
            target, family_vectors[family], l1_parameters[family],
            l1_distances[family], dual, "l1",
        )
        l1[family]["dual_orbit_seeds"] = [
            [*coordinate, *encode(value)] for coordinate, value in l1_seeds[family]
        ]
    g_l1_product = l1_parameters["G"][0] * l1_parameters["G"][1]
    check(exact_sign(g_l1_product) > 0, "L1 G optimizer not compact")
    l1["G"]["compact_sign_product_a_times_b"] = encode(g_l1_product)
    l1["G"]["optimizer_lies_in_compact_open_stratum"] = True
    l1["P"]["compact_attainment"] = "zero belongs to compact family P"
    l1["F"]["compact_attainment"] = (
        "60 is the compact-F infimum: the linear-family minimizer is zero, "
        "and compact F points scale continuously to zero"
    )

    linf_parameters = {
        "P": [Z, F5(Fraction(1, 4), Fraction(-1, 4))],
        "F": [
            F5(Fraction(-3, 40), Fraction(1, 40)),
            F5(Fraction(-1, 8), Fraction(3, 40)),
            F5(Fraction(-3, 10)),
        ],
        "G": [
            F5(Fraction(1, 20), Fraction(-1, 20)),
            F5(Fraction(-1, 4), Fraction(1, 20)),
            F5(Fraction(-3, 10), Fraction(1, 10)),
        ],
    }
    linf_distances = {
        "P": F5(Fraction(1, 2)),
        "F": S / 5,
        "G": F5(Fraction(1, 2), Fraction(-1, 10)),
    }
    linf_entries = {
        "P": [
            ((0, 1, 8), F5(Fraction(1, 2))),
            ((0, 2, 8), F5(Fraction(1, 2))),
        ],
        "F": [
            ((0, 1, 5), F5(Fraction(1, 2), Fraction(-1, 10))),
            ((0, 1, 8), F5(Fraction(1, 4), Fraction(-1, 20))),
            ((0, 1, 9), F5(Fraction(1, 4), Fraction(-1, 20))),
            ((0, 2, 8), S / 5),
        ],
        "G": [
            ((0, 1, 4), F5(Fraction(1, 4), Fraction(-1, 20))),
            ((0, 1, 7), S / 5),
            ((0, 2, 6), F5(Fraction(1, 4), Fraction(-1, 20))),
            ((0, 2, 8), F5(Fraction(1, 2), Fraction(-1, 10))),
        ],
    }
    linf = {}
    for family in ("P", "F", "G"):
        dual = sparse_covector(linf_entries[family])
        linf[family] = lp_record(
            target, family_vectors[family], linf_parameters[family],
            linf_distances[family], dual, "linfinity",
        )
        linf[family]["dual_entries"] = [
            [*coordinate, *encode(value)] for coordinate, value in linf_entries[family]
        ]
    f_linf_product = linf_parameters["F"][0] * linf_parameters["F"][1]
    g_linf_product = linf_parameters["G"][0] * linf_parameters["G"][1]
    check(exact_sign(f_linf_product) < 0, "Linfinity F optimizer not compact")
    check(exact_sign(g_linf_product) > 0, "Linfinity G optimizer not compact")
    linf["F"]["compact_sign_product_a_times_b"] = encode(f_linf_product)
    linf["F"]["optimizer_lies_in_compact_open_stratum"] = True
    linf["G"]["compact_sign_product_a_times_b"] = encode(g_linf_product)
    linf["G"]["optimizer_lies_in_compact_open_stratum"] = True

    return {
        "face_count": len(faces),
        "nonzero": sum(bool(x) for x in target),
        "norm": encode(norm),
        "channels": {name: encode(value) for name, value in zip(CHANNELS, channel_values)},
        "jacobi": jacobi,
        "families": records,
        "distances": raw_distances,
        "l1": l1,
        "l1_distances": l1_distances,
        "linfinity": linf,
        "linfinity_distances": linf_distances,
    }


def verify_certificate(path: Path = DEFAULT):
    cert = json.loads(path.read_text())
    check(cert.get("schema") == "oph.b14.oriented_face_bracket_selector.v1", "schema mismatch")
    body = dict(cert)
    stored_hash = body.pop("certificate_sha256", None)
    check(stored_hash == obj_sha(body), "certificate self-hash mismatch")
    check(cert["implementation_pins"]["producer_sha256"] == file_sha(PRODUCER), "producer pin mismatch")
    check(cert["implementation_pins"]["verifier_sha256"] == file_sha(Path(__file__).resolve()), "verifier pin mismatch")
    upstream = cert["upstream"]
    for field, input_path in (
        ("source_manifest_sha256", MANIFEST), ("stage1_basis_sha256", BASIS),
        ("stage1_receipt_sha256", RECEIPT), ("stage2_system_sha256", SYSTEM),
        ("compact_locus_file_sha256", COMPACT),
    ):
        check(upstream[field] == file_sha(input_path), f"{field} pin mismatch")
    compact = json.loads(COMPACT.read_text())
    check(upstream["compact_locus_self_sha256"] == compact["certificate_sha256"], "compact self pin mismatch")

    replayed = replay()
    source = cert["source_face_bracket"]
    check(source["oriented_face_count"] == replayed["face_count"] == 20, "face count mismatch")
    check(source["nonzero_upper_structure_constants"] == replayed["nonzero"] == 60, "support mismatch")
    check(source["upper_coordinate_squared_norm"] == replayed["norm"] == encode(F5(60)), "norm mismatch")
    check(source["reynolds_coordinates"] == [[13, *encode(F5(60))]], "R13 coordinate mismatch")
    check(cert["channel_coordinates"] == replayed["channels"], "channel replay mismatch")
    check(cert["jacobi_failure"] == replayed["jacobi"], "Jacobi replay mismatch")
    check(replayed["jacobi"]["value_histogram"] == {"-1": 120, "1": 120}, "Jacobi census mismatch")

    discriminator = cert["orthogonal_compact_locus_discriminator"]
    check(discriminator["families"] == replayed["families"], "projection replay mismatch")
    check(discriminator["metric_is_source_derived"] is False, "metric boundary weakened")
    check(discriminator["minimum_hs_or_jacobi_repair_is_source_derived"] is False, "repair boundary weakened")
    d = replayed["distances"]
    check(exact_sign(d["F"] - d["G"]) > 0 and exact_sign(d["P"] - d["G"]) > 0,
          "G is not uniquely nearest")
    check(discriminator["unique_nearest_family"] == "G", "nearest-family label mismatch")
    check(discriminator["exact_gaps"] == {
        "distance_F_minus_distance_G": encode(d["F"] - d["G"]),
        "distance_P_minus_distance_G": encode(d["P"] - d["G"]),
    }, "distance gaps mismatch")

    robustness = cert["endpoint_norm_robustness"]
    check(robustness["l1"]["families"] == replayed["l1"], "L1 primal/dual replay mismatch")
    check(robustness["linfinity"]["families"] == replayed["linfinity"],
          "Linfinity primal/dual replay mismatch")
    l1d = replayed["l1_distances"]
    linfd = replayed["linfinity_distances"]
    check(exact_sign(l1d["P"] - l1d["G"]) > 0 and
          exact_sign(l1d["F"] - l1d["G"]) > 0,
          "L1 does not uniquely distinguish G")
    check(exact_sign(linfd["P"] - linfd["G"]) > 0 and
          exact_sign(linfd["F"] - linfd["G"]) > 0,
          "Linfinity does not uniquely distinguish G")
    check(robustness["l1"]["unique_nearest_compact_family_by_minimum_or_infimum"] == "G",
          "L1 nearest-family label mismatch")
    check(robustness["linfinity"]["unique_nearest_compact_family"] == "G",
          "Linfinity nearest-family label mismatch")
    check(robustness["l1"]["exact_gap_P_minus_G"] == encode(l1d["P"] - l1d["G"]),
          "L1 P/G gap mismatch")
    check(robustness["l1"]["exact_gap_F_minus_G"] == encode(l1d["F"] - l1d["G"]),
          "L1 F/G gap mismatch")
    check(robustness["linfinity"]["exact_gap_P_minus_G"] ==
          encode(linfd["P"] - linfd["G"]), "Linfinity P/G gap mismatch")
    check(robustness["linfinity"]["exact_gap_F_minus_G"] ==
          encode(linfd["F"] - linfd["G"]), "Linfinity F/G gap mismatch")
    check(robustness["repair_norm_or_minimization_rule_is_source_derived"] is False,
          "norm-robustness boundary weakened")
    check(robustness["conclusion_status"] ==
          "EXACT_THREE_NORM_ROBUST_CONDITIONAL_DISCRIMINATOR__NOT_SOURCE_SELECTION",
          "norm-robustness status mismatch")
    return {
        "verified": str(path),
        "identity": "B_face=60*R13",
        "jacobi_nonzero": replayed["jacobi"]["nonzero_count"],
        "unique_nearest": "G",
        "robust_norms": ["L1", "L2", "Linfinity"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT)
    args = parser.parse_args()
    print(json.dumps(verify_certificate(args.certificate), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
