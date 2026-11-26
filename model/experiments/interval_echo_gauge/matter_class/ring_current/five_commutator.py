"""Compare the exact marked B14 witness to three existing native operator maps."""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from sympy import Matrix

from .apparatus import cell, geometry

ROOT = Path(__file__).resolve().parents[5]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def frozen_response(weights):
    """Contract the upstream Reynolds tensors, with independent rational and radical parts."""
    base = ROOT / "code"
    cert = json.loads(
        (base / "b14_jacobi/b14_compact_locus.certificate.json").read_text()
    )
    basis = json.loads(
        (
            base
            / "a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
        ).read_text()
    )["basis"]
    reports = {}
    x, y = weights
    cases = {
        "P1": ("PLANE", 3),
        "P2": ("PLANE", 6),
        "F2": ("SU3PF", 11),
        "G2": ("SU3MF", 11),
    }
    for label, (family, dim) in cases.items():
        point = next(
            r
            for r in cert["samples"]
            if r["family"] == family and r["compact"] and r["derived_dimension"] == dim
        )
        out = [[Fraction(), Fraction()] for _ in range(12)]
        for index, a, ad, b, bd in point["x_vector"]:
            for k, i, j, n, d in basis[index]["entries"]:
                wedge = int(x[i] * y[j] - x[j] * y[i])
                out[k][0] += Fraction(n * a * wedge, d * ad)
                out[k][1] += Fraction(n * b * wedge, d * bd)
        # Divide by 12^2 for normalized five projections. The complete physical
        # scale of the channel is still free. Use the sample's b for the unit map.
        coeff = point["channel_values"].get(
            "t_ff_to_p", point["channel_values"].get("t_ff_to_m", [1, 1, 0, 1])
        )
        scale = 144 * Fraction(coeff[0], coeff[1])
        require(coeff[2] == 0 and scale != 0, "Unexpected sample normalization")
        reports[label] = [[c / scale for c in row] for row in out]
    return reports


def certificate():
    rings, anti, _ = geometry()
    identity = np.eye(12, dtype=int)
    J = identity[list(anti)]
    odd = identity - J
    A = 0 * identity
    for a, b in cell().seams:
        A[a, b] = A[b, a] = 1
    require(
        np.array_equal(J @ J, identity)
        and np.trace(odd) == 12
        and np.array_equal(A @ A @ odd, 5 * odd)
        and np.trace(A @ odd) == 0,
        "Odd space does not split into the two marked triplets",
    )
    five12 = 6 * (identity + J) - np.ones_like(identity)
    require(
        np.array_equal(A @ five12, -five12)
        and np.array_equal(five12 @ five12, 12 * five12),
        "Wrong five projection",
    )
    # The B14 chart interchanges native 9/10. Marks 0,1 and their antipodes
    # are unchanged. Check eigenvalue markings in the actual native chart.
    chart = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 11]
    pair = five12[[0, 1]]
    predicted = frozen_response(pair[:, chart])
    for name, raw in predicted.items():
        v = np.array(raw, dtype=object)[chart, :]
        if name in ("P1", "P2"):
            require(not v.any(), "Plane is active on the five")
        else:
            sign = 1 if name == "F2" else -1
            require(
                list(v[4]) == [-1, 0]
                and list(v[8]) == [Fraction(-1, 2), Fraction(sign, 2)],
                "Two-port channel scale or polarity changed",
            )
            rw = v[:, [1, 0]].copy()
            rw[:, 0] *= 5
            require(
                v.any() and np.array_equal(A @ v, sign * rw),
                "Wrong nonzero B14 channel",
            )
    # K numbers are integer seam sums. i*K/5 is the observed signed mean.
    bare = []
    for ring in rings:
        K = 0 * identity
        for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
            K[b, a] = 1
            K[a, b] = -1
        bare.append(K)
    upper = np.triu_indices(12, 1)

    def rank(ops):
        rows = np.stack([k[upper] for k in ops]).tolist()
        return Matrix(rows).to_DM().rank()

    def comm(a, b):
        return a @ b - b @ a

    def weighted(ops, weights):
        return sum(int(v) * k for v, k in zip(weights, ops, strict=True))

    parsed = {}
    for name, ops in (("bare", bare), ("odd", [odd @ k @ odd for k in bare])):
        brackets = [comm(a, b) for a, b in combinations(ops, 2)]
        summary = [
            rank(ops),
            rank(brackets),
            rank(ops + brackets),
            rank(ops + brackets[:1]),
        ]
        target = comm(*(weighted(ops, w) for w in pair))
        require(
            not A[4, 8] and all(k[4, 8] == 0 for k in bare), "Raw witness is a seam"
        )
        if name == "odd":
            require(
                all(
                    np.array_equal(k @ odd, 2 * k) and np.array_equal(k @ A, A @ k)
                    for k in ops
                ),
                "An odd current mixes native triplets",
            )
            require(
                all(np.array_equal(ops[anti[q]], -k) for q, k in enumerate(ops)),
                "Even controls are visible on odd carrier",
            )
        parsed[name] = {
            "ranks": summary,
            "five_pair_entry_4_8": int(target[4, 8]),
            "five_pair_max": int(abs(target).max()),
        }
    require(
        parsed["bare"]
        == {"ranks": [11, 44, 44, 12], "five_pair_entry_4_8": -36, "five_pair_max": 72},
        "Bare algebra changed",
    )
    require(
        parsed["odd"]
        == {"ranks": [6, 6, 6, 6], "five_pair_entry_4_8": 0, "five_pair_max": 0},
        "Sorted algebra changed",
    )
    require(
        not comm(np.diag(pair[0]), np.diag(pair[1])).any(),
        "Occupation controls fail to commute",
    )
    vectors = {k: [[str(x), str(y)] for x, y in v] for k, v in predicted.items()}
    return {
        "B14_unit_channel_in_B14_chart": vectors,
        "native_operator_checks": parsed,
        "occupation_pair_max": 0,
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=1))
