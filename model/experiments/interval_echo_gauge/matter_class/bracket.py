"""The P1 source sample and the frozen graded-module ballot."""

from fractions import Fraction
from functools import cache
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
FREEZE = ROOT / "code/b15_matter_freeze/matter_class_freeze_v1.json"


@cache
def frozen():
    data = json.loads(FREEZE.read_text())
    migrations = json.loads(
        Path(__file__).with_name("source_metadata.json").read_text()
    )
    for source in data["provenance"]["pinned_files"]:
        name = source["path"]
        raw = (ROOT / name).read_bytes()
        digest = sha256(raw).hexdigest()
        if digest == source["sha256"]:
            continue
        change = migrations[name]
        if digest != change["current_sha256"]:
            raise ValueError("Unrecognized B15 source: " + name)
        old = json.loads(raw)
        for path, value in change["historical_metadata"]:
            obj = old
            for key in path[:-1]:
                obj = obj[key]
            if value is None:
                del obj[path[-1]]
            else:
                obj[path[-1]] = value
        recovered = (
            json.dumps(old, indent=change["indent"], sort_keys=True) + "\n"
        ).encode()
        if sha256(recovered).hexdigest() != source["sha256"]:
            raise ValueError("Historical B15 source mismatch: " + name)
    return data


def multiply(left, right):
    """Matrix product in Z[sqrt(5)] with arbitrary-precision integers."""
    a, b = left
    c, d = right
    return np.stack((a @ c + 5 * (b @ d), a @ d + b @ c))


@cache
def sample():
    frozen()
    cert = json.loads(
        (ROOT / "code/b14_jacobi/b14_compact_locus.certificate.json").read_text()
    )
    point = next(
        x
        for x in cert["samples"]
        if x["family"] == "PLANE" and x["parameters"] == ["1", "0"]
    )
    basis = json.loads(
        (
            ROOT
            / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
        ).read_text()
    )["basis"]
    # t[:,i,:,j] holds 60 times the column [e_i,e_j] in the port basis.
    t = np.zeros((2, 12, 12, 12), dtype=object)
    for index, a, ad, b, bd in point["x_vector"]:
        for out, i, j, n, d in basis[index]["entries"]:
            for part, num, denom in [(0, a, ad), (1, b, bd)]:
                val, rem = divmod(60 * num * n, denom * d)
                assert rem == 0
                t[part, i, out, j] += val
                t[part, j, out, i] -= val
    assert np.array_equal(t, -t.swapaxes(2, 3))
    for i, j in combinations(range(12), 2):
        comm = multiply(t[:, i], t[:, j]) - multiply(t[:, j], t[:, i])
        image = multiply(t[:, i, :, j], t.reshape(2, 12, -1)).reshape(2, 12, 12)
        assert np.array_equal(comm, image)
    assert not t.sum(axis=1).any()
    flat = t.reshape(2, 12, 144)
    gram = multiply(flat, flat.swapaxes(1, 2))
    trace = np.trace(gram, axis1=1, axis2=2).reshape(2, 1, 1)
    assert trace.any()
    # 3*G^2=trace(G)*G for this positive Gram matrix. Its rank is three.
    assert np.array_equal(
        3 * multiply(gram, gram),
        multiply(trace, gram.reshape(2, 1, -1)).reshape(gram.shape),
    )
    # The derived image is the orthogonal complement of the zero-action center.
    tensor = (t[0].astype(float) + np.sqrt(5) * t[1].astype(float)) / 60
    frame = []
    pivots = []
    for i, j in combinations(range(12), 2):
        col = tensor[i, :, j].copy()
        for v in frame:
            col -= v * (v @ col)
        if np.linalg.norm(col) > 1e-9:
            frame.append(col / np.linalg.norm(col))
            pivots.append([i, j])
        if len(frame) == 3:
            break
    q = np.array(frame).T
    adjoint = q.T @ tensor @ q
    np.testing.assert_allclose(tensor @ q, q @ adjoint, atol=2e-13)
    return tensor, q, adjoint, pivots


def grade_candidates(items, cap, validator):
    for mults in validator.multisets(items, cap):
        present = [(i, n) for i, n in enumerate(mults) if n]
        if not present:
            continue
        totals = [n for _, n in present]
        for left in product(*(range(n + 1) for n in totals)):
            right = tuple(n - j for n, j in zip(totals, left))
            if left <= right:
                yield [i for i, _ in present], left, right


def normalized_turn_trace(label):
    """Trace divided by dimension of one full turn around the fixed compact root."""
    if "A1(" in label:
        # On P2 use its first factor; on F2/G2 use their A1 factor.
        n = int(label.split("A1(", 1)[1].split(")")[0])
        return (-1) ** (n - 1)
    if label == "triv":
        return 1
    # In pure A2 the root SU(2) acts on coordinates 1,2. U(3) top
    # pattern (p+q,q,0) restricts to U(2) rows (i,j) interlacing it.
    p, q = map(int, label[3:-1].split(","))
    trace = sum(
        (-1) ** (i - j) * (i - j + 1) for i in range(q, p + q + 1) for j in range(q + 1)
    )
    return Fraction(trace, (p + 1) * (q + 1) * (p + q + 2) // 2)


@cache
def ballot():
    d = frozen()
    spec = importlib.util.spec_from_file_location(
        "b15_freeze", FREEZE.with_name("validate_matter_class_freeze.py")
    )
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    by_type = {}
    cap = d["budget"]["dimension_cap"]
    for kind, entry in d["enumeration"]["per_bracket_type"].items():
        items = validator.items_for_type(kind, cap)
        assert [x[0] for x in items] == [x["label"] for x in entry["admissible_irreps"]]
        count = 0
        matches = []
        for ids, left, right in grade_candidates(items, cap, validator):
            count += 1
            dims = [
                sum(k * items[i][1] for i, k in zip(ids, side))
                for side in [left, right]
            ]
            if sorted(dims) != [1, 2]:
                continue
            # For grades of dimension <=2 each pure A2 summand must be trivial.
            # The reference fixes only a relative sign; the opposite central
            # action would give the same channel and must also be allowed.
            signatures = sorted(
                (
                    dim,
                    sum(
                        k * items[i][1] * normalized_turn_trace(items[i][0])
                        for i, k in zip(ids, side)
                    ),
                )
                for dim, side in zip(dims, [left, right])
            )
            if signatures in ([(1, 1), (2, -2)], [(1, -1), (2, 2)]):
                matches.append(
                    [
                        [[items[i][0], k] for i, k in zip(ids, side) if k]
                        for side in [left, right]
                    ]
                )
        assert count == entry["counts"]["graded_up_to_flip"]
        by_type[kind] = dict(candidates=count, matched=matches)
    return {s["id"]: by_type[s["bracket_type"]] for s in d["enumeration"]["strata"]}
