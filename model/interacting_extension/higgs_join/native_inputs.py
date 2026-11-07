"""Independent exact replay of the retained OPH current and Higgs matrices.

The public fixtures are read as data. Every matrix identity used here is
recomputed; importing the simulator's service stack is unnecessary.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as s
from .._source import verify_sources

PIN = "3cf3c14810671a16e61b7915f69079a1773a12a2"
ROOT = Path(__file__).resolve().parents[3]


SOURCE_NAMES = (
    "model/observer_particles/results/representation_identification_operators.json",
    "model/observer_particles/results/representation_identification.json",
    "model/observer_particles/results/unified_yukawa_operators.json",
    "model/observer_particles/results/unified_yukawa_action.json",
    "model/common_source/native_gauge_object_results.json",
    "model/common_source/native_family_protocol_results.json",
    "model/common_source/native_family_protocol_results.npz",
    "model/common_source/native_repair_mass_results.npz",
    "model/common_source/native_family_record_memory.py",
    "model/common_source/native_gauge_family_records.py",
    "model/common_source/native_mass_dressing.py",
    "model/common_source/native_mass_bath.py",
    "model/common_source/slow_geometry.py",
    "model/native_quantum_gauss/algebra.py",
)


def load_sparse(record):
    return s.SparseMatrix(*record["shape"],
                          {(i, j): s.sympify(v) for i, j, v in record["entries"]})


def clean(value):
    if isinstance(value, s.MatrixBase):
        return value.applyfunc(lambda x: s.simplify(s.expand(x)))
    return s.simplify(s.expand(value))


def assert_zero(value, label):
    values = list(value) if isinstance(value, s.MatrixBase) else [value]
    if any(clean(x) != 0 for x in values):
        raise ArithmeticError(label)


def exterior_derivative(generator, subsets):
    indices = {tuple(row): j for j, row in enumerate(subsets)}
    entries = {}
    for col, subset in enumerate(subsets):
        for position, old in enumerate(subset):
            rest = [x for x in subset if x != old]
            for new in range(5):
                coefficient = generator[new, old]
                if coefficient == 0 or new in rest:
                    continue
                destination = tuple(sorted([*rest, new]))
                insertion = sum(x < new for x in rest)
                key = indices[destination], col
                entries[key] = entries.get(key, 0) + (-1) ** (position + insertion) * coefficient
    return clean(s.SparseMatrix(len(subsets), len(subsets), entries))


def native_data(repo):
    repo = Path(repo)
    base = repo / "model/observer_particles/results"
    packet = json.loads((base / "representation_identification_operators.json").read_text())
    vertices = json.loads((base / "unified_yukawa_operators.json").read_text())["vertices"]
    for filename, receiptname in (
        ("representation_identification_operators.json", "representation_identification.json"),
        ("unified_yukawa_operators.json", "unified_yukawa_action.json"),
    ):
        receipt = json.loads((base / receiptname).read_text())
        digest = hashlib.sha256((base / filename).read_bytes()).hexdigest()
        if digest != receipt["artifacts"][filename]:
            raise ArithmeticError("retained artifact checksum: " + filename)
    subsets = tuple(tuple(packet["basis_subsets"][i]) for i in packet["even_indices"])
    if len(subsets) != 15 or {len(row) for row in subsets} != {2, 4}:
        raise ArithmeticError("nonvacuum even exterior carrier")
    one_body = tuple(load_sparse(x) for x in packet["spin_lifted_one_body_generators"])
    currents = tuple(exterior_derivative(x, subsets) for x in one_body)
    operators = {k: load_sparse(v) for k, v in packet["even_identification_operators"].items()}
    for g in currents:
        assert_zero(g + g.H, "anti-Hermitian port current")
    if s.Matrix.hstack(*(g.reshape(225, 1) for g in currents)).rank() != 12:
        raise ArithmeticError("complete current rank")
    q, c3, c2 = (operators[k] for k in ("q", "color_C2", "weak_C2"))
    for op in (q, c3, c2):
        assert_zero(op - s.diag(*op.diagonal()), "diagonal charge/Casimir readout")
    blocks = []
    for charge in sorted(set(q.diagonal()), key=float):
        ids = tuple(i for i in range(15) if q[i, i] == charge)
        outside = tuple(i for i in range(15) if i not in ids)
        if any(g.extract(ids, outside) != s.zeros(len(ids), len(outside)) for g in currents):
            raise ArithmeticError("current mixes native irreducible blocks")
        color = clean(s.trace(c3.extract(ids, ids)) / len(ids))
        assert_zero(c3.extract(ids, ids) - color * s.eye(len(ids)), "scalar color Casimir")
        blocks.append({"q6Y": str(charge), "indices": ids, "dimension": len(ids),
                       "color_casimir": color,
                       "currents": tuple(g.extract(ids, ids) for g in currents)})
    casimir = clean(2 * c3 + c2 / 2 + q * q / 12)
    for g in currents:
        assert_zero(casimir * g - g * casimir, "native trace Casimir centrality")
    channel_names = {(-4, 1): "up", (1, 2): "down", (-3, 6): "charged_lepton"}
    channels = {}
    ward_checks = 0
    for vertex in vertices:
        origin = tuple(map(int, vertex["constituent_q6Y"]))
        name = channel_names[origin]
        comps = tuple(load_sparse(x) for x in vertex["coefficient_matrices"])
        conjugate = bool(vertex["uses_conjugate_of_same_Higgs_doublet"])
        for g, g5 in zip(currents, one_body, strict=True):
            h = -g5[3:, 3:].T if conjugate else g5[3:, 3:]
            for k in range(2):
                residual = g.T * comps[k] + comps[k] * g
                residual += sum((h[j, k] * comps[j] for j in range(2)), s.zeros(15))
                assert_zero(residual, "all-current Higgs Ward identity")
                ward_checks += 1
        support = {i for matrix in comps for (i, _), value in matrix.todok().items() if value != 0}
        left = tuple(i for i in sorted(support) if c2[i, i] != 0)
        right = tuple(i for i in sorted(support) if c2[i, i] == 0)
        cl, cr = casimir[left[0], left[0]], casimir[right[0], right[0]]
        if any(casimir[i, i] != cl for i in left) or any(casimir[i, i] != cr for i in right):
            raise ArithmeticError("uniform support Casimir")
        channels[name] = {"origin": origin, "components": comps,
                          "left": left, "right": right, "left_casimir": cl,
                          "right_casimir": cr}
    return {"currents": currents, "operators": operators, "blocks": blocks,
            "channels": channels, "higgs_ward_checks": ward_checks,
            "casimir": casimir, "subsets": subsets}


def manifest(repo):
    repo = Path(repo)
    verify_sources(repo, SOURCE_NAMES)
    hashes = {name: hashlib.sha256((repo / name).read_bytes()).hexdigest()
              for name in SOURCE_NAMES}
    return {"pinned_source_commit": PIN,
            "source_sha256": hashes,
            "source_bytes_match_pinned_commit": True,
            "empirical_mass_or_mixing_inputs": []}

