"""Independent current, geometry and count checks for the face-pulse experiments."""

from itertools import combinations
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm
from model.executable_model import carrier

HERE = Path(__file__).parent


def test_direction_analyser_from_microscopic_continuity():
    data = json.loads((HERE / "results.json").read_text())
    first = data["labs"][0]
    cell = carrier.cell()
    words = list(combinations(range(12), 2))
    index = {w: i for i, w in enumerate(words)}
    n = np.array([[p in w for w in words] for p in range(12)], dtype=float)
    native = {}
    total = np.zeros((len(words), len(words)))
    for a, b in cell.seams:
        swap = np.eye(len(words))[
            [index[tuple(sorted({a: b, b: a}.get(v, v) for v in w))] for w in words]
        ]
        h = (np.eye(len(words)) - swap) / 2 - 16 * np.diag(n[a] * n[b])
        native[(a, b)] = h
        total += h

    def loop(face):
        current = np.zeros_like(total, dtype=complex)
        for a, b in zip(face, (*face[1:], face[0])):
            h = native[min(a, b), max(a, b)]
            current += 1j * (h * n[b][None, :] - n[b][:, None] * h)
        return current

    ground = eigh(total)[1][:, 0]
    source = loop(first["source"])
    plus = expm(-1j * source) @ ground
    minus = expm(1j * source) @ ground
    batch = first["experiments"][2]
    face = next(r for r in batch["reads"] if r["face"] == [2, 3, 6])
    detector = loop(face["face"])
    energies, v = eigh(detector)
    for arm, state in enumerate((plus, minus)):
        psi = expm(-4j * total) @ state
        weights = abs(v.conj().T @ psi) ** 2
        actual = [
            sum(weights[energies > 0.1]),
            sum(weights[abs(energies) < 0.1]),
            sum(weights[energies < -0.1]),
        ]
        np.testing.assert_allclose(
            actual, face["probabilities"][data["times"].index(4)][arm], atol=2e-11
        )


def test_physical_reflection_and_native_dark_states():
    data = json.loads((HERE / "results.json").read_text())
    lab = data["labs"][1]
    p = lab["reflection"]
    a, b, c = lab["source"]
    assert (p[a], p[b], p[c]) in ((c, b, a), (a, c, b), (b, a, c))
    assert all(p[v] == v for v in lab["gateway_ports"])
    edges = set(carrier.cell().seams)
    assert {tuple(sorted((p[a], p[b]))) for a, b in edges} == edges
    for g in data["labs"]:
        for arm in g["experiments"][:2]:
            assert arm["source_weight"] < 1e-20
    for r in lab["experiments"][2]["reads"]:
        if r["face"] != lab["source"]:
            prob = np.array(r["probabilities"])
            np.testing.assert_allclose(prob[:, 0, :], prob[:, 1, :], atol=2e-11)


def test_all_launch_count_denominators_and_exact_controls():
    data = json.loads((HERE / "flight.json").read_text())
    assert list(data["times"]) == [20, 23, 24, 27, 28]
    for arm in data["panels"]:
        for r in arm["reads"]:
            counts = np.array(r["counts"])
            assert np.all(counts.sum(axis=1) == data["shots"])
            diff = np.diff(counts[:, :9].reshape(2, 3, 3).sum(axis=1), axis=0)[0]
            observed = (diff[2] - diff[0]) / data["shots"]
            lo, hi = r["exact_binomial_contrast_interval"]
            assert lo <= observed <= hi
            assert np.all(counts[:, [0, 2]] == 0)
            if arm["bridge_cut"] or arm["load"] == 1:
                assert lo <= 0 <= hi
            else:
                assert not lo <= 0 <= hi
    herald = data["panels"][1]
    assert 0.199 < herald["success"] < 0.201
