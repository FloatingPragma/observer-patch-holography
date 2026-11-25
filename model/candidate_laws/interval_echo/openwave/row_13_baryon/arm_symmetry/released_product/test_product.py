from itertools import combinations, permutations, product
import json

import numpy as np

from . import measure
from ...frames import EPS


def test_connectivity_projector_commutes_on_every_occupied_word():
    ref = json.loads((measure.source.HERE / "results.json").read_text())
    arms = np.array(ref["source"]["arm_ports"])
    adj = measure.carrier.cell().adjacency().toarray()
    singlet = np.outer(EPS.ravel(), EPS.ravel())
    for omit in range(4):
        rails = np.delete(arms, omit, axis=1)
        words = list(product(*rails))
        conn = np.array(
            [sum(adj[x, y] for x, y in combinations(w, 2)) >= 2 for w in words]
        )
        assert conn.sum() == (19 if omit >= 2 else 8)
        projector = np.diag(conn.astype(int))
        np.testing.assert_array_equal(projector @ singlet, singlet @ projector)
        np.testing.assert_array_equal(projector @ singlet, (omit >= 2) * singlet)


def test_spectral_fidelity_inequality_for_all_admitted_mixed_states():
    words = list(product(range(3), repeat=3))
    D = np.zeros((27, 27), complex)
    for f in measure.local_bases.bases():
        basis = np.kron(np.kron(f, f), f)
        weights = [
            sum(w[i] == w[j] for i, j in combinations(range(3), 2)) for w in words
        ]
        D += basis.conj().T @ np.diag(weights) @ basis
    eig = np.linalg.eigvalsh(D)
    np.testing.assert_allclose(eig, [0] + [3] * 16 + [6] * 10, atol=1e-12)
    np.testing.assert_allclose(D @ EPS.ravel(), 0, atol=1e-12)
    assert len(list(permutations(range(3)))) == 6
    for w, n in zip(words, weights, strict=True):
        assert n == (len(set(w)) < 3) + 2 * (len(set(w)) == 1)


def test_numerical_ambiguities_allow_both_counts_on_every_uniform_interval():
    spans = [[0, 0, 7], [0, 2, 1], [1, 3, 2], [3, 3, 5]]
    for mask in [[], [0], [2], [1, 3], [0, 1, 2, 3]]:
        all_counts = []
        for pick in product(range(3), range(1, 4), range(1, 4)):
            all_counts.append(
                7 * (0 in mask) + 5 * (3 in mask) + sum(w in mask for w in pick)
            )
        assert measure.event(spans, mask) == (min(all_counts), max(all_counts))


def test_every_panel_from_complete_archived_banks():
    doc = json.loads((measure.HERE / "records.json").read_text())
    assert len(doc["panels"]) == 19
    assert doc["words"] == [list(w) for w in measure.WORDS]
    for r in doc["panels"]:
        score = measure.decode(r["bins"], r["determinant_connected"])
        for k, v in score.items():
            np.testing.assert_allclose(v, r[k], atol=1e-9, rtol=0)
        if r["read"] == "native":
            assert r["fidelity"][0] > 0.98
            if r["omit"] in (0, 1):
                assert r["joint_per_million"] == [0, 0]
            else:
                assert r["connected_singlet"][0] > 0.995
                assert r["joint_per_million"][0] > 32
        else:
            assert r["connected_singlet"][1] < 0.95
        for b in r["bins"]:
            assert all(0 <= a <= z < 29 and n > 0 for a, z, n in b)
            assert sum(n for a, z, n in b) == 2**32
