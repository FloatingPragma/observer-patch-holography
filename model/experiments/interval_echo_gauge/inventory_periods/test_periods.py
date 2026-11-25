from itertools import product
from math import pi
import json
import numpy as np

from .measure import HERE, character, decode, geometry, local_blocks


def test_every_small_inventory_trace_from_native_matrices():
    for m in range(2, 7):
        groups = {}
        for w in product(range(3), repeat=m):
            groups.setdefault(tuple(w.count(i) for i in range(3)), []).append(w)
        for k, words in groups.items():
            dim, a, b = character(*k)
            h, q = local_blocks(words, 0, (0, 1), 8)
            np.testing.assert_allclose(
                [dim, a, b], [len(words), q.trace() / 8, h.trace()], atol=1e-10
            )


def test_all_lattice_and_two_load_index():
    for m in range(2, 96):
        assert character(0, m, 0) == (1, 1, 0)
        assert character(m - 1, 1, 0)[2] == 17
        assert character(m - 1, 0, 1)[1:] == (4, 65)
        assert character(m - 2, 2, 0)[1:] == (m - 1, 17 * (m - 2))
        assert 4 * 17 * (m - 2) - 65 * (m - 1) == 3 * m - 71


def test_complete_readbacks_and_independent_discriminators():
    data = json.loads((HERE / "records.json").read_text())
    ports, edge = geometry()
    assert len(ports) == 4 and 0 in edge
    for r in data["rows"]:
        d, a, b = character(4 - r["ones"] - r["twos"], r["ones"], r["twos"])
        u = decode(r["banks"], d)
        obs = np.linalg.det(u)
        assert abs(obs - complex(*r["det"])) < 2e-12
        if r["control"] in ("native", "reverse"):
            pred = np.exp(-2j * pi * (6 * a + 55 * b) / 59)
            assert abs(obs - pred) < r["joint_radius"]
            assert (abs(obs - 1) > r["joint_radius"]) == (r["ones"] == 1)
        elif r["control"] in ("commutator", "sham"):
            assert abs(obs - 1) < r["joint_radius"]
            if r["control"] == "commutator" and r["ones"]:
                assert r["max_transition"] > 0.49
        else:
            assert (
                abs(obs - np.exp(-2j * pi * (6 * a + 55 * b) / 59)) > r["joint_radius"]
            )
