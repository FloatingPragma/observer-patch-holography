"""Independent label projection, exact gauge and long-exposure checks."""

from itertools import permutations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import measure as ex


DATA = json.loads(Path(__file__).with_name("results.json").read_text())


def named_projection(states, native):
    # Attach permanently distinct labels; move the moving label with its load.
    # A flip of parity is inferred from where the actual label goes.
    named = [(tuple(w), head) for holes, head in states for w in permutations(holes)]
    index = {word: i for i, word in enumerate(named)}
    lift = np.zeros((len(named), len(named)))
    qplus = np.zeros((len(named), len(states)))
    qminus = qplus.copy()
    for k, (holes, head) in enumerate(states):
        for perm, tagged in enumerate(permutations(holes)):
            i = index[tagged, head]
            qplus[i, k] = 1 / np.sqrt(2)
            qminus[i, k] = (-1) ** perm / np.sqrt(2)
            lift[i, i] = native[k, k]
            for m in np.flatnonzero(native[:, k]):
                if m == k:
                    continue
                dest, newheads = states[m]
                (start,) = set(holes) - set(dest)
                (end,) = set(dest) - set(holes)
                nxt = tuple(end if v == start else v for v in tagged)
                lift[index[nxt, newheads], i] = native[m, k]
    return qplus.T @ lift @ qplus, qminus.T @ lift @ qminus


def test_full_labeled_exchange_and_independent_arrival():
    for odd, arm, parking in [
        ((), "two", True),
        ((1,), "odd", True),
        ((), "closed_spur", False),
    ]:
        states, h, g = ex.generators(2, frozenset(odd), parking)
        a, b = named_projection(states, h)
        np.testing.assert_allclose(h, a, atol=2e-13)
        np.testing.assert_allclose(g, b, atol=2e-13)
        selected = [i for i, (w, _) in enumerate(states) if w == states[0][0]]
        for key, m in [("native", a), ("fermion", b)]:
            u = expm(-1j * m / 32)
            for _ in range(14):
                u = u @ u
            bank = next(
                v for v in DATA["free"][arm]["records"][key] if v["time"] == 512
            )
            assert abs(sum(abs(u[selected, 0]) ** 2) - bank["returns"]) < 3.0e-9
    for arm in DATA["free"].values():
        for banks in arm["records"].values():
            for b in banks:
                assert sum(b["counts"]) == ex.SHOTS
    assert DATA["free"]["two"]["fixed_address_fibers"] == {"1": 82, "2": 2}
    run = DATA["free"]["two"]["records"]
    assert abs(run["native"][-2]["returns"] - run["fermion"][-2]["returns"]) > 0.30
    for run in ("one", "closed_spur"):
        d = DATA["free"][run]["records"]
        np.testing.assert_allclose(
            [r["returns"] for r in d["native"]],
            [r["returns"] for r in d["fermion"]],
            atol=1e-11,
        )


def test_two_arm_operators_and_odd_winding():
    # Inspect each real program before taking the fermionic sign.
    arms = ex.programs()
    for prog in arms:
        ring = [e for e in prog if ex.SPUR not in e]
        assert sorted(ring) == sorted(ex.RING)
    for odd in (frozenset(), frozenset(range(ex.N)), frozenset(range(0, ex.N, 3))):
        first = ex.schedule(arms[0], odd, False, False)
        second = ex.schedule(arms[1], odd, False, False)
        assert first[0] == second[0]
        assert first[0][1] ^ ex.prepared(2)[1] == (1 << ex.N) - 1
        assert first[1] == [1, 0] and second[1] == [0, 1]
        np.testing.assert_allclose(first[2], second[2], atol=2e-13)
        for p in (False, True):
            v = [ex.schedule(a, odd, True, p)[2] for a in arms]
            assert abs(v[0] + v[1]) < 2e-13
    result = ex.winding_qubit(np.random.default_rng(1))
    x = np.array([[0, 1], [1, 0]])
    for name, sign in [("native", 1), ("fermion", -1)]:
        a, b = np.array(result[name]["transfers"])
        np.testing.assert_array_equal(a, sign * b)
        np.testing.assert_array_equal(b, -x / 2)
        assert DATA["winding_qubit"][name]["intrinsic_odd_vs_even"] == [-1.0, -1.0]
    for name, run in DATA["pulse"].items():
        b, d, *refused = run["banks_bright_dark_transfer_refused_source_refused"][
            "real"
        ]
        assert b + d + sum(refused) == ex.SHOTS
        assert (b == 0 and d > 0) if "fermion" in name else (b > 0 and d == 0)
