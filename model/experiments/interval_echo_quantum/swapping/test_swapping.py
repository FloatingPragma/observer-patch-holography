from itertools import product
import json

import numpy as np
import pytest

from . import circuit, compiler, measure

LINE, ROUTES = circuit.apparatus()


def test_native_gates_and_routes():
    circuit.check_readers()
    ex, co = circuit.LAW.local_term_parts()
    for w in circuit.LAW.SRC:
        v = np.eye(circuit.LAW.DIM)[circuit.LAW.SRC.index(w)]
        swapped = np.eye(circuit.LAW.DIM)[circuit.LAW.SRC.index(w[::-1])]
        np.testing.assert_allclose(ex @ v, (v - swapped) / 2, atol=1e-14)
        np.testing.assert_allclose(co @ v, (w[0] - w[1]) ** 2 * v, atol=1e-14)
    board = circuit.carrier.federation(0)
    protected = [(12 * s, 12 * s + 1) for s in (LINE[1], LINE[3])]
    for owner in LINE:
        assert set(compiler.guards(owner * 12)).isdisjoint(set(circuit.rails(owner)[0]))
        np.testing.assert_allclose(circuit.settings_probability(owner), [0.5, 0.5])
    for route in ROUTES:
        perm = list(range(board.n))
        for a, b in route:
            e = tuple(sorted((a, b)))
            assert e in board.seams and e not in protected
            perm[a], perm[b] = perm[b], perm[a]
        assert len([i for i, p in enumerate(perm) if i != p]) == 2
    for edge in [e for route in ROUTES[4:6] for e in route]:
        assert set(edge).isdisjoint(circuit.rails(LINE[2])[0])
    # Each physically decoded Bell column yields both commuting syndrome bits.
    lr, rr = circuit.rails(LINE[2])
    words = [
        tuple(sorted((lr[1 - a], rr[1 - b]))) for a, b in product(range(2), repeat=2)
    ]
    decode = np.array(
        [
            [circuit.bob_analyzer({w: 1}, LINE[2]).get(v, 0) for w in words]
            for v in words
        ]
    )
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    for i, row in enumerate(decode):
        expected = (
            (np.eye(4) + (-1) ** (i // 2) * np.kron(x, x))
            @ (np.eye(4) + (-1) ** (i % 2) * np.kron(z, z))
            / 4
        )
        np.testing.assert_allclose(np.outer(row.conj(), row), expected, atol=3e-14)
    assert max(compiler.phase_bound().values()) < 1e-13


def closed_form(mode, menu):
    """Independent distribution from parity projectors, including each pointer."""
    prob = np.zeros((2, 2, 2, 2, 3, 3))
    ax, cz = circuit.MENU[menu]
    for b, d, x, z in product(range(2), repeat=4):
        a, c = ax[x], cz[z]
        if mode == "cut_alice":
            prob[b, d, x, z, 2, :2] = 1 / 8
            continue
        for s, t in product(range(2), repeat=2):
            if mode == "separable_bob":
                p = (1 + (-1) ** (s + b) * np.cos(a)) * (
                    1 + (-1) ** (t + d) * np.cos(c)
                )
            else:
                corr = (-1) ** d * np.cos(a) * np.cos(c)
                if mode == "native":
                    corr += (-1) ** b * np.sin(a) * np.sin(c)
                p = 1 + (-1) ** (s + t) * corr
            prob[b, d, x, z, s, t] = p / 16
    return prob


@pytest.mark.parametrize("mode,menu", product(circuit.MODES, circuit.MENU))
def test_every_pointer(mode, menu):
    target = closed_form(mode, menu)
    for phase in (circuit.phase, compiler.native_phase):
        actual = circuit.distributions(LINE, ROUTES, mode, menu, phase)
        np.testing.assert_allclose(actual, target, atol=1e-13, rtol=0)
        # Both unconditioned outer detectors are blind to the remote setting,
        # also when sliced by the earlier central syndrome.
        np.testing.assert_allclose(np.diff(actual.sum(axis=4), axis=2), 0, atol=1e-13)
        np.testing.assert_allclose(np.diff(actual.sum(axis=5), axis=3), 0, atol=1e-13)


def test_replay_and_record_only_witnesses():
    saved = json.loads((measure.HERE / "records.json").read_text())
    assert measure.run() == saved
    for mode, menu in product(circuit.MODES, circuit.MENU):
        bank = saved["banks"][f"{mode}/{menu}"]
        s = measure.score(bank["counts"], mode, menu)
        assert s == bank["score"]
        assert sum(map(sum, s["settings"])) == measure.SHOTS
        assert s["outer_quiets"] == (
            [measure.SHOTS, 0] if mode == "cut_alice" else [0, 0]
        )
        assert all(
            v["max_remote_rate_difference"] < v["joint_allowance"]
            for v in s["no_signalling_C_A"]
        )
        if menu == "chsh":
            assert abs(s["erased_S"]) < 0.02
            if mode == "native":
                assert s["lower_99"] > 2.8 and np.min(s["conditional_lower_99"]) > 2.78
            else:
                assert s["S"] + (s["S"] - s["lower_99"]) < 2
        elif mode == "native":
            assert s["lower_99"] > 1.4
        else:
            assert s["upper_99"] < 1
    counts = np.array(saved["banks"]["native/chsh"]["counts"])
    counts[0, 0, 0, 0, 0, 0] -= 1
    with pytest.raises(ValueError, match="complete"):
        measure.score(counts, "native", "chsh")
