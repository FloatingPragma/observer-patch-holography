"""Regression checks for the unknown direction, return phase and refusal budget."""

import json

import numpy as np
from scipy.linalg import expm

from model.executable_model.symmetry import rotations
from .material import chart, encode
from .measure import HERE, IDENTITY, closed
from .pulse import LAW, col, primitive


def test_actual_echo_erases_each_fixed_head_and_coherently_exports_both_bit_values():
    ex, ct = LAW.local_term_parts()
    h = ex + 8 * ct
    echo = LAW.echo()
    w = expm((-1j * np.pi / 8) * h)
    for a in (0, 1):
        g = 1 - a
        effects, amplitude = primitive(a)
        for bit in (0, 1):
            start = col(False, (a, g)[bit], a)
            goal = col(False, g, (g, a)[bit])
            np.testing.assert_allclose(w @ start, effects[bit] * goal, atol=3e-13)

        # Use independent matrix exponentials, including the private reference
        # on the actual endpoint, and an independently deleted echo.
        def hold(t):
            return np.diag([np.exp(-8j * t * p[a] ** 2) for p in LAW.SRC] * 11)

        before = hold(3 * np.pi / 16) @ expm((-1j * np.pi / 16) * h)
        after = expm((-1j * 3 * np.pi / 2) * h) @ hold(np.pi / 16)
        state = col(False, g, g)
        goal = col(True, 0, g)
        np.testing.assert_allclose(
            after @ echo @ before @ state, amplitude * goal, atol=4e-12
        )
        assert np.linalg.norm((after @ before @ state)[:9]) ** 2 < 1e-24
        assert abs(np.vdot(goal, after @ echo @ before @ state)) ** 2 > 1 - 1e-12


def test_physical_threefold_return_on_an_unknown_reference_and_real_material_pi_control():
    lab = chart()
    assert len({*lab["storage"], *lab["probes"], *lab["stock"], lab["reference"]}) == 13
    cycle = rotations()[13]
    state = np.random.default_rng(3584).normal(size=(8, 5, 2))
    state = state + 1j * np.roll(state, 1, axis=1)
    state /= np.linalg.norm(state)
    full = closed([cycle] * 3)
    private = closed([cycle] * 3, True)
    observed = np.einsum("ji,isk->jsk", full, state)
    np.testing.assert_allclose(observed, state, atol=1e-11)
    np.testing.assert_allclose(
        private[np.ix_([1, 2, 4], [1, 2, 4])], -np.eye(3), atol=1e-11
    )
    words, _ = encode(IDENTITY)
    assert len(set(w for w, _ in words)) == 8


def test_receipt_keeps_refusals_and_covers_every_perpendicular_pair():
    result = json.loads((HERE / "results.json").read_text())
    angles = result["turns"]
    assert len({tuple(p) for p in angles["perpendicular_pairs"]}) == 15
    assert angles["max_half_turn_order_operator_error"] < 1e-11
    delivery = result["free_delivery"]
    assert delivery["conditional_PT"] < -0.46 and delivery["Choi_change_norm"] < 1e-16
    normal = delivery["boundary_full_turn"]
    control = delivery["boundary_full_turn_with_material_private_phase"]
    for bank in (normal, control):
        assert sum(bank["counts"]) == 2**38
        assert abs(sum(bank["probabilities"]) - 1) < 1e-12
    assert normal["dark_fraction_99"][1] < 0.001
    assert control["dark_fraction_99"][0] > 0.999

    filtered = result["filtered_delivery"]
    assert filtered["conditional_PT"] < -0.49 and filtered["Choi_change_norm"] < 1e-18
    for key in ("boundary_full_turn", "boundary_full_turn_with_material_private_phase"):
        assert sum(filtered[key]["counts"]) == 2**38
        assert len(filtered[key]["counts"]) == 6
    assert filtered["boundary_full_turn"]["dark_fraction_99"][1] < 0.005
    assert (
        filtered["boundary_full_turn_with_material_private_phase"]["dark_fraction_99"][
            0
        ]
        > 0.995
    )
