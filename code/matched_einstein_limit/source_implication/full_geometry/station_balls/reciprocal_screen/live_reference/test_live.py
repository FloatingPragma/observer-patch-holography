import json

import numpy as np
import sympy as sp

from . import acquire, apparatus


def test_exact_detector_and_native_port_order():
    identity = sp.eye(4)
    psi = sp.Matrix([0, 1, 1, 0])
    phi = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2)
    block = psi * psi.T / 2
    z = sp.diag(1, 1, -1, -1)
    had = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
    yplus = sp.Matrix([1, sp.I]) / sp.sqrt(2)
    phase = identity + (sp.I - 1) * block
    detector = (z * phase * z * phase) * sp.kronecker_product(had, sp.eye(2))
    both = sp.kronecker_product(yplus, yplus)
    assert sp.simplify(detector.H * both * both.H * detector - phi * phi.T) == sp.zeros(
        4
    )
    for odd in (0, 1):
        # Overall phase cancels from two identical marked pulses at dwell pi.
        g = identity + ((-1) ** (8 + odd) - 1) * block
        assert z * g * z * g == (sp.diag(1, -1, -1, 1) if odd else identity)
        tag = apparatus.parity_read(odd)
        for source_bit, op in enumerate(tag):
            np.testing.assert_allclose(
                op.conj().T @ op, (source_bit == odd) * np.eye(2), atol=3e-14
            )
    for mount in apparatus.geometry():
        np.testing.assert_allclose(apparatus.shuttle(mount), np.eye(2), atol=3e-13)


def test_joint_coherence_survives_and_initial_labels_must_be_independent():
    effects = acquire.pair_effects()
    # A genuinely entangled coin state in this algebraic check ensures the
    # claimed channel isometries do not implicitly refactor either copy.
    u = np.array([1, 1j, 0.25, 1]) / np.sqrt(3.0625)
    v = np.array([0, 1, 1, 0]) / np.sqrt(2)
    x = np.stack([0.3 * np.outer(u, u.conj()), 0.7 * np.outer(v, v.conj())])
    y = np.stack([0.4 * np.outer(v, v.conj()), 0.6 * np.outer(u, u.conj())])
    k = np.broadcast_to(np.eye(4), (2, 1, 1, 4, 4))
    read = acquire.instrument(k, x, y, effects)[0, 0]
    assert abs(read[1] - sum(np.trace(a.T @ b).real for a, b in zip(y, x)) / 4) < 3e-14
    c = acquire.controls(effects)
    for name, expected in zip(
        (
            "independent_shared_marginals",
            "copied_lottery_label",
            "off_diagonal_xx_live",
            "off_diagonal_xx_marginals_only",
            "pure_left_y_zero_flight",
        ),
        (1 / 8, 1 / 4, 1 / 8, 1 / 16, 0),
        strict=True,
    ):
        assert abs(c[name] - expected) < 3e-14
    assert abs(apparatus.report()["omitted_pulses_phi_acceptance"] - 0.5) < 1e-14


def test_actual_both_directions_with_live_and_stale_reference():
    stored = json.loads((acquire.HERE / "results.json").read_text())
    replay = acquire.run()
    parent = json.loads((acquire.HERE.parent / "results.json").read_text())
    assert stored["independent_flight_and_reference_pairs"] == 63 * acquire.STARTS
    for key, c in stored["zero_flight_standard"].items():
        assert (
            c["sampled_99"][0]
            < replay["zero_flight_standard"][key]["population"]
            < c["sampled_99"][1]
        )
    for q, sets in stored["charges"].items():
        for key, row in sets.items():
            fresh = replay["charges"][q][key]
            if key == "held_stale":
                assert (
                    fresh["population_maximum_asymmetry"]
                    >= row["count_directional_asymmetry_lower"]
                    > 1e-7
                )
            else:
                assert fresh["population_maximum_asymmetry"] < 4e-14
            previous = {
                "equilibrium": "matched",
                "held_live": "held_recalibrated",
                "held_stale": "held_fixed",
            }[key]
            for side in ("forward", "reverse"):
                counts = np.asarray(row[side]["counts"])
                assert (
                    counts.shape == (5, 26)
                    and np.all(counts >= 0)
                    and np.all(counts.sum(1) == acquire.STARTS)
                )
                lo, hi = row[side]["count_area_gain_interval"]
                assert lo < fresh[side]["population_area_gain"] < hi
                assert lo < parent["charges"][q][previous][side]["determinant"] < hi
