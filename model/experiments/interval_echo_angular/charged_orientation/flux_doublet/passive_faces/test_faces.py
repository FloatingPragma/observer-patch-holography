"""Full-shell native support, measured face flux, and monitored transport."""

import json
from math import pi

import numpy as np
from ..measure import faces, links
from .measure import HERE, twice_native_crossings
from .monitored import material_face, pair


def test_native_two_hops_and_arbitrary_entangled_face_records():
    blocks = twice_native_crossings()
    rng = np.random.default_rng(1122)
    x = rng.normal(size=(11, 11, 11, 3)) + 1j * rng.normal(size=(11, 11, 11, 3))
    x /= np.linalg.norm(x)
    target = np.zeros_like(x)
    target[0, 0, 0] = x[0, 0, 0]
    for sign in (-1, 1):
        ops = [blocks[sign, k] / 2 for k in (True, False, True)]

        def lap(p):
            return np.einsum("ja,kb,lc,abcr->jklr", *ops, p, optimize=True)

        np.testing.assert_allclose(lap(lap(x)), target, atol=1e-13)


def test_stored_head_response_and_native_escape_from_its_projected_doublet():
    for sign in (-1, 1):
        chart = links(1, root=5)
        for face in faces():
            a, n, b = material_face(face, chart, sign)
            phase = np.exp(-1j * sign * pi / 10)
            np.testing.assert_allclose(
                [a, n, b], [phase / 8, 1 / 8, phase**2 / 64], atol=1e-13
            )
        for phase in [0.0, 0.304, 1.2]:
            h, p, z, leak = pair(sign, phase)
            out = (np.eye(4) - p @ p.conj().T) @ h @ p
            np.testing.assert_allclose(out.conj().T @ out, leak, atol=1e-13)
            np.testing.assert_allclose(leak, -4 * z, atol=1e-13)


def test_both_signs_twenty_faces_control_refusals_and_finite_monitor():
    r = json.loads((HERE / "results.json").read_text())
    assert [
        a["flux"] for a in r["classes"] if a["negative_half_turn_pairs"] == 15
    ] == list(range(1, 20, 2))
    assert len(r["records"]) == 120
    for row in r["records"]:
        for pair_of_settings in row["sampled_X_Y_by_lap"]:
            assert all(
                sum(counts) == r["shots_per_axis"] for counts in pair_of_settings
            )
    m = json.loads((HERE / "monitored.json").read_text())
    assert len(m["material_face_counts"]) == 40
    for sign, _, inferred, *banks in m["material_face_counts"]:
        assert inferred == -sign
        assert all(sum(c) == m["starts_per_port_bank"] for xy in banks for c in xy)
    for row in m["records"]:
        assert row["multiplicities"] == [2, 4, 6]
        assert max(row["energy_levels"]) < 0
        final = row["appointments"][-1]
        assert (
            sum(final["sharp_port_0_counts_with_refusals"]) == m["starts_per_port_bank"]
        )
        assert final["lower_doublet_yield_range"][0] > 0.98
        assert final["transfer_operator_error"] < 0.0091
