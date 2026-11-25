import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import material, measure


def test_entire_released_material_and_one_seam_resolvent():
    from model.executable_model import federation

    b = federation(0)
    ports = material.RAILS[0] + material.RAILS[1]
    ix = {x: i for i, x in enumerate(ports)}
    active = set(material.preparation()["active"])
    cuts = set(material.preparation()["cuts"])
    assert active.isdisjoint(cuts)
    assert active.union(cuts) == {e for e in b.seams if set(e).intersection(ports)}
    for k, sign in ((8, -1), (8, 1), (4, -1)):
        h = material.native_matrix(k, sign)
        independent = np.zeros(h.shape)
        for p, j in ix.items():
            if p in material.preparation()["buffers"]:
                independent[j, j] = k * ((1 + sign) ** 2 - 1)
        for a, c in active:
            i, j = ix[a], ix[c]
            independent[i, i] += k + 0.5
            independent[j, j] += k + 0.5
            independent[i, j] = independent[j, i] = -0.5
        np.testing.assert_allclose(independent, h, atol=2e-13)
        n = measure.N
        np.testing.assert_array_equal(h[:n, :n], h[n:, n:])
        assert material.exact_lines(k, sign) != 0
    # Join two arbitrary onsite records and a random symmetric internal
    # dressing on either side through exactly one seam. Several frequencies
    # test the bound independently of the matched material construction.
    rng = np.random.default_rng(495)
    a = rng.normal(size=(6, 6))
    d = rng.normal(size=(7, 7))
    z = 0.4 + 0.2j
    for rank in (1, 2):
        link = np.zeros((6, 7))
        link[0, 5] = 0.7
        if rank == 2:
            link[1, 4] = 0.9
        total = np.block([[a + a.T, link], [link.T, d + d.T]])
        u = np.linalg.inv(z * np.eye(13) - total)[6:, :6]
        assert np.linalg.matrix_rank(u, 1e-10) == rank


def test_filtered_dynamics_local_phase_and_rotations():
    logged = json.loads((Path(__file__).with_name("results.json")).read_text())
    for r in logged["cases"]:
        k, sign = r["kappa"], r["charge"]
        h = material.native_matrix(k, sign)[: measure.N, : measure.N]
        eig = r["spectral"]["level"]
        at = np.eye(measure.N, dtype=complex)[:, 1]
        for t in 2.0 ** np.arange(13) / 8:
            at = (at + np.exp(1j * t * eig) * expm(-1j * t * h) @ at) / 2
        at = expm(-1j * r["flight"] * h) @ at
        s = r["spectral"]
        assert (
            abs(abs(at[-2]) ** 2 - s["filtered_gateway_read"]["measured_mixed_arrival"])
            < 0.0002
        )
        assert s["filter_conditional_line_weight"] > 0.9999
        for path in (r, s):
            d = path.get("orders", path.get("filtered_spatial_orders"))
            assert all(c[1] == 0 for c in d["perpendicular_half_turn_orders"])
            assert all(c[1] == 0 for c in d["one_full_spatial_turn"])
            assert all(c[0] == 0 for c in d["buffer_swap_and_local_phase_orders"])
            assert all(sum(c) == measure.SHOTS for group in d.values() for c in group)
        assert abs(s["independent_energy_banks"][0]["measured_energy"] - eig) < 0.015
        w = s["two_rail_gateway_weight"]
        v = s["in_line_logical_phase"]["one_receiver_gateway"]["plus_minus_other_line"]
        np.testing.assert_allclose(
            np.array(v) / measure.SHOTS,
            [(1 - w) ** 2, w * w, 2 * w * (1 - w)],
            atol=0.001,
        )
        assert s["in_line_logical_phase"]["whole_second_rail"][
            "plus_minus_other_line"
        ] == [0, measure.SHOTS, 0]
        for case in (r["gateway"], s["filtered_gateway_read"]):
            assert (
                case["unnormalized_partial_transpose_min"]
                < -logged["pt_simultaneous_99_radius"]
            )
        np.testing.assert_allclose(
            r["one_cut"]["analytic_singular_squared"][1], 0, atol=1e-20
        )
