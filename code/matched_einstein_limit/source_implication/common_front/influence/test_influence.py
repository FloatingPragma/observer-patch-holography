"""Full-board moment reductions, degree-two channel response and physical read records."""

import json

import numpy as np
from scipy import sparse

from . import measure
from ...action_response.four_ticks import parse_interval
from ...directional_stress.bridge import tangents


def compress(h, heads, x, n):
    launch = np.zeros((h.shape[0], 4))
    launch[4 * x : 4 * x + 4] = np.eye(4)
    for _ in range(n):
        launch = h @ launch
    return np.einsum("pij,ji->p", launch.reshape(-1, 4, 4), heads)


def test_full_connected_memory_laurent_coefficient():
    states = tangents()
    s = states["reference"]
    native = measure.action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    receipt = json.loads((measure.HERE / "results.json").read_text())
    for q in (-1, 1):
        hs = [native.generator(board, route, q, odd) / 2 for odd in (False, True)]

        def coefficient(x):
            c = [
                [sum(compress(h, b, p, n) for h, b in zip(hs, x)) for n in (1, 2, 3)]
                for p in (0, 39)
            ]
            m1 = sparse.csr_matrix((board.n, board.n), dtype=complex)
            for h, b in zip(hs, x):
                co = h.tocoo()
                m1 += sparse.csr_matrix(
                    (co.data * b[co.col % 4, co.row % 4], (co.row // 4, co.col // 4)),
                    shape=m1.shape,
                )
            np.testing.assert_allclose(m1[:, 0].toarray().ravel(), c[0][0])
            matrix = (
                c[0][2][39]
                - np.vdot(c[1][0], c[0][1])
                - np.vdot(c[1][1], c[0][0])
                + np.vdot(c[1][0], m1 @ c[0][0])
            )
            return matrix

        vals = {"reference": coefficient(s)}
        for name in ("pairing", "ordinary"):
            h = 1e-6
            d = states[name]
            vals[name] = (coefficient(s + h * d) - coefficient(s - h * d)) / (2 * h)
        for k, w in vals.items():
            for label, value in (("real", w.real), ("imaginary", w.imag)):
                assert value in parse_interval(
                    receipt["charges"][str(q)]["memory_" + label][k]
                )


def test_marginal_invariance_and_source_coherence():
    native = measure.action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    d = tangents()["pairing"]
    for q in (-1, 1):
        mats = [native.generator(board, route, q, o) / 2 for o in (False, True)]
        # Check the partial trace of H and H^2 for EVERY spatial matrix
        # element, including the boundary rows of the connected board.
        for n in (1, 2):
            delta = sparse.csr_matrix((board.n, board.n), dtype=complex)
            for h, b in zip(mats, d):
                a = (h if n == 1 else h @ h).tocoo()
                delta += sparse.csr_matrix(
                    (a.data * b[a.col % 4, a.row % 4], (a.row // 4, a.col // 4)),
                    shape=delta.shape,
                )
            assert max(abs(delta).data, default=0) < 1e-10
        pilot = measure.effects(board, route, q)["dark"]
        exact = sum(
            np.trace(np.array(a, dtype=complex) @ b).real for a, b in zip(pilot, d)
        )
        assert -2.23 < exact < -2.19


def test_saved_counts_and_relative_transfer():
    data = json.loads((measure.HERE / "results.json").read_text())
    assert measure.run() == data
    traces = json.loads((measure.HERE / "counts.json").read_text())
    n = traces["shots_per_setting"]
    for q, row in traces["charges"].items():
        for time, bank in row["counts"].items():
            for controls in bank.values():
                for setting, obs in controls.items():
                    if setting.startswith("dark"):
                        assert sum(obs) == n
                    else:
                        assert np.array_equal(
                            np.sum(np.array(obs, dtype=object), axis=1),
                            np.array([n] * 4, dtype=object),
                        )
        exact = data["charges"][q]
        assert (
            parse_interval(exact["local_sharp_early_arrival"]["coherent_fraction"]).b
            < 1e-6
        )
        assert (
            parse_interval(exact["memory_real"]["collar_over_ordinary"]).b
            < parse_interval(
                exact["dark_incoherent_quadratic"]["collar_over_ordinary"]
            ).a
        )
        assert (
            parse_interval(exact["dark_incoherent_quadratic"]["collar_over_ordinary"]).b
            < parse_interval(data["entropy_collar_over_ordinary"]).a
        )
        # The count bounds apply to three finite dwells; the limiting and
        # finite quantities agree in sign on both changes, for each readout.
        read = row["three_dwell_nominal_intervals_99pct"]
        for change in ("ordinary", "pairing"):
            assert (
                parse_interval(read["memory"][change]).a > 0
                and parse_interval(read["dark"][change]).b < 0
            )


def test_finite_retarded_transmission_and_exterior_support():
    from . import retarded

    data = retarded.run()
    reference = json.loads((measure.HERE / "retarded.json").read_text())
    assert data["exterior_ports"] == 956 and data["marked_ports"] == [0, 1, 38, 39]
    for q, charge in data["charges"].items():
        assert charge["exterior_receiver"] == 5 and len(charge["settings"]) == 5
        for key, row in charge["settings"].items():
            native = reference["charges"][q]["settings"][key]
            assert row["maximum_sparse_residual_and_schur_difference"] < 1e-12
            for name in ("ordinary", "pairing"):
                np.testing.assert_allclose(
                    row[name]["connected_39_from_0"],
                    native[name]["connected_39_from_0"],
                    atol=1e-12,
                )
                assert abs(complex(*row[name]["blank_return_slope"])) > 1e-11
        finite = json.loads((measure.HERE / "finite.json").read_text())["charges"][q]
        cert = finite["derivatives"]
        blank = finite["exterior_5_return_derivative"]
        assert parse_interval(blank["ordinary"][0]).a > 0
        if q == "-1":
            assert parse_interval(blank["pairing"][0]).a > 0
        else:
            assert parse_interval(blank["pairing"][0]).b < 0
        for k in (0, 1):
            for arm in ("ordinary", "pairing"):
                assert charge["settings"]["0,1"][arm]["blank_return_slope"][
                    k
                ] in parse_interval(blank[arm][k])
        numeric = charge["settings"]["0,1"]
        for k, key in enumerate(("real", "imag")):
            for source in ("ordinary", "pairing"):
                assert numeric[source]["connected_39_from_0"][k] in parse_interval(
                    cert[key][source]
                )
