"""Recheck the formed, product and ordinary responses and the measurement ledger."""

import json

import numpy as np
from numpy.testing import assert_allclose

from ....null_operator import measure as material
from ..reciprocal_screen.live_reference import apparatus
from . import bounds, certificate, held, measure


def data(name):
    return json.loads((measure.HERE / name).read_text())


def test_marginal_matching_ordinary_energy_and_no_veto():
    records = data("results.json")
    rho = np.asarray(records["state_real"]) + 1j * np.asarray(records["state_imag"])
    source = measure.source_entropies(rho)
    assert_allclose(source["mutual_information"], 0.04172495550544, atol=1e-11)
    assert_allclose(source["at_product"], 0, atol=1e-13)
    assert_allclose(source["at_formed"], -0.293848431429, atol=1e-10)
    meter = material.parked_terms()
    assert [8 * v + t for v, t in zip(meter["contrast"], meter["exchange"])] == [
        0,
        8,
        8,
        9,
        9,
    ]
    bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
    for r in (0, 1):
        for s in (0, 1):
            assert_allclose(
                apparatus.bell_effect(r, s), np.outer(bell, bell), atol=1e-12
            )


def test_angular_pencil_basis_and_independent_charge_calibration():
    cert = data("certificate.json")["charges"]
    for charge, sign in (("-1", -1), ("1", 1)):
        table = np.asarray(cert[charge]["(1, 38)"]["probabilities_times_1e18"]) * 1e-18
        ref, half, full, _, ordinary = table
        f = measure.jacobian_derivative(full, 3 * full - 4 * half + ref)
        o = measure.jacobian_derivative(full, 2 * (ordinary - full))
        dt = 1e-5
        step = 3 * full - 4 * half + ref
        jac = measure.camera.image
        assert_allclose(
            (jac(full + dt * step)[0] - jac(full - dt * step)[0]) / (2 * dt),
            f,
            atol=2e-10,
        )
        p = np.linalg.solve(o, f)
        assert np.sign(np.linalg.det(p)) == sign
        left = np.array([[0.2, -2.4], [1.1, 0.4]])
        right = np.array([[-0.6, 0.1], [1.8, 1.1]])
        transformed = np.linalg.solve(left @ o @ right, left @ f @ right)
        assert_allclose(transformed, np.linalg.solve(right, p @ right))
        assert np.sign(np.linalg.det(-7 * transformed)) == sign
        assert_allclose(np.linalg.eigvals(np.linalg.solve(o.T, f.T)).sum(), np.trace(p))
        # Removing all parity mismatches makes the ordinary image tangent vanish.
        assert_allclose(measure.jacobian_derivative(full, -2 * full), 0, atol=1e-14)


def test_certified_populations_and_two_independent_count_runs():
    exact = data("certificate.json")
    first = data("results.json")
    second = data("held.json")
    low = data("bounds.json")
    for q, views in exact["charges"].items():
        for edge, records in views.items():
            intervals = certificate.decode(records["probabilities_times_1e18"])
            for ray, expected in (
                ("forward", intervals),
                ("reverse", intervals.transpose(0, 2, 1)),
            ):
                for raw, get_bands, starts in (
                    (
                        first["charges"][q][edge][ray]["counts"],
                        bounds.ranges,
                        measure.SHOTS,
                    ),
                    (
                        second["charges"][q][edge][ray]["counts"],
                        held.bands,
                        held.STARTS,
                    ),
                ):
                    assert (np.array(raw).sum(axis=2) == starts).all()
                    for b, e in zip(
                        get_bands(raw).flat, expected[: len(raw)].flat, strict=True
                    ):
                        assert b.a <= e.a and e.b <= b.b
            if edge == "(1, 38)":
                assert (
                    bounds.analyze(intervals, formed=True) == records["formed_source"]
                )
                for ray in ("forward", "reverse"):
                    a = second["charges"][q][edge][ray]["at_formed"]
                    first_pass = low[q][edge][ray]
                    if q == "-1":
                        assert a["relative_determinant"][1] < -1.27
                        assert a["energy_from_area"][0] > 29.5
                        assert first_pass["energy_from_area"][0] > 1.4
                    else:
                        assert a["relative_determinant"][0] > 3.8
                        assert a["energy_from_area"][1] < -1.94
                        assert first_pass["energy_from_area"][1] < -9
    for edge in ("(1, 38)", "(22, 44)"):
        assert (
            exact["charges"]["-1"][edge]["formed_source"]["relative_determinant"][1] < 0
        )
        assert (
            exact["charges"]["1"][edge]["formed_source"]["relative_determinant"][0] > 0
        )


def test_native_one_pixel_mutual_coin_and_integer_contraction():
    records = data("results.json")
    cert = data("certificate.json")["charges"]
    board = measure.camera.front.optics.pairing.formation.federation(1)
    route = measure.camera.front.optics.pairing.formation.region(board)
    start, end = measure.mounts(board, route)[1, 38]
    rho = np.asarray(records["state_real"]) + 1j * np.asarray(records["state_imag"])
    tau = np.kron(*measure.camera.front.optics.pairing.formation.marginals(rho))
    waves = measure.kernels(board, route, [start[0], end[0]], -1)
    p = measure.responses(waves, rho, tau)[:, 1, 0]
    even = waves[0, 1, 0]
    unit = certificate.DEN
    real = np.vectorize(lambda x: int(round(x * unit)), otypes=[object])
    integer = certificate.optical(
        (real(even.real), real(even.imag)), (real(rho.real), real(rho.imag))
    )
    assert_allclose(integer / (4 * unit**4), p[2], atol=3e-18)
    expected = np.asarray(cert["-1"]["(1, 38)"]["probabilities_times_1e18"])[:4, 0, 0]
    assert_allclose(p, expected * 1e-18, atol=3e-15, rtol=0)
