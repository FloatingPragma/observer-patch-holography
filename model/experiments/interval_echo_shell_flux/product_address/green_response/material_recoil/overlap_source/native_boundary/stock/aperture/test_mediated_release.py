import json
import numpy as np
from scipy.linalg import expm

from model.executable_model.carrier import Carrier
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    mediated_release as lab,
    medium_field,
    clock_to_flight,
    undisturbed_medium,
)


def test_clock_moments_match_direct_charge_preserving_two_arm_interference():
    board = Carrier("test", 5, [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)])
    mask = [(0, 1)]
    h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
    hb, _, _ = lab.cold.motion(board, mask, 1, lab.K)
    x = 1
    insert = np.zeros((len(i), board.n))
    rows = np.flatnonzero(i == x)
    insert[rows, j[rows]] = 1
    for t in (lab.CLOCK_STEP, lab.CLOCK_STEP / 2):
        exact = (
            expm(1j * t * hb.toarray())
            @ insert.T
            @ expm(-1j * t * h.toarray())
            @ insert
        )
        moments = lab.ramsey_moments(h, hb, i, j, x)
        np.testing.assert_allclose(
            np.polynomial.polynomial.polyval(t, moments),
            exact.diagonal(),
            rtol=0,
            atol=2e-13,
        )
    assert (
        medium_field.edge_certificate(board, mask, i, j, h)["integer_edge_equalities"]
        == 20
    )


def test_exterior_current_for_a_complex_coherent_bath_without_monitoring():
    board = Carrier("chain", 5, [(0, 1), (1, 2), (2, 3), (3, 4)])
    mask = [(0, 1)]
    hb, _, _ = lab.cold.motion(board, mask, 1, lab.K)
    h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
    rng = np.random.default_rng(184)
    factor = rng.normal(size=(5, 2)) + 1j * rng.normal(size=(5, 2))
    factor /= np.linalg.norm(factor)
    rho = factor @ factor.conj().T
    p = rho.diagonal().real
    field = 4 * lab.K * (board.adjacency() @ p + p)
    comm = hb @ rho - (hb @ rho).conj().T
    wave = factor[j] / 2
    hw = h @ wave
    local = 2 * (abs(hw) ** 2 - (wave.conj() * (h @ hw)).real).sum(axis=1)
    shells = lab.cold.distances(board, mask[0])
    measured = np.bincount(shells[i], weights=local).cumsum()
    for r, (u, v) in enumerate(undisturbed_medium.cuts(board, shells), 1):
        np.testing.assert_allclose(
            (sum(field[u] - field[v]) + sum(comm[u, v].real)) / 4,
            measured[r],
            atol=2e-12,
        )


def test_ordinary_clock_accuracy_and_source_only_uncertainty():
    data = json.loads((lab.HERE / "mediated_release.json").read_bytes())
    decoded = lab.decode(data)
    for d in decoded["clocks"]:
        a, b = d["frequency_change"]
        x, y = d["from_ordinary"]
        assert a <= y and x <= b
    x = medium_field.predict(data)["exterior_cumulative_change"]
    exact = (
        np.asarray(data["exact_cumulative_accelerations"]["thermal"])
        - data["exact_cumulative_accelerations"]["blank"]
    )
    assert all(
        a - 1e-12 <= v <= b + 1e-12 for (a, b), v in zip(x, exact[1:], strict=True)
    )
    assert x[0][0] > 0 and x[1][1] < 0
    flight = clock_to_flight.run()
    a, b = flight["exterior_native_minus_scalar_interval"][0]
    err = exact[1] - flight["comparisons"]["ordinary_calibrated"]["acceleration"][1]
    assert 0 < a < err < b


def test_separate_copy_source_reads_and_driven_conservation():
    from fractions import Fraction
    from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
        driven_medium,
    )

    medium = json.loads((lab.HERE / "mediated_release.json").read_bytes())
    free = json.loads((lab.HERE / "undisturbed_medium.json").read_bytes())
    read = undisturbed_medium.read(free, medium)
    np.testing.assert_allclose(
        undisturbed_medium.shell_field(medium),
        free["stationary_shell_field_from_counts"],
    )
    for kind, rows in read.items():
        actual = (
            np.array(free["arms"]["written_" + kind]["acceleration"])
            - free["arms"]["blank_" + kind]["acceleration"]
        )
        assert all(
            v["from_counts"][0] <= actual[v["radius"]] <= v["from_counts"][1]
            for v in rows
        )
        assert rows[0]["observed_enclosed_changes"][0][0] > 0
    assert (
        read["aligned"][0]["from_counts"][1] < read["equilibrium"][0]["from_counts"][0]
    )
    driven = json.loads((lab.HERE / "driven_medium.json").read_bytes())
    controls = driven_medium.decode(driven)
    assert driven["distant_affected"] == driven["distant_enrolled"] == 61
    for v, truth in zip(
        driven["edge"], driven["stationary_profile_change_per_injection"], strict=True
    ):
        lo, hi = controls[str(v)]["distant_occupation_per_injection_change"]
        assert lo < float(Fraction(truth)) < hi and lo * hi > 0
    assert (
        Fraction(driven["stationary_profile_change_per_injection"][0])
        * Fraction(driven["stationary_profile_change_per_injection"][1])
        < 0
    )
