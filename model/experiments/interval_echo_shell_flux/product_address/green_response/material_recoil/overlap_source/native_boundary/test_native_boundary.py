"""Independent checks of the complete seam, board boundary and one-stop transport."""

from fractions import Fraction as F
import json
import numpy as np
from scipy.linalg import expm

from . import budget, patch, precision, source, transport


def test_native_swap_work_and_free_spectrum():
    row = source.run_case(8)[1]
    assert row["memory_swap"] == dict(
        before=F(17, 2),
        after=F(51, 4),
        restored=F(17, 2),
        work_on_two_contacts=F(17, 4),
    )
    h, levels = source.energy_projectors(8)
    v = source.LAW.echo()[:, 1]
    rho = np.diag(abs(v) ** 2)
    reduced = rho.reshape(11, 9, 11, 9)
    scrambled = np.kron(
        reduced.trace(axis1=1, axis2=3), reduced.trace(axis1=0, axis2=2)
    )
    u = expm(-0.437j * h / 4)
    for state, name in ((rho, "retained_law"), (scrambled, "independent_law")):
        evolved = u @ state @ u.conj().T
        for e, proj, den in levels:
            np.testing.assert_allclose(
                np.trace(evolved @ (proj / den)),
                float(row["sharp"][name][e]),
                atol=1e-13,
            )


def test_boundary_for_random_noncompact_preparations():
    rng = np.random.default_rng(271129)
    n, edges = patch.board(1)
    local = source.run_case(8)[1]
    for p in [0, 0.1, 0.5, 0.7, 1]:
        occupied = rng.random(n) < p
        record = patch.weigh(occupied, edges, local)
        a, b = edges.T // 12
        count = np.count_nonzero(occupied[a] ^ occupied[b])
        assert record["loss_of_pairing_gain"] == count * F(17, 8)
        assert record["gain_relative_to_initial"] <= F(17, 260)


def test_transport_counts_and_receiver_contrast():
    data = json.loads((source.HERE / "precision.json").read_text())
    for row in data["source_counts"]:
        assert (
            sum(row[k] for k in ("coincident", "separate", "one_lost", "both_lost"))
            == data["trials_per_origin"]
        )
    for panel in data["panels"]:
        effect = panel["signed_per_expected_release"]
        error = panel["simultaneous_99_half_width"]
        assert abs(effect - panel["predicted"]) < error and effect + error < 0
    origin = (4, 4, 3)
    modes, eig = transport.spectral(10, 1 / 8)
    for alpha in transport.ALPHAS:
        p = precision.terminal(origin, alpha)
        np.testing.assert_allclose(
            p, alpha * transport.green(origin, 10, modes, eig, alpha), atol=9e-17
        )
        assert (p * p).sum() < transport.whole_space_bound(alpha, 1 / 8)


def test_aperture_units_and_multiple_layers():
    result = budget.run()
    for ratio, num in zip(
        result["favorable_lensing_ratios"],
        result["minimum_maximal_layers"],
        strict=True,
    ):
        assert (num - 1) * 17 / 260 < ratio <= num * 17 / 260
    assert result["planck_minimum_maximal_layers"] == 82


def test_native_free_curvature_and_local_persistence():
    from . import free

    _, edges, _, _ = free.layout(2)
    board = free.Carrier("closed test", 24, edges)
    record = free.run_case(board)
    assert record["enclosed_energy_second_derivative_at_zero"] == [
        "-85/32",
        "0",
        "0",
        "0",
    ]
    assert 1.90 < record["times"]["infinite_mean"]["enclosed_energy"][0] < 1.91
    # Independent small-time integration checks the double commutator.
    parts = free.operators(board)
    h, _, _, _, _, _, starts, weights = parts
    dt = 1e-4
    evolved = expm(-1j * dt * h.toarray()) @ starts
    value = free.read(evolved, weights, parts)["enclosed_energy"][0]
    np.testing.assert_allclose(2 * (value - 17 / 8) / dt**2, -85 / 32, atol=2e-5)
