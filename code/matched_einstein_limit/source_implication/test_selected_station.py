import numpy as np
from . import selected_station, full_station


def test_binary_station_can_be_tuned_and_both_constrained_curvatures_have_opposite_signs():
    r = selected_station.run()
    assert r["endpoint_mismatches"][0] < 0 < r["endpoint_mismatches"][1]
    assert r["calibrated"]["ordinary_eta"] > 0
    np.testing.assert_allclose(r["calibrated"]["fraction"], 1, atol=1e-9)
    assert r["constrained"]["fixed_volume_gradient_norm"] < 1e-17
    low, high = r["constrained"]["eigenvalues"]
    assert low < -1e-9 and high > 5e-12
    assert r["held_at_frozen_calibration"][0]["transfer_at_frozen_eta"] > 1.5


def test_complete_head_maximum_and_transfer_at_separately_calibrated_gain():
    r = full_station.run()
    for row in r["stations"]:
        lower, upper = row["bracket_residuals"]
        assert lower * upper < 0 and row["eta"] > 0
        assert row["gradient_norm"] < 1e-12 and max(row["hessian_eigenvalues"]) < -0.05
        np.testing.assert_allclose(
            row["hessian_eigenvalues"], row["hessian_smaller_step"], rtol=1e-7
        )
    held = r["held_reference_and_eta"]
    assert all(row["pairing_fraction"] < 0.12 for row in held[:-1])
    for point in r["selected_finite_neighbours"]:
        assert (
            abs(point["volume_error"]) < 1e-13
            and point["change_generalized_entropy"] < -1e-7
        )
    prior = r["independent_ordinary_calibration"]
    np.testing.assert_allclose(prior["reference_rates"], [31 / 64, 33 / 64])
    assert min(prior["selected_rates"]) > 0 and max(prior["selected_rates"]) < 1
    assert max(prior["stability"]["hessian_eigenvalues"]) < -0.05
    for point in prior["held_finite_differences"]:
        assert abs(point["fraction"] - 1) < 0.001


def test_right_material_rail_distinguishes_selected_state_from_free_source():
    import json
    from pathlib import Path
    from model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry import (
        measure as phase,
    )

    path = Path(__file__).with_name("full_algebra") / "results.json"
    packet = json.loads(path.read_text())
    blocks = np.array(packet["x_real"]) + 1j * np.array(packet["x_imag"])
    product, change = full_station.head_source()
    # Only the left reversal parity is dephased, so the native right transfer is
    # the sign-even contact in every arm. Sum over the unread left subsystem.
    gate = phase.marked_gate(1, 8)[:, :, 0, :]

    def right_pointer(state):
        reduced = np.trace(state.reshape(2, 2, 2, 2), axis1=0, axis2=2)
        return float(np.trace(gate[0] @ reduced @ gate[0].conj().T).real)

    np.testing.assert_allclose(
        right_pointer(product), right_pointer(product + change), atol=2e-13
    )
    np.testing.assert_allclose(right_pointer(product + change), 0.02010697, atol=5e-9)
    np.testing.assert_allclose(right_pointer(blocks.sum(axis=0)), 0.50797378, atol=5e-9)
    assert right_pointer(blocks.sum(axis=0)) - right_pointer(product) > 0.4878


def test_entropy_increase_and_area_deficit_on_both_held_source_transfers():
    from .signed_transfers import run

    data = run()
    ordinary, collar = data["tests"]
    assert ordinary["matched"][3] > 0.29 and abs(collar["matched"][3]) < 1e-13
    for arm in (ordinary, collar):
        a, v, s, _ = arm["matched"]
        assert a < 0 and abs(v) < 1e-15 and s > 0.23
        np.testing.assert_allclose(-s / a, data["eta"], rtol=1e-9)
        residuals = []
        for sample in arm["finite"]:
            assert abs(sample["volume_difference"]) < 4e-14
            assert (
                sample["right_marginal_error"] < 1e-13
                and sample["transfer_gate_error"] < 1e-13
            )
            residuals.append(abs(sample["residual"]))
        assert residuals[1] < 0.52 * residuals[0] and residuals[2] < 0.52 * residuals[1]
