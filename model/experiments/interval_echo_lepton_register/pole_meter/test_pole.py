import json
from math import pi

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_lepton_register import measure as marked
from . import measure as meter


def routed_split(psi, trail, exposure):
    """Move the two endpoints together, pulse the actual 99-column restriction, undo."""
    swaps = list(zip(trail[:-2], trail[1:-1]))
    for a, b in swaps:
        psi[[a, b]] = psi[[b, a]]
    psi[trail[-2:]] = expm(-1j * exposure * marked.seam_block(1)) @ psi[trail[-2:]]
    for a, b in reversed(swaps):
        psi[[a, b]] = psi[[b, a]]


def test_material_ramsey_and_failed_return_on_the_full_board():
    board, _, _, modes, ref, trails = meter.setup()
    h, e0, _, _, _ = meter.spectrum(1)
    cells = np.flatnonzero(np.linalg.norm(ref.reshape(-1, 12), axis=1) > 1e-10)
    assert len(cells) == 1
    block = slice(12 * cells[0], 12 * cells[0] + 12)
    sharp = trails[-1][0]
    # Synthesize the remote mode from scratch on an equal material chart.
    compact = meter.compact_preparation(np.real_if_close(ref[block]).real, sharp)
    analyzer = marked.reader(marked.addresses((0, 1)), 2 / 9, 1)
    for k in range(4):
        source = trails[k][0]
        psi = np.eye(1, board.n, k=source, dtype=complex).ravel()
        routed_split(psi, trails[k], pi / 2)
        emit = analyzer.conj().T if k < 3 else compact
        psi[:12] = emit @ psi[:12]
        psi[block] = compact @ psi[block]
        assert abs(np.vdot(modes[:, k], psi)) ** 2 > 0.499999999
        time = 9 / 16
        evolved = expm_multiply(-1j * time * h, psi)
        evolved[:12] = emit.conj().T @ evolved[:12]
        evolved[block] = compact.conj().T @ evolved[block]
        analytic = np.vdot(
            modes[:, k], expm_multiply(-1j * time * h, modes[:, k])
        ) * np.exp(1j * time * e0)
        for phi in [0, pi / 2, 0.773]:
            trial = evolved.copy()
            marked.phase(trial[:12], source, phi, 1)
            routed_split(trial, trails[k], 3 * pi / 2)
            probs = abs(trial) ** 2
            first, last = trails[k][0], trails[k][-1]
            terminal = [
                probs[first],
                probs[last],
                probs.sum() - probs[first] - probs[last],
            ]
            np.testing.assert_allclose(
                terminal, meter.outcomes(analytic, phi), atol=2e-12
            )


def test_spectral_gap_character_weight_and_phase_scan():
    h, ref_energy, energies, _, weights = meter.spectrum(1)
    board, projection, denominator, _, _, _ = meter.setup()
    native = meter.material_reference.anchor.defect_generator(board, -1, marked.KAPPA, 0)
    np.testing.assert_allclose(native.toarray(), h, atol=2e-13)
    cert = meter.geometry_certificate(h, projection, denominator)
    assert cert["minimum_integer_dominance_margin"] > 0
    assert cert["compact_dimension"] == 40
    assert abs(ref_energy - 35) < 2e-12
    assert np.count_nonzero(abs(energies - 35) < 1e-9) == 40
    outside = energies[abs(energies - 35) > 1e-9]
    assert min(abs(outside - 35)) > cert["proved_separation"]
    a, _, plus, minus = marked.addresses((0, 1))
    for theta in [0, 2 / 9, pi / 3, 0.712, -0.19]:
        modes = (
            marked.reader(marked.addresses((0, 1)), theta, 1)
            .conj()
            .T[:, [a, minus, plus]]
        )
        projected = projection[:12, :12] @ modes / denominator
        values = np.sum(abs(projected) ** 2, axis=0)
        expected = (2 - np.cos(2 * theta + 4 * pi * np.arange(3) / 3)) / 9
        np.testing.assert_allclose(values, expected, atol=3e-13)
        assert np.min(values) >= 1 / 9 - 1e-12
        if theta == 2 / 9:
            np.testing.assert_allclose(values, weights[:3])
    # A different free-contact coefficient changes the pole and centroids.
    changed_h, changed_ref, _, _, _ = meter.spectrum(3)
    assert abs(changed_ref - 41) < 2e-12
    np.testing.assert_allclose(changed_h - h, board.laplacian().toarray())


def test_acquisitions_with_certified_resolution_and_shared_reference():
    stored = json.loads((meter.HERE / "results.json").read_text())
    fresh = meter.measure()
    np.testing.assert_allclose(
        stored["reference_pole_weights"], fresh["reference_pole_weights"], atol=1e-12
    )
    cal = fresh["independent_six_seam_store_calibration"]
    assert (
        abs(cal["frequency_per_clock_gap"] - 35) < 2 * cal["frequency_standard_error"]
    )
    exact_weight = np.array(stored["reference_pole_weights"])
    shot = stored["simultaneous_99_percent_sampling_radius"]
    tail = stored["off_pole_filter_bound"]
    for idx, (bank, replay) in enumerate(
        zip(
            stored["randomized_wait_banks"], fresh["randomized_wait_banks"], strict=True
        )
    ):
        for j, (channel, expected) in enumerate(
            zip(bank["source_channels"], replay["source_channels"], strict=True)
        ):
            assert sum(channel["port_counts"]) == meter.SHOTS
            np.testing.assert_allclose(
                channel["mean_from_law"], expected["mean_from_law"], atol=2e-12
            )
            assert abs(channel["mean_from_counts"] - expected["mean_from_law"]) < shot
            if idx == 0:
                assert (
                    exact_weight[j] - 1e-11
                    <= expected["mean_from_law"]
                    <= exact_weight[j] + tail
                )
            else:
                assert abs(channel["mean_from_counts"]) < 2 * shot
    assert stored["ideal_energy_reader"]["rigorous_common_pole_ceiling"] < 8 / 9
    assert 0.57 < fresh["ideal_energy_reader"]["equal_prior_guess"] < 0.58
    for setting in stored["finite_wait_moments"]:
        assert np.all(
            np.array(setting["short_wait_port_counts"]).sum(axis=1) == meter.SHOTS
        )
        moment = setting["reference_energy"] + np.array(
            setting["relative_slope_from_nominal_wait"]
        )
        assert np.all(
            abs(moment - setting["generator_centroid"])
            < shot / meter.SMALL_TIME + np.asarray(setting["finite_wait_bias_bounds"])
        )
