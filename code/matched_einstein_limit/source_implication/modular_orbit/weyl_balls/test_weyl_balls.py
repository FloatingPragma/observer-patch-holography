import json
from pathlib import Path

import numpy as np
from scipy import linalg

from .axial_response import audit as axial_audit
from .boosts import hopping, tilted_reading
from .measure import make_board, read_ball, regional_covariances, sites_in_ball
from .packet import audit as packet_audit


def test_ground_projector_and_calibrated_hopping_on_all_modes():
    # Full-space diagonalization checks the Fourier implementation and the twist sign.
    n, length, floor = 4, 5.0, 0.01
    board = make_board(n, length, 2.5, floor)
    sites = np.indices((n,) * 3).reshape(3, -1).T
    cov, delta = regional_covariances(board, sites)
    action = hopping(n, board.spacing).toarray()
    energies, modes = linalg.eigh(action)
    measured = modes[:, energies < 0]
    assert np.linalg.norm(cov - (floor * np.eye(len(cov)) + (1 - 2 * floor) * measured @ measured.conj().T)) < 1e-10
    assert np.linalg.norm(action @ delta - delta @ action) < 1e-10
    assert abs(np.trace(action @ delta) / length**3 - board.energy) < 1e-12
    assert linalg.eigvalsh(cov + delta)[0] >= floor - 1e-12


def test_shifted_ball_has_the_same_pair_marginals_and_entropy_bound():
    board = make_board(16, 8, 2.5, 0.002)
    outer = sites_in_ball(board.spacing, 1.21, (0, 0, 0))
    inner = sites_in_ball(board.spacing, 0.87, (0.14, 0, 0))
    lookup = {tuple(site): j for j, site in enumerate(outer)}
    modes = np.array([2 * lookup[tuple(site)] + i for site in inner for i in (0, 1)])
    for full, direct in zip(regional_covariances(board, outer), regional_covariances(board, inner)):
        assert np.array_equal(full[np.ix_(modes, modes)], direct)
    row = read_ball(board, 0.87, (0.14, 0, 0))
    assert 0 < row["relative_entropy"] <= row["relative_entropy_upper_bound"]
    assert row["largest_onsite_source_change"] < 1e-15
    assert row["weighted_hopping_energy"] > 0
    # A tracial product reference has zero first-law response to the identical delta.
    assert row["product_reference_entropy_change"] < row["product_reference_modular_change"] == 0
    assert row["entropy_change"] > 0


def test_boosts_restrict_the_same_covariance_and_resolve_signed_momentum():
    board = make_board(8, 4, 2.5, 0.002)
    h = hopping(8, board.spacing)
    radius = 0.72
    rest = tilted_reading(board, h, radius, (1, 0, 0), 0.0)
    fixed = read_ball(board, radius)
    for key in ["weighted_hopping_energy", "entropy_change", "modular_change"]:
        assert np.isclose(rest[key], fixed[key], rtol=1e-9, atol=1e-12)
    plus, minus = [tilted_reading(board, h, radius, (s, 0, 0), 0.2) for s in (1, -1)]
    for key in ["weighted_hopping_energy", "entropy_change", "continuum_radiation_prediction"]:
        assert plus[key] > minus[key]
    for row in (rest, plus, minus):
        assert 0 <= row["relative_entropy"] <= row["relative_entropy_upper_bound"]


def test_archived_frames_span_symmetric_tensors():
    data = json.loads(Path(__file__).with_name("boosts.json").read_text())
    contractions = []
    for frame in data["frames"]:
        n, eta = np.array(frame["direction"]), frame["rapidity"]
        velocity = np.r_[np.cosh(eta), -n * np.sinh(eta)]
        contractions.append(np.outer(velocity, velocity)[np.triu_indices(4)])
    assert np.linalg.matrix_rank(contractions) == 10


def test_free_charged_packet_arrives_with_a_measured_energy_response():
    launch_first, launch_last, receive_first, receive_last, *central = packet_audit()["readings"]
    for key in ["occupation_excess", "weighted_hopping_energy", "modular_change", "entropy_change"]:
        assert launch_first[key] > launch_last[key] > 0
        assert receive_last[key] > receive_first[key] > 0
    for row in [launch_first, launch_last, receive_first, receive_last, *central]:
        assert 0 <= row["relative_entropy"] <= row["relative_entropy_upper_bound"]


def test_equal_rest_energy_fails_axial_geometry_and_the_free_signal_is_invariant():
    data = axial_audit()
    matter = data["weyl_action"]
    rods = data["sheet_einstein_over_coupling_at_equal_rest_energy"]
    assert np.isclose(matter["rest"], data["cutoff"]**4 / (96 * np.pi**2))
    assert np.isclose(rods["rest"], matter["rest"], rtol=1e-14, atol=0)
    assert rods["axial_null"] == [0, 0] and min(matter["axial_null"]) > 0
    assert np.isclose(np.mean(matter["axial_null"]), 4 * matter["rest"] / 3)
    assert np.isclose(matter["trace"], 0, atol=1e-14)
    assert np.isclose(rods["trace"], -2 * matter["rest"])
    assert np.isclose(rods["mean_axial_boost_over_rest"], 1)
    assert np.isclose(matter["mean_axial_boost_over_rest"], 1 + 4 * np.sinh(data["rapidity"])**2 / 3)
    pulse = data["finite_signal_on_two_preparations"]
    assert pulse["region_entropy_change"] > 0 and pulse["maximum_source_change_in_lesser_signal"] > 0.01
    assert max(pulse["maximum_vacuum_sum_error"], pulse["maximum_source_sum_error"]) < 1e-14
