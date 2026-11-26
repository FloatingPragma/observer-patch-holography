"""Regression checks for the recoded generator, its recoil and its decay controls."""

from functools import cache

import numpy as np
from numpy.testing import assert_allclose

from .measure import field_symbol, run
from .reconstruct import reconstruct


@cache
def receipt():
    return run()


def test_recoded_family_and_virtual_gap():
    result = receipt()
    locality = result["locality"]
    assert locality["analyzer_census"]["disjoint_tilings"] == 240
    assert locality["analyzer_census"]["largest_number_of_central_seams"] == 0
    assert locality["face_mass_seams"] == [[4, 5], [4, 11], [5, 11]]
    mass = result["spectrum"]
    assert max(mass["virtual_tau_at_allowed_mu"]) < mass["lab_windows"][1][0]
    assert_allclose(
        result["field"]["records"][0]["positive_energies"][:3],
        mass["rest_energies"],
        atol=2e-14,
    )


def test_wilson_corners_leave_common_rest_recoil():
    masses = np.array([0.01, 0.25, 0.6])
    pole_zero = np.linalg.eigvalsh(field_symbol(np.diag(masses), np.zeros(3), 1)[0])[
        6::2
    ]
    for corner in ((np.pi, 0, 0), (np.pi, np.pi, 0), (np.pi, np.pi, np.pi)):
        h, _, _ = field_symbol(np.diag(masses), np.array(corner), 1)
        assert_allclose(
            np.linalg.eigvalsh(h)[6::2] - pole_zero, 4 * np.count_nonzero(corner) ** 2
        )
    momentum = 1 / 512
    h, _, _ = field_symbol(np.diag(masses), np.array([momentum, 0, 0]), 1)
    recoil = (np.linalg.eigvalsh(h)[6::2] ** 2 - pole_zero**2) / momentum**2
    assert_allclose(recoil, 1, atol=2e-6, rtol=0)


def test_weak_first_pair_and_port_control():
    data = receipt()
    weak = data["weak"]
    assert np.argwhere(
        np.array(weak["widths_divided_by_coupling_squared"]) > 0
    ).tolist() == [[1, 0], [2, 0], [2, 1]]
    assert weak["cut_daughter_probability"] == 0
    last = weak["original_start_probabilities"]["banks"][-1]
    assert last["daughter_probabilities"][0] > 8e-6
    assert weak["maximum_bin_refinement_difference"] < 1e-9
    controls = data["locality"]["weak_current_controls"]
    assert_allclose(controls["independent_wrong_fraction"], 2 / 3)
    assert controls["coherent_wrong_fraction"] < 1e-14


def test_reconstruct_from_independently_numbered_outputs():
    data = receipt()
    observers = data["observer_records"]
    field = data["field"]
    assert (
        max(row["lower_bound"] for row in observers["native_frequency_check"]) > 0.009
    )
    geometry = observers["geometry"]
    permuted = dict(geometry, counts=[geometry["counts"][i] for i in (1, 2, 0, 3)])
    outcome = reconstruct(
        permuted,
        field["rest_interferometer"],
        field["recoil_interferometers"],
        observers["native_free_return"],
    )
    assert outcome["mapping"]["geometry_channel_for_each_formation_tag"] == (0, 1, 2)
    assert (
        outcome["common_speed_squared_interval"][0]
        < 1
        < outcome["common_speed_squared_interval"][1]
    )
