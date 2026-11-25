import json
from fractions import Fraction

import numpy as np

from . import angular
from .measure import run


def test_triplet_projection_is_preserved_by_the_native_waits_and_rotations():
    cell = angular.cap.cell()
    distance = angular.shortest_path(cell.adjacency(), directed=False, unweighted=True)
    # 4Q = A + B/sqrt(5). Test each rational coefficient separately.
    a = (distance == 0).astype(int) - (distance == 3)
    b = (distance == 1).astype(int) - (distance == 2)
    np.testing.assert_array_equal(5 * (a @ a) + b @ b, 20 * a)
    np.testing.assert_array_equal(a @ b + b @ a, 4 * b)
    assert np.trace(a) == 12 and np.trace(b) == 0
    h = angular.cap.native_energy(
        [angular.cap.seam_indices(e) for e in cell.seams], np.ones(30)
    )[0]
    twice = np.rint(2 * h[1:13, 1:13]).astype(int)
    np.testing.assert_array_equal(twice, 2 * h[1:13, 1:13])
    for coeff in (a, b):
        np.testing.assert_array_equal(coeff @ twice, twice @ coeff)
        for permutation in angular.modular.rotations():
            np.testing.assert_array_equal(
                coeff[np.ix_(permutation, permutation)], coeff
            )


def test_native_nonlinear_remainder_and_stationary_pole():
    data = run()
    native = data["written"]
    endpoints = [json.loads(value) for value in native["entropy_slopes"]]
    assert (
        endpoints[0][0] > 0.00330
        and endpoints[1] == [0.0, 0.0]
        and endpoints[2][1] < -0.00294
    )
    assert (
        Fraction(native["fisher_lower"])
        < Fraction(1, 324)
        < Fraction(native["fisher_upper"])
    )
    transfer = native["reversed_to_shared"]
    actual = json.loads(transfer["entropy_change"])[0]
    predicted = json.loads(transfer["linear_entropy_prediction"])[0]
    error = json.loads(transfer["relative_entropy_remainder"])[0]
    assert 0 < actual < 0.02 * predicted and abs(actual - (predicted - error)) < 4e-18
    for row in native["resolved_annuli"]:
        for arm in row["information_cost_reverse_fresh_shared"]:
            assert (
                abs(arm["mean_nats"] - arm["measured_nats"])
                <= 5 * arm["standard_error"]
            )
    cell = data["angular"]
    assert cell["mode_meter"]["rotations"] == 30
    assert cell["mode_meter"]["worst_rotation_error"] < 5e-13
    for arm in cell["sampled_classification"]:
        assert sum(arm["axial_transverse_other"]) == arm["trials"]
    regions = cell["nested_regions"]
    fractions = [r["full_information_fraction_home"] for r in regions]
    assert fractions == sorted(fractions) and abs(fractions[-1] - 1) < 1e-11
    assert (
        regions[0]["home"]["linear_entropy_slope"]
        > 0
        > regions[1]["home"]["linear_entropy_slope"]
    )
    assert abs(regions[-1]["home"]["linear_entropy_slope"]) < 1e-17
