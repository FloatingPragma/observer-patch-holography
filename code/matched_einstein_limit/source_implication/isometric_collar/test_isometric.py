import numpy as np
import pytest

from . import certify, measure


def test_exact_affine_witness_and_sharp_preparation():
    data = certify.certificate()
    assert data["exact_rank"] == 12
    assert data["certified_entropy_slope_floor"] > .113
    assert data["exact_kernel_tangent_l1_correction"] < 1.1e-6
    assert data["reachable_entropy_slope_floor"] > 2e-14
    assert data["two_menu_left_inverse_error"] < 1e-5
    assert max(data["known_left_inverse_error"],data["known_right_inverse_error"]) < 1e-5
    assert len(set(data["span_address_ids"])) == 16
    assert data["held_even_entropy_slope_floor"] > .01
    assert data["held_never_dephased_slope_floor"] > .0024
    assert data["calibrated_total_floor"] > .001


def test_optimized_unconditional_source_and_all_six_menus():
    data = measure.run()
    source = data["operational_half_tick_six"]
    assert abs(max(map(abs,source["biases"]))) < 1
    assert source["entropy_derivative"] == pytest.approx(2.535e-9,rel=.02)
    assert source["entropy_centered"] == pytest.approx(source["entropy_derivative"],rel=.001)
    assert source["max_range_shift"] < 2e-15
    assert source["cap_four_max_range_shift"] > 1e-13
    assert len(source["addresses"]) == 50
    assert min(source["marginal_norms"]) > 1e-9
    assert source["min_output_eigenvalue"] > 6e-9
    assert data["held_selected_source"]["same_direction_even_slope"] == pytest.approx(.01703385, abs=1e-7)
    for row in data["menus"]:
        assert row["rank"] == (9 if row["fixed"] else 15 if row["menu"]=="all" else 12)


def test_certified_fixed_point_representation_and_weights():
    from scipy.sparse import csr_matrix

    h = csr_matrix([[4,-2],[-2,1]],dtype=np.int64)
    z = np.array([[1,-3],[-4,7]],dtype=object)*certify.UNIT+5
    high = np.asarray(z//certify.BASE,dtype=np.int64)
    low = np.asarray(z%certify.BASE,dtype=np.int64)
    a,b = certify.divided_product(h,high,low,5)
    np.testing.assert_array_equal(a.astype(object)*certify.BASE+b, h.toarray().astype(object)@z//32)
