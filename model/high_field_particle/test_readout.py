from fractions import Fraction as F

import pytest

from .readout import creation_certificate, speed_bound, finite_preparation_example


def test_complete_local_charge_creation_is_even_in_full_fock_space():
    value = creation_certificate()
    assert value['full_Fock_completeness_and_grading_exact']
    assert value['disjoint_creators_commute_exactly']
    assert set(value['successful_charge_by_sublattice']) == {-1, 1}
    assert value['failure_and_out_of_code_records_retained']
    assert not value['success_renormalized']
    assert not value['single_physical_fermion_under_grading_preserving_dictionary']


@pytest.mark.parametrize('args', [(F(0), F(1), F(1)),
                                 (F(1), F(-1), F(1)),
                                 (1.0, F(1), F(1)),
                                 (F(1), F(1), True)])
def test_invalid_refinement_inputs_rejected(args):
    with pytest.raises(ValueError):
        speed_bound(*args)


def test_bounded_rest_energy_cannot_hide_a_clock_conversion():
    # Independently form a source J,B and clock factor. The source speed
    # obeys the claimed bound whenever E_min>=B/8 and J/B<=lambda_max.
    a, j, b, clock = F(1, 64), F(3), F(120), F(2, 5)
    energy = clock*b/8
    assert a*clock*72*j == speed_bound(a, j/b, energy)
    assert speed_bound(a/2, j/b, energy) == speed_bound(a, j/b, energy)/2


def test_ground_plan_does_not_certify_an_unknown_infinite_volume_threshold():
    value = finite_preparation_example()
    assert value['finite_ground_bound_applies']
    assert not value['value_1_over_100_certified_below_Yarotsky_threshold']
    assert not value['native_tape_executed']
